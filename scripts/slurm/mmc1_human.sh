#!/usr/bin/env bash
# Submit and execute the human-only MMC1 download, subsampling, and PCoA DAG.
# Run this orchestrator from the repository root:
#   bash scripts/slurm/mmc1_human.sh

set -euo pipefail

SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${SCRIPT_PATH}")/../.." && pwd)}}"
MODE="${MODE:-submit}"
METADATA_DIR="${METADATA_DIR:-data/mmc1_16s_cleaning}"
PROJECTS_FILE="${PROJECTS_FILE:-${METADATA_DIR}/projects.txt}"
RUN_MANIFEST="${RUN_MANIFEST:-${METADATA_DIR}/approved_runs.tsv}"
METADATA="${METADATA:-${METADATA_DIR}/cleaned_metadata.csv}"
RUN_ROOT="${RUN_ROOT:-runs/mmc1_human}"
GG2_DIR="${GG2_DIR:-data/gg2}"
PCOA_CONDA_ENV="${PCOA_CONDA_ENV:-qiime2-amplicon-2024.10}"
SUBSAMPLE_CONDA_ENV="${SUBSAMPLE_CONDA_ENV:-cli-tools}"
SEED="${SEED:-11}"
TRIM_LENGTH="${TRIM_LENGTH:-120}"
SAMPLING_DEPTH="${SAMPLING_DEPTH:-1000}"

activate_environment() {
  local environment_name="$1"
  local conda_base
  conda_base="$(conda info --base)"
  # shellcheck source=/dev/null
  source "${conda_base}/etc/profile.d/conda.sh"
  conda activate "${environment_name}"
}

project_for_task() {
  local task_id="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
  sed -n "$((task_id + 1))p" "${PROJECTS_FILE}"
}

preflight() {
  cd "${REPO_ROOT}"
  for command_name in python conda sbatch; do
    command -v "${command_name}" >/dev/null 2>&1 || {
      echo "Missing required command: ${command_name}" >&2
      return 1
    }
  done
  for required_file in \
    "${PROJECTS_FILE}" \
    "${RUN_MANIFEST}" \
    "${METADATA}" \
    src/download_ena_amplicon.py \
    src/subsample_fastq.py \
    src/run_pipeline.py \
    "${GG2_DIR}/2024.09.backbone.full-length.fna.qza" \
    "${GG2_DIR}/2024.09.phylogeny.id.nwk.qza"; do
    [[ -s "${required_file}" ]] || {
      echo "Missing or empty required file: ${required_file}" >&2
      return 1
    }
  done
  conda env list | awk '{print $1}' | grep -Fx "${PCOA_CONDA_ENV}" >/dev/null || {
    echo "Missing conda environment: ${PCOA_CONDA_ENV}" >&2
    return 1
  }
  conda env list | awk '{print $1}' | grep -Fx "${SUBSAMPLE_CONDA_ENV}" >/dev/null || {
    echo "Missing conda environment: ${SUBSAMPLE_CONDA_ENV}" >&2
    return 1
  }

  python - "${PROJECTS_FILE}" "${RUN_MANIFEST}" "${METADATA}" <<'PY'
import csv
import sys
from collections import Counter
from pathlib import Path

projects_path, approved_path, metadata_path = map(Path, sys.argv[1:])
projects = [line.strip() for line in projects_path.read_text().splitlines() if line.strip()]
if not projects or len(projects) != len(set(projects)):
    raise SystemExit("projects.txt is empty or contains duplicates")
with approved_path.open(newline="", encoding="utf-8-sig") as handle:
    approved = list(csv.DictReader(handle, delimiter="\t"))
with metadata_path.open(newline="", encoding="utf-8-sig") as handle:
    metadata = list(csv.DictReader(handle))
runs = [row["run_accession"] for row in approved]
sample_ids = [row["sample-id"] for row in metadata]
if len(runs) != len(set(runs)) or len(sample_ids) != len(set(sample_ids)):
    raise SystemExit("approved runs or metadata sample IDs are not unique")
if set(runs) != set(sample_ids):
    raise SystemExit("approved runs and cleaned metadata sample IDs differ")
counts = Counter(row["project_accession"] for row in approved)
if set(projects) != set(counts) or any(counts[p] == 0 for p in projects):
    raise SystemExit("projects.txt and approved run projects differ")
if any(row.get("study_hmr") != "Human" or row.get("host_species_resolved") != "Human" for row in metadata):
    raise SystemExit("non-human rows survived cleaned metadata")
print(f"Metadata preflight: {len(projects)} projects, {len(runs)} unique approved runs")
PY
  df -h "${REPO_ROOT}"
}

submit_jobs() {
  [[ -z "${SLURM_JOB_ID:-}" ]] || {
    echo "Run MODE=submit on the login node, not through sbatch." >&2
    exit 2
  }
  preflight
  mkdir -p "${RUN_ROOT}/logs"
  local project_count array_range
  project_count="$(awk 'NF {count++} END {print count+0}' "${PROJECTS_FILE}")"
  array_range="0-$((project_count - 1))%5"
  export REPO_ROOT METADATA_DIR PROJECTS_FILE RUN_MANIFEST METADATA RUN_ROOT
  export GG2_DIR PCOA_CONDA_ENV SUBSAMPLE_CONDA_ENV SEED TRIM_LENGTH SAMPLING_DEPTH

  local download_job subsample_job full_job sub50_job sub25_job sub10_job
  download_job="$(sbatch --parsable --chdir="${REPO_ROOT}" \
    --partition=short --time=23:59:00 --cpus-per-task=4 --mem=8G \
    --array="${array_range}" --job-name=mmc1-download \
    --output="${RUN_ROOT}/logs/download-%A_%a.out" \
    --export=ALL,MODE=download "${SCRIPT_PATH}")"
  download_job="${download_job%%;*}"

  full_job="$(sbatch --parsable --chdir="${REPO_ROOT}" \
    --partition=highmem --time=14-00:00:00 --cpus-per-task=32 --mem=256G \
    --dependency="afterok:${download_job}" --job-name=mmc1-pcoa-full \
    --output="${RUN_ROOT}/logs/pcoa-full-%j.out" \
    --export=ALL,MODE=pcoa,LEVEL=full "${SCRIPT_PATH}")"
  full_job="${full_job%%;*}"

  subsample_job="$(sbatch --parsable --chdir="${REPO_ROOT}" \
    --partition=short --time=23:59:00 --cpus-per-task=8 --mem=16G \
    --array="${array_range}" --dependency="afterok:${download_job}" \
    --job-name=mmc1-subsample \
    --output="${RUN_ROOT}/logs/subsample-%A_%a.out" \
    --export=ALL,MODE=subsample "${SCRIPT_PATH}")"
  subsample_job="${subsample_job%%;*}"

  sub50_job="$(submit_pcoa sub50 32 192G "${subsample_job}")"
  sub25_job="$(submit_pcoa sub25 24 128G "${subsample_job}")"
  sub10_job="$(submit_pcoa sub10 16 96G "${subsample_job}")"

  python - "${RUN_ROOT}/submission.json" "${download_job}" "${subsample_job}" \
    "${full_job}" "${sub50_job}" "${sub25_job}" "${sub10_job}" \
    "${project_count}" <<'PY'
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

path = Path(sys.argv[1])
download, subsample, full, sub50, sub25, sub10, project_count = sys.argv[2:]
payload = {
    "submitted_utc": datetime.now(timezone.utc).isoformat(),
    "project_count": int(project_count),
    "jobs": {
        "download_array": {"id": download, "dependency": None},
        "full": {"id": full, "dependency": f"afterok:{download}"},
        "subsample_array": {"id": subsample, "dependency": f"afterok:{download}"},
        "sub50": {"id": sub50, "dependency": f"afterok:{subsample}"},
        "sub25": {"id": sub25, "dependency": f"afterok:{subsample}"},
        "sub10": {"id": sub10, "dependency": f"afterok:{subsample}"},
    },
}
path.write_text(json.dumps(payload, indent=2) + "\n")
PY
  echo "Submitted MMC1 human workflow:"
  echo "  download array: ${download_job} (${array_range})"
  echo "  full PCoA:      ${full_job} (afterok:${download_job})"
  echo "  subsample array:${subsample_job} (afterok:${download_job})"
  echo "  50% PCoA:       ${sub50_job} (afterok:${subsample_job})"
  echo "  25% PCoA:       ${sub25_job} (afterok:${subsample_job})"
  echo "  10% PCoA:       ${sub10_job} (afterok:${subsample_job})"
  echo "  record:         ${RUN_ROOT}/submission.json"
}

submit_pcoa() {
  local level="$1" cpus="$2" memory="$3" dependency="$4"
  local submission
  submission="$(sbatch --parsable --chdir="${REPO_ROOT}" \
    --partition=highmem --time=14-00:00:00 --cpus-per-task="${cpus}" --mem="${memory}" \
    --dependency="afterok:${dependency}" --job-name="mmc1-pcoa-${level}" \
    --output="${RUN_ROOT}/logs/pcoa-${level}-%j.out" \
    --export="ALL,MODE=pcoa,LEVEL=${level}" "${SCRIPT_PATH}")"
  echo "${submission%%;*}"
}

run_download() {
  cd "${REPO_ROOT}"
  local project stagger output_dir
  project="$(project_for_task)"
  [[ -n "${project}" ]] || { echo "No project for array task" >&2; exit 2; }
  stagger="$((SLURM_ARRAY_TASK_ID % 5 * 5))"
  output_dir="data/fastq_data/${project}/full"
  echo "Project: ${project}; stagger: ${stagger}s; started: $(date -Is)"
  sleep "${stagger}"
  python src/download_ena_amplicon.py "${project}" \
    --output-dir "${output_dir}" \
    --run-manifest "${RUN_MANIFEST}" \
    --workers 4 --retries 5 --timeout 120 \
    --retry-until-complete --max-runtime 84000
  [[ -s "${output_dir}/.ena_download_complete.json" ]] || {
    echo "Downloader returned without a completion marker: ${project}" >&2
    exit 1
  }
  echo "Completed and validated: ${project}; finished: $(date -Is)"
}

validate_subset() {
  local project="$1" output_dir="$2" marker="$3" percent="$4"
  python - "${RUN_MANIFEST}" "${project}" "${output_dir}" "${marker}" "${percent}" "${SEED}" <<'PY'
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

allowlist, project, output_dir, marker, percent, seed = sys.argv[1:]
with Path(allowlist).open(newline="", encoding="utf-8-sig") as handle:
    expected = {
        f"{row['run_accession']}_1.fastq.gz"
        for row in csv.DictReader(handle, delimiter="\t")
        if row["project_accession"] == project
    }
output = Path(output_dir)
actual = {path.name for path in output.glob("*.fastq.gz") if path.is_file()}
empty = sorted(path.name for path in output.glob("*.fastq.gz") if path.stat().st_size == 0)
if actual != expected or empty:
    raise SystemExit(
        f"subset validation failed for {project} {percent}%: "
        f"missing={len(expected-actual)}, extra={len(actual-expected)}, empty={len(empty)}"
    )
Path(marker).write_text(json.dumps({
    "project_accession": project,
    "percent": int(percent),
    "seed": int(seed),
    "fastq_count": len(actual),
    "validated_utc": datetime.now(timezone.utc).isoformat(),
}, indent=2) + "\n")
PY
}

run_subsample() {
  cd "${REPO_ROOT}"
  activate_environment "${SUBSAMPLE_CONDA_ENV}"
  command -v seqkit >/dev/null 2>&1 || { echo "seqkit is unavailable" >&2; exit 1; }
  local project input_dir level percent output_dir marker
  project="$(project_for_task)"
  input_dir="data/fastq_data/${project}/full"
  mkdir -p "${RUN_ROOT}/subsampling_timings"
  [[ -s "${input_dir}/.ena_download_complete.json" ]] || {
    echo "Validated full download marker missing for ${project}" >&2
    exit 1
  }
  for percent in 50 25 10; do
    level="sub${percent}"
    output_dir="data/fastq_data/${project}/${level}"
    marker="${output_dir}/.subsampling_complete.json"
    if [[ -s "${marker}" ]] && validate_subset "${project}" "${output_dir}" "${marker}" "${percent}"; then
      echo "Validated existing ${project} ${level}; skipping regeneration"
      continue
    fi
    python src/subsample_fastq.py \
      --input-dir "${input_dir}" --output-dir "${output_dir}" \
      --percent "${percent}" --seed "${SEED}" \
      --timings-tsv "${RUN_ROOT}/subsampling_timings/${project}_${level}.tsv"
    validate_subset "${project}" "${output_dir}" "${marker}" "${percent}"
    echo "Completed ${project} ${level}: $(date -Is)"
  done
}

run_pcoa() {
  cd "${REPO_ROOT}"
  activate_environment "${PCOA_CONDA_ENV}"
  local level="${LEVEL:?LEVEL is required}" run_dir="${RUN_ROOT}/${LEVEL}"
  local project fastq_dir expected actual marker
  local -a pipeline_args=()
  while IFS= read -r project; do
    [[ -n "${project}" ]] || continue
    fastq_dir="data/fastq_data/${project}/${level}"
    if [[ "${level}" == "full" ]]; then
      marker="${fastq_dir}/.ena_download_complete.json"
    else
      marker="${fastq_dir}/.subsampling_complete.json"
    fi
    [[ -s "${marker}" ]] || { echo "Missing completion marker: ${marker}" >&2; exit 1; }
    expected="$(awk -F '\t' -v p="${project}" 'NR>1 && $1==p {n++} END {print n+0}' "${RUN_MANIFEST}")"
    actual="$(find "${fastq_dir}" -maxdepth 1 -type f -name '*_1.fastq.gz' | wc -l | tr -d ' ')"
    [[ "${expected}" == "${actual}" ]] || {
      echo "FASTQ count mismatch for ${project}/${level}: expected ${expected}, found ${actual}" >&2
      exit 1
    }
    pipeline_args+=(--study "${project}=${fastq_dir}")
  done < "${PROJECTS_FILE}"

  pipeline_args+=(
    --run-dir "${run_dir}"
    --metadata "${METADATA}"
    --color-by environment_harmonized
    --color-by disease_harmonized
    --color-by source_study
    --trim-length "${TRIM_LENGTH}"
    --min-reads 0
    --sampling-depth "${SAMPLING_DEPTH}"
    --gg2-dir "${GG2_DIR}"
    --threads "${SLURM_CPUS_PER_TASK}"
  )
  if [[ -s "${run_dir}/run_manifest.json" && -s "${run_dir}/run_state.json" ]]; then
    pipeline_args+=(--resume)
  elif [[ -e "${run_dir}/run_manifest.json" || -e "${run_dir}/run_state.json" ]]; then
    echo "Incomplete resume metadata under ${run_dir}; refusing an ambiguous restart" >&2
    exit 1
  fi

  echo "Starting ${level} PCoA: $(date -Is); CPUs=${SLURM_CPUS_PER_TASK}; depth=${SAMPLING_DEPTH}"
  python src/run_pipeline.py "${pipeline_args[@]}"
  for required_output in \
    distance_matrix_unweighted_unifrac.tsv \
    pcoa_coordinates_unweighted_unifrac.txt \
    pcoa_environment_harmonized.png \
    pcoa_disease_harmonized.png \
    pcoa_source_study.png; do
    [[ -s "${run_dir}/results/${required_output}" ]] || {
      echo "Missing final output: ${run_dir}/results/${required_output}" >&2
      exit 1
    }
  done
  echo "Completed ${level} PCoA: $(date -Is); results=${run_dir}/results"
}

cd "${REPO_ROOT}"
case "${MODE}" in
  submit) submit_jobs ;;
  preflight) preflight ;;
  download) run_download ;;
  subsample) run_subsample ;;
  pcoa) run_pcoa ;;
  *) echo "MODE must be submit, preflight, download, subsample, or pcoa" >&2; exit 2 ;;
esac
