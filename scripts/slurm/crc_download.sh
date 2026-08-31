#!/usr/bin/env bash
#SBATCH --job-name=crc-full-fastq
#SBATCH --partition=short
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=slurm-crc-full-fastq-%j.out

set -euo pipefail

# Submit from the repository root:
#   sbatch scripts/slurm/crc_download.sh
#
# The defaults download every ENA run in ERP005534 and PRJEB46665. Paired-end
# runs retain both mates, but downstream pipeline scripts select only
# *_1.fastq.gz. The job is resumable and retries incomplete runs.
#
# Optional overrides:
#   DOWNLOAD_WORKERS=2 sbatch scripts/slurm/crc_download.sh
#   ACCESSIONS="ERP005534 PRJEB46665" sbatch scripts/slurm/crc_download.sh
#   CONDA_ENV=qiime2-amplicon-2024.10 sbatch scripts/slurm/crc_download.sh

REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
CONDA_ENV="${CONDA_ENV:-qiime2-amplicon-2024.10}"
ACCESSIONS="${ACCESSIONS:-ERP005534 PRJEB46665}"
OUT_ROOT="${OUT_ROOT:-${REPO_ROOT}/data/fastq_data/CRC_cross_study}"
DOWNLOAD_WORKERS="${DOWNLOAD_WORKERS:-${SLURM_CPUS_PER_TASK:-4}}"
FASTQ_DL_ATTEMPTS="${FASTQ_DL_ATTEMPTS:-5}"
RUN_ATTEMPTS="${RUN_ATTEMPTS:-2}"
RETRY_SLEEP_SECONDS="${RETRY_SLEEP_SECONDS:-120}"

STATE_ROOT="${OUT_ROOT}/download_state"
MANIFEST_ROOT="${STATE_ROOT}/manifests"
RUN_STATE_ROOT="${STATE_ROOT}/runs"
STATUS_ROOT="${STATE_ROOT}/status"
RUN_MANIFEST="${STATE_ROOT}/run_manifest.tsv"
STATUS_TSV="${OUT_ROOT}/download_status.tsv"
SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"

export REPO_ROOT CONDA_ENV OUT_ROOT DOWNLOAD_WORKERS FASTQ_DL_ATTEMPTS
export RUN_ATTEMPTS RETRY_SLEEP_SECONDS STATE_ROOT MANIFEST_ROOT
export RUN_STATE_ROOT STATUS_ROOT RUN_MANIFEST STATUS_TSV SCRIPT_PATH

timestamp() {
  date "+%Y-%m-%dT%H:%M:%S%z"
}

require_command() {
  local name="$1"
  if ! command -v "${name}" >/dev/null 2>&1; then
    echo "Missing required command in PATH: ${name}" >&2
    exit 1
  fi
}

activate_conda_environment() {
  if command -v conda >/dev/null 2>&1; then
    local conda_base
    conda_base="$(conda info --base)"
    # shellcheck source=/dev/null
    source "${conda_base}/etc/profile.d/conda.sh"
    conda activate "${CONDA_ENV}"
  fi
}

sanitize_status_message() {
  printf "%s" "$1" | tr '\t\r\n' '   '
}

write_run_status() {
  local project="$1"
  local run_accession="$2"
  local status="$3"
  local attempts="$4"
  local start_stamp="$5"
  local end_stamp="$6"
  local elapsed_seconds="$7"
  local file_names="$8"
  local total_bytes="$9"
  local message="${10}"
  local status_file="${STATUS_ROOT}/${project}__${run_accession}.tsv"
  local status_tmp="${status_file}.tmp.$$"

  mkdir -p "${STATUS_ROOT}"
  message="$(sanitize_status_message "${message}")"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${project}" "${run_accession}" "${status}" "${attempts}" \
    "${start_stamp}" "${end_stamp}" "${elapsed_seconds}" \
    "${file_names}" "${total_bytes}" "${message}" > "${status_tmp}"
  mv "${status_tmp}" "${status_file}"
}

completed_marker_is_valid() {
  local marker="$1"
  local output_dir="$2"
  local file_name expected_bytes actual_bytes
  local found=0

  [[ -s "${marker}" ]] || return 1
  while IFS=$'\t' read -r file_name expected_bytes || [[ -n "${file_name}" ]]; do
    [[ -n "${file_name}" ]] || continue
    found=1
    [[ "${expected_bytes}" =~ ^[1-9][0-9]*$ ]] || return 1
    [[ -s "${output_dir}/${file_name}" ]] || return 1
    actual_bytes="$(wc -c < "${output_dir}/${file_name}" | tr -d '[:space:]')"
    [[ "${actual_bytes}" -eq "${expected_bytes}" ]] || return 1
  done < "${marker}"
  [[ "${found}" -eq 1 ]]
}

download_run() {
  local project="$1"
  local run_accession="$2"
  local output_dir="${OUT_ROOT}/${project}/full"
  local run_state="${RUN_STATE_ROOT}/${project}/${run_accession}"
  local marker="${run_state}/complete.txt"
  local log_file="${run_state}/fastq-dl.log"
  local start_epoch start_stamp end_epoch end_stamp elapsed
  local attempt=0
  local exit_code=1
  local file_path file_name file_names file_bytes total_bytes
  local -a fastq_files=()

  mkdir -p "${output_dir}" "${run_state}"
  start_epoch="$(date +%s)"
  start_stamp="$(timestamp)"

  if completed_marker_is_valid "${marker}" "${output_dir}"; then
    file_names=""
    total_bytes=0
    while IFS=$'\t' read -r file_name file_bytes || [[ -n "${file_name}" ]]; do
      [[ -n "${file_name}" ]] || continue
      if [[ -n "${file_names}" ]]; then
        file_names+=";"
      fi
      file_names+="${file_name}"
      total_bytes="$((total_bytes + file_bytes))"
    done < "${marker}"
    end_epoch="$(date +%s)"
    end_stamp="$(timestamp)"
    elapsed="$((end_epoch - start_epoch))"
    write_run_status "${project}" "${run_accession}" "skipped" "0" \
      "${start_stamp}" "${end_stamp}" "${elapsed}" "${file_names}" \
      "${total_bytes}" "Previously completed and files are present"
    echo "[${project}/${run_accession}] already complete; skipping"
    return 0
  fi

  if [[ -e "${marker}" ]]; then
    mv "${marker}" "${marker}.stale.$(date +%s)"
  fi

  printf "[%s/%s] started at %s\n" \
    "${project}" "${run_accession}" "${start_stamp}" | tee -a "${log_file}"

  while [[ "${attempt}" -lt "${RUN_ATTEMPTS}" ]]; do
    attempt="$((attempt + 1))"
    printf "[%s/%s] wrapper attempt %s of %s\n" \
      "${project}" "${run_accession}" "${attempt}" "${RUN_ATTEMPTS}" \
      | tee -a "${log_file}"

    set +e
    (
      cd "${run_state}"
      fastq-dl \
        --accession "${run_accession}" \
        --provider ena \
        --only-provider \
        --protocol https \
        --max-attempts "${FASTQ_DL_ATTEMPTS}" \
        --sleep 30 \
        --cpus 1 \
        --outdir "${output_dir}" \
        --prefix "${run_accession}"
    ) >> "${log_file}" 2>&1
    exit_code="$?"
    set -e

    if [[ -f "${output_dir}/${run_accession}-run-info.tsv" ]]; then
      mv "${output_dir}/${run_accession}-run-info.tsv" \
        "${run_state}/${run_accession}-run-info.tsv"
    fi

    if [[ "${exit_code}" -eq 0 ]]; then
      break
    fi

    echo "[${project}/${run_accession}] attempt ${attempt} failed with exit ${exit_code}" \
      | tee -a "${log_file}" >&2
    if [[ "${attempt}" -lt "${RUN_ATTEMPTS}" ]]; then
      echo "[${project}/${run_accession}] retrying in ${RETRY_SLEEP_SECONDS}s" \
        | tee -a "${log_file}"
      sleep "${RETRY_SLEEP_SECONDS}"
    fi
  done

  end_epoch="$(date +%s)"
  end_stamp="$(timestamp)"
  elapsed="$((end_epoch - start_epoch))"

  if [[ "${exit_code}" -ne 0 ]]; then
    write_run_status "${project}" "${run_accession}" "failed" "${attempt}" \
      "${start_stamp}" "${end_stamp}" "${elapsed}" "" "0" \
      "fastq-dl exited ${exit_code}; see ${log_file}"
    return 1
  fi

  # Normalize ENA single-end naming so downstream forward-read discovery finds it.
  if [[ -s "${output_dir}/${run_accession}.fastq.gz" ]]; then
    mv -f "${output_dir}/${run_accession}.fastq.gz" \
      "${output_dir}/${run_accession}_1.fastq.gz"
  fi

  for file_path in \
    "${output_dir}/${run_accession}_1.fastq.gz" \
    "${output_dir}/${run_accession}_2.fastq.gz"; do
    if [[ -s "${file_path}" ]]; then
      fastq_files+=("${file_path}")
    fi
  done

  if [[ "${#fastq_files[@]}" -eq 0 ]]; then
    write_run_status "${project}" "${run_accession}" "failed" "${attempt}" \
      "${start_stamp}" "${end_stamp}" "${elapsed}" "" "0" \
      "fastq-dl exited successfully but no non-empty FASTQ was found"
    return 1
  fi

  file_names=""
  total_bytes=0
  local marker_tmp="${marker}.tmp.$$"
  : > "${marker_tmp}"
  for file_path in "${fastq_files[@]}"; do
    file_name="$(basename "${file_path}")"
    file_bytes="$(wc -c < "${file_path}" | tr -d '[:space:]')"
    printf "%s\t%s\n" "${file_name}" "${file_bytes}" >> "${marker_tmp}"
    if [[ -n "${file_names}" ]]; then
      file_names+=";"
    fi
    file_names+="${file_name}"
    total_bytes="$((total_bytes + file_bytes))"
  done
  mv "${marker_tmp}" "${marker}"

  write_run_status "${project}" "${run_accession}" "downloaded" "${attempt}" \
    "${start_stamp}" "${end_stamp}" "${elapsed}" "${file_names}" \
    "${total_bytes}" "ENA download and MD5 validation completed"
  echo "[${project}/${run_accession}] complete (${file_names})"
}

if [[ "${1:-}" == "--download-run" ]]; then
  if [[ "$#" -ne 3 ]]; then
    echo "Internal usage: $0 --download-run PROJECT RUN_ACCESSION" >&2
    exit 2
  fi
  require_command fastq-dl
  download_run "$2" "$3"
  exit $?
fi

cd "${REPO_ROOT}"

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Repo: ${REPO_ROOT}"
echo "Conda env: ${CONDA_ENV}"
echo "Accessions: ${ACCESSIONS}"
echo "Output root: ${OUT_ROOT}"
echo "Concurrent downloads: ${DOWNLOAD_WORKERS}"
echo "Fastq-DL attempts per invocation: ${FASTQ_DL_ATTEMPTS}"
echo "Wrapper attempts per run: ${RUN_ATTEMPTS}"

if ! [[ "${DOWNLOAD_WORKERS}" =~ ^[1-9][0-9]*$ ]]; then
  echo "DOWNLOAD_WORKERS must be a positive integer: ${DOWNLOAD_WORKERS}" >&2
  exit 2
fi
if ! [[ "${FASTQ_DL_ATTEMPTS}" =~ ^[1-9][0-9]*$ ]]; then
  echo "FASTQ_DL_ATTEMPTS must be a positive integer: ${FASTQ_DL_ATTEMPTS}" >&2
  exit 2
fi
if ! [[ "${RUN_ATTEMPTS}" =~ ^[1-9][0-9]*$ ]]; then
  echo "RUN_ATTEMPTS must be a positive integer: ${RUN_ATTEMPTS}" >&2
  exit 2
fi

activate_conda_environment
require_command fastq-dl
require_command python
require_command wget
require_command xargs

echo "Fastq-DL: $(command -v fastq-dl)"
fastq-dl --version

mkdir -p "${OUT_ROOT}" "${MANIFEST_ROOT}" "${RUN_STATE_ROOT}" "${STATUS_ROOT}"
printf "project\trun_accession\n" > "${RUN_MANIFEST}.tmp"

metadata_failures=0
for project in ${ACCESSIONS}; do
  manifest_dir="${MANIFEST_ROOT}/${project}"
  manifest="${manifest_dir}/${project}-run-info.tsv"
  metadata_log="${manifest_dir}/fastq-dl-metadata.log"
  mkdir -p "${manifest_dir}" "${OUT_ROOT}/${project}/full"

  echo "[${project}] fetching run metadata"
  set +e
  (
    cd "${manifest_dir}"
    fastq-dl \
      --accession "${project}" \
      --provider ena \
      --only-provider \
      --protocol https \
      --max-attempts "${FASTQ_DL_ATTEMPTS}" \
      --sleep 30 \
      --only-download-metadata \
      --outdir "${manifest_dir}" \
      --prefix "${project}"
  ) >> "${metadata_log}" 2>&1
  metadata_exit="$?"
  set -e

  if [[ "${metadata_exit}" -ne 0 && -s "${manifest}" ]]; then
    echo "[${project}] metadata refresh failed; using cached ${manifest}" >&2
  elif [[ "${metadata_exit}" -ne 0 || ! -s "${manifest}" ]]; then
    echo "[${project}] metadata retrieval failed; see ${metadata_log}" >&2
    metadata_failures="$((metadata_failures + 1))"
    continue
  fi

  project_rows="${manifest_dir}/run-manifest-rows.tmp.$$"
  set +e
  PROJECT="${project}" MANIFEST="${manifest}" python - <<'PY' \
    > "${project_rows}"
import csv
import os
import sys
from collections import Counter
from pathlib import Path

project = os.environ["PROJECT"]
manifest = Path(os.environ["MANIFEST"])
with manifest.open(newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    if "run_accession" not in (reader.fieldnames or []):
        raise SystemExit(f"Missing run_accession column in {manifest}")
    rows = list(reader)

runs = sorted({row["run_accession"].strip() for row in rows})

runs = [run for run in runs if run]
if not runs:
    raise SystemExit(f"No run accessions found in {manifest}")
for run in runs:
    print(f"{project}\t{run}")

layouts = Counter(row.get("library_layout", "UNKNOWN") or "UNKNOWN" for row in rows)
total_bytes = sum(
    sum(int(value) for value in row.get("fastq_bytes", "").split(";") if value)
    for row in rows
)
layout_summary = ", ".join(f"{key}={layouts[key]}" for key in sorted(layouts))
print(
    f"[{project}] discovered {len(runs)} run(s); {layout_summary}; "
    f"advertised ENA payload={total_bytes / 10**9:.1f} GB",
    file=sys.stderr,
)
PY
  parse_exit="$?"
  set -e
  if [[ "${parse_exit}" -ne 0 ]]; then
    echo "[${project}] could not parse run metadata: ${manifest}" >&2
    metadata_failures="$((metadata_failures + 1))"
    mv "${project_rows}" "${project_rows}.failed"
    continue
  fi
  cat "${project_rows}" >> "${RUN_MANIFEST}.tmp"
  mv "${project_rows}" "${manifest_dir}/run-manifest-rows.tsv"
done

mv "${RUN_MANIFEST}.tmp" "${RUN_MANIFEST}"
run_count="$(( $(wc -l < "${RUN_MANIFEST}") - 1 ))"
echo "Total runs queued: ${run_count}"
echo "Run manifest: ${RUN_MANIFEST}"

worker_exit=0
if [[ "${run_count}" -gt 0 ]]; then
  set +e
  tail -n +2 "${RUN_MANIFEST}" \
    | xargs -P "${DOWNLOAD_WORKERS}" -n 2 \
      bash "${SCRIPT_PATH}" --download-run
  worker_exit="$?"
  set -e
fi

STATUS_ROOT="${STATUS_ROOT}" STATUS_TSV="${STATUS_TSV}" \
  RUN_MANIFEST="${RUN_MANIFEST}" python - <<'PY'
import csv
import os
from pathlib import Path

status_root = Path(os.environ["STATUS_ROOT"])
status_tsv = Path(os.environ["STATUS_TSV"])
run_manifest = Path(os.environ["RUN_MANIFEST"])
columns = [
    "project",
    "run_accession",
    "status",
    "wrapper_attempts",
    "start",
    "end",
    "elapsed_seconds",
    "files",
    "total_bytes",
    "message",
]

with run_manifest.open(newline="", encoding="utf-8") as handle:
    expected = {
        (row["project"], row["run_accession"])
        for row in csv.DictReader(handle, delimiter="\t")
    }

rows_by_run = {}
for path in status_root.glob("*.tsv"):
    with path.open(newline="", encoding="utf-8") as handle:
        values = next(csv.reader(handle, delimiter="\t"), None)
    if values is None or len(values) != len(columns):
        continue
    row = dict(zip(columns, values))
    key = (row["project"], row["run_accession"])
    if key in expected:
        rows_by_run[key] = row

for project, run_accession in expected - rows_by_run.keys():
    rows_by_run[(project, run_accession)] = {
        "project": project,
        "run_accession": run_accession,
        "status": "failed",
        "wrapper_attempts": "0",
        "start": "",
        "end": "",
        "elapsed_seconds": "0",
        "files": "",
        "total_bytes": "0",
        "message": "Worker exited without writing a status record",
    }

rows = list(rows_by_run.values())
rows.sort(key=lambda row: (row["project"], row["run_accession"]))
status_tsv.parent.mkdir(parents=True, exist_ok=True)
with status_tsv.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(
        handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)

counts = {}
for row in rows:
    counts[row["status"]] = counts.get(row["status"], 0) + 1
print(f"Status report: {status_tsv}")
print("Status counts: " + ", ".join(f"{key}={counts[key]}" for key in sorted(counts)))
PY

failed_runs="$(awk -F '\t' 'NR > 1 && $3 == "failed" {count++} END {print count+0}' "${STATUS_TSV}")"
echo "Metadata failures: ${metadata_failures}"
echo "Failed runs: ${failed_runs}"
echo "Job finished: $(date)"

if [[ "${metadata_failures}" -ne 0 || "${worker_exit}" -ne 0 || "${failed_runs}" -ne 0 ]]; then
  echo "Download job completed with failures. Rerun the same sbatch command to resume." >&2
  exit 1
fi

echo "All discovered runs are complete."
