#!/usr/bin/env bash
#SBATCH --job-name=crc-subsample
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=slurm-crc-subsample-%j.out

set -euo pipefail

# Rebuild the existing sub10, sub25, and sub50 directories for the complete
# metadata-selected 16S amplicon cohort (129 ERP005534 + 109 PRJEB46665 runs).
# Submit from the repository root:
#   sbatch scripts/slurm/crc_subsample.sh
# Optional environment override:
#   CONDA_ENV=qiime2-amplicon-2024.10 sbatch scripts/slurm/crc_subsample.sh

REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
cd "${REPO_ROOT}"

CONDA_ENV="${CONDA_ENV:-qiime2-amplicon-2024.10}"
DATA_ROOT="${DATA_ROOT:-data/fastq_data/CRC_cross_study}"
METADATA="${METADATA:-${DATA_ROOT}/crc_cross_study_metadata.tsv}"
SEED="${SEED:-11}"
JOB_TAG="${SLURM_JOB_ID:-manual-$(date +%Y%m%d%H%M%S)}"
TIMINGS_TSV="${DATA_ROOT}/crc_subsample_${JOB_TAG}_timings.tsv"
DETAIL_TIMINGS_DIR="${DATA_ROOT}/subsample_timings/${JOB_TAG}"
STAGING_PARENT="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
STAGING_ROOT=""

format_seconds() {
  local seconds="$1"
  printf "%02d:%02d:%02d" \
    "$((seconds / 3600))" \
    "$(((seconds % 3600) / 60))" \
    "$((seconds % 60))"
}

run_timed_step() {
  local step_name="$1"
  shift

  local start_epoch end_epoch elapsed exit_code
  local start_stamp end_stamp elapsed_hms
  start_epoch="$(date +%s)"
  start_stamp="$(date -Is)"
  echo "${step_name} started: ${start_stamp}"

  set +e
  "$@"
  exit_code="$?"
  set -e

  end_epoch="$(date +%s)"
  end_stamp="$(date -Is)"
  elapsed="$((end_epoch - start_epoch))"
  elapsed_hms="$(format_seconds "${elapsed}")"
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${step_name}" "${start_stamp}" "${end_stamp}" "${elapsed}" \
    "${elapsed_hms}" "${exit_code}" >> "${TIMINGS_TSV}"

  echo "${step_name} finished: ${end_stamp} (elapsed ${elapsed_hms}, exit ${exit_code})"
  if [[ "${exit_code}" -ne 0 ]]; then
    echo "Failed step: ${step_name}. Timing record: ${TIMINGS_TSV}" >&2
    exit "${exit_code}"
  fi
}

require_file() {
  local path="$1"
  if [[ ! -s "${path}" ]]; then
    echo "Missing or empty required file: ${path}" >&2
    exit 1
  fi
}

require_dir() {
  local path="$1"
  if [[ ! -d "${path}" ]]; then
    echo "Missing required directory: ${path}" >&2
    exit 1
  fi
}

require_command() {
  local name="$1"
  if ! command -v "${name}" >/dev/null 2>&1; then
    echo "Missing required command in PATH: ${name}" >&2
    exit 1
  fi
}

cleanup_staging() {
  if [[ -n "${STAGING_ROOT}" && -d "${STAGING_ROOT}" ]]; then
    rm -rf -- "${STAGING_ROOT}"
  fi
}
trap cleanup_staging EXIT

stage_selected_inputs() {
  STAGING_ROOT="$(mktemp -d "${STAGING_PARENT%/}/crc-amplicon-inputs.XXXXXX")"

  python - "${METADATA}" "${DATA_ROOT}" "${STAGING_ROOT}" <<'PY'
import csv
import sys
from collections import Counter
from pathlib import Path

metadata_path = Path(sys.argv[1])
data_root = Path(sys.argv[2])
staging_root = Path(sys.argv[3])
expected_counts = {"ERP005534": 129, "PRJEB46665": 109}
required_columns = {
    "sample-id",
    "study",
    "disease_status",
    "diagnosis_detail",
    "library_strategy",
    "fastq_filename",
}

with metadata_path.open(newline="", encoding="utf-8-sig") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    fields = set(reader.fieldnames or [])
    missing_columns = sorted(required_columns - fields)
    if missing_columns:
        raise SystemExit(f"Metadata missing required columns: {missing_columns}")
    rows = list(reader)

sample_ids = [row["sample-id"].strip() for row in rows]
if not sample_ids or any(not sample_id for sample_id in sample_ids):
    raise SystemExit("Metadata contains an empty sample-id")
duplicates = sorted(sample_id for sample_id, count in Counter(sample_ids).items() if count > 1)
if duplicates:
    raise SystemExit(f"Duplicate metadata sample-id values: {duplicates}")

observed_counts = Counter(row["study"].strip() for row in rows)
if dict(observed_counts) != expected_counts:
    raise SystemExit(
        f"Unexpected metadata study counts: {dict(observed_counts)}; "
        f"expected {expected_counts}"
    )

missing_sources = []
empty_sources = []
for row in rows:
    study = row["study"].strip()
    sample_id = row["sample-id"].strip()
    filename = row["fastq_filename"].strip()
    if study not in expected_counts:
        raise SystemExit(f"Unexpected study {study!r} for {sample_id}")
    if row["library_strategy"].strip() != "AMPLICON":
        raise SystemExit(f"Non-amplicon metadata row for {sample_id}")
    if filename != f"{sample_id}_1.fastq.gz":
        raise SystemExit(
            f"FASTQ filename mismatch for {sample_id}: {filename!r}"
        )

    source = data_root / study / "full" / filename
    if not source.is_file():
        missing_sources.append(str(source))
        continue
    if source.stat().st_size == 0:
        empty_sources.append(str(source))
        continue

    destination_dir = staging_root / study
    destination_dir.mkdir(parents=True, exist_ok=True)
    (destination_dir / filename).symlink_to(source.resolve())

if missing_sources:
    preview = "\n".join(missing_sources[:20])
    suffix = "\n..." if len(missing_sources) > 20 else ""
    raise SystemExit(
        f"Missing {len(missing_sources)} selected full FASTQ(s):\n{preview}{suffix}"
    )
if empty_sources:
    preview = "\n".join(empty_sources[:20])
    suffix = "\n..." if len(empty_sources) > 20 else ""
    raise SystemExit(
        f"Found {len(empty_sources)} empty selected full FASTQ(s):\n{preview}{suffix}"
    )

for study, expected in expected_counts.items():
    observed = len(list((staging_root / study).glob("*_1.fastq.gz")))
    print(f"{study}: staged {observed}/{expected} selected amplicon FASTQs")
    if observed != expected:
        raise SystemExit(f"Staging count mismatch for {study}")
PY
}

clear_previous_subsamples() {
  local study subset output_dir nested_count
  for study in ERP005534 PRJEB46665; do
    for subset in sub10 sub25 sub50; do
      output_dir="${DATA_ROOT}/${study}/${subset}"
      if [[ -L "${output_dir}" ]]; then
        echo "Refusing to clear symlinked output directory: ${output_dir}" >&2
        exit 1
      fi
      mkdir -p "${output_dir}"
      nested_count="$(find "${output_dir}" -mindepth 1 -maxdepth 1 -type d | wc -l | tr -d ' ')"
      if [[ "${nested_count}" != "0" ]]; then
        echo "Refusing to clear ${output_dir}; it contains nested directories." >&2
        exit 1
      fi
      find "${output_dir}" -mindepth 1 -maxdepth 1 \( -type f -o -type l \) -delete
      echo "Cleared previous files from ${output_dir}"
    done
  done
}

subsample_study() {
  local study="$1"
  local percent="$2"
  local subset="sub${percent}"
  local detail_tsv="${DETAIL_TIMINGS_DIR}/${study}_${subset}.tsv"

  python src/subsample_fastq.py \
    --input-dir "${STAGING_ROOT}/${study}" \
    --output-dir "${DATA_ROOT}/${study}/${subset}" \
    --percent "${percent}" \
    --seed "${SEED}" \
    --timings-tsv "${detail_tsv}"
}

verify_outputs() {
  python - "${METADATA}" "${DATA_ROOT}" <<'PY'
import csv
import sys
from collections import defaultdict
from pathlib import Path

metadata_path = Path(sys.argv[1])
data_root = Path(sys.argv[2])
expected = defaultdict(set)
with metadata_path.open(newline="", encoding="utf-8-sig") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        expected[row["study"].strip()].add(row["fastq_filename"].strip())

for study in ("ERP005534", "PRJEB46665"):
    for subset in ("sub10", "sub25", "sub50"):
        output_dir = data_root / study / subset
        observed = {path.name for path in output_dir.glob("*_1.fastq.gz")}
        missing = sorted(expected[study] - observed)
        extra = sorted(observed - expected[study])
        empty = sorted(path.name for path in output_dir.glob("*_1.fastq.gz") if path.stat().st_size == 0)
        print(
            f"{study}/{subset}: {len(observed)}/{len(expected[study])} FASTQs; "
            f"missing={len(missing)}, extra={len(extra)}, empty={len(empty)}"
        )
        if missing or extra or empty:
            raise SystemExit(
                f"Output validation failed for {study}/{subset}: "
                f"missing={missing[:10]}, extra={extra[:10]}, empty={empty[:10]}"
            )
PY
}

mkdir -p "${DETAIL_TIMINGS_DIR}"
printf "step\tstart\tend\telapsed_seconds\telapsed_hms\texit_code\n" > "${TIMINGS_TSV}"

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Repo: ${REPO_ROOT}"
echo "Conda env: ${CONDA_ENV}"
echo "Metadata: ${METADATA}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-4}"
echo "Memory: ${SLURM_MEM_PER_NODE:-24G}"
echo "Summary timings: ${TIMINGS_TSV}"
echo "Detailed timings: ${DETAIL_TIMINGS_DIR}"

if command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base)"
  # shellcheck source=/dev/null
  source "${CONDA_BASE}/etc/profile.d/conda.sh"
  conda activate "${CONDA_ENV}"
fi

require_command python
require_command seqkit
require_file "${METADATA}"
require_file "src/subsample_fastq.py"
require_dir "${DATA_ROOT}/ERP005534/full"
require_dir "${DATA_ROOT}/PRJEB46665/full"
require_dir "${STAGING_PARENT}"

echo "Python: $(command -v python)"
echo "SeqKit: $(command -v seqkit)"
seqkit version

# The destructive step occurs only after every selected full FASTQ has passed
# metadata, existence, and non-empty-file validation.
run_timed_step "1_stage_and_validate_selected_inputs" stage_selected_inputs
run_timed_step "2_clear_previous_subsamples" clear_previous_subsamples
run_timed_step "3_subsample_erp005534_10" subsample_study ERP005534 10
run_timed_step "4_subsample_prjeb46665_10" subsample_study PRJEB46665 10
run_timed_step "5_subsample_erp005534_25" subsample_study ERP005534 25
run_timed_step "6_subsample_prjeb46665_25" subsample_study PRJEB46665 25
run_timed_step "7_subsample_erp005534_50" subsample_study ERP005534 50
run_timed_step "8_subsample_prjeb46665_50" subsample_study PRJEB46665 50
run_timed_step "9_verify_subsample_outputs" verify_outputs

echo "Timings:"
cat "${TIMINGS_TSV}"
echo "Job finished: $(date)"
