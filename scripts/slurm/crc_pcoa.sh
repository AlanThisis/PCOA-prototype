#!/usr/bin/env bash
#SBATCH --job-name=crc-pcoa
#SBATCH --partition=short
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --output=slurm-crc-pcoa-%j.out

set -euo pipefail

# Run one CRC cross-study PCoA level with fixed resources.
# Required submission variable: LEVEL=full, sub10, sub25, or sub50.
# Example:
#   sbatch --job-name=crc-pcoa-sub10 --export=ALL,LEVEL=sub10 \
#     scripts/slurm/crc_pcoa.sh

REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
cd "${REPO_ROOT}"

LEVEL="${LEVEL:-}"
CONDA_ENV="${CONDA_ENV:-qiime2-amplicon-2024.10}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
DATA_ROOT="${DATA_ROOT:-data/fastq_data/CRC_cross_study}"
METADATA="${METADATA:-${DATA_ROOT}/crc_cross_study_metadata.tsv}"
GG2_DIR="${GG2_DIR:-data/gg2}"
TRIM_LENGTH="${TRIM_LENGTH:-120}"
MIN_READS="${MIN_READS:-0}"
SAMPLING_DEPTH="${SAMPLING_DEPTH:-1000}"
JOB_TAG="${SLURM_JOB_ID:-manual-$(date +%Y%m%d%H%M%S)}"
RUN_DIR="${RUN_DIR:-runs/crc-cross-study/${LEVEL}-${JOB_TAG}}"
RESUME="${RESUME:-0}"
FULL_STAGE_ROOT="${FULL_STAGE_ROOT:-work/crc-cross-study-full-amplicon-inputs}"
START_EPOCH="$(date +%s)"

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

validate_level() {
  case "${LEVEL}" in
    full|sub10|sub25|sub50) ;;
    *)
      echo "LEVEL must be one of: full, sub10, sub25, sub50" >&2
      exit 2
      ;;
  esac
  if [[ "${RESUME}" != "0" && "${RESUME}" != "1" ]]; then
    echo "RESUME must be 0 or 1" >&2
    exit 2
  fi
}

validate_selected_inputs() {
  python - "${METADATA}" "${DATA_ROOT}" "${LEVEL}" <<'PY'
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

metadata_path = Path(sys.argv[1])
data_root = Path(sys.argv[2])
level = sys.argv[3]
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
duplicates = sorted(sample_id for sample_id, count in Counter(sample_ids).items() if count > 1)
if not sample_ids or any(not sample_id for sample_id in sample_ids):
    raise SystemExit("Metadata contains an empty sample-id")
if duplicates:
    raise SystemExit(f"Duplicate metadata sample-id values: {duplicates}")

observed_counts = Counter(row["study"].strip() for row in rows)
if dict(observed_counts) != expected_counts:
    raise SystemExit(
        f"Unexpected metadata study counts: {dict(observed_counts)}; "
        f"expected {expected_counts}"
    )

expected_files = defaultdict(set)
for row in rows:
    study = row["study"].strip()
    sample_id = row["sample-id"].strip()
    filename = row["fastq_filename"].strip()
    if row["library_strategy"].strip() != "AMPLICON":
        raise SystemExit(f"Non-amplicon metadata row for {sample_id}")
    if filename != f"{sample_id}_1.fastq.gz":
        raise SystemExit(f"FASTQ filename mismatch for {sample_id}: {filename!r}")
    expected_files[study].add(filename)

missing = []
empty = []
for study, filenames in expected_files.items():
    source_dir = data_root / study / ("full" if level == "full" else level)
    if not source_dir.is_dir():
        raise SystemExit(f"Missing input directory: {source_dir}")
    for filename in filenames:
        path = source_dir / filename
        if not path.is_file():
            missing.append(str(path))
        elif path.stat().st_size == 0:
            empty.append(str(path))

    observed = {path.name for path in source_dir.glob("*_1.fastq.gz")}
    extras = sorted(observed - filenames)
    # The raw full download intentionally contains non-amplicon runs. The
    # subsampled directories, however, must contain exactly this cohort.
    if level != "full" and extras:
        raise SystemExit(
            f"Found {len(extras)} unexpected FASTQs in {source_dir}: {extras[:20]}"
        )
    print(
        f"{study}/{level}: selected={len(filenames)}, present={len(filenames) - len([p for p in missing if f'/{study}/' in p])}, "
        f"all_forward_files={len(observed)}, ignored_extras={len(extras) if level == 'full' else 0}"
    )

if missing:
    raise SystemExit(f"Missing {len(missing)} selected FASTQs: {missing[:20]}")
if empty:
    raise SystemExit(f"Found {len(empty)} empty selected FASTQs: {empty[:20]}")
PY
}

stage_full_inputs() {
  local study stage_dir regular_count nested_count
  for study in ERP005534 PRJEB46665; do
    stage_dir="${FULL_STAGE_ROOT}/${study}"
    if [[ -L "${stage_dir}" ]]; then
      echo "Refusing to use symlinked staging directory: ${stage_dir}" >&2
      exit 1
    fi
    mkdir -p "${stage_dir}"
    regular_count="$(find "${stage_dir}" -mindepth 1 -maxdepth 1 -type f | wc -l | tr -d ' ')"
    nested_count="$(find "${stage_dir}" -mindepth 1 -maxdepth 1 -type d | wc -l | tr -d ' ')"
    if [[ "${regular_count}" != "0" || "${nested_count}" != "0" ]]; then
      echo "Refusing to rebuild staging directory containing regular files/directories: ${stage_dir}" >&2
      exit 1
    fi
    find "${stage_dir}" -mindepth 1 -maxdepth 1 -type l -delete
  done

  python - "${METADATA}" "${DATA_ROOT}" "${FULL_STAGE_ROOT}" <<'PY'
import csv
import sys
from pathlib import Path

metadata_path = Path(sys.argv[1])
data_root = Path(sys.argv[2])
stage_root = Path(sys.argv[3])
with metadata_path.open(newline="", encoding="utf-8-sig") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

for row in rows:
    study = row["study"].strip()
    filename = row["fastq_filename"].strip()
    source = (data_root / study / "full" / filename).resolve(strict=True)
    destination = stage_root / study / filename
    destination.symlink_to(source)

for study, expected in (("ERP005534", 129), ("PRJEB46665", 109)):
    observed = len(list((stage_root / study).glob("*_1.fastq.gz")))
    print(f"{study}/full staging: {observed}/{expected} selected amplicon FASTQs")
    if observed != expected:
        raise SystemExit(f"Full-input staging count mismatch for {study}")
PY
}

format_seconds() {
  local seconds="$1"
  printf "%02d:%02d:%02d" \
    "$((seconds / 3600))" \
    "$(((seconds % 3600) / 60))" \
    "$((seconds % 60))"
}

validate_level

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID:-manual}"
echo "Level: ${LEVEL}"
echo "Repo: ${REPO_ROOT}"
echo "Run directory: ${RUN_DIR}"
echo "Conda env: ${CONDA_ENV}"
echo "Partition: ${SLURM_JOB_PARTITION:-short}"
echo "Threads: ${THREADS}"
echo "Memory: ${SLURM_MEM_PER_NODE:-96G}"
echo "Resume: ${RESUME}"

if command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base)"
  # shellcheck source=/dev/null
  source "${CONDA_BASE}/etc/profile.d/conda.sh"
  conda activate "${CONDA_ENV}"
fi

require_command python
require_command qiime
require_command deblur
require_file "${METADATA}"
require_file "${GG2_DIR}/2024.09.backbone.full-length.fna.qza"
require_file "${GG2_DIR}/2024.09.phylogeny.id.nwk.qza"
require_file "src/run_pipeline.py"
require_dir "${DATA_ROOT}/ERP005534/full"
require_dir "${DATA_ROOT}/PRJEB46665/full"

echo "Python: $(command -v python)"
echo "QIIME: $(command -v qiime)"
echo "Deblur: $(command -v deblur)"
echo "Validating QIIME Greengenes2 plugin..."
qiime greengenes2 --help >/dev/null

echo "Validating metadata-selected ${LEVEL} inputs..."
validate_selected_inputs

if [[ "${LEVEL}" == "full" ]]; then
  echo "Staging only metadata-selected amplicon FASTQs from the mixed full download..."
  stage_full_inputs
  ERP_INPUT="${FULL_STAGE_ROOT}/ERP005534"
  PRJEB_INPUT="${FULL_STAGE_ROOT}/PRJEB46665"
else
  ERP_INPUT="${DATA_ROOT}/ERP005534/${LEVEL}"
  PRJEB_INPUT="${DATA_ROOT}/PRJEB46665/${LEVEL}"
fi

PIPELINE_ARGS=(
  --study "ERP005534=${ERP_INPUT}"
  --study "PRJEB46665=${PRJEB_INPUT}"
  --metadata "${METADATA}"
  --color-by disease_status
  --color-by diagnosis_detail
  --color-by study
  --run-dir "${RUN_DIR}"
  --trim-length "${TRIM_LENGTH}"
  --min-reads "${MIN_READS}"
  --sampling-depth "${SAMPLING_DEPTH}"
  --gg2-dir "${GG2_DIR}"
  --threads "${THREADS}"
)
if [[ "${RESUME}" == "1" ]]; then
  PIPELINE_ARGS+=(--resume)
fi

python src/run_pipeline.py "${PIPELINE_ARGS[@]}"

END_EPOCH="$(date +%s)"
ELAPSED="$((END_EPOCH - START_EPOCH))"
echo "Job finished: $(date)"
echo "Total elapsed: $(format_seconds "${ELAPSED}") (${ELAPSED} seconds)"
echo "Results: ${RUN_DIR}/results"
echo "Timings: ${RUN_DIR}/timings"
