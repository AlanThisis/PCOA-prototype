#!/usr/bin/env bash
#SBATCH --job-name=mmc20-download
#SBATCH --partition=short
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=0-19%5
#SBATCH --output=slurm-mmc20-download-%A_%a.out

set -euo pipefail

# Submit this file from the repository root:
#   sbatch scripts/slurm/mmc20_download.sh
#
# Each array task handles one project. At most five projects run concurrently,
# and each project downloads up to four forward FASTQs concurrently.

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
PYTHON_BIN="${PYTHON_BIN:-python}"
WORKERS="${SLURM_CPUS_PER_TASK:-4}"

ACCESSIONS=(
  PRJNA950484
  PRJNA646468
  PRJDB9293
  PRJNA227062
  PRJNA734525
  PRJEB4927
  PRJNA321051
  PRJNA632856
  PRJNA780512
  PRJNA647539
  PRJNA726058
  PRJNA495633
  PRJNA533177
  PRJNA1148754
  PRJNA674379
  PRJNA1062161
  PRJNA997108
  PRJNA505133
  PRJEB13307
  PRJEB17696
)

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "This script must be submitted with sbatch as a SLURM array." >&2
  exit 2
fi

if ((SLURM_ARRAY_TASK_ID < 0 || SLURM_ARRAY_TASK_ID >= ${#ACCESSIONS[@]})); then
  echo "Invalid array index: ${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

cd "${REPO_ROOT}"

ACCESSION="${ACCESSIONS[SLURM_ARRAY_TASK_ID]}"
OUTPUT_DIR="data/fastq_data/${ACCESSION}/full"
START_EPOCH="$(date +%s)"
STAGGER_SECONDS="$(((SLURM_ARRAY_TASK_ID % 5) * 5))"

echo "Accession: ${ACCESSION}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Host: $(hostname)"
echo "Started: $(date -Is)"
echo "Output: ${OUTPUT_DIR}"
echo "Workers: ${WORKERS}"
echo "Metadata request stagger: ${STAGGER_SECONDS} seconds"

# Avoid making five ENA metadata requests in the same second when a new wave of
# array tasks starts. The Python client also retries transient ENA failures.
sleep "${STAGGER_SECONDS}"

"${PYTHON_BIN}" src/download_ena_amplicon.py "${ACCESSION}" \
  --output-dir "${OUTPUT_DIR}" \
  --workers "${WORKERS}"

END_EPOCH="$(date +%s)"
ELAPSED_SECONDS="$((END_EPOCH - START_EPOCH))"

echo "Finished: $(date -Is)"
echo "SLURM task elapsed seconds: ${ELAPSED_SECONDS}"
echo "Persistent timing summary: ${OUTPUT_DIR}/ena_amplicon_forward_summary.tsv"
echo "Per-run timing/status: ${OUTPUT_DIR}/ena_amplicon_forward_status.tsv"
