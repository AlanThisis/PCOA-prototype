#!/usr/bin/env bash
set -euo pipefail

# Run this file with bash from the repository root. It submits five jobs:
#
#   full PCoA -------------------------------> runs independently
#
#   subsample ---> sub50 PCoA
#             |--> sub25 PCoA                 (three PCoAs run in parallel)
#             `--> sub10 PCoA
#
# Usage:
#   bash scripts/slurm/mmc20_analysis.sh
#
# Optional overrides, for example:
#   PCOA_CONDA_ENV=my-qiime-env SUBSAMPLE_CONDA_ENV=my-cli-env \
#     bash scripts/slurm/mmc20_analysis.sh

SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${SCRIPT_PATH}")/../.." && pwd)}"
PCOA_CONDA_ENV="${PCOA_CONDA_ENV:-qiime2-amplicon-2024.10}"
SUBSAMPLE_CONDA_ENV="${SUBSAMPLE_CONDA_ENV:-cli-tools}"
MODE="${MODE:-submit}"

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

# Allow recovery jobs to operate on a comma-separated subset without changing
# the default 20-project analysis set used by the PCoA jobs.
if [[ -n "${ACCESSIONS_CSV:-}" ]]; then
  IFS=',' read -r -a ACCESSIONS <<< "${ACCESSIONS_CSV}"
fi

activate_environment() {
  local environment_name="$1"
  if command -v conda >/dev/null 2>&1; then
    local conda_base
    conda_base="$(conda info --base)"
    # shellcheck source=/dev/null
    source "${conda_base}/etc/profile.d/conda.sh"
    conda activate "${environment_name}"
  fi
}

submit_jobs() {
  if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    echo "MODE=submit must be run on the login node with bash, not through sbatch." >&2
    exit 2
  fi
  if ! command -v sbatch >/dev/null 2>&1; then
    echo "sbatch is not available." >&2
    exit 1
  fi

  export REPO_ROOT PCOA_CONDA_ENV SUBSAMPLE_CONDA_ENV

  local full_submission full_job
  local subsample_submission subsample_job
  full_submission="$(
    sbatch --parsable \
      --chdir="${REPO_ROOT}" \
      --job-name=mmc20-pcoa-full \
      --cpus-per-task=16 \
      --mem=64G \
      --output=slurm-mmc20-pcoa-full-%j.out \
      --export=ALL,MODE=pcoa,LEVEL=full \
      "${SCRIPT_PATH}"
  )"
  full_job="${full_submission%%;*}"

  subsample_submission="$(
    sbatch --parsable \
      --chdir="${REPO_ROOT}" \
      --job-name=mmc20-subsample \
      --cpus-per-task=8 \
      --mem=32G \
      --output=slurm-mmc20-subsample-%j.out \
      --export=ALL,MODE=subsample \
      "${SCRIPT_PATH}"
  )"
  subsample_job="${subsample_submission%%;*}"

  local level submission job_id
  local sub50_job=""
  local sub25_job=""
  local sub10_job=""
  for level in sub50 sub25 sub10; do
    submission="$(
      sbatch --parsable \
        --chdir="${REPO_ROOT}" \
        --job-name="mmc20-pcoa-${level}" \
        --cpus-per-task=16 \
        --mem=64G \
        --output="slurm-mmc20-pcoa-${level}-%j.out" \
        --dependency="afterok:${subsample_job}" \
        --export="ALL,MODE=pcoa,LEVEL=${level}" \
        "${SCRIPT_PATH}"
    )"
    job_id="${submission%%;*}"
    case "${level}" in
      sub50) sub50_job="${job_id}" ;;
      sub25) sub25_job="${job_id}" ;;
      sub10) sub10_job="${job_id}" ;;
    esac
  done

  echo "Submitted five MMC20 jobs:"
  echo "  full PCoA:  ${full_job} (starts independently)"
  echo "  subsample:  ${subsample_job} (starts independently)"
  echo "  sub50 PCoA: ${sub50_job} (afterok:${subsample_job})"
  echo "  sub25 PCoA: ${sub25_job} (afterok:${subsample_job})"
  echo "  sub10 PCoA: ${sub10_job} (afterok:${subsample_job})"
}

run_subsampling() {
  local python_bin="${PYTHON_BIN:-python}"
  local max_parallel="${MAX_PARALLEL:-4}"
  local seed="${SEED:-11}"
  local timing_root="${TIMING_ROOT:-runs/mmc20/subsampling_timings}"

  if ! [[ "${max_parallel}" =~ ^[1-9][0-9]*$ ]]; then
    echo "MAX_PARALLEL must be a positive integer." >&2
    exit 2
  fi

  activate_environment "${SUBSAMPLE_CONDA_ENV}"
  for command_name in "${python_bin}" seqkit xargs; do
    if ! command -v "${command_name}" >/dev/null 2>&1; then
      echo "Missing required command: ${command_name}" >&2
      exit 1
    fi
  done
  if [[ ! -s src/subsample_fastq.py ]]; then
    echo "Missing src/subsample_fastq.py" >&2
    exit 1
  fi

  mkdir -p "${timing_root}"

  run_subsample_task() {
    set -euo pipefail
    local token="$1"
    local accession="${token%%:*}"
    local percent="${token##*:}"
    local subset="sub${percent}"
    local input_dir="data/fastq_data/${accession}/full"
    local output_dir="data/fastq_data/${accession}/${subset}"
    local timings="${timing_root}/${accession}_${subset}.tsv"

    echo "Starting ${accession} ${subset}: $(date -Is)"
    "${python_bin}" src/subsample_fastq.py \
      --input-dir "${input_dir}" \
      --output-dir "${output_dir}" \
      --percent "${percent}" \
      --seed "${seed}" \
      --timings-tsv "${timings}"
    echo "Finished ${accession} ${subset}: $(date -Is)"
  }
  export -f run_subsample_task
  export python_bin seed timing_root

  local accession percent input_dir
  local -a tasks=()
  for percent in 50 25 10; do
    for accession in "${ACCESSIONS[@]}"; do
      input_dir="data/fastq_data/${accession}/full"
      if [[ ! -d "${input_dir}" ]]; then
        echo "Missing full FASTQ directory: ${input_dir}" >&2
        exit 1
      fi
      tasks+=("${accession}:${percent}")
    done
  done

  echo "Started: $(date -Is)"
  echo "Tasks: ${#tasks[@]}"
  echo "Concurrent project/level tasks: ${max_parallel}"
  echo "Conda environment: ${SUBSAMPLE_CONDA_ENV}"
  printf '%s\n' "${tasks[@]}" \
    | xargs -n 1 -P "${max_parallel}" bash -c 'run_subsample_task "$1"' _
  echo "Finished: $(date -Is)"
  echo "Timing files: ${timing_root}"
}

run_pcoa() {
  local level="${LEVEL:-}"
  local threads="${SLURM_CPUS_PER_TASK:-16}"
  local trim_length="${TRIM_LENGTH:-120}"
  local sampling_depth="${SAMPLING_DEPTH:-1000}"
  local metadata="${METADATA:-data/metadata/mmc20_16s_metadata.csv}"
  local gg2_dir="${GG2_DIR:-data/gg2}"
  local run_dir="${RUN_DIR:-runs/mmc20/${level}}"
  local resume="${RESUME:-0}"

  case "${level}" in
    full|sub50|sub25|sub10) ;;
    *)
      echo "LEVEL must be one of: full, sub50, sub25, sub10." >&2
      exit 2
      ;;
  esac
  if [[ "${resume}" != "0" && "${resume}" != "1" ]]; then
    echo "RESUME must be 0 or 1." >&2
    exit 2
  fi

  activate_environment "${PCOA_CONDA_ENV}"
  for command_name in python qiime deblur; do
    if ! command -v "${command_name}" >/dev/null 2>&1; then
      echo "Missing required command: ${command_name}" >&2
      exit 1
    fi
  done
  for required_file in \
    "${metadata}" \
    "${gg2_dir}/2024.09.backbone.full-length.fna.qza" \
    "${gg2_dir}/2024.09.phylogeny.id.nwk.qza" \
    src/run_pipeline.py; do
    if [[ ! -s "${required_file}" ]]; then
      echo "Missing or empty required file: ${required_file}" >&2
      exit 1
    fi
  done

  local accession fastq_dir fastq_count
  local total_fastqs=0
  local -a pipeline_args=()
  for accession in "${ACCESSIONS[@]}"; do
    fastq_dir="data/fastq_data/${accession}/${level}"
    if [[ ! -d "${fastq_dir}" ]]; then
      echo "Missing FASTQ directory: ${fastq_dir}" >&2
      exit 1
    fi
    fastq_count="$(find "${fastq_dir}" -maxdepth 1 -type f -name '*_1.fastq.gz' | wc -l | tr -d ' ')"
    if [[ "${fastq_count}" == "0" ]]; then
      echo "No forward FASTQs found in ${fastq_dir}" >&2
      exit 1
    fi
    echo "${accession}/${level}: ${fastq_count} forward FASTQs"
    total_fastqs="$((total_fastqs + fastq_count))"
    pipeline_args+=(--study "${accession}=${fastq_dir}")
  done

  pipeline_args+=(
    --run-dir "${run_dir}"
    --metadata "${metadata}"
    --color-by environment_harmonized
    --color-by disease_harmonized
    --color-by source_study
    --trim-length "${trim_length}"
    --min-reads 0
    --sampling-depth "${sampling_depth}"
    --gg2-dir "${gg2_dir}"
    --threads "${threads}"
  )
  if [[ "${resume}" == "1" ]]; then
    pipeline_args+=(--resume)
  fi

  echo "Started: $(date -Is)"
  echo "Level: ${level}"
  echo "FASTQs: ${total_fastqs}"
  echo "Run directory: ${run_dir}"
  echo "Conda environment: ${PCOA_CONDA_ENV}"
  echo "Threads: ${threads}"
  echo "Trim length: ${trim_length}"
  echo "Sampling depth: ${sampling_depth}"
  python src/run_pipeline.py "${pipeline_args[@]}"
  echo "Finished: $(date -Is)"
  echo "Results: ${run_dir}/results"
  echo "Timings: ${run_dir}/timings"
}

cd "${REPO_ROOT}"

case "${MODE}" in
  submit) submit_jobs ;;
  subsample) run_subsampling ;;
  pcoa) run_pcoa ;;
  *)
    echo "MODE must be one of: submit, subsample, pcoa." >&2
    exit 2
    ;;
esac
