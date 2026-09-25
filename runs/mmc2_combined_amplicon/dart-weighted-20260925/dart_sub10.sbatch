#!/usr/bin/env bash
#SBATCH --job-name=wuf-sub10
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=/home/yal090/PCOA-prototype/runs/mmc2_combined_amplicon/dart-weighted-20260925/logs/%j-wuf-sub10.out
#SBATCH --error=/home/yal090/PCOA-prototype/runs/mmc2_combined_amplicon/dart-weighted-20260925/logs/%j-wuf-sub10.err

set -euo pipefail
source /opt/software/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10
cd /home/yal090/PCOA-prototype

dart=/home/yal090/.conda/envs/dartunifrac-0.3.0/bin/dartunifrac
test -x "$dart"

python -u runs/mmc2_combined_amplicon/dart-weighted-20260925/run_weighted_dart.py sub10 \
  --dart "$dart" \
  --threads "$SLURM_CPUS_PER_TASK"
