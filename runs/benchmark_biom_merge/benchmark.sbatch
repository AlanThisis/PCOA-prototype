#!/usr/bin/env bash
#SBATCH --job-name=biom-merge-bench
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/home/yal090/PCOA-prototype/runs/benchmark_biom_merge/logs/%j.out
#SBATCH --error=/home/yal090/PCOA-prototype/runs/benchmark_biom_merge/logs/%j.err

set -euo pipefail
source /opt/software/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10
cd /home/yal090/PCOA-prototype

mkdir -p runs/benchmark_biom_merge/logs

# MMC20: 20 per-bioproject deblur outputs
deblur_dirs=(runs/mmc20/full/work/deblur/*/workflow)

echo "=== biom merge() vs concat() benchmark ==="
echo "Found ${#deblur_dirs[@]} deblur directories"
echo ""

python runs/benchmark_biom_merge/benchmark.py \
  --deblur-dirs "${deblur_dirs[@]}" \
  --repeats 3
