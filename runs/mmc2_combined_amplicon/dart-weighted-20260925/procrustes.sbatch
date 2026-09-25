#!/usr/bin/env bash
#SBATCH --job-name=wuf-procrustes
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/home/yal090/PCOA-prototype/runs/mmc2_combined_amplicon/dart-weighted-20260925/logs/%j-procrustes.out
#SBATCH --error=/home/yal090/PCOA-prototype/runs/mmc2_combined_amplicon/dart-weighted-20260925/logs/%j-procrustes.err

set -euo pipefail
source /opt/software/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10
cd /home/yal090/PCOA-prototype

run_root="runs/mmc2_combined_amplicon/dart-weighted-20260925"

python -u src/compare_subsampling_procrustes.py \
  --full   "$run_root/full/results/pcoa_coordinates_weighted_unifrac.txt" \
  --sub50  "$run_root/sub50/results/pcoa_coordinates_weighted_unifrac.txt" \
  --sub25  "$run_root/sub25/results/pcoa_coordinates_weighted_unifrac.txt" \
  --sub10  "$run_root/sub10/results/pcoa_coordinates_weighted_unifrac.txt" \
  --out-dir "$run_root/procrustes/"
