#!/usr/bin/env bash
#SBATCH --job-name=perm-sub25
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --output=/home/yal090/PCOA-prototype/runs/mmc2_combined_amplicon/dart-weighted-20260925/logs/%j-perm-sub25.out
#SBATCH --error=/home/yal090/PCOA-prototype/runs/mmc2_combined_amplicon/dart-weighted-20260925/logs/%j-perm-sub25.err

set -euo pipefail
source /opt/software/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10
cd /home/yal090/PCOA-prototype

python -u src/run_permanova.py \
  --distance-matrix runs/mmc2_combined_amplicon/dart-20260921/cached/sub25/work/dart/unweighted_unifrac_distance_matrix.tsv.zst \
  --metadata runs/mmc2_combined_amplicon/dart-20260921/full/metadata/illumina_pipeline_metadata.tsv \
  --grouping-column environment_harmonized \
  --permutations 999 \
  --max-samples 100000 \
  --out-dir runs/mmc2_combined_amplicon/dart-weighted-20260925/permanova/unweighted/sub25/
