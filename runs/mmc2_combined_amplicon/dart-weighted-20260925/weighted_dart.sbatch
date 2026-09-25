#!/usr/bin/env bash
#SBATCH --job-name=weighted-unifrac
#SBATCH --partition=long
#SBATCH --time=36-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --output=/home/yal090/PCOA-prototype/runs/mmc2_combined_amplicon/dart-weighted-20260925/logs/%j-%x.out
#SBATCH --error=/home/yal090/PCOA-prototype/runs/mmc2_combined_amplicon/dart-weighted-20260925/logs/%j-%x.err

set -euo pipefail
source /opt/software/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10
cd /home/yal090/PCOA-prototype

dart=/home/yal090/.conda/envs/dartunifrac-0.3.0/bin/dartunifrac
test -x "$dart"

run_root="runs/mmc2_combined_amplicon/dart-weighted-20260925"
mkdir -p "$run_root/logs"

# ── Weighted DART on all 4 levels ──
for level in full sub50 sub25 sub10; do
  echo "=== Weighted UniFrac: $level ==="
  python -u "$run_root/run_weighted_dart.py" "$level" \
    --dart "$dart" \
    --threads "$SLURM_CPUS_PER_TASK"
  echo ""
done

# ── Procrustes across subsampling levels ──
echo "=== Procrustes: weighted UniFrac across subsampling levels ==="
mkdir -p "$run_root/procrustes"
python -u src/compare_subsampling_procrustes.py \
  --full   "$run_root/full/results/pcoa_coordinates_weighted_unifrac.txt" \
  --sub50  "$run_root/sub50/results/pcoa_coordinates_weighted_unifrac.txt" \
  --sub25  "$run_root/sub25/results/pcoa_coordinates_weighted_unifrac.txt" \
  --sub10  "$run_root/sub10/results/pcoa_coordinates_weighted_unifrac.txt" \
  --out-dir "$run_root/procrustes/"
echo ""

# ── PERMANOVA on sub10 (environment_harmonized) ──
echo "=== PERMANOVA: sub10, environment_harmonized ==="
python -u src/run_permanova.py \
  --distance-matrix "$run_root/sub10/work/dart/weighted_unifrac_distance_matrix.tsv.zst" \
  --metadata runs/mmc2_combined_amplicon/dart-20260921/full/metadata/illumina_pipeline_metadata.tsv \
  --grouping-column environment_harmonized \
  --permutations 999 \
  --max-samples 20000 \
  --out-dir "$run_root/permanova/sub10/"
echo ""

# ── PERMANOVA on sub25 (environment_harmonized) ──
echo "=== PERMANOVA: sub25, environment_harmonized ==="
python -u src/run_permanova.py \
  --distance-matrix "$run_root/sub25/work/dart/weighted_unifrac_distance_matrix.tsv.zst" \
  --metadata runs/mmc2_combined_amplicon/dart-20260921/full/metadata/illumina_pipeline_metadata.tsv \
  --grouping-column environment_harmonized \
  --permutations 999 \
  --max-samples 50000 \
  --out-dir "$run_root/permanova/sub25/"

echo ""
echo "=== All done ==="
