#!/usr/bin/env bash
set -euo pipefail
cd /home/yal090/PCOA-prototype
dir="runs/mmc2_combined_amplicon/dart-weighted-20260925"
mkdir -p "$dir/logs"

# ── 4 weighted DART jobs (independent) ──
dart_full=$(sbatch --parsable "$dir/dart_full.sbatch")
dart_sub50=$(sbatch --parsable "$dir/dart_sub50.sbatch")
dart_sub25=$(sbatch --parsable "$dir/dart_sub25.sbatch")
dart_sub10=$(sbatch --parsable "$dir/dart_sub10.sbatch")
echo "Weighted DART:  full=$dart_full  sub50=$dart_sub50  sub25=$dart_sub25  sub10=$dart_sub10"

# ── 4 PERMANOVA jobs (independent, uses existing unweighted DMs) ──
perm_sub10=$(sbatch --parsable "$dir/permanova_sub10.sbatch")
perm_sub25=$(sbatch --parsable "$dir/permanova_sub25.sbatch")
perm_sub50=$(sbatch --parsable "$dir/permanova_sub50.sbatch")
perm_full=$(sbatch --parsable "$dir/permanova_full.sbatch")
echo "PERMANOVA:      sub10=$perm_sub10  sub25=$perm_sub25  sub50=$perm_sub50  full=$perm_full"

# ── Procrustes (depends on all 4 weighted DART jobs) ──
procrustes=$(sbatch --parsable --dependency=afterok:${dart_full}:${dart_sub50}:${dart_sub25}:${dart_sub10} "$dir/procrustes.sbatch")
echo "Procrustes:     $procrustes (after DART: $dart_full,$dart_sub50,$dart_sub25,$dart_sub10)"

echo ""
echo "Total: 9 jobs submitted"
