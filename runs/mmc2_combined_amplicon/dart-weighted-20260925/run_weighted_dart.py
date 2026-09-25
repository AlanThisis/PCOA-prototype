#!/usr/bin/env python3
"""Run weighted DART/fPCoA on existing rarefied GG2-mapped QIIME tables."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "src"))

from dart_unifrac import dart_version, prepare_dart_inputs, run_dart_unifrac
from pipeline_lib import TimingRecorder


CACHED_RAREFIED_TABLES = {
    "full": REPO / "runs/mmc2_combined_amplicon/dart-20260921/full/work/qiime2/rarefied-backbone-mapped-table.qza",
    "sub50": REPO / "runs/mmc2_combined_amplicon/repair-20260917/sub50/work/qiime2/rarefied-backbone-mapped-table.qza",
    "sub25": REPO / "runs/mmc2_combined_amplicon/repair-20260917/sub25/work/qiime2/rarefied-backbone-mapped-table.qza",
    "sub10": REPO / "runs/mmc2_combined_amplicon/sub10/work/qiime2/rarefied-backbone-mapped-table.qza",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("level", choices=("full", "sub10", "sub25", "sub50"))
    parser.add_argument("--dart", type=Path, required=True)
    parser.add_argument("--threads", type=int, required=True)
    args = parser.parse_args()

    cached_table = CACHED_RAREFIED_TABLES[args.level]
    if not cached_table.is_file() or cached_table.stat().st_size == 0:
        raise FileNotFoundError(f"Missing cached rarefied table: {cached_table}")

    output = REPO / "runs/mmc2_combined_amplicon/dart-weighted-20260925" / args.level
    work = output / "work/dart"
    results = output / "results"
    results.mkdir(parents=True, exist_ok=True)
    timing = TimingRecorder(output / "timings/dart.tsv", component="dart-weighted")

    biom_fp, tree_fp, stats = prepare_dart_inputs(
        cached_table,
        REPO / "data/gg2/2024.09.phylogeny.id.nwk.qza",
        work,
    )
    result = run_dart_unifrac(
        executable=str(args.dart),
        biom_fp=biom_fp,
        tree_fp=tree_fp,
        work_dir=work,
        threads=args.threads,
        sketch_size=2048,
        seed=1337,
        bbits=16,
        compress=True,
        timing=timing,
        weighted=True,
    )
    if result["sample_count"] != stats["sample_count"]:
        raise RuntimeError("DART sample count differs from the cached rarefied table")

    coordinates = results / "pcoa_coordinates_weighted_unifrac.txt"
    shutil.copy2(result["ordination_path"], coordinates)
    summary = {
        "unifrac_engine": "dart",
        "metric": "weighted_unifrac",
        "sampling_depth": 1000,
        "rarefied_samples": stats["sample_count"],
        "rarefied_features": stats["feature_count"],
        "cached_rarefied_table": str(cached_table),
        "cached_rarefied_table_sha256": hashlib.sha256(cached_table.read_bytes()).hexdigest(),
        "dart": {
            "version": dart_version(str(args.dart)),
            "method": "dmh",
            "weighted": True,
            "sketch_size": 2048,
            "seed": 1337,
            "bbits": 16,
            "threads": args.threads,
        },
        "pcoa": {"selected_method": "dart-fpcoa", "dimensions": result["dimensions"]},
        "distance_tsv": {
            "exported": True,
            "backend_path": str(result["distance_path"]),
        },
    }
    (results / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Weighted DART/fPCoA ready: {args.level}, {stats['sample_count']} samples", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
