#!/usr/bin/env python3
"""Run PERMANOVA on a UniFrac distance matrix grouped by a metadata column.

Usage:
  python src/run_permanova.py \
    --distance-matrix dm.tsv.zst \
    --metadata metadata.tsv \
    --grouping-column environment_harmonized \
    --out-dir results/permanova/
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
import skbio


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--distance-matrix", type=Path, required=True,
                        help="Distance matrix TSV (or .tsv.zst for zstd-compressed).")
    parser.add_argument("--metadata", type=Path, required=True,
                        help="Sample metadata TSV with a sample-id column.")
    parser.add_argument("--grouping-column", type=str, default="environment_harmonized",
                        help="Metadata column to group by (default: environment_harmonized).")
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--max-samples", type=int, default=50000,
                        help="Safety limit: error if distance matrix exceeds this many samples.")
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args()


def load_distance_matrix(path: Path, max_samples: int) -> skbio.DistanceMatrix:
    if path.suffix == ".zst" or str(path).endswith(".tsv.zst"):
        print(f"Decompressing {path} via zstd...", flush=True)
        with tempfile.NamedTemporaryFile(suffix=".tsv", delete=True) as tmp:
            subprocess.run(
                ["zstd", "-dc", "-o", tmp.name, str(path)],
                check=True,
            )
            print(f"Reading decompressed distance matrix from {tmp.name}...", flush=True)
            dm = skbio.DistanceMatrix.read(tmp.name)
    else:
        dm = skbio.DistanceMatrix.read(str(path))
    n = dm.shape[0]
    print(f"Distance matrix: {n} samples", flush=True)
    if n > max_samples:
        raise ValueError(
            f"Distance matrix has {n} samples, exceeding --max-samples {max_samples}. "
            "Use a subsampled level or increase the limit."
        )
    return dm


def load_metadata(path: Path, grouping_column: str) -> pd.Series:
    md = pd.read_csv(path, sep="\t", dtype=str)
    id_col = "sample-id" if "sample-id" in md.columns else md.columns[0]
    if grouping_column not in md.columns:
        raise ValueError(
            f"Column {grouping_column!r} not found in metadata. "
            f"Available: {list(md.columns)}"
        )
    md = md[[id_col, grouping_column]].dropna()
    md = md.set_index(id_col)[grouping_column]
    return md


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    dm = load_distance_matrix(args.distance_matrix, args.max_samples)
    grouping = load_metadata(args.metadata, args.grouping_column)

    shared = sorted(set(dm.ids) & set(grouping.index))
    print(f"Samples in distance matrix: {len(dm.ids)}", flush=True)
    print(f"Samples in metadata with {args.grouping_column}: {len(grouping)}", flush=True)
    print(f"Shared samples: {len(shared)}", flush=True)

    if len(shared) < 3:
        raise ValueError(f"Need at least 3 shared samples, found {len(shared)}")

    dm_filtered = dm.filter(shared)
    grouping_filtered = grouping.loc[shared]

    group_counts = grouping_filtered.value_counts()
    print(f"\nGroup counts ({args.grouping_column}):", flush=True)
    for group, count in group_counts.items():
        print(f"  {group}: {count}", flush=True)

    if len(group_counts) < 2:
        raise ValueError("PERMANOVA requires at least 2 groups")

    print(f"\nRunning PERMANOVA with {args.permutations} permutations...", flush=True)
    result = skbio.stats.distance.permanova(
        dm_filtered, grouping_filtered, permutations=args.permutations
    )
    print(f"\n{result}", flush=True)

    summary = {
        "test": "PERMANOVA",
        "grouping_column": args.grouping_column,
        "permutations": args.permutations,
        "sample_count": len(shared),
        "group_count": len(group_counts),
        "groups": {str(k): int(v) for k, v in group_counts.items()},
        "test_statistic": float(result["test statistic"]),
        "p_value": float(result["p-value"]),
        "sample_size": int(result["sample size"]),
        "number_of_groups": int(result["number of groups"]),
        "distance_matrix": str(args.distance_matrix),
        "metadata": str(args.metadata),
    }

    out_json = args.out_dir / "permanova_results.json"
    out_json.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\nWritten: {out_json}", flush=True)

    out_tsv = args.out_dir / "permanova_results.tsv"
    pd.DataFrame([summary]).to_csv(out_tsv, sep="\t", index=False)
    print(f"Written: {out_tsv}", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
