#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

DEFAULT_BIOM_FP = Path("work/deblur/workflow/all.biom")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a feature table from Deblur workflow BIOM output."
    )
    parser.add_argument("--biom-fp", type=Path, default=DEFAULT_BIOM_FP)
    parser.add_argument("--results-dir", type=Path, default=Path("results/forward_only"))
    return parser.parse_args()


def load_biom_table(biom_fp: Path):
    try:
        import biom
    except ImportError as exc:
        raise RuntimeError(
            "BIOM input requires the 'biom' package. Install biom-format in this environment."
        ) from exc

    return biom.load_table(str(biom_fp))


def build_feature_table(biom_table) -> pd.DataFrame:
    table = biom_table.to_dataframe(dense=True).transpose()
    table.index.name = "sample_id"
    table.columns.name = "feature_sequence"
    return table.astype(int)


def main() -> int:
    args = parse_args()
    args.biom_fp = args.biom_fp.resolve()
    args.results_dir = args.results_dir.resolve()
    args.results_dir.mkdir(parents=True, exist_ok=True)

    if not args.biom_fp.exists():
        raise FileNotFoundError(f"Deblur BIOM output not found: {args.biom_fp}")

    biom_table = load_biom_table(args.biom_fp)
    feature_table = build_feature_table(biom_table)
    feature_fp = args.results_dir / "feature_table.tsv"
    feature_table.to_csv(feature_fp, sep="\t")

    print(f"Feature table shape: {feature_table.shape}")
    print("Samples in feature table:", ", ".join(feature_table.index.tolist()))
    print(f"Finished. Feature table written to: {feature_fp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
