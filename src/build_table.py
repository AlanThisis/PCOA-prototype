#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


SAMPLE_SUFFIX_PATTERN = re.compile(r"_[0-9]+_L[0-9]{3}_R[12]_[0-9]{3}$")
ENA_READ_SUFFIX_PATTERN = re.compile(r"^(?P<accession>[EDS]RR[0-9]+)_[12]$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a feature table from Deblur BIOM output."
    )
    parser.add_argument("--work-dir", type=Path, default=Path("work/deblur"))
    parser.add_argument("--results-dir", type=Path, default=Path("results/forward_only"))
    parser.add_argument("--biom-fp", type=Path, default=None)
    return parser.parse_args()


def normalize_sample_id(sample_id: str) -> str:
    demux_normalized = SAMPLE_SUFFIX_PATTERN.sub("", sample_id)
    ena_match = ENA_READ_SUFFIX_PATTERN.match(demux_normalized)
    if ena_match:
        return ena_match.group("accession")
    return demux_normalized


def read_feature_table_from_biom(biom_fp: Path) -> pd.DataFrame:
    from biom import load_table

    biom_table = load_table(str(biom_fp))
    return biom_table.to_dataframe(dense=True).T


def build_feature_table_from_biom(biom_fp: Path) -> pd.DataFrame:
    feature_table = read_feature_table_from_biom(biom_fp)

    feature_table.index = feature_table.index.to_series().map(normalize_sample_id)
    feature_table = feature_table.groupby(level=0, sort=True).sum()
    feature_table = feature_table.reindex(sorted(feature_table.columns), axis=1)
    feature_table.index.name = "sample_id"
    feature_table.columns.name = "feature_id"
    return feature_table


def main() -> int:
    args = parse_args()
    args.work_dir = args.work_dir.resolve()
    args.results_dir = args.results_dir.resolve()
    args.results_dir.mkdir(parents=True, exist_ok=True)
    biom_fp = (
        args.biom_fp.resolve()
        if args.biom_fp is not None
        else args.work_dir / "workflow" / "all.biom"
    )
    if not biom_fp.exists():
        raise FileNotFoundError(f"Deblur BIOM file not found: {biom_fp}")

    feature_table = build_feature_table_from_biom(biom_fp)
    feature_fp = args.results_dir / "feature_table.tsv"
    feature_table.to_csv(feature_fp, sep="\t")

    print(f"Feature table shape: {feature_table.shape}")
    print("Samples in feature table:", ", ".join(feature_table.index.tolist()))
    print(f"Finished. Feature table written to: {feature_fp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
