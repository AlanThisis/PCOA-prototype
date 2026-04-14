#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import tempfile
from pathlib import Path

import pandas as pd

from pipeline_lib import resolve_executable, run_command


SAMPLE_SUFFIX_PATTERN = re.compile(r"_[0-9]+_L[0-9]{3}_R[12]_[0-9]{3}$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a feature table from Deblur BIOM output."
    )
    parser.add_argument("--work-dir", type=Path, default=Path("work/deblur"))
    parser.add_argument("--results-dir", type=Path, default=Path("results/forward_only"))
    parser.add_argument("--biom-fp", type=Path, default=None)
    return parser.parse_args()


def normalize_sample_id(sample_id: str) -> str:
    return SAMPLE_SUFFIX_PATTERN.sub("", sample_id)


def convert_biom_to_tsv(biom_fp: Path, tsv_fp: Path, biom_executable_path: str) -> None:
    run_command(
        [
            biom_executable_path,
            "convert",
            "-i",
            str(biom_fp),
            "-o",
            str(tsv_fp),
            "--to-tsv",
        ]
    )


def read_biom_tsv(tsv_fp: Path) -> pd.DataFrame:
    table = pd.read_csv(tsv_fp, sep="\t", skiprows=1)
    first_col = table.columns[0]
    table = table.rename(columns={first_col: "feature_id"}).set_index("feature_id")
    return table.T


def build_feature_table_from_biom(
    biom_fp: Path, biom_executable_path: str
) -> pd.DataFrame:
    with tempfile.TemporaryDirectory(prefix="biom-convert-") as tmp_dir:
        tsv_fp = Path(tmp_dir) / "feature_table.tsv"
        convert_biom_to_tsv(biom_fp, tsv_fp, biom_executable_path)
        feature_table = read_biom_tsv(tsv_fp)

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

    biom_executable_path = resolve_executable("biom")
    feature_table = build_feature_table_from_biom(biom_fp, biom_executable_path)
    feature_fp = args.results_dir / "feature_table.tsv"
    feature_table.to_csv(feature_fp, sep="\t")

    print(f"Feature table shape: {feature_table.shape}")
    print("Samples in feature table:", ", ".join(feature_table.index.tolist()))
    print(f"Finished. Feature table written to: {feature_fp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
