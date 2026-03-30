#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

from pipeline_lib import parse_deblur_clean_fasta


CLEAN_PATTERN = re.compile(r"^(?P<sample>.+)\.trim\d+\.fasta\.derep\.clean$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a feature table from Deblur .clean files."
    )
    parser.add_argument("--work-dir", type=Path, default=Path("work/deblur"))
    parser.add_argument("--results-dir", type=Path, default=Path("results/forward_only"))
    return parser.parse_args()


def sample_id_from_clean(clean_fp: Path) -> str:
    name = clean_fp.name
    match = CLEAN_PATTERN.match(name)
    if match:
        return match.group("sample")
    return name.split(".", 1)[0]


def discover_clean_files(work_dir: Path) -> list[Path]:
    clean_files = sorted(work_dir.rglob("*.clean"))
    if not clean_files:
        raise FileNotFoundError(f"No Deblur .clean files found under {work_dir}")
    return clean_files


def build_feature_table(clean_files: list[Path]) -> pd.DataFrame:
    by_sample = {}
    for clean_fp in clean_files:
        by_sample[sample_id_from_clean(clean_fp)] = parse_deblur_clean_fasta(clean_fp)

    all_features = sorted({seq for counts in by_sample.values() for seq in counts})
    table = pd.DataFrame(0, index=sorted(by_sample.keys()), columns=all_features, dtype=int)
    for sample_id, counts in by_sample.items():
        for seq, n in counts.items():
            table.at[sample_id, seq] = n

    table.index.name = "sample_id"
    table.columns.name = "feature_sequence"
    return table


def main() -> int:
    args = parse_args()
    args.work_dir = args.work_dir.resolve()
    args.results_dir = args.results_dir.resolve()
    args.results_dir.mkdir(parents=True, exist_ok=True)

    clean_files = discover_clean_files(args.work_dir)
    feature_table = build_feature_table(clean_files)
    feature_fp = args.results_dir / "feature_table.tsv"
    feature_table.to_csv(feature_fp, sep="\t")

    print(f"Feature table shape: {feature_table.shape}")
    print("Samples in feature table:", ", ".join(feature_table.index.tolist()))
    print(f"Finished. Feature table written to: {feature_fp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
