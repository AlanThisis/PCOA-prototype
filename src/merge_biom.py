#!/usr/bin/env python3
"""Merge two or more Deblur all.biom + all.seqs.fa outputs into one.

Usage:
  python src/merge_biom.py \
    --deblur-dirs work/deblur_study1/workflow work/deblur_study2/workflow \
    --out-dir work/deblur_merged
"""
from __future__ import annotations

import argparse
from pathlib import Path

import biom
import biom.util
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Merge Deblur outputs from multiple studies.")
    parser.add_argument("--deblur-dirs", nargs="+", type=Path, required=True,
                        help="Deblur workflow output dirs (each must have all.biom and all.seqs.fa).")
    parser.add_argument("--out-dir", type=Path, required=True,
                        help="Output directory for merged all.biom and all.seqs.fa.")
    return parser.parse_args()


def load_seqs(fa_fp: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    current_id: str | None = None
    with fa_fp.open() as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                current_id = line[1:].split()[0]
            elif current_id is not None:
                seqs[current_id] = seqs.get(current_id, "") + line
    return seqs


def main() -> int:
    args = parse_args()
    args.out_dir = args.out_dir.resolve()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    tables = []
    all_seqs: dict[str, str] = {}

    for d in args.deblur_dirs:
        d = d.resolve()
        biom_fp = d / "all.biom"
        fa_fp = d / "all.seqs.fa"
        if not biom_fp.exists():
            raise FileNotFoundError(f"Missing: {biom_fp}")
        if not fa_fp.exists():
            raise FileNotFoundError(f"Missing: {fa_fp}")
        t = biom.load_table(str(biom_fp))
        print(f"  {d.name}: {t.shape[1]} samples, {t.shape[0]} features")
        tables.append(t)
        all_seqs.update(load_seqs(fa_fp))

    merged = tables[0]
    for t in tables[1:]:
        merged = merged.merge(t)

    print(f"Merged: {merged.shape[1]} samples, {merged.shape[0]} features")

    out_biom = args.out_dir / "all.biom"
    with biom.util.biom_open(str(out_biom), "w") as f:
        merged.to_hdf5(f, "merge_biom.py")
    print(f"Written: {out_biom}")

    out_fa = args.out_dir / "all.seqs.fa"
    feature_ids = set(merged.ids(axis="observation"))
    with out_fa.open("w") as f:
        for seq_id, seq in all_seqs.items():
            if seq_id in feature_ids:
                f.write(f">{seq_id}\n{seq}\n")
    print(f"Written: {out_fa}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
