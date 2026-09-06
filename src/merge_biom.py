#!/usr/bin/env python3
"""Merge two or more Deblur all.biom + all.seqs.fa outputs into one.

Usage:
  python src/merge_biom.py \
    --deblur-dirs work/deblur_study1/workflow work/deblur_study2/workflow \
    --out-dir work/deblur_merged
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import biom
import biom.util

from pipeline_lib import TimingRecorder, add_timing_argument, run_timed_main


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Merge Deblur outputs from multiple studies.")
    parser.add_argument("--deblur-dirs", nargs="+", type=Path, required=True,
                        help="Deblur workflow output dirs (each must have all.biom and all.seqs.fa).")
    parser.add_argument("--out-dir", type=Path, required=True,
                        help="Output directory for merged all.biom and all.seqs.fa.")
    parser.add_argument(
        "--skip-empty",
        action="store_true",
        help="Skip valid Deblur outputs containing zero samples and zero features.",
    )
    parser.add_argument(
        "--run-state",
        type=Path,
        help="Pipeline run_state.json used to exclude non-completed Deblur studies.",
    )
    add_timing_argument(parser)
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


def run(args: argparse.Namespace, timing: TimingRecorder) -> int:
    args.out_dir = args.out_dir.resolve()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    tables = []
    all_seqs: dict[str, str] = {}
    stage_states = None
    if args.run_state is not None:
        with args.run_state.open(encoding="utf-8") as handle:
            stage_states = json.load(handle).get("stages", {})

    for d in args.deblur_dirs:
        d = d.resolve()
        if stage_states is not None:
            stage_status = stage_states.get(f"deblur:{d.parent.name}", {}).get("status")
            if stage_status != "completed":
                print(f"  {d.parent.name}: excluded with stage status {stage_status}")
                timing.skipped(
                    "load_deblur_output",
                    item=str(d),
                    message=f"stage status {stage_status}",
                )
                continue
        biom_fp = d / "all.biom"
        fa_fp = d / "all.seqs.fa"
        if not biom_fp.exists():
            raise FileNotFoundError(f"Missing: {biom_fp}")
        if not fa_fp.exists():
            raise FileNotFoundError(f"Missing: {fa_fp}")
        with timing.step("load_deblur_output", item=str(d)):
            t = biom.load_table(str(biom_fp))
            loaded_seqs = load_seqs(fa_fp)
        if t.shape == (0, 0) and not loaded_seqs:
            if not args.skip_empty:
                raise ValueError(f"Empty Deblur output: {d}")
            print(f"  {d.name}: excluded empty Deblur output")
            timing.skipped(
                "load_deblur_output",
                item=str(d),
                message="zero samples and zero features",
            )
            continue
        print(f"  {d.name}: {t.shape[1]} samples, {t.shape[0]} features")
        tables.append(t)
        all_seqs.update(loaded_seqs)

    if not tables:
        raise ValueError("No non-empty Deblur outputs are available to merge")

    with timing.step("merge_feature_tables", item=f"{len(tables)} tables"):
        merged = tables[0]
        for t in tables[1:]:
            merged = merged.merge(t)

    print(f"Merged: {merged.shape[1]} samples, {merged.shape[0]} features")

    out_biom = args.out_dir / "all.biom"
    with timing.step("write_merged_biom", item=str(out_biom)):
        with biom.util.biom_open(str(out_biom), "w") as f:
            merged.to_hdf5(f, "merge_biom.py")
    print(f"Written: {out_biom}")

    out_fa = args.out_dir / "all.seqs.fa"
    with timing.step("write_representative_sequences", item=str(out_fa)):
        feature_ids = set(merged.ids(axis="observation"))
        with out_fa.open("w") as f:
            for seq_id, seq in all_seqs.items():
                if seq_id in feature_ids:
                    f.write(f">{seq_id}\n{seq}\n")
    print(f"Written: {out_fa}")

    return 0


def main() -> int:
    args = parse_args()
    timing = TimingRecorder(args.timings_tsv, component="merge_biom")
    return run_timed_main(timing, lambda: run(args, timing))


if __name__ == "__main__":
    raise SystemExit(main())
