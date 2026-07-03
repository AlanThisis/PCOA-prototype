#!/usr/bin/env python3
"""Plot a PCoA ordination colored by a metadata column.

Examples
--------
python src/plot_pcoa_crossstudy.py \\
  --pcoa results/prjeb44533_sub10/pcoa_coordinates_unweighted_unifrac.txt \\
  --metadata data/PRJEB44533/metadata.csv \\
  --color-by severity_groups \\
  --out results/prjeb44533_sub10/pcoa_severity.png

python src/plot_pcoa_crossstudy.py \\
  --pcoa results/crc_cross_study/pcoa_coordinates_unweighted_unifrac.txt \\
  --metadata data/fastq_data/CRC_cross_study/crc_cross_study_metadata.tsv \\
  --color-by disease_status \\
  --out results/crc_cross_study/pcoa_disease.png \\
  --title "CRC cross-study — disease status"

The metadata file can be:
  - CSV from get_ENA_metadata.py  (join key: run_accessions column)
  - TSV with a sample-id column   (join key: sample-id column)
PCoA sample IDs are matched after stripping trailing _1 / _2 suffixes.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import skbio


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot PCoA coordinates colored by a metadata column."
    )
    parser.add_argument("--pcoa", type=Path, required=True,
                        help="OrdinationResults file (from unifrac.py).")
    parser.add_argument("--metadata", type=Path, required=True,
                        help="Metadata CSV/TSV. Must have either a run_accessions "
                             "column (get_ENA_metadata.py output) or a sample-id column.")
    parser.add_argument("--color-by", required=True,
                        help="Column name to color points by.")
    parser.add_argument("--out", type=Path, required=True,
                        help="Output PNG path.")
    parser.add_argument("--pc", type=int, nargs=2, default=[1, 2], metavar=("X", "Y"),
                        help="Which PCs to plot (1-based, default: 1 2).")
    parser.add_argument("--title", default=None,
                        help="Plot title. Defaults to the output filename stem.")
    return parser.parse_args()


def strip_read_suffix(sample_id: str) -> str:
    return re.sub(r"(_R?[12](_001)?)$", "", sample_id)


def load_id_to_label(metadata_path: Path, color_by: str) -> dict[str, str]:
    """Return {sample_or_run_id: label}.

    Supports two metadata formats:
    - run_accessions column (get_ENA_metadata.py): maps each ERR to the label
    - sample-id column (legacy TSV): maps sample-id directly
    """
    sep = "\t" if metadata_path.suffix == ".tsv" else ","
    id_to_label: dict[str, str] = {}

    with metadata_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=sep)
        fields = reader.fieldnames or []
        if color_by not in fields:
            raise ValueError(
                f"Column {color_by!r} not in metadata. Available: {', '.join(fields)}"
            )

        use_run_accessions = "run_accessions" in fields
        use_sample_id = "sample-id" in fields

        for row in reader:
            label = row[color_by].strip()
            if use_run_accessions:
                for run in row["run_accessions"].split(";"):
                    run = run.strip()
                    if run:
                        id_to_label[run] = label
            elif use_sample_id:
                sid = row["sample-id"].strip()
                if sid:
                    id_to_label[sid] = label

    return id_to_label


def main() -> int:
    args = parse_args()

    pcoa = skbio.OrdinationResults.read(str(args.pcoa))
    px, py = args.pc[0] - 1, args.pc[1] - 1
    coords = pcoa.samples.iloc[:, [px, py]].copy()
    coords.columns = ["x", "y"]

    eigvals = pcoa.eigvals
    pct = eigvals / eigvals.sum() * 100

    try:
        id_to_label = load_id_to_label(args.metadata, args.color_by)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    labels = [id_to_label.get(strip_read_suffix(sid), "Unknown") for sid in coords.index]
    unmapped = labels.count("Unknown")
    if unmapped:
        print(f"Warning: {unmapped}/{len(labels)} samples not found in metadata (labeled 'Unknown').",
              file=sys.stderr)

    unique_labels = sorted(set(labels))
    cmap = plt.get_cmap("tab10" if len(unique_labels) <= 10 else "tab20")
    color_map = {label: cmap(i) for i, label in enumerate(unique_labels)}

    fig, ax = plt.subplots(figsize=(8, 6))
    for label in unique_labels:
        idx = [i for i, l in enumerate(labels) if l == label]
        ax.scatter(
            coords["x"].iloc[idx],
            coords["y"].iloc[idx],
            label=f"{label} (n={len(idx)})",
            color=color_map[label],
            alpha=0.7,
            s=40,
            edgecolors="none",
        )

    ax.set_xlabel(f"PC{args.pc[0]} ({pct.iloc[px]:.1f}%)")
    ax.set_ylabel(f"PC{args.pc[1]} ({pct.iloc[py]:.1f}%)")
    ax.set_title(args.title or args.out.stem)
    ax.legend(title=args.color_by, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
    plt.tight_layout()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.out, dpi=150)
    print(f"Saved to {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
