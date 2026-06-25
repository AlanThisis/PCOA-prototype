#!/usr/bin/env python3
"""Plot cross-study PCoA colored by disease status and study.

Usage:
  uv run --with scikit-bio --with matplotlib --with pandas \
    src/plot_pcoa_crossstudy.py \
    --pcoa results/crc_cross_study/pcoa_coordinates_unweighted_unifrac.txt \
    --metadata data/fastq_data/CRC_cross_study/crc_cross_study_metadata.tsv \
    --out-dir results/crc_cross_study
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import skbio


DISEASE_COLORS = {"cancer": "#d62728", "control": "#2ca02c"}
STUDY_MARKERS = {"PRJEB46665": "o", "ERP005534": "s"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pcoa", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    ord_res = skbio.io.read(
        str(args.pcoa),
        format="ordination",
        into=skbio.stats.ordination.OrdinationResults,
    )
    coords = ord_res.samples
    prop = ord_res.proportion_explained

    meta = pd.read_csv(args.metadata, sep="\t", index_col="sample-id")

    # Align — sample IDs in PCoA may have _1 suffix from ENA
    def strip_suffix(sid: str) -> str:
        if sid.endswith("_1") and sid[:-2] in meta.index:
            return sid[:-2]
        return sid

    coords.index = [strip_suffix(s) for s in coords.index]
    coords = coords[coords.index.isin(meta.index)]

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # Plot 1: colored by disease status
    ax = axes[0]
    for disease, group in meta.groupby("disease_status"):
        idx = [s for s in group.index if s in coords.index]
        ax.scatter(
            coords.loc[idx, coords.columns[0]],
            coords.loc[idx, coords.columns[1]],
            c=DISEASE_COLORS.get(disease, "gray"),
            label=disease,
            s=60,
            alpha=0.85,
        )
    ax.set_xlabel(f"PC1 ({prop.iloc[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({prop.iloc[1]*100:.1f}%)")
    ax.set_title("Unweighted UniFrac — colored by disease")
    ax.legend(title="disease_status")

    # Plot 2: colored by study
    ax = axes[1]
    study_colors = {"PRJEB46665": "#1f77b4", "ERP005534": "#ff7f0e"}
    for study, group in meta.groupby("study"):
        idx = [s for s in group.index if s in coords.index]
        ax.scatter(
            coords.loc[idx, coords.columns[0]],
            coords.loc[idx, coords.columns[1]],
            c=study_colors.get(study, "gray"),
            marker=STUDY_MARKERS.get(study, "o"),
            label=study,
            s=60,
            alpha=0.85,
        )
    ax.set_xlabel(f"PC1 ({prop.iloc[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({prop.iloc[1]*100:.1f}%)")
    ax.set_title("Unweighted UniFrac — colored by study")
    ax.legend(title="study")

    fig.suptitle(
        "Cross-study CRC PCoA (PRJEB46665 + ERP005534, GG2 full-length backbone)",
        fontsize=13,
    )
    fig.tight_layout()
    out_fp = args.out_dir / "pcoa_cross_study.png"
    fig.savefig(out_fp, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_fp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
