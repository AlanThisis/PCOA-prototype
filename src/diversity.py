#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute beta diversity and PCoA from a feature table."
    )
    parser.add_argument("--results-dir", type=Path, default=Path("results/forward_only"))
    parser.add_argument("--metric", default="braycurtis")
    return parser.parse_args()


def compute_pcoa(feature_table: pd.DataFrame, metric: str):
    ids = feature_table.index.tolist()
    matrix = feature_table.to_numpy(dtype=float)
    dm = beta_diversity(metric=metric, counts=matrix, ids=ids)
    ord_res = pcoa(dm)
    return dm, ord_res


def main() -> int:
    args = parse_args()
    args.results_dir = args.results_dir.resolve()
    feature_fp = args.results_dir / "feature_table.tsv"
    if not feature_fp.exists():
        raise FileNotFoundError(f"Feature table not found: {feature_fp}")

    feature_table = pd.read_csv(feature_fp, sep="\t", index_col=0)
    dm, ord_res = compute_pcoa(feature_table, metric=args.metric)

    dm_fp = args.results_dir / f"distance_matrix_{args.metric}.tsv"
    coords_fp = args.results_dir / "pcoa_coordinates.tsv"
    plot_fp = args.results_dir / "pcoa_plot.png"

    dm_df = pd.DataFrame(dm.data, index=dm.ids, columns=dm.ids)
    dm_df.index.name = "sample_id"
    dm_df.to_csv(dm_fp, sep="\t")

    coords = ord_res.samples.copy()
    coords.index.name = "sample_id"
    coords.to_csv(coords_fp, sep="\t")

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(coords["PC1"], coords["PC2"])
    for sample_id, xv, yv in zip(coords.index, coords["PC1"], coords["PC2"]):
        ax.annotate(sample_id, (xv, yv), fontsize=8)
    ax.set_xlabel(f"PC1 ({ord_res.proportion_explained['PC1'] * 100:.2f}%)")
    ax.set_ylabel(f"PC2 ({ord_res.proportion_explained['PC2'] * 100:.2f}%)")
    ax.set_title("PCoA (Bray-Curtis from Deblur feature table)")
    fig.tight_layout()
    fig.savefig(plot_fp, dpi=150)

    print(f"Finished. Outputs written under: {args.results_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
