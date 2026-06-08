#!/usr/bin/env python3
"""Run UniFrac phylogenetic diversity analysis using Greengenes2 and QIIME2.

Requires:
  - QIIME2 amplicon environment active (qiime2-amplicon-2026.1)
  - q2-greengenes2 plugin installed in that environment
  - GG2 tree artifact downloaded to --gg2-dir (see below)

Download GG2 tree (one-time, ~308 MB):
  mkdir -p data/gg2
  wget -P data/gg2 https://ftp.microbio.me/greengenes_release/current/2024.09.phylogeny.asv.nwk.qza
"""
from __future__ import annotations

import argparse
import os
import shutil
from pathlib import Path

import biom
import matplotlib.pyplot as plt
import pandas as pd
import skbio

from pipeline_lib import resolve_executable, run_command


GG2_TREE_FILENAME = "2024.09.phylogeny.asv.nwk.qza"
GG2_FTP_BASE = "https://ftp.microbio.me/greengenes_release/current/"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute UniFrac PCoA using Greengenes2 and QIIME2.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--deblur-dir",
        type=Path,
        required=True,
        help="Deblur workflow output dir containing all.biom and all.seqs.fa.",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        required=True,
        help="Directory to write final TSV and PNG outputs.",
    )
    parser.add_argument(
        "--gg2-dir",
        type=Path,
        default=Path("data/gg2"),
        help="Directory containing GG2 .qza artifacts (default: data/gg2).",
    )
    parser.add_argument(
        "--sampling-depth",
        type=int,
        default=None,
        help="Rarefaction depth (default: minimum sample depth in table).",
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=None,
        help="Directory for intermediate QIIME2 artifacts (default: <results-dir>/../<name>_unifrac_work).",
    )
    return parser.parse_args()


def auto_sampling_depth(biom_fp: Path) -> int:
    table = biom.load_table(str(biom_fp))
    return int(table.sum(axis="sample").min())


def write_minimal_metadata(biom_fp: Path, metadata_fp: Path) -> None:
    table = biom.load_table(str(biom_fp))
    pd.DataFrame({"sample-id": list(table.ids())}).to_csv(
        metadata_fp, sep="\t", index=False
    )


def import_artifact(qiime: str, input_path: Path, artifact_type: str, output_qza: Path, input_format: str | None = None) -> None:
    if output_qza.exists():
        print(f"  skipping (exists): {output_qza.name}")
        return
    cmd = [
        qiime, "tools", "import",
        "--type", artifact_type,
        "--input-path", str(input_path),
        "--output-path", str(output_qza),
    ]
    if input_format:
        cmd += ["--input-format", input_format]
    run_command(cmd)


def export_artifact(qiime: str, artifact_fp: Path, export_dir: Path) -> None:
    export_dir.mkdir(parents=True, exist_ok=True)
    run_command([
        qiime, "tools", "export",
        "--input-path", str(artifact_fp),
        "--output-path", str(export_dir),
    ])


def run_qiime2_unifrac(
    biom_fp: Path,
    seqs_fp: Path,
    gg2_tree_fp: Path,
    sampling_depth: int,
    work_dir: Path,
    qiime: str,
) -> Path:
    """Run QIIME2 import → filter → core-metrics. Returns core-metrics output dir."""
    work_dir.mkdir(parents=True, exist_ok=True)

    table_qza = work_dir / "table.qza"
    rep_seqs_qza = work_dir / "rep-seqs.qza"
    filtered_table_qza = work_dir / "gg2-filtered-table.qza"
    metadata_tsv = work_dir / "metadata.tsv"
    core_metrics_dir = work_dir / "core-metrics"

    print("Step 1/4: Importing feature table...")
    import_artifact(qiime, biom_fp, "FeatureTable[Frequency]", table_qza, input_format="BIOMV210Format")

    print("Step 2/4: Importing representative sequences...")
    import_artifact(qiime, seqs_fp, "FeatureData[Sequence]", rep_seqs_qza)

    print("Step 3/4: Filtering features to GG2 tree coverage...")
    if filtered_table_qza.exists():
        print(f"  skipping (exists): {filtered_table_qza.name}")
    else:
        run_command([
            qiime, "greengenes2", "filter-features",
            "--i-feature-table", str(table_qza),
            "--i-reference", str(gg2_tree_fp),
            "--o-filtered-feature-table", str(filtered_table_qza),
        ])

    write_minimal_metadata(biom_fp, metadata_tsv)

    print(f"Step 4/4: Running core-metrics-phylogenetic (depth={sampling_depth})...")
    if core_metrics_dir.exists():
        shutil.rmtree(core_metrics_dir)
    run_command([
        qiime, "diversity", "core-metrics-phylogenetic",
        "--i-phylogeny", str(gg2_tree_fp),
        "--i-table", str(filtered_table_qza),
        "--p-sampling-depth", str(sampling_depth),
        "--m-metadata-file", str(metadata_tsv),
        "--output-dir", str(core_metrics_dir),
    ])

    return core_metrics_dir


def plot_pcoa(ordination_fp: Path, plot_fp: Path, title: str) -> None:
    ord_res = skbio.io.read(str(ordination_fp), format="ordination", into=skbio.stats.ordination.OrdinationResults)
    coords = ord_res.samples
    prop = ord_res.proportion_explained

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(coords.iloc[:, 0], coords.iloc[:, 1], alpha=0.7, s=30)
    for sample_id, row in coords.iterrows():
        ax.annotate(sample_id, (row.iloc[0], row.iloc[1]), fontsize=5, alpha=0.7)
    ax.set_xlabel(f"PC1 ({prop.iloc[0] * 100:.2f}%)")
    ax.set_ylabel(f"PC2 ({prop.iloc[1] * 100:.2f}%)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(plot_fp, dpi=150)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    args.deblur_dir = args.deblur_dir.resolve()
    args.results_dir = args.results_dir.resolve()
    args.gg2_dir = args.gg2_dir.resolve()

    biom_fp = args.deblur_dir / "all.biom"
    seqs_fp = args.deblur_dir / "all.seqs.fa"
    gg2_tree_fp = args.gg2_dir / GG2_TREE_FILENAME

    for fp in [biom_fp, seqs_fp]:
        if not fp.exists():
            raise FileNotFoundError(f"Required input not found: {fp}")

    if not gg2_tree_fp.exists():
        raise FileNotFoundError(
            f"GG2 tree artifact not found: {gg2_tree_fp}\n"
            f"Download it with:\n"
            f"  mkdir -p {args.gg2_dir}\n"
            f"  wget -P {args.gg2_dir} {GG2_FTP_BASE}{GG2_TREE_FILENAME}"
        )

    qiime = resolve_executable("qiime")

    # Ensure R_HOME is set so rpy2 (used by q2-composition) can find R libs.
    # When the qiime2 env is activated, the conda activation script sets this;
    # derive it from the qiime executable location as a fallback.
    if not os.environ.get("R_HOME"):
        r_home = Path(qiime).resolve().parent.parent / "lib" / "R"
        if r_home.is_dir():
            os.environ["R_HOME"] = str(r_home)

    sampling_depth = args.sampling_depth if args.sampling_depth is not None else auto_sampling_depth(biom_fp)
    print(f"Sampling depth: {sampling_depth}")

    if args.work_dir is not None:
        work_dir = args.work_dir.resolve()
    else:
        work_dir = (args.results_dir.parent / f"{args.results_dir.name}_unifrac_work").resolve()

    core_metrics_dir = run_qiime2_unifrac(
        biom_fp=biom_fp,
        seqs_fp=seqs_fp,
        gg2_tree_fp=gg2_tree_fp,
        sampling_depth=sampling_depth,
        work_dir=work_dir,
        qiime=qiime,
    )

    args.results_dir.mkdir(parents=True, exist_ok=True)

    print("Exporting UniFrac distance matrix...")
    dm_export_dir = work_dir / "unweighted_unifrac_dm_export"
    export_artifact(qiime, core_metrics_dir / "unweighted_unifrac_distance_matrix.qza", dm_export_dir)
    shutil.copy(
        dm_export_dir / "distance-matrix.tsv",
        args.results_dir / "distance_matrix_unweighted_unifrac.tsv",
    )

    print("Exporting UniFrac PCoA...")
    pcoa_export_dir = work_dir / "unweighted_unifrac_pcoa_export"
    export_artifact(qiime, core_metrics_dir / "unweighted_unifrac_pcoa_results.qza", pcoa_export_dir)
    ordination_fp = pcoa_export_dir / "ordination.txt"
    shutil.copy(ordination_fp, args.results_dir / "pcoa_coordinates_unweighted_unifrac.txt")

    print("Generating PCoA plot...")
    plot_pcoa(
        ordination_fp=ordination_fp,
        plot_fp=args.results_dir / "pcoa_plot_unweighted_unifrac.png",
        title="PCoA — Unweighted UniFrac (Greengenes2)",
    )

    print(f"\nFinished. Outputs written under: {args.results_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
