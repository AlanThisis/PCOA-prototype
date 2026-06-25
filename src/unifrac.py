#!/usr/bin/env python3
"""Run UniFrac phylogenetic diversity analysis using Greengenes2 and QIIME2.

Requires:
  - QIIME2 amplicon environment active (qiime2-amplicon-2026.1)
  - q2-greengenes2 plugin installed in that environment

One-time GG2 artifact downloads (~175 MB total):
  mkdir -p data/gg2
  wget -P data/gg2 https://ftp.microbio.me/greengenes_release/current/2024.09.backbone.v4.fna.qza
  wget -P data/gg2 https://ftp.microbio.me/greengenes_release/current/2024.09.phylogeny.id.nwk.qza

Approach:
  Because this pipeline processes ENA studies not deposited in Qiita, the ASVs
  are absent from the pre-placed GG2 ASV tree. Instead, this script uses GG2's
  non-v4-16s action (closed-reference OTU picking via vsearch at 99% identity)
  to map de novo Deblur ASVs against the GG2 backbone (~10K reference genomes).
  The backbone-mapped table is then used with the GG2 ID phylogeny for UniFrac.
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


GG2_BACKBONE_FILENAME = "2024.09.backbone.full-length.fna.qza"
GG2_ID_TREE_FILENAME = "2024.09.phylogeny.id.nwk.qza"
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
        help="Rarefaction depth (default: minimum sample depth after backbone mapping).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Threads for vsearch and core-metrics (default: 4).",
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


def import_artifact(
    qiime: str,
    input_path: Path,
    artifact_type: str,
    output_qza: Path,
    input_format: str | None = None,
) -> None:
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
    gg2_backbone_fp: Path,
    gg2_tree_fp: Path,
    sampling_depth: int | None,
    threads: int,
    work_dir: Path,
    qiime: str,
) -> tuple[Path, int]:
    """Run QIIME2 import → backbone mapping → core-metrics.

    Returns (core-metrics output dir, sampling depth used).
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    table_qza = work_dir / "table.qza"
    rep_seqs_qza = work_dir / "rep-seqs.qza"
    backbone_table_qza = work_dir / "backbone-mapped-table.qza"
    backbone_reps_qza = work_dir / "backbone-representatives.qza"
    metadata_tsv = work_dir / "metadata.tsv"
    core_metrics_dir = work_dir / "core-metrics"

    print("Step 1/4: Importing feature table...")
    import_artifact(qiime, biom_fp, "FeatureTable[Frequency]", table_qza, input_format="BIOMV210Format")

    print("Step 2/4: Importing representative sequences...")
    import_artifact(qiime, seqs_fp, "FeatureData[Sequence]", rep_seqs_qza)

    print("Step 3/4: Mapping ASVs to GG2 backbone via closed-reference OTU picking...")
    if backbone_table_qza.exists() and backbone_reps_qza.exists():
        print(f"  skipping (exists): {backbone_table_qza.name}")
    else:
        run_command([
            qiime, "greengenes2", "non-v4-16s",
            "--i-table", str(table_qza),
            "--i-sequences", str(rep_seqs_qza),
            "--i-backbone", str(gg2_backbone_fp),
            "--p-threads", str(threads),
            "--o-mapped-table", str(backbone_table_qza),
            "--o-representatives", str(backbone_reps_qza),
        ])

    # Determine sampling depth from backbone-mapped table (may differ from raw table)
    if sampling_depth is None:
        sampling_depth = auto_sampling_depth_from_qza(qiime, backbone_table_qza, work_dir)
    print(f"  Sampling depth: {sampling_depth}")

    write_minimal_metadata(biom_fp, metadata_tsv)

    print(f"Step 4/4: Running core-metrics-phylogenetic (depth={sampling_depth})...")
    if core_metrics_dir.exists():
        shutil.rmtree(core_metrics_dir)
    run_command([
        qiime, "diversity", "core-metrics-phylogenetic",
        "--i-phylogeny", str(gg2_tree_fp),
        "--i-table", str(backbone_table_qza),
        "--p-sampling-depth", str(sampling_depth),
        "--m-metadata-file", str(metadata_tsv),
        "--p-n-jobs-or-threads", str(threads),
        "--output-dir", str(core_metrics_dir),
    ])

    return core_metrics_dir, sampling_depth


def auto_sampling_depth_from_qza(qiime: str, table_qza: Path, work_dir: Path) -> int:
    export_dir = work_dir / "_depth_check_export"
    export_artifact(qiime, table_qza, export_dir)
    return auto_sampling_depth(export_dir / "feature-table.biom")


def plot_pcoa(ordination_fp: Path, plot_fp: Path, title: str) -> None:
    ord_res = skbio.io.read(
        str(ordination_fp),
        format="ordination",
        into=skbio.stats.ordination.OrdinationResults,
    )
    coords = ord_res.samples
    prop = ord_res.proportion_explained

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(coords.iloc[:, 0], coords.iloc[:, 1], alpha=0.7, s=30)
    for sample_id, row in coords.iterrows():
        ax.annotate(sample_id, (row.iloc[0], row.iloc[1]), fontsize=5, alpha=0.6)
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
    gg2_backbone_fp = args.gg2_dir / GG2_BACKBONE_FILENAME
    gg2_tree_fp = args.gg2_dir / GG2_ID_TREE_FILENAME

    for fp in [biom_fp, seqs_fp]:
        if not fp.exists():
            raise FileNotFoundError(f"Required input not found: {fp}")

    missing = [fp for fp in [gg2_backbone_fp, gg2_tree_fp] if not fp.exists()]
    if missing:
        raise FileNotFoundError(
            f"GG2 artifact(s) not found: {[str(f) for f in missing]}\n"
            f"Download with:\n"
            f"  mkdir -p {args.gg2_dir}\n"
            f"  wget -P {args.gg2_dir} {GG2_FTP_BASE}{GG2_BACKBONE_FILENAME}\n"
            f"  wget -P {args.gg2_dir} {GG2_FTP_BASE}{GG2_ID_TREE_FILENAME}"
        )

    qiime = resolve_executable("qiime")

    # Ensure R_HOME is set so rpy2 (loaded by q2-composition at startup) can
    # find R shared libs. The conda activation script sets this on env activate;
    # derive it from the qiime executable as a fallback for non-interactive use.
    if not os.environ.get("R_HOME"):
        r_home = Path(qiime).resolve().parent.parent / "lib" / "R"
        if r_home.is_dir():
            os.environ["R_HOME"] = str(r_home)

    # Also ensure the env bin is on PATH so QIIME2 subprocesses (vsearch, etc.)
    # can be found without requiring the user to have activated the env.
    env_bin = str(Path(qiime).resolve().parent)
    if env_bin not in os.environ.get("PATH", ""):
        os.environ["PATH"] = env_bin + os.pathsep + os.environ.get("PATH", "")

    if args.work_dir is not None:
        work_dir = args.work_dir.resolve()
    else:
        work_dir = (args.results_dir.parent / f"{args.results_dir.name}_unifrac_work").resolve()

    core_metrics_dir, sampling_depth = run_qiime2_unifrac(
        biom_fp=biom_fp,
        seqs_fp=seqs_fp,
        gg2_backbone_fp=gg2_backbone_fp,
        gg2_tree_fp=gg2_tree_fp,
        sampling_depth=args.sampling_depth,
        threads=args.threads,
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
        title=f"PCoA — Unweighted UniFrac (Greengenes2 backbone, depth={sampling_depth})",
    )

    print(f"\nFinished. Outputs written under: {args.results_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
