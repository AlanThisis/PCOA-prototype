#!/usr/bin/env python3
"""Run UniFrac analysis using Greengenes2 plus DART or QIIME 2.

Requires:
  - A QIIME2 amplicon environment with q2-greengenes2 installed
    (see https://docs.qiime2.org and https://github.com/biocore/q2-greengenes2)

One-time GG2 artifact downloads (~175 MB total):
  mkdir -p data/gg2
  wget -P data/gg2 https://ftp.microbio.me/greengenes_release/current/2024.09.backbone.full-length.fna.qza
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
import json
import os
import shutil
import tempfile
from pathlib import Path

DEFAULT_CACHE_DIR = Path(tempfile.gettempdir()) / "pcoa-prototype-cache"
for env_name, dirname in [
    ("MPLCONFIGDIR", "matplotlib"),
    ("NUMBA_CACHE_DIR", "numba"),
    ("XDG_CACHE_HOME", "xdg"),
]:
    env_path = DEFAULT_CACHE_DIR / dirname
    env_path.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault(env_name, str(env_path))

import biom
import matplotlib.pyplot as plt
import skbio

from dart_unifrac import (
    DART_METHOD,
    DART_PCOA_DIMENSIONS,
    dart_version,
    prepare_dart_inputs,
    run_dart_unifrac,
    validate_dart_parameters,
)

from pipeline_lib import (
    TimingRecorder,
    add_timing_argument,
    resolve_executable,
    run_command,
    run_timed_main,
)


GG2_BACKBONE_FILENAME = "2024.09.backbone.full-length.fna.qza"
GG2_ID_TREE_FILENAME = "2024.09.phylogeny.id.nwk.qza"
GG2_FTP_BASE = "https://ftp.microbio.me/greengenes_release/current/"
# DART --weighted computes normalized weighted UniFrac, so the QIIME backend
# must use weighted_normalized_unifrac for the two engines to agree.
QIIME_UNIFRAC_METRICS = {
    "unweighted": "unweighted_unifrac",
    "weighted": "weighted_normalized_unifrac",
}


def metric_label(metric: str) -> str:
    """Filename stem for a UniFrac metric, e.g. 'weighted_unifrac'."""
    return f"{metric}_unifrac"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute GG2-mapped UniFrac PCoA using DART or QIIME2.",
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
        default=1000,
        help="Rarefaction depth after backbone mapping (default: 1000).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Threads for vsearch and UniFrac beta diversity (default: 4).",
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=None,
        help="Directory for intermediate QIIME2 artifacts (default: <results-dir>/../<name>_unifrac_work).",
    )
    parser.add_argument(
        "--refresh-input-artifacts",
        action="store_true",
        help=(
            "Re-import the feature table and representative sequences and rerun "
            "GG2 mapping. Used when upstream Deblur or merge outputs were regenerated."
        ),
    )
    parser.add_argument(
        "--unifrac-engine",
        choices=("dart", "qiime"),
        default="dart",
        help="UniFrac/PCoA backend after QIIME GG2 mapping and rarefaction (default: dart).",
    )
    parser.add_argument(
        "--metric",
        choices=tuple(QIIME_UNIFRAC_METRICS),
        default="unweighted",
        help=(
            "UniFrac variant (default: unweighted). 'weighted' is normalized "
            "weighted UniFrac on both engines."
        ),
    )
    parser.add_argument(
        "--dartunifrac-executable",
        default="dartunifrac",
        help="DartUniFrac executable name or path (default: dartunifrac).",
    )
    parser.add_argument(
        "--dart-sketch-size",
        type=int,
        default=2048,
        help="DART DMH sketch size (default: 2048).",
    )
    parser.add_argument(
        "--dart-seed",
        type=int,
        default=1337,
        help="DART sketch/fPCoA random seed (default: 1337).",
    )
    parser.add_argument(
        "--dart-bbits",
        type=int,
        choices=(16, 32, 64),
        default=16,
        help="DART sketch register width (default: 16).",
    )
    parser.add_argument(
        "--pcoa-method",
        choices=("auto", "eigh", "fsvd"),
        default="auto",
        help="PCoA eigensolver. Auto uses the memory budget to choose exact or FSVD.",
    )
    parser.add_argument(
        "--pcoa-dimensions",
        type=int,
        default=10,
        help="Axes requested from FSVD (default: 10).",
    )
    parser.add_argument(
        "--pcoa-memory-budget-gb",
        type=float,
        help="Memory budget for auto PCoA selection; defaults to SLURM allocation or available RAM.",
    )
    parser.add_argument(
        "--export-distance-tsv",
        choices=("auto", "always", "never"),
        default="auto",
        help="Export the distance QZA to TSV; auto skips projected exports over 20 GiB.",
    )
    add_timing_argument(parser)
    return parser.parse_args()


def auto_sampling_depth(biom_fp: Path) -> int:
    table = biom.load_table(str(biom_fp))
    sample_depths = table.sum(axis="sample")
    positive_depths = sample_depths[sample_depths > 0]
    if positive_depths.size == 0:
        raise ValueError(
            "Cannot choose a rarefaction depth: no samples retained reads after "
            "Greengenes2 backbone mapping."
        )
    return int(positive_depths.min())


def system_available_memory_gb() -> float:
    slurm_mb = os.environ.get("SLURM_MEM_PER_NODE")
    if slurm_mb:
        try:
            return float(slurm_mb) / 1024
        except ValueError:
            pass
    slurm_per_cpu_mb = os.environ.get("SLURM_MEM_PER_CPU")
    slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm_per_cpu_mb and slurm_cpus:
        try:
            return float(slurm_per_cpu_mb) * int(slurm_cpus) / 1024
        except ValueError:
            pass
    try:
        page_size = os.sysconf("SC_PAGE_SIZE")
        available_pages = os.sysconf("SC_AVPHYS_PAGES")
        return page_size * available_pages / (1024**3)
    except (AttributeError, OSError, ValueError):
        return 4.0


def choose_pcoa_method(
    requested: str,
    sample_count: int,
    memory_budget_gb: float | None,
) -> dict[str, float | int | str]:
    budget_gb = memory_budget_gb or system_available_memory_gb()
    estimate_bytes = 72 * sample_count * sample_count
    estimate_gb = estimate_bytes / (1024**3)
    method = requested
    if requested == "auto":
        method = "eigh" if estimate_gb <= 0.8 * budget_gb else "fsvd"
    return {
        "requested_method": requested,
        "selected_method": method,
        "sample_count": sample_count,
        "estimated_exact_memory_gb": estimate_gb,
        "memory_budget_gb": budget_gb,
        "auto_budget_fraction": 0.8,
    }


def should_export_distance_tsv(policy: str, sample_count: int) -> tuple[bool, float]:
    projected_gib = (16 * sample_count * sample_count) / (1024**3)
    if policy == "always":
        return True, projected_gib
    if policy == "never":
        return False, projected_gib
    return projected_gib <= 20, projected_gib


def import_artifact(
    qiime: str,
    input_path: Path,
    artifact_type: str,
    output_qza: Path,
    input_format: str | None = None,
    timing: TimingRecorder | None = None,
    step: str = "import_artifact",
) -> None:
    if output_qza.exists():
        print(f"  skipping (exists): {output_qza.name}")
        if timing is not None:
            timing.skipped(step, item=str(output_qza), message="artifact exists")
        return
    cmd = [
        qiime, "tools", "import",
        "--type", artifact_type,
        "--input-path", str(input_path),
        "--output-path", str(output_qza),
    ]
    if input_format:
        cmd += ["--input-format", input_format]
    run_command(
        cmd,
        timing=timing,
        step=step,
        item=str(output_qza),
    )


def export_artifact(
    qiime: str,
    artifact_fp: Path,
    export_dir: Path,
    *,
    timing: TimingRecorder | None = None,
    step: str = "export_artifact",
) -> None:
    if export_dir.exists():
        shutil.rmtree(export_dir)
    export_dir.parent.mkdir(parents=True, exist_ok=True)
    run_command(
        [
            qiime,
            "tools",
            "export",
            "--input-path",
            str(artifact_fp),
            "--output-path",
            str(export_dir),
        ],
        timing=timing,
        step=step,
        item=str(artifact_fp),
    )


def run_qiime2_unifrac(
    biom_fp: Path,
    seqs_fp: Path,
    gg2_backbone_fp: Path,
    gg2_tree_fp: Path,
    sampling_depth: int | None,
    threads: int,
    work_dir: Path,
    qiime: str,
    timing: TimingRecorder,
    pcoa_method: str = "auto",
    pcoa_dimensions: int = 10,
    pcoa_memory_budget_gb: float | None = None,
    metric: str = "unweighted",
) -> tuple[Path, int, dict[str, float | int | str]]:
    """Run QIIME2 import → backbone mapping → rarefied UniFrac PCoA.

    Returns (work dir containing UniFrac artifacts, sampling depth used).
    """
    sampling_depth, rarefied_stats = prepare_qiime2_rarefied_inputs(
        biom_fp=biom_fp,
        seqs_fp=seqs_fp,
        gg2_backbone_fp=gg2_backbone_fp,
        sampling_depth=sampling_depth,
        threads=threads,
        work_dir=work_dir,
        qiime=qiime,
        timing=timing,
    )
    label = metric_label(metric)
    distance_matrix_qza = work_dir / f"{label}_distance_matrix.qza"
    pcoa_qza = work_dir / f"{label}_pcoa_results.qza"

    retained_sample_count = rarefied_stats["sample_count"]
    pcoa_policy = choose_pcoa_method(
        pcoa_method, retained_sample_count, pcoa_memory_budget_gb
    )
    effective_pcoa_dimensions = min(pcoa_dimensions, retained_sample_count)
    pcoa_policy["requested_dimensions"] = pcoa_dimensions
    pcoa_policy["dimensions"] = effective_pcoa_dimensions
    selected_pcoa_method = str(pcoa_policy["selected_method"])
    print(
        f"  PCoA method: {selected_pcoa_method} "
        f"(exact estimate {pcoa_policy['estimated_exact_memory_gb']:.2f} GiB; "
        f"budget {pcoa_policy['memory_budget_gb']:.2f} GiB)",
        flush=True,
    )

    print(f"Step 5/5: Computing {metric} UniFrac distance matrix and PCoA...")
    if distance_matrix_qza.exists():
        distance_matrix_qza.unlink()
    run_command(
        [
            qiime,
            "diversity",
            "beta-phylogenetic",
            "--i-phylogeny",
            str(gg2_tree_fp),
            "--i-table",
            str(work_dir / "rarefied-backbone-mapped-table.qza"),
            "--p-metric",
            QIIME_UNIFRAC_METRICS[metric],
            "--p-threads",
            str(threads),
            "--o-distance-matrix",
            str(distance_matrix_qza),
        ],
        timing=timing,
        step=label,
        item=str(distance_matrix_qza),
    )

    if pcoa_qza.exists():
        pcoa_qza.unlink()
    pcoa_command = [
        qiime,
        "diversity",
        "pcoa",
        "--i-distance-matrix",
        str(distance_matrix_qza),
        "--o-pcoa",
        str(pcoa_qza),
    ]
    if selected_pcoa_method == "fsvd":
        pcoa_command.extend(
            ("--p-number-of-dimensions", str(effective_pcoa_dimensions))
        )
    run_command(
        pcoa_command,
        timing=timing,
        step="pcoa",
        item=str(pcoa_qza),
    )

    return work_dir, sampling_depth, pcoa_policy


def prepare_qiime2_rarefied_inputs(
    biom_fp: Path,
    seqs_fp: Path,
    gg2_backbone_fp: Path,
    sampling_depth: int | None,
    threads: int,
    work_dir: Path,
    qiime: str,
    timing: TimingRecorder,
) -> tuple[int, dict[str, int]]:
    """Import, map to the GG2 backbone, and rarefy a Deblur feature table."""
    work_dir.mkdir(parents=True, exist_ok=True)

    table_qza = work_dir / "table.qza"
    rep_seqs_qza = work_dir / "rep-seqs.qza"
    backbone_table_qza = work_dir / "backbone-mapped-table.qza"
    backbone_reps_qza = work_dir / "backbone-representatives.qza"
    rarefied_table_qza = work_dir / "rarefied-backbone-mapped-table.qza"

    print("Step 1/5: Importing feature table...")
    import_artifact(
        qiime,
        biom_fp,
        "FeatureTable[Frequency]",
        table_qza,
        input_format="BIOMV210Format",
        timing=timing,
        step="import_feature_table",
    )

    print("Step 2/5: Importing representative sequences...")
    import_artifact(
        qiime,
        seqs_fp,
        "FeatureData[Sequence]",
        rep_seqs_qza,
        timing=timing,
        step="import_representative_sequences",
    )

    print("Step 3/5: Mapping ASVs to GG2 backbone via closed-reference OTU picking...")
    if backbone_table_qza.exists() and backbone_reps_qza.exists():
        print(f"  skipping (exists): {backbone_table_qza.name}")
        timing.skipped(
            "gg2_non_v4_16s_mapping",
            item=str(backbone_table_qza),
            message="mapped table and representatives exist",
        )
    else:
        run_command(
            [
                qiime,
                "greengenes2",
                "non-v4-16s",
                "--i-table",
                str(table_qza),
                "--i-sequences",
                str(rep_seqs_qza),
                "--i-backbone",
                str(gg2_backbone_fp),
                "--p-threads",
                str(threads),
                "--o-mapped-table",
                str(backbone_table_qza),
                "--o-representatives",
                str(backbone_reps_qza),
            ],
            timing=timing,
            step="gg2_non_v4_16s_mapping",
            item=str(backbone_table_qza),
        )

    # Determine sampling depth from backbone-mapped table (may differ from raw table)
    if sampling_depth is None:
        sampling_depth = auto_sampling_depth_from_qza(
            qiime, backbone_table_qza, work_dir, timing
        )
    else:
        timing.skipped(
            "determine_sampling_depth",
            item=str(sampling_depth),
            message="sampling depth supplied by --sampling-depth",
        )
    print(f"  Sampling depth: {sampling_depth}")

    print(f"Step 4/5: Rarefying mapped table (depth={sampling_depth})...")
    if rarefied_table_qza.exists():
        rarefied_table_qza.unlink()
    run_command(
        [
            qiime,
            "feature-table",
            "rarefy",
            "--i-table",
            str(backbone_table_qza),
            "--p-sampling-depth",
            str(sampling_depth),
            "--o-rarefied-table",
            str(rarefied_table_qza),
        ],
        timing=timing,
        step="rarefy_feature_table",
        item=f"depth={sampling_depth}",
    )

    rarefied_stats = table_stats_from_qza(
        qiime,
        rarefied_table_qza,
        work_dir / "_rarefied_table_export",
        timing,
        "export_rarefied_table_for_policy",
    )
    retained_sample_count = rarefied_stats["sample_count"]
    if retained_sample_count < 2:
        raise ValueError(
            "PCoA requires at least two samples after rarefaction; "
            f"found {retained_sample_count}"
        )
    return sampling_depth, rarefied_stats


def table_stats_from_qza(
    qiime: str,
    table_qza: Path,
    export_dir: Path,
    timing: TimingRecorder,
    step: str,
) -> dict[str, int]:
    export_artifact(qiime, table_qza, export_dir, timing=timing, step=step)
    table = biom.load_table(str(export_dir / "feature-table.biom"))
    return {"sample_count": table.shape[1], "feature_count": table.shape[0]}


def auto_sampling_depth_from_qza(
    qiime: str,
    table_qza: Path,
    work_dir: Path,
    timing: TimingRecorder,
) -> int:
    export_dir = work_dir / "_depth_check_export"
    export_artifact(
        qiime,
        table_qza,
        export_dir,
        timing=timing,
        step="export_mapped_table_for_depth",
    )
    with timing.step(
        "calculate_sampling_depth", item=str(export_dir / "feature-table.biom")
    ):
        return auto_sampling_depth(export_dir / "feature-table.biom")


def configure_writable_caches(work_dir: Path) -> None:
    cache_dir = work_dir / "_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    for env_name, dirname in [
        ("NUMBA_CACHE_DIR", "numba"),
        ("MPLCONFIGDIR", "matplotlib"),
        ("XDG_CACHE_HOME", "xdg"),
    ]:
        env_path = cache_dir / dirname
        env_path.mkdir(parents=True, exist_ok=True)
        os.environ[env_name] = str(env_path)


def clear_input_artifact_cache(work_dir: Path) -> None:
    for filename in (
        "table.qza",
        "rep-seqs.qza",
        "backbone-mapped-table.qza",
        "backbone-representatives.qza",
    ):
        (work_dir / filename).unlink(missing_ok=True)
    shutil.rmtree(work_dir / "_depth_check_export", ignore_errors=True)


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


def run(args: argparse.Namespace, timing: TimingRecorder) -> int:
    if args.pcoa_dimensions < 2:
        raise ValueError("--pcoa-dimensions must be at least two")
    if args.pcoa_memory_budget_gb is not None and args.pcoa_memory_budget_gb <= 0:
        raise ValueError("--pcoa-memory-budget-gb must be greater than zero")
    if args.unifrac_engine == "dart":
        validate_dart_parameters(
            sketch_size=args.dart_sketch_size,
            seed=args.dart_seed,
            bbits=args.dart_bbits,
            pcoa_dimensions=args.pcoa_dimensions,
        )
        if args.pcoa_method != "auto":
            raise ValueError(
                "--pcoa-method applies to the QIIME backend; use 'auto' with DART"
            )
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

    configure_writable_caches(work_dir)
    if args.refresh_input_artifacts:
        clear_input_artifact_cache(work_dir)

    args.results_dir.mkdir(parents=True, exist_ok=True)
    label = metric_label(args.metric)
    distance_result = args.results_dir / f"distance_matrix_{label}.tsv"
    if args.unifrac_engine == "qiime":
        unifrac_work_dir, sampling_depth, pcoa_policy = run_qiime2_unifrac(
            biom_fp=biom_fp,
            seqs_fp=seqs_fp,
            gg2_backbone_fp=gg2_backbone_fp,
            gg2_tree_fp=gg2_tree_fp,
            sampling_depth=args.sampling_depth,
            threads=args.threads,
            work_dir=work_dir,
            qiime=qiime,
            timing=timing,
            pcoa_method=args.pcoa_method,
            pcoa_dimensions=args.pcoa_dimensions,
            pcoa_memory_budget_gb=args.pcoa_memory_budget_gb,
            metric=args.metric,
        )
        export_distance, projected_distance_gib = should_export_distance_tsv(
            args.export_distance_tsv, int(pcoa_policy["sample_count"])
        )
        if export_distance:
            print("Exporting UniFrac distance matrix...")
            dm_export_dir = work_dir / f"{label}_dm_export"
            export_artifact(
                qiime,
                unifrac_work_dir / f"{label}_distance_matrix.qza",
                dm_export_dir,
                timing=timing,
                step="export_distance_matrix",
            )
            with timing.step("copy_distance_matrix"):
                shutil.copy(dm_export_dir / "distance-matrix.tsv", distance_result)
        else:
            distance_result.unlink(missing_ok=True)
            timing.skipped(
                "export_distance_matrix",
                item=str(distance_result),
                message=(
                    f"policy={args.export_distance_tsv}; "
                    f"projected={projected_distance_gib:.2f} GiB"
                ),
            )

        print("Exporting UniFrac PCoA...")
        pcoa_export_dir = work_dir / f"{label}_pcoa_export"
        export_artifact(
            qiime,
            unifrac_work_dir / f"{label}_pcoa_results.qza",
            pcoa_export_dir,
            timing=timing,
            step="export_pcoa",
        )
        ordination_fp = pcoa_export_dir / "ordination.txt"
        distance_backend_path = work_dir / f"{label}_distance_matrix.qza"
    else:
        sampling_depth, rarefied_stats = prepare_qiime2_rarefied_inputs(
            biom_fp=biom_fp,
            seqs_fp=seqs_fp,
            gg2_backbone_fp=gg2_backbone_fp,
            sampling_depth=args.sampling_depth,
            threads=args.threads,
            work_dir=work_dir,
            qiime=qiime,
            timing=timing,
        )
        export_distance, projected_distance_gib = should_export_distance_tsv(
            args.export_distance_tsv, rarefied_stats["sample_count"]
        )
        dart_work_dir = work_dir / "dart"
        dart_biom, dart_tree, extracted_stats = prepare_dart_inputs(
            work_dir / "rarefied-backbone-mapped-table.qza",
            gg2_tree_fp,
            dart_work_dir,
        )
        if extracted_stats != rarefied_stats:
            raise RuntimeError(
                "Extracted DART BIOM dimensions do not match the rarefied QIIME table"
            )
        dart_executable = resolve_executable(args.dartunifrac_executable)
        detected_dart_version = dart_version(dart_executable)
        dart_result = run_dart_unifrac(
            executable=dart_executable,
            biom_fp=dart_biom,
            tree_fp=dart_tree,
            work_dir=dart_work_dir,
            threads=args.threads,
            sketch_size=args.dart_sketch_size,
            seed=args.dart_seed,
            bbits=args.dart_bbits,
            compress=not export_distance,
            timing=timing,
            weighted=args.metric == "weighted",
        )
        if int(dart_result["sample_count"]) != rarefied_stats["sample_count"]:
            raise RuntimeError(
                "DartUniFrac ordination sample count does not match the rarefied table"
            )
        ordination_fp = Path(str(dart_result["ordination_path"]))
        distance_backend_path = Path(str(dart_result["distance_path"]))
        if export_distance:
            with timing.step("copy_distance_matrix"):
                shutil.copy(distance_backend_path, distance_result)
        else:
            distance_result.unlink(missing_ok=True)
            timing.skipped(
                "export_distance_matrix",
                item=str(distance_result),
                message=(
                    f"policy={args.export_distance_tsv}; compressed DART matrix "
                    f"retained at {distance_backend_path}"
                ),
            )
        pcoa_policy = {
            "requested_method": "dart-fpcoa",
            "selected_method": "dart-fpcoa",
            "sample_count": rarefied_stats["sample_count"],
            "requested_dimensions": DART_PCOA_DIMENSIONS,
            "dimensions": int(dart_result["dimensions"]),
        }

    with timing.step("copy_pcoa_coordinates"):
        shutil.copy(
            ordination_fp,
            args.results_dir / f"pcoa_coordinates_{label}.txt",
        )

    print("Generating PCoA plot...")
    with timing.step("plot_unlabeled_pcoa"):
        plot_pcoa(
            ordination_fp=ordination_fp,
            plot_fp=args.results_dir / f"pcoa_plot_{label}.png",
            title=(
                f"PCoA — {args.metric.capitalize()} UniFrac "
                f"(Greengenes2 backbone, {args.unifrac_engine}, depth={sampling_depth})"
            ),
        )

    mapped_stats = table_stats_from_qza(
        qiime,
        work_dir / "backbone-mapped-table.qza",
        work_dir / "_mapped_table_summary_export",
        timing,
        "export_mapped_table_for_summary",
    )
    analysis_summary = {
        "unifrac_engine": args.unifrac_engine,
        "metric": QIIME_UNIFRAC_METRICS[args.metric],
        "sampling_depth": sampling_depth,
        "mapped_samples": mapped_stats["sample_count"],
        "mapped_features": mapped_stats["feature_count"],
        "rarefied_samples": pcoa_policy["sample_count"],
        "pcoa": pcoa_policy,
        "pcoa_dimensions": args.pcoa_dimensions,
        "distance_tsv": {
            "policy": args.export_distance_tsv,
            "exported": export_distance,
            "projected_size_gib": projected_distance_gib,
            "backend_path": str(distance_backend_path),
        },
    }
    if args.unifrac_engine == "dart":
        analysis_summary["dart"] = {
            "version": detected_dart_version,
            "method": DART_METHOD,
            "weighted": args.metric == "weighted",
            "sketch_size": args.dart_sketch_size,
            "seed": args.dart_seed,
            "bbits": args.dart_bbits,
            "threads": args.threads,
        }
    with (args.results_dir / "analysis_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(analysis_summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(f"\nFinished. Outputs written under: {args.results_dir}")
    return 0


def main() -> int:
    args = parse_args()
    timing = TimingRecorder(args.timings_tsv, component="unifrac")
    return run_timed_main(timing, lambda: run(args, timing))


if __name__ == "__main__":
    raise SystemExit(main())
