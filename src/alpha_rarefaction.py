#!/usr/bin/env python3
"""Generate a Faith's PD alpha-rarefaction visualization from a mapped table.

The input table should be the Greengenes2 backbone-mapped table produced by
``src/unifrac.py``.  Running on the mapped table makes the diagnostic describe
the same feature space used by the downstream UniFrac analysis.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import tempfile
import zipfile
from pathlib import Path

import biom

DEFAULT_CACHE_DIR = Path(tempfile.gettempdir()) / "pcoa-prototype-cache"
for env_name, dirname in (
    ("MPLCONFIGDIR", "matplotlib"),
    ("NUMBA_CACHE_DIR", "numba"),
    ("XDG_CACHE_HOME", "xdg"),
):
    env_path = DEFAULT_CACHE_DIR / dirname
    env_path.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault(env_name, str(env_path))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from pipeline_lib import (
    TimingRecorder,
    add_timing_argument,
    resolve_executable,
    run_command,
    run_timed_main,
)


DEPTH_COLUMN = re.compile(r"^depth-(\d+)_iter-(\d+)$")
FORWARD_READ_SUFFIX = re.compile(r"(_R?1(_001)?)$")
LEVEL_COLORS = {
    "Full": "#2563EB",
    "50%": "#0D9488",
    "25%": "#EA580C",
    "10%": "#9333EA",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    run_parser = subparsers.add_parser(
        "run", help="Run QIIME 2 Faith's PD alpha rarefaction and plot it."
    )
    _add_run_arguments(run_parser)
    comparison_parser = subparsers.add_parser(
        "compare", help="Combine two or more existing curve-summary TSV files."
    )
    comparison_parser.add_argument(
        "--curve",
        action="append",
        required=True,
        metavar="LABEL=PATH",
        help="Curve label and faith_pd_curve_summary.tsv path; repeat per level.",
    )
    comparison_parser.add_argument("--output", type=Path, required=True)
    comparison_parser.add_argument("--reference-depth", type=int, default=1000)
    group_parser = subparsers.add_parser(
        "group-plot",
        help="Plot an embedded metadata grouping from an existing alpha-rarefaction QZV.",
    )
    group_parser.add_argument("--qzv", type=Path, required=True)
    group_parser.add_argument("--metadata-column", required=True)
    group_parser.add_argument("--output", type=Path, required=True)
    group_parser.add_argument("--reference-depth", type=int, default=1000)
    return parser.parse_args()


def _add_run_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--mapped-table",
        type=Path,
        required=True,
        help="Greengenes2 backbone-mapped FeatureTable[Frequency] QZA.",
    )
    parser.add_argument(
        "--phylogeny",
        type=Path,
        required=True,
        help="Rooted GG2 phylogeny QZA used for Faith's PD.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for the QZV, depth tables, and run summary.",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        help="Optional QIIME 2-compatible sample metadata TSV.",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        help=(
            "Maximum rarefaction depth. By default, use the nearest-rank "
            "percentile selected by --max-depth-percentile, while ensuring "
            "the reference depth is included when the data permit."
        ),
    )
    parser.add_argument(
        "--max-depth-percentile",
        type=float,
        default=90.0,
        help="Mapped-depth percentile used for automatic max depth (default: 90).",
    )
    parser.add_argument(
        "--reference-depth",
        type=int,
        default=1000,
        help="Pipeline depth highlighted in summary tables (default: 1000).",
    )
    parser.add_argument(
        "--min-depth",
        type=int,
        default=1,
        help="Minimum alpha-rarefaction depth (default: 1).",
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=10,
        help="Number of depths evaluated by QIIME 2 (default: 10).",
    )
    parser.add_argument(
        "--iterations",
        type=int,
        default=10,
        help="Repeated rarefactions at each depth (default: 10).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing Faith's PD QZV in the output directory.",
    )
    add_timing_argument(parser)


def nearest_rank_percentile(values: list[int], percentile: float) -> int:
    """Return a deterministic nearest-rank percentile for positive depths."""
    if not values:
        raise ValueError("Cannot calculate a percentile from no values")
    if not 0 < percentile <= 100:
        raise ValueError("--max-depth-percentile must be in (0, 100]")
    ordered = sorted(values)
    rank = max(1, math.ceil((percentile / 100) * len(ordered)))
    return ordered[rank - 1]


def choose_max_depth(
    depths: list[int],
    *,
    requested: int | None,
    percentile: float,
    reference_depth: int,
) -> tuple[int, str]:
    positive = [depth for depth in depths if depth > 0]
    if not positive:
        raise ValueError("Mapped table contains no samples with positive read depth")
    maximum = max(positive)
    if requested is not None:
        if requested <= 0:
            raise ValueError("--max-depth must be greater than zero")
        if requested > maximum:
            raise ValueError(
                f"--max-depth ({requested}) exceeds the maximum mapped sample "
                f"depth ({maximum})"
            )
        return requested, "explicit"
    percentile_depth = nearest_rank_percentile(positive, percentile)
    selected = min(maximum, max(percentile_depth, reference_depth))
    return selected, f"nearest-rank-p{percentile:g}"


def evenly_spaced_depths(min_depth: int, max_depth: int, steps: int) -> list[int]:
    if min_depth < 1:
        raise ValueError("--min-depth must be at least one")
    if max_depth <= min_depth:
        raise ValueError("Maximum depth must be greater than --min-depth")
    if steps < 2:
        raise ValueError("--steps must be at least two")
    if max_depth - min_depth < steps:
        raise ValueError(
            f"--steps ({steps}) exceeds the possible steps between "
            f"{min_depth} and {max_depth}"
        )
    interval = (max_depth - min_depth) / (steps - 1)
    # q2-diversity uses ``np.linspace(..., dtype=int)``, which truncates each
    # intermediate value toward zero. Mirror that behavior so the retention
    # table describes the exact depths represented in the QZV.
    return sorted({int(min_depth + interval * index) for index in range(steps)})


def write_depth_outputs(
    sample_depths: list[tuple[str, int]],
    candidate_depths: list[int],
    reference_depth: int,
    output_dir: Path,
) -> None:
    depth_path = output_dir / "mapped_sample_depths.tsv"
    with depth_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("sample-id", "mapped_reads", "eligible_at_reference_depth"))
        for sample_id, depth in sorted(sample_depths):
            writer.writerow((sample_id, depth, str(depth >= reference_depth).lower()))

    positive_count = sum(depth > 0 for _, depth in sample_depths)
    retention_path = output_dir / "sample_retention_by_depth.tsv"
    with retention_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ("depth", "retained_samples", "positive_samples", "retained_percent")
        )
        for depth in candidate_depths:
            retained = sum(sample_depth >= depth for _, sample_depth in sample_depths)
            percent = 100 * retained / positive_count if positive_count else 0
            writer.writerow((depth, retained, positive_count, f"{percent:.6f}"))


def strip_forward_read_suffix(sample_id: str) -> str:
    return FORWARD_READ_SUFFIX.sub("", sample_id)


def prepare_qiime_metadata(
    metadata_path: Path, table_sample_ids: list[str], output_path: Path
) -> Path:
    """Write metadata whose IDs exactly match the mapped feature table."""
    delimiter = "\t" if metadata_path.suffix.lower() in {".tsv", ".txt"} else ","
    with metadata_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        fields = reader.fieldnames or []
        id_field = next(
            (
                field
                for field in ("sample-id", "sample_id", "run_accession", "run_accessions")
                if field in fields
            ),
            None,
        )
        if id_field is None:
            raise ValueError(
                "Metadata requires sample-id, sample_id, run_accession, or "
                "run_accessions"
            )
        rows = list(reader)

    rows_by_id: dict[str, dict[str, str]] = {}
    for row in rows:
        raw_ids = row.get(id_field, "")
        values = raw_ids.split(";") if id_field == "run_accessions" else [raw_ids]
        for value in values:
            identifier = value.strip()
            if not identifier:
                continue
            if identifier in rows_by_id:
                raise ValueError(f"Duplicate metadata identifier: {identifier}")
            rows_by_id[identifier] = row

    matched: list[tuple[str, dict[str, str]]] = []
    missing: list[str] = []
    for table_id in table_sample_ids:
        row = rows_by_id.get(table_id)
        if row is None:
            row = rows_by_id.get(strip_forward_read_suffix(table_id))
        if row is None:
            missing.append(table_id)
        else:
            matched.append((table_id, row))
    if missing:
        examples = ", ".join(missing[:10])
        raise ValueError(
            f"Metadata is missing {len(missing):,} of {len(table_sample_ids):,} "
            f"mapped-table sample IDs. Examples: {examples}"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_fields = ["sample-id", *[field for field in fields if field != id_field]]
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=output_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for table_id, row in matched:
            output_row = {field: row.get(field, "") for field in output_fields}
            output_row["sample-id"] = table_id
            writer.writerow(output_row)
    return output_path


def find_faith_pd_csv(export_dir: Path) -> Path:
    matches = sorted(export_dir.rglob("faith_pd.csv"))
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one faith_pd.csv in exported QZV, found {len(matches)}"
        )
    return matches[0]


def summarize_faith_pd(csv_path: Path, positive_samples: int) -> pd.DataFrame:
    """Collapse QIIME iterations per sample, then summarize samples by depth."""
    frame = pd.read_csv(csv_path)
    columns_by_depth: dict[int, list[str]] = {}
    for column in frame.columns:
        match = DEPTH_COLUMN.match(str(column))
        if match:
            columns_by_depth.setdefault(int(match.group(1)), []).append(str(column))
    if not columns_by_depth:
        raise ValueError(f"No QIIME alpha-rarefaction depth columns found in {csv_path}")

    rows: list[dict[str, float | int]] = []
    for depth, columns in sorted(columns_by_depth.items()):
        values = frame[columns].apply(pd.to_numeric, errors="coerce")
        # First average repeated rarefactions within each sample. This ensures
        # every retained sample contributes once to the cohort summary.
        sample_means = values.mean(axis=1, skipna=True).dropna()
        retained = int(sample_means.shape[0])
        rows.append(
            {
                "depth": depth,
                "retained_samples": retained,
                "positive_samples": positive_samples,
                "retained_percent": (
                    100 * retained / positive_samples if positive_samples else 0
                ),
                "faith_pd_median": float(sample_means.median()),
                "faith_pd_q25": float(sample_means.quantile(0.25)),
                "faith_pd_q75": float(sample_means.quantile(0.75)),
            }
        )
    return pd.DataFrame(rows)


def _reference_line(ax: plt.Axes, reference_depth: int, depths: pd.Series) -> None:
    if float(depths.min()) <= reference_depth <= float(depths.max()):
        ax.axvline(
            reference_depth,
            color="#374151",
            linestyle="--",
            linewidth=1.2,
            alpha=0.8,
        )


def plot_curve_summary(
    curve: pd.DataFrame, output: Path, reference_depth: int
) -> None:
    fig, (alpha_ax, retention_ax) = plt.subplots(
        2,
        1,
        figsize=(9, 7.5),
        sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1]},
    )
    color = "#0F766E"
    x = curve["depth"].to_numpy()
    median = curve["faith_pd_median"].to_numpy()
    q25 = curve["faith_pd_q25"].to_numpy()
    q75 = curve["faith_pd_q75"].to_numpy()
    alpha_ax.plot(x, median, color=color, linewidth=2.4, marker="o", markersize=4)
    alpha_ax.fill_between(x, q25, q75, color=color, alpha=0.2, linewidth=0)
    alpha_ax.set_ylabel("Faith's PD")
    alpha_ax.set_title("Faith's PD alpha rarefaction")
    alpha_ax.text(
        0.01,
        0.98,
        "Median across samples; band = interquartile range",
        transform=alpha_ax.transAxes,
        va="top",
        fontsize=9,
        color="#4B5563",
    )

    retention_ax.plot(
        x,
        curve["retained_samples"].to_numpy(),
        color="#2563EB",
        linewidth=2.2,
        marker="o",
        markersize=4,
    )
    retention_ax.set_xlabel("Rarefaction depth (mapped reads per sample)")
    retention_ax.set_ylabel("Samples retained")
    positive_samples = int(curve["positive_samples"].iloc[0])
    if positive_samples:
        secondary = retention_ax.secondary_yaxis(
            "right",
            functions=(
                lambda count: 100 * count / positive_samples,
                lambda percent: percent * positive_samples / 100,
            ),
        )
        secondary.set_ylabel("Samples retained (%)")

    for ax in (alpha_ax, retention_ax):
        _reference_line(ax, reference_depth, curve["depth"])
        ax.grid(True, color="#D1D5DB", linewidth=0.7, alpha=0.7)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def parse_curve_spec(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise ValueError(f"Curve must use LABEL=PATH syntax: {value}")
    label, raw_path = value.split("=", 1)
    if not label.strip() or not raw_path.strip():
        raise ValueError(f"Curve must use LABEL=PATH syntax: {value}")
    path = Path(raw_path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Curve summary not found: {path}")
    return label.strip(), path


def plot_curve_comparison(
    curves: list[tuple[str, pd.DataFrame]],
    output: Path,
    reference_depth: int,
) -> None:
    if len(curves) < 2:
        raise ValueError("Comparison requires at least two --curve inputs")
    expected_depths = curves[0][1]["depth"].tolist()
    for label, curve in curves[1:]:
        if curve["depth"].tolist() != expected_depths:
            raise ValueError(
                f"{label} uses a different depth grid; rerun all levels with "
                "the same --min-depth, --max-depth, and --steps"
            )

    fig, (alpha_ax, retention_ax) = plt.subplots(
        2,
        1,
        figsize=(9.5, 8),
        sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1]},
    )
    fallback_colors = plt.get_cmap("tab10")
    for index, (label, curve) in enumerate(curves):
        color = LEVEL_COLORS.get(label, fallback_colors(index % 10))
        alpha_ax.plot(
            curve["depth"],
            curve["faith_pd_median"],
            label=label,
            color=color,
            linewidth=2.3,
            marker="o",
            markersize=3.5,
        )
        retention_ax.plot(
            curve["depth"],
            curve["retained_percent"],
            label=label,
            color=color,
            linewidth=2.1,
            marker="o",
            markersize=3.5,
        )

    alpha_ax.set_ylabel("Median Faith's PD")
    alpha_ax.set_title("Faith's PD alpha rarefaction across FASTQ levels")
    alpha_ax.legend(frameon=False, ncol=min(4, len(curves)))
    retention_ax.set_xlabel("Rarefaction depth (mapped reads per sample)")
    retention_ax.set_ylabel("Samples retained (%)")
    retention_ax.set_ylim(0, 105)
    for ax in (alpha_ax, retention_ax):
        _reference_line(ax, reference_depth, curves[0][1]["depth"])
        ax.grid(True, color="#D1D5DB", linewidth=0.7, alpha=0.7)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def load_group_summary_from_qzv(qzv: Path, metadata_column: str) -> pd.DataFrame:
    """Load QIIME's precomputed grouped rarefaction summary from a QZV."""
    suffix = f"/data/faith_pd-{metadata_column}.jsonp"
    with zipfile.ZipFile(qzv) as archive:
        matches = [name for name in archive.namelist() if name.endswith(suffix)]
        if len(matches) != 1:
            raise ValueError(
                f"Expected one embedded Faith's PD grouping for {metadata_column!r}; "
                f"found {len(matches)}"
            )
        payload = archive.read(matches[0]).decode("utf-8")
    match = re.fullmatch(
        r"load_data\([^,]+,\s*[^,]+,\s*(.*)\);?\s*", payload, flags=re.DOTALL
    )
    if match is None:
        raise ValueError(f"Unexpected grouped alpha-rarefaction payload in {qzv}")
    decoded = json.loads(match.group(1))
    frame = pd.DataFrame(decoded["data"], columns=decoded["columns"])
    required = {
        metadata_column,
        "_alpha_rarefaction_depth_column_",
        "25%",
        "50%",
        "75%",
        "count",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Grouped alpha-rarefaction data lacks: {sorted(missing)}")
    return frame.rename(
        columns={
            metadata_column: "group",
            "_alpha_rarefaction_depth_column_": "depth",
            "25%": "faith_pd_q25",
            "50%": "faith_pd_median",
            "75%": "faith_pd_q75",
            "count": "retained_samples",
        }
    )


def plot_group_summary(
    frame: pd.DataFrame,
    output: Path,
    metadata_column: str,
    reference_depth: int,
) -> None:
    groups = sorted(frame["group"].astype(str).unique())
    colors = plt.get_cmap("tab10")
    fig, (alpha_ax, retention_ax) = plt.subplots(
        2,
        1,
        figsize=(10.5, 8.5),
        sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1]},
    )
    for index, group in enumerate(groups):
        curve = frame.loc[frame["group"].astype(str) == group].sort_values("depth")
        color = colors(index % 10)
        x = pd.to_numeric(curve["depth"]).to_numpy()
        median = pd.to_numeric(curve["faith_pd_median"]).to_numpy()
        q25 = pd.to_numeric(curve["faith_pd_q25"]).to_numpy()
        q75 = pd.to_numeric(curve["faith_pd_q75"]).to_numpy()
        retained = pd.to_numeric(curve["retained_samples"]).to_numpy()
        retained_percent = 100 * retained / retained[0] if retained[0] else retained
        alpha_ax.plot(x, median, label=group, color=color, linewidth=2.2)
        alpha_ax.fill_between(x, q25, q75, color=color, alpha=0.12, linewidth=0)
        retention_ax.plot(
            x, retained_percent, label=group, color=color, linewidth=2.0
        )

    label = (
        "sampling site"
        if metadata_column == "environment_harmonized"
        else metadata_column.replace("_", " ")
    )
    alpha_ax.set_title(f"Faith's PD alpha rarefaction by {label}")
    alpha_ax.set_ylabel("Faith's PD")
    alpha_ax.legend(frameon=False, ncol=min(4, len(groups)), fontsize=9)
    retention_ax.set_xlabel("Rarefaction depth (mapped reads per sample)")
    retention_ax.set_ylabel("Samples retained (%)")
    retention_ax.set_ylim(0, 105)
    for ax in (alpha_ax, retention_ax):
        _reference_line(ax, reference_depth, frame["depth"])
        ax.grid(True, color="#D1D5DB", linewidth=0.7, alpha=0.7)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def group_plot(args: argparse.Namespace) -> int:
    qzv = args.qzv.expanduser().resolve()
    if not qzv.is_file():
        raise FileNotFoundError(f"QZV not found: {qzv}")
    frame = load_group_summary_from_qzv(qzv, args.metadata_column)
    output = args.output.expanduser().resolve()
    plot_group_summary(frame, output, args.metadata_column, args.reference_depth)
    print(f"Grouped alpha-rarefaction plot: {output}")
    return 0


def compare(args: argparse.Namespace) -> int:
    curves = [
        (label, pd.read_csv(path, sep="\t"))
        for label, path in (parse_curve_spec(value) for value in args.curve)
    ]
    plot_curve_comparison(curves, args.output.resolve(), args.reference_depth)
    print(f"Combined alpha-rarefaction plot: {args.output.resolve()}")
    return 0


def configure_qiime_environment(qiime: str, output_dir: Path) -> None:
    env_bin = str(Path(qiime).resolve().parent)
    if env_bin not in os.environ.get("PATH", "").split(os.pathsep):
        os.environ["PATH"] = env_bin + os.pathsep + os.environ.get("PATH", "")
    if not os.environ.get("R_HOME"):
        r_home = Path(qiime).resolve().parent.parent / "lib" / "R"
        if r_home.is_dir():
            os.environ["R_HOME"] = str(r_home)
    cache_root = output_dir / "_cache"
    for env_name, child in (
        ("NUMBA_CACHE_DIR", "numba"),
        ("MPLCONFIGDIR", "matplotlib"),
        ("XDG_CACHE_HOME", "xdg"),
    ):
        path = cache_root / child
        path.mkdir(parents=True, exist_ok=True)
        os.environ[env_name] = str(path)


def run(args: argparse.Namespace, timing: TimingRecorder) -> int:
    mapped_table = args.mapped_table.resolve()
    phylogeny = args.phylogeny.resolve()
    output_dir = args.output_dir.resolve()
    metadata = args.metadata.resolve() if args.metadata else None
    for path in (mapped_table, phylogeny):
        if not path.is_file():
            raise FileNotFoundError(f"Required input not found: {path}")
    if metadata is not None and not metadata.is_file():
        raise FileNotFoundError(f"Metadata not found: {metadata}")
    if args.reference_depth < 1:
        raise ValueError("--reference-depth must be at least one")
    if args.iterations < 1:
        raise ValueError("--iterations must be at least one")

    output_dir.mkdir(parents=True, exist_ok=True)
    qiime = resolve_executable("qiime")
    configure_qiime_environment(qiime, output_dir)

    output_qzv = output_dir / "faith_pd_alpha_rarefaction.qzv"
    if output_qzv.exists():
        if not args.overwrite:
            raise FileExistsError(
                f"Output already exists: {output_qzv}; pass --overwrite to replace it"
            )
        output_qzv.unlink()

    table_export = output_dir / "_mapped_table_export"
    if table_export.exists():
        shutil.rmtree(table_export)
    run_command(
        [
            qiime,
            "tools",
            "export",
            "--input-path",
            str(mapped_table),
            "--output-path",
            str(table_export),
        ],
        timing=timing,
        step="export_mapped_table",
        item=str(mapped_table),
    )
    biom_path = table_export / "feature-table.biom"
    table = biom.load_table(str(biom_path))
    table_sample_ids = [str(sample_id) for sample_id in table.ids(axis="sample")]
    sample_depths = [
        (str(sample_id), int(table.data(sample_id, axis="sample").sum()))
        for sample_id in table_sample_ids
    ]
    max_depth, max_depth_policy = choose_max_depth(
        [depth for _, depth in sample_depths],
        requested=args.max_depth,
        percentile=args.max_depth_percentile,
        reference_depth=args.reference_depth,
    )
    candidate_depths = evenly_spaced_depths(args.min_depth, max_depth, args.steps)
    write_depth_outputs(
        sample_depths, candidate_depths, args.reference_depth, output_dir
    )

    command = [
        qiime,
        "diversity",
        "alpha-rarefaction",
        "--i-table",
        str(mapped_table),
        "--i-phylogeny",
        str(phylogeny),
        "--p-metrics",
        "faith_pd",
        "--p-min-depth",
        str(args.min_depth),
        "--p-max-depth",
        str(max_depth),
        "--p-steps",
        str(args.steps),
        "--p-iterations",
        str(args.iterations),
        "--o-visualization",
        str(output_qzv),
    ]
    qiime_metadata = None
    if metadata is not None:
        qiime_metadata = prepare_qiime_metadata(
            metadata, table_sample_ids, output_dir / "qiime_metadata.tsv"
        )
        command.extend(("--m-metadata-file", str(qiime_metadata)))
    run_command(
        command,
        timing=timing,
        step="faith_pd_alpha_rarefaction",
        item=str(output_qzv),
    )

    positive_depths = [depth for _, depth in sample_depths if depth > 0]
    visualization_export = output_dir / "faith_pd_alpha_rarefaction_export"
    if visualization_export.exists():
        shutil.rmtree(visualization_export)
    run_command(
        [
            qiime,
            "tools",
            "export",
            "--input-path",
            str(output_qzv),
            "--output-path",
            str(visualization_export),
        ],
        timing=timing,
        step="export_alpha_rarefaction_visualization",
        item=str(output_qzv),
    )
    with timing.step("summarize_faith_pd"):
        curve = summarize_faith_pd(
            find_faith_pd_csv(visualization_export), len(positive_depths)
        )
        curve_summary = output_dir / "faith_pd_curve_summary.tsv"
        curve.to_csv(curve_summary, sep="\t", index=False)
    curve_plot = output_dir / "faith_pd_alpha_rarefaction.png"
    with timing.step("plot_faith_pd_alpha_rarefaction", item=str(curve_plot)):
        plot_curve_summary(curve, curve_plot, args.reference_depth)

    summary = {
        "metric": "faith_pd",
        "mapped_table": str(mapped_table),
        "phylogeny": str(phylogeny),
        "metadata": str(metadata) if metadata is not None else None,
        "qiime_metadata": str(qiime_metadata) if qiime_metadata is not None else None,
        "sample_count": len(sample_depths),
        "positive_sample_count": len(positive_depths),
        "reference_depth": args.reference_depth,
        "samples_retained_at_reference_depth": sum(
            depth >= args.reference_depth for depth in positive_depths
        ),
        "min_depth": args.min_depth,
        "max_depth": max_depth,
        "max_depth_policy": max_depth_policy,
        "steps": args.steps,
        "iterations": args.iterations,
        "evaluated_depths": candidate_depths,
        "output_qzv": str(output_qzv),
        "output_png": str(curve_plot),
        "curve_summary_tsv": str(curve_summary),
    }
    with (output_dir / "alpha_rarefaction_summary.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(
        f"Faith's PD alpha rarefaction complete: {output_qzv}\n"
        f"Static plot: {curve_plot}\n"
        f"Mapped samples: {len(sample_depths):,}; max depth: {max_depth:,} "
        f"({max_depth_policy})",
        flush=True,
    )
    return 0


def main() -> int:
    args = parse_args()
    if args.command == "compare":
        return compare(args)
    if args.command == "group-plot":
        return group_plot(args)
    timing = TimingRecorder(args.timings_tsv, component="alpha_rarefaction")
    return run_timed_main(timing, lambda: run(args, timing))


if __name__ == "__main__":
    raise SystemExit(main())
