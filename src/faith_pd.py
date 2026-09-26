#!/usr/bin/env python3
"""Faith's PD alpha diversity at the pipeline rarefaction depth.

``run`` computes Faith's PD on the same rarefied, GG2-mapped table that
``src/unifrac.py`` used, so alpha and beta diversity share one random draw.
It writes one value per sample and, optionally, a box plot plus a
Kruskal-Wallis test for each metadata grouping column.

``compare`` takes ``faith_pd.tsv`` files from several runs (for example
full / 50% / 25% / 10% read subsampling) and shows how Faith's PD shifts
between them on a common set of samples.

Examples:
  python src/faith_pd.py run \\
    --rarefied-table runs/x/work/qiime2/rarefied-backbone-mapped-table.qza \\
    --phylogeny data/gg2/2024.09.phylogeny.id.nwk.qza \\
    --output-dir runs/x/results --metadata meta.tsv --group-by body_site

  python src/faith_pd.py compare --run full=a/faith_pd.tsv --run sub10=b/faith_pd.tsv \\
    --metadata meta.tsv --group-by body_site --output-dir comparison/
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
from pathlib import Path

from alpha_rarefaction import configure_qiime_environment

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from pipeline_lib import (  # noqa: E402
    TimingRecorder,
    add_timing_argument,
    resolve_executable,
    run_command,
    run_timed_main,
)
from plot_pcoa import UNKNOWN_LABEL, load_id_to_label, strip_read_suffix  # noqa: E402

BOX_COLOR = "#3987e5"
INK = "#1f2933"
MUTED = "#6b7280"
GRID = "#e5e7eb"
# Ordinal blue ramp (validated light->dark); the most reads gets the darkest step.
LEVEL_RAMP = ("#0d366b", "#1c5cab", "#3987e5", "#86b6ef")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = parser.add_subparsers(dest="command", required=True)

    run_parser = sub.add_parser("run", help="Compute Faith's PD for one run.")
    run_parser.add_argument("--rarefied-table", type=Path, required=True,
                            help="Rarefied GG2-mapped FeatureTable[Frequency] QZA.")
    run_parser.add_argument("--phylogeny", type=Path, required=True,
                            help="Rooted GG2 phylogeny QZA.")
    run_parser.add_argument("--output-dir", type=Path, required=True,
                            help="Directory for faith_pd.tsv and box plots.")
    run_parser.add_argument("--work-dir", type=Path,
                            help="Intermediate QIIME files (default: <output-dir>/_faith_pd_work).")
    run_parser.add_argument("--metadata", type=Path)
    run_parser.add_argument("--group-by", action="append", default=[],
                            help="Metadata column for a box plot; repeatable. Needs --metadata.")
    add_timing_argument(run_parser)

    compare_parser = sub.add_parser("compare", help="Compare faith_pd.tsv across runs.")
    compare_parser.add_argument("--run", action="append", required=True,
                                metavar="LABEL=FAITH_PD_TSV",
                                help="Repeat per run, reference (e.g. full) first.")
    compare_parser.add_argument("--metadata", type=Path, required=True)
    compare_parser.add_argument("--group-by", required=True)
    compare_parser.add_argument("--cohort-ids", type=Path,
                                help="Sample IDs to use; default: samples present in every run.")
    compare_parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def safe_name(value: str) -> str:
    # Same rule as run_pipeline.safe_output_name so the pipeline can predict filenames.
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._")
    if not safe:
        raise ValueError(f"Metadata column cannot form a safe output filename: {value!r}")
    return safe


def read_faith_pd(path: Path) -> pd.Series:
    frame = pd.read_csv(path, sep="\t", dtype={"sample-id": str})
    return frame.set_index("sample-id")["faith_pd"].astype(float)


def labelled(values: pd.Series, metadata: Path, column: str) -> tuple[pd.DataFrame, int]:
    """Attach metadata labels; drop samples without a label."""
    labels = load_id_to_label(metadata, column)
    frame = pd.DataFrame({"faith_pd": values})
    frame["group"] = [labels.get(strip_read_suffix(s), UNKNOWN_LABEL) for s in frame.index]
    unlabelled = int((frame["group"] == UNKNOWN_LABEL).sum())
    return frame[frame["group"] != UNKNOWN_LABEL], unlabelled


def kruskal(frame: pd.DataFrame) -> dict[str, float | int]:
    samples = [group["faith_pd"].to_numpy() for _, group in frame.groupby("group")]
    h, p = stats.kruskal(*samples)
    n = len(frame)
    # Epsilon-squared effect size (Tomczak & Tomczak 2014): H / (n - 1).
    return {"n": n, "groups": len(samples), "kruskal_h": float(h), "p_value": float(p),
            "epsilon2": float(h / (n - 1))}


def style_axes(ax: plt.Axes) -> None:
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=10)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def plot_groups(frame: pd.DataFrame, output: Path, column: str, test: dict) -> None:
    order = frame.groupby("group")["faith_pd"].median().sort_values(ascending=False).index
    data = [frame.loc[frame["group"] == g, "faith_pd"].to_numpy() for g in order]
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.boxplot(
        data, widths=0.6, showfliers=False, patch_artist=True,
        boxprops={"facecolor": BOX_COLOR, "edgecolor": BOX_COLOR, "alpha": 0.85},
        medianprops={"color": INK, "linewidth": 2},
        whiskerprops={"color": MUTED}, capprops={"color": MUTED},
    )
    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels([f"{g}\nn={len(d):,}" for g, d in zip(order, data)])
    ax.set_ylabel("Faith's PD", color=INK)
    ax.set_title(
        f"Faith's PD by {column}   (Kruskal-Wallis H = {test['kruskal_h']:,.0f}, "
        f"ε² = {test['epsilon2']:.3f})",
        loc="left", color=INK, fontsize=12,
    )
    style_axes(ax)
    fig.text(0.01, 0.01, "Boxes: median and IQR; whiskers: 1.5 × IQR; outliers not drawn.",
             color=MUTED, fontsize=8)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(output, dpi=200)
    plt.close(fig)


def export_faith_pd(qiime: str, table: Path, phylogeny: Path, work_dir: Path,
                    timing: TimingRecorder) -> pd.Series:
    vector_qza = work_dir / "faith_pd_vector.qza"
    vector_qza.unlink(missing_ok=True)
    run_command(
        [qiime, "diversity", "alpha-phylogenetic",
         "--i-table", str(table), "--i-phylogeny", str(phylogeny),
         "--p-metric", "faith_pd", "--o-alpha-diversity", str(vector_qza)],
        timing=timing, step="faith_pd", item=str(vector_qza),
    )
    export_dir = work_dir / "faith_pd_export"
    shutil.rmtree(export_dir, ignore_errors=True)
    run_command(
        [qiime, "tools", "export", "--input-path", str(vector_qza),
         "--output-path", str(export_dir)],
        timing=timing, step="export_faith_pd", item=str(vector_qza),
    )
    exported = pd.read_csv(export_dir / "alpha-diversity.tsv", sep="\t", index_col=0)
    values = exported.iloc[:, 0].astype(float)
    values.index = values.index.astype(str)
    values.index.name = "sample-id"
    return values.rename("faith_pd")


def run(args: argparse.Namespace, timing: TimingRecorder) -> int:
    if args.group_by and args.metadata is None:
        raise ValueError("--group-by requires --metadata")
    table = args.rarefied_table.resolve()
    phylogeny = args.phylogeny.resolve()
    for path in [table, phylogeny] + ([args.metadata] if args.metadata else []):
        if not path.is_file():
            raise FileNotFoundError(f"Required input not found: {path}")
    output_dir = args.output_dir.resolve()
    work_dir = (args.work_dir or output_dir / "_faith_pd_work").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)

    qiime = resolve_executable("qiime")
    configure_qiime_environment(qiime, work_dir)
    values = export_faith_pd(qiime, table, phylogeny, work_dir, timing)
    values.to_frame().to_csv(output_dir / "faith_pd.tsv", sep="\t")
    print(f"Faith's PD for {len(values):,} samples: {output_dir / 'faith_pd.tsv'}", flush=True)

    summary: dict[str, object] = {
        "metric": "faith_pd", "rarefied_table": str(table), "phylogeny": str(phylogeny),
        "sample_count": int(len(values)), "median": float(values.median()), "groupings": {},
    }
    for column in args.group_by:
        frame, unlabelled = labelled(values, args.metadata, column)
        with timing.step("kruskal_wallis", item=column):
            test = kruskal(frame)
        plot = output_dir / f"faith_pd_{safe_name(column)}.png"
        with timing.step("plot_faith_pd_groups", item=str(plot)):
            plot_groups(frame, plot, column, test)
        medians = frame.groupby("group")["faith_pd"].median().sort_values(ascending=False)
        summary["groupings"][column] = {  # type: ignore[index]
            **test, "unlabelled_samples_excluded": unlabelled,
            "group_medians": {str(k): float(v) for k, v in medians.items()},
            "plot": str(plot),
        }
        print(f"  {column}: H={test['kruskal_h']:.1f} eps2={test['epsilon2']:.3f} -> {plot}",
              flush=True)
    with (output_dir / "faith_pd_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return 0


def parse_runs(specs: list[str]) -> dict[str, Path]:
    runs: dict[str, Path] = {}
    for spec in specs:
        label, sep, path = spec.partition("=")
        if not sep or not label or not path:
            raise ValueError(f"--run must be LABEL=PATH, got {spec!r}")
        if label in runs:
            raise ValueError(f"Duplicate --run label {label!r}")
        runs[label] = Path(path)
    if len(runs) < 2:
        raise ValueError("compare needs at least two --run values")
    return runs


def compare_frames(runs: dict[str, Path], cohort_ids: Path | None) -> tuple[pd.DataFrame, str]:
    """Wide table (samples x runs) restricted to the comparison cohort."""
    wide = pd.concat({label: read_faith_pd(path) for label, path in runs.items()}, axis=1)
    if cohort_ids is not None:
        cohort = [s.strip() for s in cohort_ids.read_text().splitlines() if s.strip()]
        missing = sorted(set(cohort) - set(wide.dropna().index))
        if missing:
            raise ValueError(f"{len(missing)} cohort IDs lack Faith's PD in some run, "
                             f"e.g. {missing[:5]}")
        return wide.loc[cohort], f"fixed cohort from {cohort_ids.name}"
    return wide.dropna(), "samples present in every run"


def plot_levels(long: pd.DataFrame, labels: list[str], output: Path, column: str) -> None:
    reference = labels[0]
    order = (long[long["run"] == reference].groupby("group")["faith_pd"].median()
             .sort_values(ascending=False).index)
    fig, ax = plt.subplots(figsize=(12, 5.5))
    width = 0.8 / len(labels)
    for index, label in enumerate(labels):
        subset = long[long["run"] == label]
        data = [subset.loc[subset["group"] == g, "faith_pd"].to_numpy() for g in order]
        positions = np.arange(len(order)) + (index - (len(labels) - 1) / 2) * width
        color = LEVEL_RAMP[index % len(LEVEL_RAMP)]
        ax.boxplot(
            data, positions=positions, widths=width * 0.85, showfliers=False,
            patch_artist=True,
            boxprops={"facecolor": color, "edgecolor": color},
            medianprops={"color": "white", "linewidth": 1.5},
            whiskerprops={"color": color}, capprops={"color": color},
        )
        ax.plot([], [], color=color, linewidth=8, label=label)
    counts = long[long["run"] == reference].groupby("group").size()
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{g}\nn={counts[g]:,}" for g in order])
    ax.set_ylabel("Faith's PD", color=INK)
    ax.legend(frameon=False, ncol=len(labels), loc="upper right", fontsize=10)
    style_axes(ax)
    fig.text(0.01, 0.01, f"Same samples in every box group; by {column}. Outliers not drawn.",
             color=MUTED, fontsize=8)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(output, dpi=200)
    plt.close(fig)


def plot_agreement(wide: pd.DataFrame, labels: list[str], output: Path) -> None:
    reference = labels[0]
    others = labels[1:]
    fig, axes = plt.subplots(1, len(others), figsize=(4.2 * len(others), 4.2),
                             sharex=True, sharey=True, squeeze=False)
    limit = float(np.nanpercentile(wide[labels].to_numpy(), 99.5))
    for ax, label in zip(axes[0], others):
        ax.hexbin(wide[reference], wide[label], gridsize=60, bins="log", cmap="Blues",
                  extent=(0, limit, 0, limit), mincnt=1)
        ax.plot([0, limit], [0, limit], color=MUTED, linewidth=1, linestyle="--")
        rho = stats.spearmanr(wide[reference], wide[label]).statistic
        ax.set_title(f"{label} vs {reference}   ρ = {rho:.3f}", loc="left", color=INK,
                     fontsize=11)
        ax.set_xlabel(f"Faith's PD, {reference}", color=MUTED)
        style_axes(ax)
    axes[0][0].set_ylabel("Faith's PD, subsampled", color=MUTED)
    fig.tight_layout()
    fig.savefig(output, dpi=200)
    plt.close(fig)


def compare(args: argparse.Namespace) -> int:
    runs = parse_runs(args.run)
    labels = list(runs)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    wide, cohort_note = compare_frames(runs, args.cohort_ids)
    group_labels = load_id_to_label(args.metadata, args.group_by)
    wide["group"] = [group_labels.get(strip_read_suffix(s), UNKNOWN_LABEL) for s in wide.index]
    wide = wide[wide["group"] != UNKNOWN_LABEL]
    long = wide.melt(id_vars="group", value_vars=labels, var_name="run",
                     value_name="faith_pd", ignore_index=False)

    reference = labels[0]
    rows = []
    for label in labels:
        test = kruskal(long[long["run"] == label])
        change = 100 * (wide[label] - wide[reference]) / wide[reference]
        rows.append({
            "run": label, **test,
            "median_faith_pd": float(wide[label].median()),
            "spearman_vs_reference": float(stats.spearmanr(wide[reference], wide[label]).statistic),
            "median_pct_change_vs_reference": float(change.median()),
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(output_dir / "faith_pd_levels_summary.tsv", sep="\t", index=False)
    by_group = (long.groupby(["group", "run"])["faith_pd"].median().unstack("run")[labels])
    by_group.to_csv(output_dir / "faith_pd_levels_group_medians.tsv", sep="\t")

    safe = safe_name(args.group_by)
    plot_levels(long, labels, output_dir / f"faith_pd_levels_{safe}.png", args.group_by)
    plot_agreement(wide, labels, output_dir / "faith_pd_levels_agreement.png")
    with (output_dir / "faith_pd_levels_summary.json").open("w", encoding="utf-8") as handle:
        json.dump({"cohort": cohort_note, "sample_count": int(len(wide)),
                   "reference": reference, "group_by": args.group_by,
                   "runs": {label: str(path) for label, path in runs.items()},
                   "levels": rows,
                   "group_medians": json.loads(by_group.to_json(orient="index"))},
                  handle, indent=2)
        handle.write("\n")
    print(summary.to_string(index=False), flush=True)
    return 0


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.command == "compare":
        return compare(args)
    timing = TimingRecorder(args.timings_tsv, component="faith_pd")
    return run_timed_main(timing, lambda: run(args, timing))


if __name__ == "__main__":
    raise SystemExit(main())
