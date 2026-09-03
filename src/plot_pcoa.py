#!/usr/bin/env python3
"""Plot PCoA coordinates using adaptive, reproducible styling."""
from __future__ import annotations

import argparse
import colorsys
import csv
import hashlib
import json
import math
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from pipeline_lib import TimingRecorder, add_timing_argument, run_timed_main

UNKNOWN_LABEL = "Unknown"
UNKNOWN_COLOR = "#C9CDD3"


@dataclass(frozen=True)
class PlotStyle:
    point_size: float
    alpha: float
    figure_width: float
    figure_height: float
    legend: str
    legend_columns: int
    legend_font_size: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pcoa", type=Path, required=True,
                        help="scikit-bio ordination text or compact PC coordinate TSV.")
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--color-by", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--pc", type=int, nargs=2, default=[1, 2], metavar=("X", "Y"))
    parser.add_argument("--title")
    parser.add_argument("--legend-title")
    parser.add_argument("--legend", choices=("auto", "show", "hide", "separate"),
                        default="auto")
    parser.add_argument("--legend-columns", type=int)
    parser.add_argument("--point-size", type=float)
    parser.add_argument("--alpha", type=float)
    parser.add_argument("--fig-width", type=float)
    parser.add_argument("--fig-height", type=float)
    parser.add_argument("--dpi", type=int, default=220)
    parser.add_argument("--variance", type=Path,
                        help="Axis/proportion TSV used with compact coordinates.")
    parser.add_argument("--palette-file", type=Path,
                        help="CSV/TSV with label and color columns.")
    parser.add_argument("--color-key", type=Path)
    parser.add_argument("--report", type=Path)
    add_timing_argument(parser)
    return parser.parse_args()


def strip_read_suffix(sample_id: str) -> str:
    return re.sub(r"(_R?[12](_001)?)$", "", sample_id)


def display_name(value: str) -> str:
    return value.replace("_", " ").strip().capitalize()


def load_id_to_label(metadata_path: Path, color_by: str) -> dict[str, str]:
    """Return sample/run IDs mapped to labels without imposing label semantics."""
    sep = "\t" if metadata_path.suffix.lower() in {".tsv", ".txt"} else ","
    result: dict[str, str] = {}
    with metadata_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=sep)
        fields = reader.fieldnames or []
        if color_by not in fields:
            raise ValueError(
                f"Column {color_by!r} not in metadata. Available: {', '.join(fields)}"
            )
        if "run_accessions" in fields:
            id_field = "run_accessions"
        elif "sample-id" in fields:
            id_field = "sample-id"
        else:
            raise ValueError("Metadata requires a 'sample-id' or 'run_accessions' column")
        for row in reader:
            label = row.get(color_by, "").strip() or UNKNOWN_LABEL
            values = row[id_field].split(";") if id_field == "run_accessions" else [row[id_field]]
            for value in values:
                sample_id = value.strip()
                if sample_id:
                    result[sample_id] = label
    return result


def _default_variance_path(path: Path) -> Path:
    if path.name == "pcoa_coordinates_pc1_pc3.tsv":
        return path.with_name("pcoa_variance_pc1_pc3.tsv")
    return path.with_name(f"{path.stem}_variance.tsv")


def load_coordinates(
    path: Path, pcs: tuple[int, int], variance_path: Path | None = None
) -> tuple[pd.DataFrame, tuple[float | None, float | None]]:
    """Load two axes from compact TSV or scikit-bio ordination text."""
    px, py = pcs
    with path.open(encoding="utf-8") as handle:
        fields = handle.readline().rstrip("\n").split("\t")
    if "sample-id" in fields and f"PC{px}" in fields and f"PC{py}" in fields:
        frame = pd.read_csv(path, sep="\t", usecols=["sample-id", f"PC{px}", f"PC{py}"])
        frame.columns = ["sample-id", "x", "y"]
        resolved_variance = variance_path or _default_variance_path(path)
        proportions: tuple[float | None, float | None] = (None, None)
        if resolved_variance.exists():
            variance = pd.read_csv(resolved_variance, sep="\t").set_index("axis")
            proportions = (
                float(variance.loc[f"PC{px}", "proportion_explained"]),
                float(variance.loc[f"PC{py}", "proportion_explained"]),
            )
        return frame, proportions

    try:
        import skbio
    except ImportError as error:
        raise RuntimeError(
            "scikit-bio is required for full ordination files; compact coordinate "
            "TSVs can be plotted without it"
        ) from error
    pcoa = skbio.OrdinationResults.read(str(path))
    indices = [px - 1, py - 1]
    frame = pcoa.samples.iloc[:, indices].copy()
    frame.columns = ["x", "y"]
    frame.insert(0, "sample-id", frame.index.astype(str))
    total = float(pcoa.eigvals.sum())
    proportions = (
        float(pcoa.eigvals.iloc[indices[0]] / total),
        float(pcoa.eigvals.iloc[indices[1]] / total),
    )
    return frame.reset_index(drop=True), proportions


def resolve_style(sample_count: int, category_count: int, args: argparse.Namespace) -> PlotStyle:
    size = max(6.0, min(38.0, 34.0 * math.sqrt(500 / max(sample_count, 500))))
    alpha = max(0.48, min(0.78, 0.78 - 0.08 * math.log10(max(sample_count, 500) / 500)))
    if category_count <= 12:
        width, height, columns, font_size = 12.0, 8.0, 1, 9.5
    elif category_count <= 30:
        width, height, columns, font_size = 14.0, 9.0, 1, 8.5
    elif category_count <= 80:
        width, height, columns, font_size = 18.0, 10.0, 2, 7.2
    else:
        width, height, columns, font_size = 14.0, 9.0, 3, 6.8
    legend = args.legend
    if legend == "auto":
        legend = "show" if category_count <= 80 else "separate"
    return PlotStyle(
        args.point_size if args.point_size is not None else size,
        args.alpha if args.alpha is not None else alpha,
        args.fig_width if args.fig_width is not None else width,
        args.fig_height if args.fig_height is not None else height,
        legend,
        args.legend_columns or columns,
        font_size,
    )


def deterministic_color(label: str) -> str:
    if label == UNKNOWN_LABEL:
        return UNKNOWN_COLOR
    digest = hashlib.sha256(label.encode("utf-8")).digest()
    hue = int.from_bytes(digest[:2], "big") / 65535
    saturation = 0.52 + digest[2] / 255 * 0.20
    lightness = 0.43 + digest[3] / 255 * 0.17
    red, green, blue = colorsys.hls_to_rgb(hue, lightness, saturation)
    return f"#{round(red * 255):02X}{round(green * 255):02X}{round(blue * 255):02X}"


def load_palette(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    frame = pd.read_csv(path, sep=sep)
    if not {"label", "color"}.issubset(frame.columns):
        raise ValueError("Palette file requires 'label' and 'color' columns")
    return dict(zip(frame["label"].astype(str), frame["color"].astype(str)))


def axis_label(pc: int, proportion: float | None) -> str:
    return f"PC{pc}" if proportion is None else f"PC{pc} ({proportion * 100:.1f}%)"


def legend_handles(order: list[str], colors: dict[str, str], counts: dict[str, int]):
    return [plt.Line2D(
        [0], [0], marker="o", linestyle="", markersize=7,
        markerfacecolor=colors[label], markeredgewidth=0,
        label=f"{label} (n={counts[label]:,})",
    ) for label in order]


def write_color_key(path: Path, order: list[str], colors: dict[str, str], counts: dict[str, int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("label", "color", "samples"))
        writer.writerows((label, colors[label], counts[label]) for label in order)


def save_separate_legend(
    path: Path, order: list[str], colors: dict[str, str], counts: dict[str, int],
    title: str, columns: int,
) -> None:
    rows = math.ceil(len(order) / columns)
    fig = plt.figure(figsize=(max(8, columns * 4.2), max(2.5, rows * 0.27 + 1.1)))
    fig.legend(handles=legend_handles(order, colors, counts), title=title,
               loc="center", frameon=False, ncol=columns, fontsize=8, title_fontsize=11)
    fig.savefig(path, dpi=220, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def run(args: argparse.Namespace, timing: TimingRecorder) -> int:
    if args.pc[0] < 1 or args.pc[1] < 1 or args.pc[0] == args.pc[1]:
        raise ValueError("--pc values must be distinct positive integers")
    with timing.step("load_pcoa", item=str(args.pcoa)):
        coords, proportions = load_coordinates(args.pcoa, tuple(args.pc), args.variance)
    with timing.step("load_metadata", item=str(args.metadata)):
        id_to_label = load_id_to_label(args.metadata, args.color_by)
        palette = load_palette(args.palette_file)

    coords["sample-id"] = coords["sample-id"].astype(str).map(strip_read_suffix)
    coords["label"] = coords["sample-id"].map(id_to_label).fillna(UNKNOWN_LABEL)
    counts = coords["label"].value_counts().to_dict()
    labels = sorted(counts, key=lambda label: (-counts[label], label.casefold()))
    if UNKNOWN_LABEL in labels:
        labels.remove(UNKNOWN_LABEL)
        labels.insert(0, UNKNOWN_LABEL)
    colors = {label: palette.get(label, deterministic_color(label)) for label in labels}
    style = resolve_style(len(coords), len(labels), args)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    color_key = args.color_key or args.out.with_name(f"{args.out.stem}_colors.tsv")
    report_path = args.report or args.out.with_name(f"{args.out.stem}_report.json")
    write_color_key(color_key, labels, colors, counts)

    with timing.step("render_pcoa_plot", item=str(args.out)):
        fig, ax = plt.subplots(figsize=(style.figure_width, style.figure_height))
        for label in labels:
            group = coords[coords["label"] == label]
            is_unknown = label == UNKNOWN_LABEL
            ax.scatter(group["x"], group["y"], color=colors[label],
                       s=style.point_size * (0.8 if is_unknown else 1.0),
                       alpha=style.alpha * (0.6 if is_unknown else 1.0),
                       linewidths=0, rasterized=True)
        ax.set_xlabel(axis_label(args.pc[0], proportions[0]), fontsize=12)
        ax.set_ylabel(axis_label(args.pc[1], proportions[1]), fontsize=12)
        ax.set_title(args.title or display_name(args.color_by), loc="left",
                     fontsize=17, fontweight="bold", pad=18)
        ax.grid(True, color="#E2E6EC", linewidth=0.65, alpha=0.75)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
        legend_title = args.legend_title or display_name(args.color_by)
        if style.legend == "show":
            ax.legend(handles=legend_handles(labels, colors, counts), title=legend_title,
                      bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False,
                      ncol=style.legend_columns, fontsize=style.legend_font_size,
                      title_fontsize=10.5, columnspacing=1.0, handletextpad=0.45)
            right = 0.70 if style.legend_columns == 1 else 0.60
        else:
            right = 0.96
        fig.subplots_adjust(left=0.08, right=right, bottom=0.10, top=0.91)
        fig.savefig(args.out, dpi=args.dpi, facecolor="white", bbox_inches="tight")
        plt.close(fig)
        legend_path: Path | None = None
        if style.legend == "separate":
            legend_path = args.out.with_name(f"{args.out.stem}_legend.png")
            save_separate_legend(legend_path, labels, colors, counts, legend_title,
                                 style.legend_columns)

    report = {
        "schema_version": 1,
        "pcoa": str(args.pcoa), "metadata": str(args.metadata),
        "color_by": args.color_by, "samples_plotted": len(coords),
        "metadata_categories": len(labels),
        "unmatched_samples": counts.get(UNKNOWN_LABEL, 0),
        "pc": args.pc, "proportion_explained": list(proportions),
        "style": asdict(style), "output": str(args.out),
        "color_key": str(color_key),
        "legend_output": str(legend_path) if legend_path else None,
    }
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Saved {args.out}", file=sys.stderr)
    print(f"Samples={len(coords):,}; categories={len(labels)}; "
          f"unmatched={counts.get(UNKNOWN_LABEL, 0):,}; "
          f"point_size={style.point_size:.1f}; alpha={style.alpha:.2f}; "
          f"legend={style.legend}", file=sys.stderr)
    return 0


def main() -> int:
    args = parse_args()
    timing = TimingRecorder(args.timings_tsv, component="plot_pcoa")
    return run_timed_main(timing, lambda: run(args, timing))


if __name__ == "__main__":
    raise SystemExit(main())
