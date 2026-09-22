#!/usr/bin/env python3
"""Compare four PCoA ordinations on pairwise and four-way shared samples."""

from __future__ import annotations

import argparse
import csv
import json
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import procrustes

from compare_unifrac_backends import ordination_coordinates


LEVELS = ("full", "sub50", "sub25", "sub10")


def load_coordinates(path: Path) -> dict[str, np.ndarray]:
    ids, coordinates = ordination_coordinates(path, 2)
    if len(ids) != len(set(ids)):
        raise ValueError(f"Duplicate sample IDs in {path}")
    if coordinates.shape[0] < 3 or not np.isfinite(coordinates).all():
        raise ValueError(f"Too few samples or non-finite coordinates in {path}")
    return dict(zip(ids, coordinates, strict=True))


def compare_pair(
    reference: dict[str, np.ndarray],
    candidate: dict[str, np.ndarray],
    ids: set[str],
) -> tuple[float, np.ndarray, np.ndarray]:
    ordered = sorted(ids)
    if len(ordered) < 3:
        raise ValueError("Procrustes requires at least three shared samples")
    first = np.stack([reference[sample] for sample in ordered])
    second = np.stack([candidate[sample] for sample in ordered])
    if np.linalg.norm(first - first.mean(axis=0)) == 0:
        raise ValueError("Reference coordinates have no variation")
    if np.linalg.norm(second - second.mean(axis=0)) == 0:
        raise ValueError("Candidate coordinates have no variation")
    aligned_first, aligned_second, m2 = procrustes(first, second)
    return float(m2), aligned_first, aligned_second


def make_figure(
    panels: list[tuple[str, str, int, float, np.ndarray, np.ndarray]],
    output: Path,
    title: str,
) -> None:
    figure, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    for axis, (first, second, count, m2, reference, candidate) in zip(
        axes.flat, panels, strict=True
    ):
        axis.scatter(reference[:, 0], reference[:, 1], s=0.4, alpha=0.25,
                     c="#536477", rasterized=True, label=first)
        axis.scatter(candidate[:, 0], candidate[:, 1], s=0.4, alpha=0.25,
                     c="#d5733a", rasterized=True, label=second)
        axis.set_title(f"{first} vs {second}   n={count:,}   M²={m2:.4f}")
        axis.set_aspect("equal", adjustable="datalim")
        axis.set_xlabel("Aligned axis 1")
        axis.set_ylabel("Aligned axis 2")
    # Individual panel titles identify which color corresponds to which level.
    figure.suptitle(f"{title}\nGray = first level; orange = second level", fontsize=15)
    figure.savefig(output, dpi=180)
    plt.close(figure)


def run(paths: dict[str, Path], out_dir: Path) -> dict[str, object]:
    if set(paths) != set(LEVELS):
        raise ValueError(f"Expected paths for {', '.join(LEVELS)}")
    data = {level: load_coordinates(paths[level]) for level in LEVELS}
    shared_all = set.intersection(*(set(samples) for samples in data.values()))
    if len(shared_all) < 3:
        raise ValueError("Fewer than three samples occur in all four ordinations")

    out_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    panels: dict[str, list[tuple[str, str, int, float, np.ndarray, np.ndarray]]] = {
        "pairwise": [], "four_way": []
    }
    for first, second in combinations(LEVELS, 2):
        pair_ids = set(data[first]) & set(data[second])
        for cohort, ids in (("pairwise", pair_ids), ("four_way", shared_all)):
            m2, reference, candidate = compare_pair(data[first], data[second], ids)
            rows.append({"level_a": first, "level_b": second, "cohort": cohort,
                         "sample_count": len(ids), "dimensions": 2, "m2": m2})
            panels[cohort].append((first, second, len(ids), m2, reference, candidate))

    with (out_dir / "procrustes_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "ordination_files": {level: str(paths[level].resolve()) for level in LEVELS},
        "samples_per_level": {level: len(data[level]) for level in LEVELS},
        "four_way_shared_samples": len(shared_all),
        "comparisons": rows,
    }
    (out_dir / "procrustes_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    make_figure(panels["pairwise"], out_dir / "procrustes_pairwise.png",
                "Pairwise shared samples — two-axis Procrustes")
    make_figure(panels["four_way"], out_dir / "procrustes_four_way.png",
                "Same four-way shared cohort — two-axis Procrustes")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for level in LEVELS:
        parser.add_argument(f"--{level}", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    summary = run({level: getattr(args, level) for level in LEVELS}, args.out_dir)
    print(f"Compared six pairs; four-way overlap: {summary['four_way_shared_samples']:,} samples")
    print(f"Results: {args.out_dir}")


if __name__ == "__main__":
    main()
