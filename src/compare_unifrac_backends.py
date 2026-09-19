#!/usr/bin/env python3
"""Compare exact QIIME UniFrac/PCoA outputs with DartUniFrac outputs."""
from __future__ import annotations

import argparse
import csv
import json
import random
import subprocess
from contextlib import contextmanager
from itertools import chain
from pathlib import Path
from typing import Iterator, TextIO

import matplotlib.pyplot as plt
import numpy as np
import skbio
from scipy.spatial import procrustes
from scipy.stats import pearsonr, spearmanr


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qiime-distance", type=Path, required=True)
    parser.add_argument("--dart-distance", type=Path, required=True)
    parser.add_argument("--qiime-ordination", type=Path, required=True)
    parser.add_argument("--dart-ordination", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--dimensions", type=int, default=10)
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--max-distance-samples",
        type=int,
        default=2500,
        help="Deterministic sample limit for distance comparisons (default: 2500).",
    )
    return parser.parse_args()


@contextmanager
def open_distance_text(path: Path) -> Iterator[TextIO]:
    if path.suffix != ".zst":
        with path.open(encoding="utf-8", newline="") as handle:
            yield handle
        return
    process = subprocess.Popen(
        ["zstd", "-dc", str(path)],
        stdout=subprocess.PIPE,
        text=True,
        encoding="utf-8",
    )
    assert process.stdout is not None
    try:
        yield process.stdout
    finally:
        process.stdout.close()
        exit_code = process.wait()
        if exit_code != 0:
            raise RuntimeError(f"zstd failed with exit code {exit_code}: {path}")


def normalized_distance_header(header: list[str], first_row: list[str]) -> list[str]:
    if len(header) == len(first_row):
        return header[1:]
    if len(header) == len(first_row) - 1:
        return header
    raise ValueError(
        "Distance matrix header/row widths are inconsistent: "
        f"header={len(header)}, row={len(first_row)}"
    )


def distance_ids(path: Path) -> list[str]:
    with open_distance_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        first_row = next(reader)
    return normalized_distance_header(header, first_row)


def read_distance_subset(path: Path, selected_ids: list[str]) -> np.ndarray:
    wanted = set(selected_ids)
    with open_distance_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        first_row = next(reader)
        column_ids = normalized_distance_header(header, first_row)
        column_index = {sample_id: index for index, sample_id in enumerate(column_ids)}
        missing = wanted - column_index.keys()
        if missing:
            raise ValueError(f"Distance matrix lacks selected IDs: {sorted(missing)[:5]}")
        selected_columns = [column_index[sample_id] for sample_id in selected_ids]
        rows: dict[str, np.ndarray] = {}
        for fields in chain((first_row,), reader):
            if not fields:
                continue
            sample_id = fields[0]
            if sample_id in wanted:
                values = fields[1:]
                rows[sample_id] = np.asarray(
                    [float(values[index]) for index in selected_columns], dtype=float
                )
    missing_rows = wanted - rows.keys()
    if missing_rows:
        raise ValueError(f"Distance matrix lacks selected rows: {sorted(missing_rows)[:5]}")
    return np.vstack([rows[sample_id] for sample_id in selected_ids])


def deterministic_common_ids(
    first: Path, second: Path, maximum: int, seed: int
) -> tuple[list[str], dict[str, int]]:
    first_ids = distance_ids(first)
    second_ids = distance_ids(second)
    if set(first_ids) != set(second_ids):
        raise ValueError(
            "Distance matrix sample IDs differ: "
            f"QIIME-only={len(set(first_ids) - set(second_ids))}, "
            f"DART-only={len(set(second_ids) - set(first_ids))}"
        )
    selected = sorted(first_ids)
    if len(selected) > maximum:
        selected = sorted(random.Random(seed).sample(selected, maximum))
    return selected, {"total": len(first_ids), "compared": len(selected)}


def validate_distance(matrix: np.ndarray, label: str) -> dict[str, float]:
    symmetry = float(np.max(np.abs(matrix - matrix.T)))
    diagonal = float(np.max(np.abs(np.diag(matrix))))
    if symmetry > 1e-8 or diagonal > 1e-8:
        raise ValueError(
            f"{label} distance matrix is invalid: "
            f"max symmetry error={symmetry}, max diagonal={diagonal}"
        )
    return {"max_symmetry_error": symmetry, "max_abs_diagonal": diagonal}


def ordination_coordinates(path: Path, ids: list[str], dimensions: int) -> np.ndarray:
    result = skbio.io.read(
        str(path),
        format="ordination",
        into=skbio.stats.ordination.OrdinationResults,
    )
    missing = set(ids) - set(result.samples.index)
    if missing:
        raise ValueError(f"Ordination lacks sample IDs: {sorted(missing)[:5]}")
    usable = min(dimensions, result.samples.shape[1])
    if usable < 2:
        raise ValueError(f"Ordination has fewer than two axes: {path}")
    return result.samples.loc[ids].iloc[:, :usable].to_numpy(dtype=float)


def exact_pcoa(matrix: np.ndarray, ids: list[str], dimensions: int) -> np.ndarray:
    distance = skbio.DistanceMatrix(matrix, ids=ids)
    result = skbio.stats.ordination.pcoa(
        distance,
        method="eigh",
        number_of_dimensions=min(dimensions, len(ids) - 1),
    )
    return result.samples.to_numpy(dtype=float)


def fsvd_pcoa(
    matrix: np.ndarray, ids: list[str], dimensions: int, seed: int
) -> np.ndarray:
    """Run the seeded scikit-bio FSVD reference used by DART's validation."""
    state = np.random.get_state()
    np.random.seed(seed)
    try:
        result = skbio.stats.ordination.pcoa(
            skbio.DistanceMatrix(matrix, ids=ids),
            method="fsvd",
            number_of_dimensions=min(dimensions, len(ids) - 1),
        )
    finally:
        np.random.set_state(state)
    return result.samples.to_numpy(dtype=float)


def procrustes_test(
    reference: np.ndarray,
    candidate: np.ndarray,
    *,
    permutations: int,
    seed: int,
) -> dict[str, float | int]:
    dimensions = min(reference.shape[1], candidate.shape[1])
    reference = reference[:, :dimensions]
    candidate = candidate[:, :dimensions]
    _, transformed, disparity = procrustes(reference, candidate)
    rng = np.random.default_rng(seed)
    exceedances = 0
    for _ in range(permutations):
        _, _, permuted_disparity = procrustes(
            reference, candidate[rng.permutation(candidate.shape[0])]
        )
        if permuted_disparity <= disparity:
            exceedances += 1
    return {
        "dimensions": dimensions,
        "m2": float(disparity),
        "similarity_1_minus_m2": float(1 - disparity),
        "permutations": permutations,
        "p_value": float((exceedances + 1) / (permutations + 1)),
        "transformed_rms": float(np.sqrt(np.mean(transformed**2))),
    }


def run(args: argparse.Namespace) -> int:
    if args.dimensions < 2 or args.permutations < 0 or args.max_distance_samples < 2:
        raise ValueError("dimensions and sample limit must be >=2; permutations >=0")
    ids, counts = deterministic_common_ids(
        args.qiime_distance,
        args.dart_distance,
        args.max_distance_samples,
        args.seed,
    )
    qiime_dm = read_distance_subset(args.qiime_distance, ids)
    dart_dm = read_distance_subset(args.dart_distance, ids)
    validation = {
        "qiime": validate_distance(qiime_dm, "QIIME"),
        "dart": validate_distance(dart_dm, "DART"),
    }
    triangle = np.triu_indices(len(ids), k=1)
    qiime_values = qiime_dm[triangle]
    dart_values = dart_dm[triangle]
    residuals = dart_values - qiime_values
    distance_metrics = {
        "pearson": float(pearsonr(qiime_values, dart_values).statistic),
        "spearman": float(spearmanr(qiime_values, dart_values).statistic),
        "rmse": float(np.sqrt(np.mean(residuals**2))),
        "mae": float(np.mean(np.abs(residuals))),
        "maximum_absolute_error": float(np.max(np.abs(residuals))),
    }

    exact_qiime = exact_pcoa(qiime_dm, ids, args.dimensions)
    exact_dart = exact_pcoa(dart_dm, ids, args.dimensions)
    fsvd_dart = fsvd_pcoa(dart_dm, ids, args.dimensions, args.seed)
    qiime_coords = ordination_coordinates(args.qiime_ordination, ids, args.dimensions)
    dart_coords = ordination_coordinates(args.dart_ordination, ids, args.dimensions)
    full_matrix_compared = counts["total"] == counts["compared"]
    procrustes_metrics = {
        "distance_only": procrustes_test(
            exact_qiime,
            exact_dart,
            permutations=args.permutations,
            seed=args.seed,
        ),
        "dart_fpcoa_only": (
            procrustes_test(
                fsvd_dart,
                dart_coords,
                permutations=args.permutations,
                seed=args.seed + 1,
            )
            if full_matrix_compared
            else None
        ),
        "end_to_end": procrustes_test(
            qiime_coords,
            dart_coords,
            permutations=args.permutations,
            seed=args.seed + 2,
        ),
    }
    thresholds = {
        "pearson_minimum": 0.98,
        "spearman_minimum": 0.98,
        "rmse_maximum": 0.03,
        "distance_procrustes_m2_maximum": 0.02,
        "fpcoa_procrustes_m2_maximum": 0.001,
    }
    gates = {
        "pearson": distance_metrics["pearson"] >= thresholds["pearson_minimum"],
        "spearman": distance_metrics["spearman"] >= thresholds["spearman_minimum"],
        "rmse": distance_metrics["rmse"] <= thresholds["rmse_maximum"],
        "distance_procrustes": (
            procrustes_metrics["distance_only"]["m2"]
            <= thresholds["distance_procrustes_m2_maximum"]
        ),
        "fpcoa_procrustes": (
            procrustes_metrics["dart_fpcoa_only"]["m2"]
            <= thresholds["fpcoa_procrustes_m2_maximum"]
            if procrustes_metrics["dart_fpcoa_only"] is not None
            else None
        ),
    }
    summary = {
        "sample_counts": counts,
        "validation": validation,
        "distance_metrics": distance_metrics,
        "procrustes": procrustes_metrics,
        "thresholds": thresholds,
        "gates": gates,
        "passed": all(value for value in gates.values() if value is not None),
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    with (args.out_dir / "comparison_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    with (args.out_dir / "distance_residuals.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("qiime_distance", "dart_distance", "residual"))
        writer.writerows(zip(qiime_values, dart_values, residuals, strict=True))

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].hexbin(qiime_values, dart_values, gridsize=60, mincnt=1)
    limits = [
        min(float(qiime_values.min()), float(dart_values.min())),
        max(float(qiime_values.max()), float(dart_values.max())),
    ]
    axes[0].plot(limits, limits, color="black", linewidth=1)
    axes[0].set(xlabel="QIIME exact", ylabel="DART DMH", title="Pairwise distances")
    axes[1].hist(residuals, bins=50)
    axes[1].set(xlabel="DART - QIIME", ylabel="Pairs", title="Distance residuals")
    fig.tight_layout()
    fig.savefig(args.out_dir / "distance_comparison.png", dpi=180)
    plt.close(fig)
    return 0 if summary["passed"] else 1


def main() -> int:
    return run(parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
