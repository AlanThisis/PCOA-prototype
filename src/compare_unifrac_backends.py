#!/usr/bin/env python3
"""Compare exact QIIME UniFrac/PCoA outputs with DartUniFrac outputs."""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import procrustes
from scipy.stats import spearmanr


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qiime-distance", type=Path, required=True)
    parser.add_argument("--dart-distance", type=Path, required=True)
    parser.add_argument("--qiime-ordination", type=Path, required=True)
    parser.add_argument(
        "--dart-exact-ordination",
        type=Path,
        required=True,
        help="Exact eigh PCoA computed from the DART distance matrix.",
    )
    parser.add_argument("--dart-ordination", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--dimensions", type=int, default=10)
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--spearman-pairs",
        type=int,
        default=10_000_000,
        help=(
            "Maximum uniformly sampled distance pairs for Spearman correlation. "
            "Pearson and RMSE always use every pair; use 0 for every pair."
        ),
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


def distance_header(handle: TextIO, path: Path) -> list[str]:
    header = handle.readline().rstrip("\n\r").split("\t")
    if header and header[0] == "":
        header = header[1:]
    if len(header) < 2:
        raise ValueError(f"Distance matrix has an invalid header: {path}")
    return header


def sample_pairs(
    sample_count: int, pair_count: int, seed: int
) -> tuple[np.ndarray, np.ndarray]:
    """Uniformly sample off-diagonal unordered pairs, with replacement."""
    rng = np.random.default_rng(seed)
    first_parts: list[np.ndarray] = []
    second_parts: list[np.ndarray] = []
    remaining = pair_count
    while remaining:
        draw = max(remaining + remaining // 20, 1024)
        first = rng.integers(0, sample_count, size=draw, dtype=np.int64)
        second = rng.integers(0, sample_count, size=draw, dtype=np.int64)
        keep = first != second
        first = first[keep][:remaining]
        second = second[keep][:remaining]
        first_parts.append(np.minimum(first, second))
        second_parts.append(np.maximum(first, second))
        remaining -= len(first)
    return np.concatenate(first_parts), np.concatenate(second_parts)


def compare_distance_matrices(
    qiime_path: Path,
    dart_path: Path,
    *,
    spearman_pairs: int,
    seed: int,
    scratch_dir: Path,
) -> tuple[
    list[str], dict[str, float | int | str], np.ndarray, np.ndarray, np.ndarray
]:
    """Compare all distance pairs while bounding only Spearman rank storage."""
    with open_distance_text(qiime_path) as qiime_handle, open_distance_text(
        dart_path
    ) as dart_handle:
        ids = distance_header(qiime_handle, qiime_path)
        dart_ids = distance_header(dart_handle, dart_path)
        if ids != dart_ids:
            if set(ids) == set(dart_ids):
                raise ValueError(
                    "Distance matrix IDs match but their order differs; reorder the "
                    "matrices before full streaming comparison"
                )
            raise ValueError(
                "Distance matrix sample IDs differ: "
                f"QIIME-only={len(set(ids) - set(dart_ids))}, "
                f"DART-only={len(set(dart_ids) - set(ids))}"
            )

        sample_count = len(ids)
        total_pairs = sample_count * (sample_count - 1) // 2
        exact_spearman = spearman_pairs == 0 or spearman_pairs >= total_pairs
        retained_pairs = total_pairs if exact_spearman else spearman_pairs
        if exact_spearman:
            sampled_rows = sampled_columns = None
            sampled_qiime_parts: list[np.ndarray] = []
            sampled_dart_parts: list[np.ndarray] = []
        else:
            sampled_rows, sampled_columns = sample_pairs(
                sample_count, retained_pairs, seed
            )
            order = np.argsort(sampled_rows, kind="stable")
            sampled_rows = sampled_rows[order]
            sampled_columns = sampled_columns[order]
            sampled_qiime = np.empty(retained_pairs, dtype=np.float64)
            sampled_dart = np.empty(retained_pairs, dtype=np.float64)

        count = 0
        sum_qiime = sum_dart = 0.0
        sum_qiime_sq = sum_dart_sq = sum_product = 0.0
        sum_squared_error = sum_absolute_error = 0.0
        maximum_absolute_error = 0.0
        sample_start = 0
        max_abs_diagonal = {"qiime": 0.0, "dart": 0.0}
        max_symmetry_error = {"qiime": 0.0, "dart": 0.0}

        scratch_dir.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix=".distance-symmetry-", dir=scratch_dir
        ) as temporary:
            temporary_path = Path(temporary)
            qiime_upper_store = np.memmap(
                temporary_path / "qiime-upper.bin",
                dtype=np.float64,
                mode="w+",
                shape=(total_pairs,),
            )
            dart_upper_store = np.memmap(
                temporary_path / "dart-upper.bin",
                dtype=np.float64,
                mode="w+",
                shape=(total_pairs,),
            )

            for row_index, (qiime_line, dart_line) in enumerate(
                zip(qiime_handle, dart_handle, strict=True)
            ):
                qiime_id, separator, qiime_text = qiime_line.partition("\t")
                dart_id, dart_separator, dart_text = dart_line.partition("\t")
                if not separator or not dart_separator or qiime_id != ids[row_index]:
                    raise ValueError(f"Invalid QIIME row {row_index} in {qiime_path}")
                if dart_id != qiime_id:
                    raise ValueError(
                        f"Distance matrix row IDs differ at row {row_index}: "
                        f"{qiime_id!r} != {dart_id!r}"
                    )
                qiime_row = np.fromstring(qiime_text, sep="\t", dtype=np.float64)
                dart_row = np.fromstring(dart_text, sep="\t", dtype=np.float64)
                if len(qiime_row) != sample_count or len(dart_row) != sample_count:
                    raise ValueError(f"Invalid distance row width for {qiime_id}")
                max_abs_diagonal["qiime"] = max(
                    max_abs_diagonal["qiime"], abs(float(qiime_row[row_index]))
                )
                max_abs_diagonal["dart"] = max(
                    max_abs_diagonal["dart"], abs(float(dart_row[row_index]))
                )

                if row_index:
                    columns = np.arange(row_index, dtype=np.int64)
                    prior_indices = (
                        columns * sample_count
                        - columns * (columns + 1) // 2
                        + row_index
                        - columns
                        - 1
                    )
                    max_symmetry_error["qiime"] = max(
                        max_symmetry_error["qiime"],
                        float(
                            np.max(
                                np.abs(
                                    qiime_row[:row_index]
                                    - qiime_upper_store[prior_indices]
                                )
                            )
                        ),
                    )
                    max_symmetry_error["dart"] = max(
                        max_symmetry_error["dart"],
                        float(
                            np.max(
                                np.abs(
                                    dart_row[:row_index]
                                    - dart_upper_store[prior_indices]
                                )
                            )
                        ),
                    )

                qiime_upper = qiime_row[row_index + 1 :]
                dart_upper = dart_row[row_index + 1 :]
                pair_count = len(qiime_upper)
                upper_start = (
                    row_index * sample_count - row_index * (row_index + 1) // 2
                )
                upper_end = upper_start + pair_count
                qiime_upper_store[upper_start:upper_end] = qiime_upper
                dart_upper_store[upper_start:upper_end] = dart_upper

                residuals = dart_upper - qiime_upper
                count += pair_count
                sum_qiime += float(qiime_upper.sum(dtype=np.float64))
                sum_dart += float(dart_upper.sum(dtype=np.float64))
                sum_qiime_sq += float(np.dot(qiime_upper, qiime_upper))
                sum_dart_sq += float(np.dot(dart_upper, dart_upper))
                sum_product += float(np.dot(qiime_upper, dart_upper))
                sum_squared_error += float(np.dot(residuals, residuals))
                sum_absolute_error += float(np.abs(residuals).sum(dtype=np.float64))
                if pair_count:
                    maximum_absolute_error = max(
                        maximum_absolute_error, float(np.abs(residuals).max())
                    )

                if exact_spearman:
                    sampled_qiime_parts.append(qiime_upper.copy())
                    sampled_dart_parts.append(dart_upper.copy())
                else:
                    assert sampled_rows is not None and sampled_columns is not None
                    sample_end = int(
                        np.searchsorted(sampled_rows, row_index, side="right")
                    )
                    columns = sampled_columns[sample_start:sample_end]
                    sampled_qiime[sample_start:sample_end] = qiime_row[columns]
                    sampled_dart[sample_start:sample_end] = dart_row[columns]
                    sample_start = sample_end

            qiime_upper_store.flush()
            dart_upper_store.flush()
            del qiime_upper_store, dart_upper_store

        if count != total_pairs:
            raise ValueError(f"Expected {total_pairs} distance pairs, observed {count}")
        if max(max_abs_diagonal.values()) > 1e-8:
            raise ValueError(
                "Distance matrix diagonal is nonzero: "
                f"QIIME={max_abs_diagonal['qiime']}, "
                f"DART={max_abs_diagonal['dart']}"
            )
        if max(max_symmetry_error.values()) > 1e-8:
            raise ValueError(
                "Distance matrix is asymmetric: "
                f"QIIME={max_symmetry_error['qiime']}, "
                f"DART={max_symmetry_error['dart']}"
            )
        if exact_spearman:
            sampled_qiime = np.concatenate(sampled_qiime_parts)
            sampled_dart = np.concatenate(sampled_dart_parts)

    numerator = count * sum_product - sum_qiime * sum_dart
    denominator = np.sqrt(
        (count * sum_qiime_sq - sum_qiime**2)
        * (count * sum_dart_sq - sum_dart**2)
    )
    pearson = numerator / denominator
    spearman = spearmanr(sampled_qiime, sampled_dart).statistic
    residual_sample = sampled_dart - sampled_qiime
    metrics: dict[str, float | int | str] = {
        "sample_count": sample_count,
        "distance_pair_count": count,
        "pearson": float(pearson),
        "pearson_scope": "all_pairs",
        "spearman": float(spearman),
        "spearman_scope": "all_pairs" if exact_spearman else "uniform_pair_sample",
        "spearman_pair_count": len(sampled_qiime),
        "spearman_seed": seed,
        "rmse": float(np.sqrt(sum_squared_error / count)),
        "rmse_scope": "all_pairs",
        "mae": float(sum_absolute_error / count),
        "maximum_absolute_error": maximum_absolute_error,
        "qiime_max_abs_diagonal": max_abs_diagonal["qiime"],
        "dart_max_abs_diagonal": max_abs_diagonal["dart"],
        "qiime_max_symmetry_error": max_symmetry_error["qiime"],
        "dart_max_symmetry_error": max_symmetry_error["dart"],
    }
    return ids, metrics, sampled_qiime, sampled_dart, residual_sample


def ordination_coordinates(path: Path, dimensions: int) -> tuple[list[str], np.ndarray]:
    ids: list[str] = []
    rows: list[list[float]] = []
    in_sites = False
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not in_sites:
                if line.startswith("Site\t"):
                    in_sites = True
                continue
            if not line.strip():
                break
            fields = line.rstrip("\n").split("\t", dimensions + 1)
            if len(fields) < dimensions + 1:
                raise ValueError(f"Ordination has fewer than {dimensions} axes: {path}")
            ids.append(fields[0])
            rows.append([float(value) for value in fields[1 : dimensions + 1]])
    if not ids:
        raise ValueError(f"Ordination contains no sample coordinates: {path}")
    return ids, np.asarray(rows, dtype=np.float64)


def ordination_for_ids(path: Path, ids: list[str], dimensions: int) -> np.ndarray:
    candidate_ids, candidate = ordination_coordinates(path, dimensions)
    candidate_index = {sample_id: index for index, sample_id in enumerate(candidate_ids)}
    missing = set(ids) - candidate_index.keys()
    extra = set(candidate_ids) - set(ids)
    if missing or extra:
        raise ValueError(
            f"Ordination sample IDs differ from distance matrices for {path}: "
            f"distance-only={len(missing)}, ordination-only={len(extra)}"
        )
    return candidate[[candidate_index[sample_id] for sample_id in ids]]


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
    if args.dimensions < 2 or args.permutations < 0 or args.spearman_pairs < 0:
        raise ValueError("dimensions must be >=2; permutations and pairs >=0")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    distance_ids, distance_metrics, qiime_values, dart_values, residuals = (
        compare_distance_matrices(
            args.qiime_distance,
            args.dart_distance,
            spearman_pairs=args.spearman_pairs,
            seed=args.seed,
            scratch_dir=args.out_dir,
        )
    )
    exact_qiime = ordination_for_ids(
        args.qiime_ordination, distance_ids, args.dimensions
    )
    exact_dart = ordination_for_ids(
        args.dart_exact_ordination, distance_ids, args.dimensions
    )
    dart_fpcoa = ordination_for_ids(
        args.dart_ordination, distance_ids, args.dimensions
    )
    procrustes_metrics = {
        "distance_only": procrustes_test(
            exact_qiime,
            exact_dart,
            permutations=args.permutations,
            seed=args.seed,
        ),
        "dart_fpcoa_only": procrustes_test(
            exact_dart,
            dart_fpcoa,
            permutations=args.permutations,
            seed=args.seed + 1,
        ),
        "end_to_end": procrustes_test(
            exact_qiime,
            dart_fpcoa,
            permutations=args.permutations,
            seed=args.seed + 2,
        ),
    }
    summary = {
        "sample_counts": {
            "total": distance_metrics["sample_count"],
            "ordination": len(distance_ids),
        },
        "distance_metrics": distance_metrics,
        "procrustes": procrustes_metrics,
    }
    with (args.out_dir / "comparison_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    with (args.out_dir / "distance_residuals.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("qiime_distance", "dart_distance", "residual"))
        output_count = min(len(qiime_values), 100_000)
        writer.writerows(
            zip(
                qiime_values[:output_count],
                dart_values[:output_count],
                residuals[:output_count],
                strict=True,
            )
        )

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
    return 0


def main() -> int:
    return run(parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
