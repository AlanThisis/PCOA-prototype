from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import compare_subsampling_procrustes as compare


def write_ordination(path: Path, ids: list[str], values: np.ndarray) -> None:
    with path.open("w") as handle:
        handle.write("Eigvals\t2\n1\t0.5\n\nSite\t2\n")
        for sample_id, (first, second) in zip(ids, values, strict=True):
            handle.write(f"{sample_id}\t{first}\t{second}\n")
        handle.write("\n")


def test_pairwise_and_four_way_cohorts(tmp_path: Path) -> None:
    points = {
        "a": (0.0, 0.0), "b": (1.0, 0.0), "c": (0.0, 1.0),
        "d": (1.0, 1.0), "e": (2.0, 1.0), "f": (1.0, 2.0),
    }
    membership = {
        "full": ["a", "b", "c", "d", "e", "f"],
        "sub50": ["f", "e", "d", "c", "b"],
        "sub25": ["b", "c", "d", "e"],
        "sub10": ["e", "d", "c", "b"],
    }
    paths = {}
    for level, ids in membership.items():
        paths[level] = tmp_path / f"{level}.txt"
        write_ordination(paths[level], ids, np.asarray([points[i] for i in ids]))

    out = tmp_path / "out"
    summary = compare.run(paths, out)
    assert summary["four_way_shared_samples"] == 4
    assert len(summary["comparisons"]) == 12
    assert all(row["sample_count"] == 4 for row in summary["comparisons"]
               if row["cohort"] == "four_way")
    assert next(row for row in summary["comparisons"]
                if row["level_a"] == "full" and row["level_b"] == "sub50"
                and row["cohort"] == "pairwise")["sample_count"] == 5
    assert max(row["m2"] for row in summary["comparisons"]) < 1e-12
    assert (out / "procrustes_pairwise.png").is_file()
    assert (out / "procrustes_four_way.png").is_file()


def test_duplicate_and_nonfinite_ids_rejected(tmp_path: Path) -> None:
    path = tmp_path / "ordination.txt"
    write_ordination(path, ["a", "a", "b"], np.asarray([[0, 0], [1, 0], [0, 1]]))
    with pytest.raises(ValueError, match="Duplicate"):
        compare.load_coordinates(path)
    write_ordination(path, ["a", "b", "c"],
                     np.asarray([[0, 0], [float("nan"), 0], [0, 1]]))
    with pytest.raises(ValueError, match="non-finite"):
        compare.load_coordinates(path)
