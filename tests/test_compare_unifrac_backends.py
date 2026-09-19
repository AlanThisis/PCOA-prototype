from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import skbio

import compare_unifrac_backends as compare


def write_distance(path: Path, ids: list[str], matrix: np.ndarray) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t" + "\t".join(ids) + "\n")
        for sample_id, row in zip(ids, matrix, strict=True):
            handle.write(sample_id + "\t" + "\t".join(map(str, row)) + "\n")


def write_ordination(path: Path, ids: list[str], matrix: np.ndarray) -> None:
    result = skbio.stats.ordination.pcoa(
        skbio.DistanceMatrix(matrix, ids=ids),
        number_of_dimensions=min(3, len(ids) - 1),
    )
    result.write(str(path), format="ordination")


def test_identical_backends_pass_all_acceptance_gates(tmp_path: Path) -> None:
    ids = ["a", "b", "c", "d"]
    matrix = np.array(
        [
            [0.0, 0.1, 0.2, 0.4],
            [0.1, 0.0, 0.25, 0.35],
            [0.2, 0.25, 0.0, 0.3],
            [0.4, 0.35, 0.3, 0.0],
        ]
    )
    qiime_dm = tmp_path / "qiime.tsv"
    dart_dm = tmp_path / "dart.tsv"
    write_distance(qiime_dm, ids, matrix)
    write_distance(dart_dm, ids, matrix)
    qiime_ord = tmp_path / "qiime-ordination.txt"
    dart_ord = tmp_path / "dart-ordination.txt"
    write_ordination(qiime_ord, ids, matrix)
    write_ordination(dart_ord, ids, matrix)
    out_dir = tmp_path / "comparison"
    args = argparse.Namespace(
        qiime_distance=qiime_dm,
        dart_distance=dart_dm,
        qiime_ordination=qiime_ord,
        dart_ordination=dart_ord,
        out_dir=out_dir,
        dimensions=3,
        permutations=9,
        seed=0,
        max_distance_samples=4,
    )

    assert compare.run(args) == 0
    summary = json.loads((out_dir / "comparison_summary.json").read_text())
    assert summary["passed"] is True
    assert summary["distance_metrics"]["rmse"] == 0
    assert summary["procrustes"]["distance_only"]["m2"] < 1e-12


def test_distance_comparison_rejects_mismatched_ids(tmp_path: Path) -> None:
    first = tmp_path / "first.tsv"
    second = tmp_path / "second.tsv"
    write_distance(first, ["a", "b"], np.array([[0, 1], [1, 0]]))
    write_distance(second, ["a", "c"], np.array([[0, 1], [1, 0]]))

    try:
        compare.deterministic_common_ids(first, second, 10, 0)
    except ValueError as exc:
        assert "sample IDs differ" in str(exc)
    else:
        raise AssertionError("mismatched IDs were accepted")
