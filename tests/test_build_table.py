from __future__ import annotations

from pathlib import Path
import sys
from types import SimpleNamespace

import pandas as pd
import pytest

import build_table
from build_table import build_feature_table


class FakeBiomTable:
    def __init__(self, feature_by_sample: pd.DataFrame) -> None:
        self._feature_by_sample = feature_by_sample

    def to_dataframe(self, dense: bool = True) -> pd.DataFrame:
        assert dense is True
        return self._feature_by_sample.transpose()


def test_build_feature_table_from_biom_with_overlapping_and_unique_features() -> None:
    biom_table = FakeBiomTable(
        pd.DataFrame(
            {
                "AAAA": [3, 2],
                "CCCC": [1, 0],
                "GGGG": [0, 5],
            },
            index=["sampleA", "sampleB"],
        )
    )
    table = build_feature_table(biom_table)

    assert table.shape == (2, 3)
    assert table.index.name == "sample_id"
    assert table.columns.name == "feature_sequence"
    assert list(table.index) == ["sampleA", "sampleB"]
    assert list(table.columns) == ["AAAA", "CCCC", "GGGG"]
    assert table.loc["sampleA", "AAAA"] == 3
    assert table.loc["sampleA", "CCCC"] == 1
    assert table.loc["sampleA", "GGGG"] == 0
    assert table.loc["sampleB", "AAAA"] == 2
    assert table.loc["sampleB", "CCCC"] == 0
    assert table.loc["sampleB", "GGGG"] == 5


def test_load_biom_table_raises_helpful_error_when_biom_missing(tmp_path: Path) -> None:
    biom_fp = tmp_path / "all.biom"
    biom_fp.write_text("placeholder")

    with pytest.raises(RuntimeError, match="requires the 'biom' package"):
        build_table.load_biom_table(biom_fp)


def test_load_biom_table_calls_biom_load_table(monkeypatch: pytest.MonkeyPatch) -> None:
    called = {}

    def fake_load_table(path: str) -> str:
        called["path"] = path
        return "fake-table"

    monkeypatch.setitem(
        sys.modules,
        "biom",
        SimpleNamespace(load_table=fake_load_table),
    )

    table = build_table.load_biom_table(Path("work/deblur/workflow/all.biom"))
    assert table == "fake-table"
    assert called["path"].endswith("work/deblur/workflow/all.biom")
