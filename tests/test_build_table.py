from __future__ import annotations

import sys
import types
from pathlib import Path

import pandas as pd

from build_table import (
    build_feature_table_from_biom,
    normalize_sample_id,
    read_feature_table_from_biom,
)


def test_normalize_sample_id_strips_demux_suffix() -> None:
    assert normalize_sample_id("L5S222_17_L001_R1_001") == "L5S222"
    assert normalize_sample_id("SRR27336825_1") == "SRR27336825"
    assert normalize_sample_id("ERR10317641_2") == "ERR10317641"
    assert normalize_sample_id("custom_sample_1") == "custom_sample_1"


def test_build_feature_table_from_biom_groups_normalized_samples(
    tmp_path: Path, monkeypatch
) -> None:
    biom_fp = tmp_path / "all.biom"
    biom_fp.write_text("placeholder\n")

    source_table = pd.DataFrame(
        {
            "L5S222_17_L001_R1_001": {"feat_b": 2, "feat_a": 3},
            "L5S222_18_L001_R1_001": {"feat_b": 1, "feat_a": 4},
            "L2S309_33_L001_R1_001": {"feat_b": 0, "feat_a": 5},
        }
    )

    def fake_read_feature_table_from_biom(_biom_fp: Path) -> pd.DataFrame:
        return source_table.T

    monkeypatch.setattr(
        "build_table.read_feature_table_from_biom",
        fake_read_feature_table_from_biom,
    )

    table = build_feature_table_from_biom(biom_fp)

    assert table.index.name == "sample_id"
    assert table.columns.name == "feature_id"
    assert list(table.index) == ["L2S309", "L5S222"]
    assert list(table.columns) == ["feat_a", "feat_b"]
    assert table.loc["L5S222", "feat_a"] == 7
    assert table.loc["L5S222", "feat_b"] == 3
    assert table.loc["L2S309", "feat_a"] == 5
    assert table.loc["L2S309", "feat_b"] == 0


def test_read_feature_table_from_biom_uses_load_table(monkeypatch) -> None:
    class FakeBiomTable:
        def to_dataframe(self, dense: bool = False) -> pd.DataFrame:
            assert dense is True
            return pd.DataFrame(
                {
                    "sample_a": {"feat_x": 1, "feat_y": 2},
                    "sample_b": {"feat_x": 0, "feat_y": 3},
                }
            )

    def fake_load_table(path: str):
        assert path == "/tmp/mock.biom"
        return FakeBiomTable()

    fake_biom_module = types.SimpleNamespace(load_table=fake_load_table)
    monkeypatch.setitem(sys.modules, "biom", fake_biom_module)

    table = read_feature_table_from_biom(Path("/tmp/mock.biom"))

    assert list(table.index) == ["sample_a", "sample_b"]
    assert list(table.columns) == ["feat_x", "feat_y"]
    assert table.loc["sample_a", "feat_y"] == 2
