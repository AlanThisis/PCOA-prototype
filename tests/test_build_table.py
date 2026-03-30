from __future__ import annotations

from collections import Counter
from pathlib import Path

from build_table import build_feature_table


def test_build_feature_table_with_overlapping_and_unique_features(
    tmp_path: Path, monkeypatch
) -> None:
    clean_a = tmp_path / "sampleA.trim250.fasta.derep.clean"
    clean_b = tmp_path / "sampleB.trim250.fasta.derep.clean"
    clean_a.write_text("placeholder\n")
    clean_b.write_text("placeholder\n")

    counters = {
        clean_a: Counter({"AAAA": 3, "CCCC": 1}),
        clean_b: Counter({"AAAA": 2, "GGGG": 5}),
    }

    def fake_parse(clean_fp: Path) -> Counter[str]:
        return counters[clean_fp]

    monkeypatch.setattr("build_table.parse_deblur_clean_fasta", fake_parse)

    table = build_feature_table([clean_a, clean_b])

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
