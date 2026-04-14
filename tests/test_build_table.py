from __future__ import annotations

from pathlib import Path

from build_table import build_feature_table_from_biom, normalize_sample_id


def test_normalize_sample_id_strips_demux_suffix() -> None:
    assert normalize_sample_id("L5S222_17_L001_R1_001") == "L5S222"
    assert normalize_sample_id("SRR27336825_1") == "SRR27336825_1"


def test_build_feature_table_from_biom_groups_normalized_samples(
    tmp_path: Path, monkeypatch
) -> None:
    biom_fp = tmp_path / "all.biom"
    biom_fp.write_text("placeholder\n")

    def fake_convert_biom_to_tsv(
        _biom_fp: Path, tsv_fp: Path, _biom_executable_path: str
    ) -> None:
        tsv_fp.write_text(
            "# Constructed from biom file\n"
            "#OTU ID\tL5S222_17_L001_R1_001\tL5S222_18_L001_R1_001\tL2S309_33_L001_R1_001\n"
            "feat_b\t2\t1\t0\n"
            "feat_a\t3\t4\t5\n"
        )

    monkeypatch.setattr("build_table.convert_biom_to_tsv", fake_convert_biom_to_tsv)

    table = build_feature_table_from_biom(biom_fp, "/env/bin/biom")

    assert table.index.name == "sample_id"
    assert table.columns.name == "feature_id"
    assert list(table.index) == ["L2S309", "L5S222"]
    assert list(table.columns) == ["feat_a", "feat_b"]
    assert table.loc["L5S222", "feat_a"] == 7
    assert table.loc["L5S222", "feat_b"] == 3
    assert table.loc["L2S309", "feat_a"] == 5
    assert table.loc["L2S309", "feat_b"] == 0
