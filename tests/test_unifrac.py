from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

import unifrac
from pipeline_lib import TimingRecorder


def test_auto_sampling_depth_ignores_zero_count_samples(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    table = unifrac.biom.Table(
        np.array([[0, 5, 12]]),
        observation_ids=["feature-1"],
        sample_ids=["zero", "minimum-positive", "larger"],
        input_is_dense=True,
    )
    monkeypatch.setattr(unifrac.biom, "load_table", lambda _: table)

    assert unifrac.auto_sampling_depth(Path("mapped.biom")) == 5


def test_auto_sampling_depth_rejects_an_all_zero_table(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    table = unifrac.biom.Table(
        np.array([[0, 0]]),
        observation_ids=["feature-1"],
        sample_ids=["zero-1", "zero-2"],
        input_is_dense=True,
    )
    monkeypatch.setattr(unifrac.biom, "load_table", lambda _: table)

    with pytest.raises(ValueError, match="no samples retained reads"):
        unifrac.auto_sampling_depth(Path("mapped.biom"))


def test_qiime_workflow_uses_distinct_scientific_timing_steps(
    tmp_path: Path, monkeypatch
) -> None:
    calls: list[tuple[list[str], str]] = []

    def fake_run_command(
        args: list[str], cwd: Path | None = None, **kwargs: object
    ) -> None:
        del cwd
        calls.append((args, str(kwargs["step"])))

    monkeypatch.setattr(unifrac, "run_command", fake_run_command)
    monkeypatch.setattr(
        unifrac,
        "table_stats_from_qza",
        lambda *args, **kwargs: {"sample_count": 100, "feature_count": 5},
    )
    timing = TimingRecorder(None, component="unifrac")

    unifrac.run_qiime2_unifrac(
        biom_fp=tmp_path / "all.biom",
        seqs_fp=tmp_path / "all.seqs.fa",
        gg2_backbone_fp=tmp_path / "backbone.qza",
        gg2_tree_fp=tmp_path / "tree.qza",
        sampling_depth=100,
        threads=4,
        work_dir=tmp_path / "work",
        qiime="/env/bin/qiime",
        timing=timing,
    )

    commands_by_step = {step: command for command, step in calls}
    assert commands_by_step["gg2_non_v4_16s_mapping"][1:3] == [
        "greengenes2",
        "non-v4-16s",
    ]
    assert commands_by_step["rarefy_feature_table"][1:3] == [
        "feature-table",
        "rarefy",
    ]
    assert commands_by_step["unweighted_unifrac"][1:3] == [
        "diversity",
        "beta-phylogenetic",
    ]
    assert commands_by_step["pcoa"][1:3] == ["diversity", "pcoa"]
    assert "--p-number-of-dimensions" not in commands_by_step["pcoa"]


def test_qiime_workflow_records_cached_mapping_as_skipped(
    tmp_path: Path, monkeypatch
) -> None:
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    for filename in [
        "table.qza",
        "rep-seqs.qza",
        "backbone-mapped-table.qza",
        "backbone-representatives.qza",
    ]:
        (work_dir / filename).write_text("cached")

    monkeypatch.setattr(unifrac, "run_command", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        unifrac,
        "table_stats_from_qza",
        lambda *args, **kwargs: {"sample_count": 100, "feature_count": 5},
    )
    timings_path = tmp_path / "timings.tsv"
    timing = TimingRecorder(timings_path, component="unifrac")

    unifrac.run_qiime2_unifrac(
        biom_fp=tmp_path / "all.biom",
        seqs_fp=tmp_path / "all.seqs.fa",
        gg2_backbone_fp=tmp_path / "backbone.qza",
        gg2_tree_fp=tmp_path / "tree.qza",
        sampling_depth=100,
        threads=4,
        work_dir=work_dir,
        qiime="/env/bin/qiime",
        timing=timing,
    )

    with timings_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    statuses = {row["step"]: row["status"] for row in rows}
    assert statuses["import_feature_table"] == "skipped"
    assert statuses["import_representative_sequences"] == "skipped"
    assert statuses["gg2_non_v4_16s_mapping"] == "skipped"
    assert statuses["determine_sampling_depth"] == "skipped"


def test_qiime_workflow_requests_dimensions_for_fsvd(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[tuple[list[str], str]] = []
    monkeypatch.setattr(
        unifrac,
        "run_command",
        lambda args, **kwargs: calls.append((args, str(kwargs["step"]))),
    )
    monkeypatch.setattr(
        unifrac,
        "table_stats_from_qza",
        lambda *args, **kwargs: {"sample_count": 60_000, "feature_count": 5},
    )

    _, _, policy = unifrac.run_qiime2_unifrac(
        biom_fp=tmp_path / "all.biom",
        seqs_fp=tmp_path / "all.seqs.fa",
        gg2_backbone_fp=tmp_path / "backbone.qza",
        gg2_tree_fp=tmp_path / "tree.qza",
        sampling_depth=1000,
        threads=4,
        work_dir=tmp_path / "work",
        qiime="/env/bin/qiime",
        timing=TimingRecorder(None, component="unifrac"),
        pcoa_method="auto",
        pcoa_dimensions=10,
        pcoa_memory_budget_gb=128,
    )

    pcoa_command = next(command for command, step in calls if step == "pcoa")
    assert policy["selected_method"] == "fsvd"
    assert pcoa_command[pcoa_command.index("--p-number-of-dimensions") + 1] == "10"


def test_qiime_workflow_clamps_fsvd_dimensions_to_retained_samples(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[tuple[list[str], str]] = []
    monkeypatch.setattr(
        unifrac,
        "run_command",
        lambda args, **kwargs: calls.append((args, str(kwargs["step"]))),
    )
    monkeypatch.setattr(
        unifrac,
        "table_stats_from_qza",
        lambda *args, **kwargs: {"sample_count": 4, "feature_count": 5},
    )

    _, _, policy = unifrac.run_qiime2_unifrac(
        biom_fp=tmp_path / "all.biom",
        seqs_fp=tmp_path / "all.seqs.fa",
        gg2_backbone_fp=tmp_path / "backbone.qza",
        gg2_tree_fp=tmp_path / "tree.qza",
        sampling_depth=1000,
        threads=4,
        work_dir=tmp_path / "work",
        qiime="/env/bin/qiime",
        timing=TimingRecorder(None, component="unifrac"),
        pcoa_method="fsvd",
        pcoa_dimensions=10,
    )

    pcoa_command = next(command for command, step in calls if step == "pcoa")
    assert pcoa_command[pcoa_command.index("--p-number-of-dimensions") + 1] == "4"
    assert policy["requested_dimensions"] == 10
    assert policy["dimensions"] == 4


def test_qiime_workflow_rejects_fewer_than_two_retained_samples(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(unifrac, "run_command", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        unifrac,
        "table_stats_from_qza",
        lambda *args, **kwargs: {"sample_count": 1, "feature_count": 5},
    )

    with pytest.raises(ValueError, match="at least two samples"):
        unifrac.run_qiime2_unifrac(
            biom_fp=tmp_path / "all.biom",
            seqs_fp=tmp_path / "all.seqs.fa",
            gg2_backbone_fp=tmp_path / "backbone.qza",
            gg2_tree_fp=tmp_path / "tree.qza",
            sampling_depth=1000,
            threads=4,
            work_dir=tmp_path / "work",
            qiime="/env/bin/qiime",
            timing=TimingRecorder(None, component="unifrac"),
            pcoa_method="fsvd",
        )


def test_clear_input_artifact_cache_removes_only_input_derived_files(
    tmp_path: Path,
) -> None:
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    for filename in (
        "table.qza",
        "rep-seqs.qza",
        "backbone-mapped-table.qza",
        "backbone-representatives.qza",
    ):
        (work_dir / filename).write_text("cached")
    depth_export = work_dir / "_depth_check_export"
    depth_export.mkdir()
    (depth_export / "feature-table.biom").write_text("cached")
    retained = work_dir / "unweighted_unifrac_distance_matrix.qza"
    retained.write_text("retained")

    unifrac.clear_input_artifact_cache(work_dir)

    assert not (work_dir / "table.qza").exists()
    assert not (work_dir / "rep-seqs.qza").exists()
    assert not (work_dir / "backbone-mapped-table.qza").exists()
    assert not (work_dir / "backbone-representatives.qza").exists()
    assert not depth_export.exists()
    assert retained.is_file()


def test_auto_pcoa_uses_exact_when_estimate_fits_budget() -> None:
    policy = unifrac.choose_pcoa_method("auto", 10_000, 64)

    assert policy["selected_method"] == "eigh"
    assert policy["estimated_exact_memory_gb"] < 64 * 0.8


def test_auto_pcoa_uses_fsvd_when_exact_exceeds_budget() -> None:
    policy = unifrac.choose_pcoa_method("auto", 60_000, 256)

    assert policy["selected_method"] == "fsvd"
    assert policy["estimated_exact_memory_gb"] > 256 * 0.8


def test_explicit_pcoa_method_overrides_memory_estimate() -> None:
    assert unifrac.choose_pcoa_method("eigh", 100_000, 1)["selected_method"] == "eigh"
    assert unifrac.choose_pcoa_method("fsvd", 10, 100)["selected_method"] == "fsvd"


def test_distance_export_auto_skips_projected_huge_tsv() -> None:
    export_small, small_size = unifrac.should_export_distance_tsv("auto", 10_000)
    export_large, large_size = unifrac.should_export_distance_tsv("auto", 60_000)

    assert export_small is True
    assert export_large is False
    assert small_size < 20 < large_size
