from __future__ import annotations

import csv
from pathlib import Path

import unifrac
from pipeline_lib import TimingRecorder


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
