from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

import biom
import biom.util
import numpy as np
import pytest

import deblur_scheduler
from pipeline_lib import TimingRecorder


def make_entry(tmp_path: Path, sample_id: str, size: int, shard_id: int = 0):
    fastq = tmp_path / f"{sample_id}_1.fastq.gz"
    fastq.write_bytes(b"x" * size)
    return deblur_scheduler.FastqEntry(
        sample_id, "study", fastq, size, size + 10, shard_id
    )


def write_manifest(path: Path, entries: list[deblur_scheduler.FastqEntry]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            (
                "sample_id",
                "study",
                "fastq_path",
                "size_bytes",
                "weight_bytes",
                "shard_id",
                "excluded_assay",
            )
        )
        for entry in entries:
            writer.writerow(
                (
                    entry.sample_id,
                    entry.study,
                    entry.fastq_path,
                    entry.size_bytes,
                    entry.weight_bytes,
                    entry.shard_id,
                    int(entry.excluded_assay),
                )
            )


def write_workflow(work_dir: Path, sample_ids: list[str]) -> None:
    workflow = work_dir / "workflow"
    workflow.mkdir(parents=True, exist_ok=True)
    table = biom.Table(
        np.ones((1, len(sample_ids))),
        observation_ids=["ACGT"],
        sample_ids=[f"{sample_id}_1" for sample_id in sample_ids],
        input_is_dense=True,
    )
    with biom.util.biom_open(str(workflow / "all.biom"), "w") as handle:
        table.to_hdf5(handle, "test")
    (workflow / "all.seqs.fa").write_text(">ACGT\nACGT\n")


def test_split_entries_balances_cumulative_work_deterministically(tmp_path: Path) -> None:
    entries = tuple(
        make_entry(tmp_path, sample_id, size)
        for sample_id, size in (("A", 90), ("B", 80), ("C", 30), ("D", 20))
    )

    left, right = deblur_scheduler.split_entries(entries)

    assert [entry.sample_id for entry in left] == ["A", "D"]
    assert [entry.sample_id for entry in right] == ["B", "C"]


def test_scheduler_isolates_one_bad_sample_and_merges_the_rest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    entries = [make_entry(tmp_path, sample_id, 10) for sample_id in ("GOOD1", "BAD", "GOOD2")]
    manifest = tmp_path / "shards.tsv"
    write_manifest(manifest, entries)
    calls: list[tuple[str, ...]] = []

    def fake_run_command(command: list[str], **_: object) -> None:
        input_manifest = Path(command[command.index("--input-manifest") + 1])
        sample_ids = [Path(line).name.removesuffix("_1.fastq.gz") for line in input_manifest.read_text().splitlines()]
        calls.append(tuple(sample_ids))
        if "BAD" in sample_ids:
            raise RuntimeError("bad sample")
        write_workflow(Path(command[command.index("--work-dir") + 1]), sample_ids)

    monkeypatch.setattr(deblur_scheduler, "run_command", fake_run_command)
    args = deblur_scheduler.parse_args(
        [
            "--shard-manifest",
            str(manifest),
            "--work-dir",
            str(tmp_path / "work"),
            "--results-dir",
            str(tmp_path / "results"),
            "--max-failed-samples",
            "1",
        ]
    )

    assert deblur_scheduler.run(args, TimingRecorder(None, "test")) == 0

    summary = json.loads((tmp_path / "results" / "deblur_processing_summary.json").read_text())
    assert summary["counts"]["failed"] == 1
    assert summary["counts"]["completed"] == 2
    assert summary["status"] == "completed_with_exclusions"
    merged = biom.load_table(str(tmp_path / "work" / "workflow" / "all.biom"))
    assert set(merged.ids(axis="sample")) == {"GOOD1_1", "GOOD2_1"}
    assert any(call == ("BAD",) for call in calls)


def test_scheduler_blocks_when_failed_sample_tolerance_is_exceeded(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    entries = [make_entry(tmp_path, sample_id, 10) for sample_id in ("GOOD", "BAD")]
    manifest = tmp_path / "shards.tsv"
    write_manifest(manifest, entries)

    def fake_run_command(command: list[str], **_: object) -> None:
        input_manifest = Path(command[command.index("--input-manifest") + 1])
        sample_ids = [Path(line).name.removesuffix("_1.fastq.gz") for line in input_manifest.read_text().splitlines()]
        if "BAD" in sample_ids:
            raise RuntimeError("bad sample")
        write_workflow(Path(command[command.index("--work-dir") + 1]), sample_ids)

    monkeypatch.setattr(deblur_scheduler, "run_command", fake_run_command)
    args = deblur_scheduler.parse_args(
        [
            "--shard-manifest",
            str(manifest),
            "--work-dir",
            str(tmp_path / "work"),
            "--results-dir",
            str(tmp_path / "results"),
            "--max-failed-fraction",
            "0",
        ]
    )

    with pytest.raises(RuntimeError, match="exceeding allowed limit 0"):
        deblur_scheduler.run(args, TimingRecorder(None, "test"))

    assert (tmp_path / "results" / "sample_processing_status.tsv").is_file()


def test_split_marker_avoids_repeating_failed_parent_on_resume(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    entries = tuple(make_entry(tmp_path, sample_id, 10) for sample_id in ("GOOD", "BAD"))
    args = deblur_scheduler.parse_args(
        [
            "--shard-manifest",
            str(tmp_path / "unused.tsv"),
            "--work-dir",
            str(tmp_path / "work"),
            "--results-dir",
            str(tmp_path / "results"),
        ]
    )
    parent_calls = 0

    def fake_run_command(command: list[str], **_: object) -> None:
        nonlocal parent_calls
        input_manifest = Path(command[command.index("--input-manifest") + 1])
        sample_ids = [Path(line).name.removesuffix("_1.fastq.gz") for line in input_manifest.read_text().splitlines()]
        if len(sample_ids) == 2:
            parent_calls += 1
        if "BAD" in sample_ids:
            raise RuntimeError("bad sample")
        write_workflow(Path(command[command.index("--work-dir") + 1]), sample_ids)

    monkeypatch.setattr(deblur_scheduler, "run_command", fake_run_command)
    timing = TimingRecorder(None, "test")
    deblur_scheduler.run_node("shard-0000", entries, args, timing)
    deblur_scheduler.run_node("shard-0000", entries, args, timing)

    assert parent_calls == 1


def test_cache_is_rejected_when_it_does_not_cover_current_membership(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    first = make_entry(tmp_path, "A", 10)
    second = make_entry(tmp_path, "B", 10)
    args = deblur_scheduler.parse_args(
        [
            "--shard-manifest",
            str(tmp_path / "unused.tsv"),
            "--work-dir",
            str(tmp_path / "work"),
            "--results-dir",
            str(tmp_path / "results"),
        ]
    )
    node_dir = tmp_path / "work" / "nodes" / "shard-0000"
    write_workflow(node_dir, ["A"])
    deblur_scheduler.write_fastq_manifest(node_dir / "fastqs.txt", (first,))
    calls = 0

    def fake_run_command(command: list[str], **_: object) -> None:
        nonlocal calls
        calls += 1
        input_manifest = Path(command[command.index("--input-manifest") + 1])
        sample_ids = [
            Path(line).name.removesuffix("_1.fastq.gz")
            for line in input_manifest.read_text().splitlines()
        ]
        write_workflow(Path(command[command.index("--work-dir") + 1]), sample_ids)

    monkeypatch.setattr(deblur_scheduler, "run_command", fake_run_command)

    deblur_scheduler.run_node("shard-0000", (first,), args, TimingRecorder(None, "test"))
    assert calls == 0

    deblur_scheduler.run_node(
        "shard-0000", (first, second), args, TimingRecorder(None, "test")
    )
    assert calls == 1


def test_cached_superset_is_reused_and_filtered_to_current_membership(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    first = make_entry(tmp_path, "A", 10)
    excluded = make_entry(tmp_path, "EXCLUDED", 10)
    args = deblur_scheduler.parse_args(
        [
            "--shard-manifest",
            str(tmp_path / "unused.tsv"),
            "--work-dir",
            str(tmp_path / "work"),
            "--results-dir",
            str(tmp_path / "results"),
        ]
    )
    node_dir = tmp_path / "work" / "nodes" / "shard-0000"
    write_workflow(node_dir, ["A", "EXCLUDED"])
    deblur_scheduler.write_fastq_manifest(
        node_dir / "fastqs.txt", (first, excluded)
    )

    def unexpected_run(*_: object, **__: object) -> None:
        raise AssertionError("covered cache should not rerun Deblur")

    monkeypatch.setattr(deblur_scheduler, "run_command", unexpected_run)
    leaves = deblur_scheduler.run_node(
        "shard-0000", (first,), args, TimingRecorder(None, "test")
    )
    observed = deblur_scheduler.merge_workflows(leaves, tmp_path / "merged")

    assert observed == {"A_1"}
    merged = biom.load_table(str(tmp_path / "merged" / "all.biom"))
    assert list(merged.ids(axis="sample")) == ["A_1"]


def test_timeout_splits_shard_and_records_singleton_timeout(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    entries = [make_entry(tmp_path, sample_id, 10) for sample_id in ("GOOD", "SLOW")]
    manifest = tmp_path / "shards.tsv"
    write_manifest(manifest, entries)

    def fake_run_command(command: list[str], **_: object) -> None:
        input_manifest = Path(command[command.index("--input-manifest") + 1])
        sample_ids = [
            Path(line).name.removesuffix("_1.fastq.gz")
            for line in input_manifest.read_text().splitlines()
        ]
        if len(sample_ids) > 1 or "SLOW" in sample_ids:
            raise subprocess.TimeoutExpired(command, 1)
        write_workflow(Path(command[command.index("--work-dir") + 1]), sample_ids)

    monkeypatch.setattr(deblur_scheduler, "run_command", fake_run_command)
    args = deblur_scheduler.parse_args(
        [
            "--shard-manifest",
            str(manifest),
            "--work-dir",
            str(tmp_path / "work"),
            "--results-dir",
            str(tmp_path / "results"),
            "--max-failed-samples",
            "1",
            "--shard-timeout-seconds",
            "2",
            "--singleton-timeout-seconds",
            "1",
        ]
    )

    assert deblur_scheduler.run(args, TimingRecorder(None, "test")) == 0
    summary = json.loads(
        (tmp_path / "results" / "deblur_processing_summary.json").read_text()
    )
    assert summary["counts"]["completed"] == 1
    assert summary["counts"]["timed_out"] == 1
    assert summary["status"] == "completed_with_exclusions"
    status_rows = list(
        csv.DictReader(
            (tmp_path / "results" / "sample_processing_status.tsv").open(),
            delimiter="\t",
        )
    )
    assert {row["sample_id"]: row["status"] for row in status_rows} == {
        "GOOD": "completed",
        "SLOW": "timed_out",
    }


def test_scheduler_fails_fast_on_infrastructure_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    entries = tuple(make_entry(tmp_path, sample_id, 10) for sample_id in ("A", "B"))
    args = deblur_scheduler.parse_args(
        [
            "--shard-manifest",
            str(tmp_path / "unused.tsv"),
            "--work-dir",
            str(tmp_path / "work"),
            "--results-dir",
            str(tmp_path / "results"),
        ]
    )
    calls = 0

    def fake_run_command(*_: object, **kwargs: object) -> None:
        nonlocal calls
        calls += 1
        log_handle = kwargs["stdout"]
        log_handle.write("OSError: No space left on device\n")
        log_handle.flush()
        raise subprocess.CalledProcessError(1, ["run_deblur.py"])

    monkeypatch.setattr(deblur_scheduler, "run_command", fake_run_command)

    with pytest.raises(deblur_scheduler.InfrastructureFailure, match="refusing"):
        deblur_scheduler.run_node("shard-0000", entries, args, TimingRecorder(None, "test"))

    assert calls == 1
    assert not (tmp_path / "work" / "nodes" / "shard-0000" / "split.json").exists()


def test_scheduler_reports_assay_exclusion_without_processing_it(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    retained = make_entry(tmp_path, "SIXTEEN_S", 10)
    excluded_base = make_entry(tmp_path, "ITS", 10, -1)
    excluded = deblur_scheduler.FastqEntry(
        excluded_base.sample_id,
        excluded_base.study,
        excluded_base.fastq_path,
        excluded_base.size_bytes,
        excluded_base.weight_bytes,
        -1,
        True,
    )
    manifest = tmp_path / "shards.tsv"
    write_manifest(manifest, [retained, excluded])
    processed: list[str] = []

    def fake_run_command(command: list[str], **_: object) -> None:
        input_manifest = Path(command[command.index("--input-manifest") + 1])
        sample_ids = [
            Path(line).name.removesuffix("_1.fastq.gz")
            for line in input_manifest.read_text().splitlines()
        ]
        processed.extend(sample_ids)
        write_workflow(Path(command[command.index("--work-dir") + 1]), sample_ids)

    monkeypatch.setattr(deblur_scheduler, "run_command", fake_run_command)
    args = deblur_scheduler.parse_args(
        [
            "--shard-manifest",
            str(manifest),
            "--work-dir",
            str(tmp_path / "work"),
            "--results-dir",
            str(tmp_path / "results"),
        ]
    )

    deblur_scheduler.run(args, TimingRecorder(None, "test"))

    assert processed == ["SIXTEEN_S"]
    summary = json.loads((tmp_path / "results" / "deblur_processing_summary.json").read_text())
    assert summary["counts"]["excluded_assay"] == 1
