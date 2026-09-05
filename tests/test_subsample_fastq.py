from __future__ import annotations

import threading
import time
from pathlib import Path

import pytest

import subsample_fastq
from subsample_fastq import subsample_fastqs, validate_percent, validate_workers


def test_validate_percent_accepts_expected_range() -> None:
    assert validate_percent(50) == 0.5
    assert validate_percent(25) == 0.25
    assert validate_percent(100) == 1.0


@pytest.mark.parametrize("bad_percent", [0, -1, 101, 250])
def test_validate_percent_rejects_invalid_values(bad_percent: float) -> None:
    with pytest.raises(ValueError, match=r"\(0, 100\]"):
        validate_percent(bad_percent)


@pytest.mark.parametrize("bad_workers", [0, -1, -10])
def test_validate_workers_rejects_nonpositive_values(bad_workers: int) -> None:
    with pytest.raises(ValueError, match="at least 1"):
        validate_workers(bad_workers)


def test_subsample_fastqs_discovers_forward_reads_only_and_builds_sample2_commands(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    (input_dir / "SRR100_1.fastq.gz").write_bytes(b"non-empty")
    (input_dir / "L2S309_33_L001_R1_001.fastq.gz").write_bytes(b"non-empty")
    (input_dir / "SRR100_2.fastq.gz").write_bytes(b"non-empty")

    output_dir = tmp_path / "out"
    commands: list[list[str]] = []

    def fake_run_command(
        args: list[str], cwd: Path | None = None, **_: object
    ) -> None:
        del cwd
        commands.append(args)

    monkeypatch.setattr(subsample_fastq, "resolve_executable", lambda _: "/env/bin/seqkit")
    monkeypatch.setattr(subsample_fastq, "run_command", fake_run_command)

    sampled_paths = subsample_fastqs(input_dir, output_dir, percent=50, seed=11)

    assert [path.name for path in sampled_paths] == [
        "L2S309_33_L001_R1_001.fastq.gz",
        "SRR100_1.fastq.gz",
    ]
    assert len(commands) == 2
    for command, sampled_path in zip(commands, sampled_paths):
        assert command == [
            "/env/bin/seqkit",
            "sample2",
            "-p",
            "0.5",
            "-s",
            "11",
            "-2",
            str(input_dir.resolve() / sampled_path.name),
            "-o",
            str(sampled_path),
        ]


def test_subsample_fastqs_preserves_relative_paths_to_avoid_name_collisions(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_dir = tmp_path / "in"
    (input_dir / "lane_a").mkdir(parents=True)
    (input_dir / "lane_b").mkdir(parents=True)
    (input_dir / "lane_a" / "sample_1.fastq.gz").write_bytes(b"non-empty")
    (input_dir / "lane_b" / "sample_1.fastq.gz").write_bytes(b"non-empty")
    output_dir = tmp_path / "out"
    commands: list[list[str]] = []

    def fake_run_command(
        args: list[str], cwd: Path | None = None, **_: object
    ) -> None:
        del cwd
        commands.append(args)

    monkeypatch.setattr(subsample_fastq, "resolve_executable", lambda _: "/env/bin/seqkit")
    monkeypatch.setattr(subsample_fastq, "run_command", fake_run_command)

    sampled_paths = subsample_fastqs(input_dir, output_dir, percent=25, seed=11)

    assert sorted(path.relative_to(output_dir).as_posix() for path in sampled_paths) == [
        "lane_a/sample_1.fastq.gz",
        "lane_b/sample_1.fastq.gz",
    ]
    assert len(commands) == 2
    assert commands[0][-1].endswith("lane_a/sample_1.fastq.gz")
    assert commands[1][-1].endswith("lane_b/sample_1.fastq.gz")


def test_subsample_fastqs_excludes_existing_output_files_under_input_tree(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    output_dir = input_dir / "subsample_50"
    output_dir.mkdir()

    source_fastq = input_dir / "ERR100_1.fastq.gz"
    prior_output_fastq = output_dir / "ERR999_1.fastq.gz"
    source_fastq.write_bytes(b"non-empty")
    prior_output_fastq.write_bytes(b"non-empty")

    commands: list[list[str]] = []

    def fake_run_command(
        args: list[str], cwd: Path | None = None, **_: object
    ) -> None:
        del cwd
        commands.append(args)

    monkeypatch.setattr(subsample_fastq, "resolve_executable", lambda _: "/env/bin/seqkit")
    monkeypatch.setattr(subsample_fastq, "run_command", fake_run_command)

    sampled_paths = subsample_fastqs(input_dir, output_dir, percent=50, seed=11)

    assert sampled_paths == [output_dir / "ERR100_1.fastq.gz"]
    assert len(commands) == 1
    assert commands[0][7] == str(source_fastq.resolve())


def test_subsample_fastqs_runs_multiple_seqkit_processes_concurrently(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    for index in range(4):
        (input_dir / f"SRR{index}_1.fastq.gz").write_bytes(b"non-empty")
    output_dir = tmp_path / "out"
    commands: list[list[str]] = []
    active = 0
    maximum_active = 0
    lock = threading.Lock()

    def fake_run_command(args: list[str], **_: object) -> None:
        nonlocal active, maximum_active
        with lock:
            commands.append(args)
            active += 1
            maximum_active = max(maximum_active, active)
        time.sleep(0.02)
        with lock:
            active -= 1

    monkeypatch.setattr(subsample_fastq, "resolve_executable", lambda _: "/env/bin/seqkit")
    monkeypatch.setattr(subsample_fastq, "run_command", fake_run_command)

    sampled_paths = subsample_fastqs(
        input_dir, output_dir, percent=10, seed=11, workers=2
    )

    assert maximum_active == 2
    assert len(commands) == 4
    assert all(command[1:3] == ["--threads", "1"] for command in commands)
    assert [path.name for path in sampled_paths] == [
        "SRR0_1.fastq.gz",
        "SRR1_1.fastq.gz",
        "SRR2_1.fastq.gz",
        "SRR3_1.fastq.gz",
    ]


def test_subsample_fastqs_propagates_parallel_worker_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    (input_dir / "ERR100_1.fastq.gz").write_bytes(b"non-empty")
    (input_dir / "ERR101_1.fastq.gz").write_bytes(b"non-empty")

    def fake_run_command(args: list[str], **_: object) -> None:
        if "ERR100_1.fastq.gz" in " ".join(args):
            raise RuntimeError("seqkit failed")

    monkeypatch.setattr(subsample_fastq, "resolve_executable", lambda _: "/env/bin/seqkit")
    monkeypatch.setattr(subsample_fastq, "run_command", fake_run_command)

    with pytest.raises(RuntimeError, match="seqkit failed"):
        subsample_fastqs(
            input_dir, tmp_path / "out", percent=50, seed=11, workers=2
        )
