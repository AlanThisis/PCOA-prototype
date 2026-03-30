from __future__ import annotations

import sys
from pathlib import Path

import pytest

import run_deblur
from run_deblur import (
    discover_inputs,
    parse_args,
    run_deblur_workflow,
    run_fastp_for_sample,
)


def test_discover_inputs_finds_all_matching_forward_fastqs(tmp_path: Path) -> None:
    for name in ("SRR100_1.fastq.gz", "SRR101_1.fastq.gz", "SRR999_2.fastq.gz"):
        (tmp_path / name).write_bytes(b"non-empty")

    found = discover_inputs(tmp_path)

    assert [fastq_path.name for fastq_path in found] == [
        "SRR100_1.fastq.gz",
        "SRR101_1.fastq.gz",
    ]


def test_discover_inputs_raises_for_empty_matching_file(tmp_path: Path) -> None:
    (tmp_path / "SRR100_1.fastq.gz").write_bytes(b"")

    with pytest.raises(RuntimeError, match="empty"):
        discover_inputs(tmp_path)


def test_discover_inputs_raises_when_no_matching_files(tmp_path: Path) -> None:
    (tmp_path / "SRR100_2.fastq.gz").write_bytes(b"non-empty")

    with pytest.raises(FileNotFoundError, match=r"\*_1.fastq.gz"):
        discover_inputs(tmp_path)


def test_parse_args_parses_jobs_to_start(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        ["run_deblur.py", "--jobs-to-start", "4", "--trim-length", "200"],
    )

    args = parse_args()

    assert args.jobs_to_start == 4
    assert args.trim_length == 200


def test_run_fastp_for_sample_builds_expected_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fastq_path = tmp_path / "SRR100_1.fastq.gz"
    fastq_path.write_bytes(b"non-empty")
    fastp_output_dir = tmp_path / "fastp"
    fastp_output_dir.mkdir()
    commands: list[tuple[list[str], Path | None]] = []

    def fake_run_command(args: list[str], cwd: Path | None = None) -> None:
        commands.append((args, cwd))

    monkeypatch.setattr(run_deblur, "run_command", fake_run_command)

    filtered_fastq_path = run_fastp_for_sample(
        fastq_path,
        fastp_output_dir,
        "/env/bin/fastp",
    )

    assert filtered_fastq_path == fastp_output_dir / "SRR100_1.fastq.gz"
    assert commands == [
        (
            [
                "/env/bin/fastp",
                "--in1",
                str(fastq_path),
                "--out1",
                str(filtered_fastq_path),
            ],
            None,
        )
    ]


def test_run_deblur_workflow_builds_expected_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fastp_output_dir = tmp_path / "fastp"
    fastp_output_dir.mkdir()
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    workflow_output_dir = work_dir / "workflow"
    commands: list[tuple[list[str], Path | None]] = []

    def fake_run_command(args: list[str], cwd: Path | None = None) -> None:
        workflow_output_dir.mkdir()
        commands.append((args, cwd))

    monkeypatch.setattr(run_deblur, "run_command", fake_run_command)

    returned_output_dir = run_deblur_workflow(
        fastp_output_dir,
        work_dir,
        250,
        "1,0.06,0.02",
        3,
        "/env/bin/deblur",
    )

    assert returned_output_dir == workflow_output_dir
    assert commands == [
        (
            [
                "/env/bin/deblur",
                "workflow",
                "--seqs-fp",
                str(fastp_output_dir),
                "--output-dir",
                str(workflow_output_dir),
                "--trim-length",
                "250",
                "--error-dist",
                "1,0.06,0.02",
                "--jobs-to-start",
                "3",
                "--keep-tmp-files",
                "--overwrite",
            ],
            None,
        )
    ]
