from __future__ import annotations

import sys
from pathlib import Path

import pytest

import run_deblur
from run_deblur import (
    discover_inputs,
    parse_args,
    run_deblur_workflow,
    stage_inputs_for_deblur,
)


def test_discover_inputs_finds_all_matching_forward_fastqs(tmp_path: Path) -> None:
    for name in (
        "SRR100_1.fastq.gz",
        "L2S309_33_L001_R1_001.fastq.gz",
        "SRR999_2.fastq.gz",
    ):
        (tmp_path / name).write_bytes(b"non-empty")

    found = discover_inputs(tmp_path)

    assert [fastq_path.name for fastq_path in found] == [
        "L2S309_33_L001_R1_001.fastq.gz",
        "SRR100_1.fastq.gz",
    ]


def test_discover_inputs_raises_for_empty_matching_file(tmp_path: Path) -> None:
    (tmp_path / "SRR100_1.fastq.gz").write_bytes(b"")

    with pytest.raises(RuntimeError, match="empty"):
        discover_inputs(tmp_path)


def test_discover_inputs_raises_when_no_matching_files(tmp_path: Path) -> None:
    (tmp_path / "SRR100_2.fastq.gz").write_bytes(b"non-empty")

    with pytest.raises(
        FileNotFoundError, match=r"\*_1.fastq.gz.*\*_R1_001.fastq.gz"
    ):
        discover_inputs(tmp_path)


def test_parse_args_parses_jobs_to_start(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_deblur.py",
            "--jobs-to-start",
            "4",
            "--trim-length",
            "200",
            "--min-reads",
            "0",
        ],
    )

    args = parse_args()

    assert args.jobs_to_start == 4
    assert args.trim_length == 200
    assert args.min_reads == 0


def test_run_deblur_workflow_builds_expected_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    seqs_dir = tmp_path / "inputs"
    seqs_dir.mkdir()
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    workflow_output_dir = work_dir / "workflow"
    commands: list[tuple[list[str], Path | None]] = []

    def fake_run_command(args: list[str], cwd: Path | None = None) -> None:
        workflow_output_dir.mkdir()
        commands.append((args, cwd))

    monkeypatch.setattr(run_deblur, "run_command", fake_run_command)

    returned_output_dir = run_deblur_workflow(
        seqs_dir,
        work_dir,
        250,
        "1,0.06,0.02",
        0,
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
                str(seqs_dir),
                "--output-dir",
                str(workflow_output_dir),
                "--trim-length",
                "250",
                "--error-dist",
                "1,0.06,0.02",
                "--min-reads",
                "0",
                "--jobs-to-start",
                "3",
                "--keep-tmp-files",
                "--overwrite",
            ],
            None,
        )
    ]


def test_stage_inputs_for_deblur_stages_only_discovered_fastqs(tmp_path: Path) -> None:
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    fastq_a = source_dir / "SRR100_1.fastq.gz"
    fastq_b = source_dir / "L2S309_33_L001_R1_001.fastq.gz"
    fastq_a.write_bytes(b"sample-a")
    fastq_b.write_bytes(b"sample-b")

    staging_dir = tmp_path / "staging"
    stage_inputs_for_deblur([fastq_b, fastq_a], staging_dir)

    staged_names = sorted(path.name for path in staging_dir.iterdir())
    assert staged_names == [fastq_b.name, fastq_a.name]


def test_stage_inputs_for_deblur_raises_on_duplicate_basenames(tmp_path: Path) -> None:
    run_a = tmp_path / "run_a"
    run_b = tmp_path / "run_b"
    run_a.mkdir()
    run_b.mkdir()
    fastq_a = run_a / "sample_1.fastq.gz"
    fastq_b = run_b / "sample_1.fastq.gz"
    fastq_a.write_bytes(b"sample-a")
    fastq_b.write_bytes(b"sample-b")

    with pytest.raises(RuntimeError, match="Duplicate FASTQ basenames"):
        stage_inputs_for_deblur([fastq_a, fastq_b], tmp_path / "staging")
