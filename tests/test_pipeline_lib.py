from __future__ import annotations

import gzip
import csv
import subprocess
import sys
from pathlib import Path

import pytest

from pipeline_lib import (
    TimingRecorder,
    fastq_to_trimmed_fasta,
    parse_deblur_clean_fasta,
    run_command,
    run_timed_main,
)


def test_fastq_to_trimmed_fasta_keeps_only_reads_at_trim_length(tmp_path: Path) -> None:
    fastq_fp = tmp_path / "sample_1.fastq.gz"
    fasta_fp = tmp_path / "trimmed.fasta"
    trim_length = 5

    with gzip.open(fastq_fp, "wt") as fout:
        fout.write(
            "@read_short\n"
            "ACG\n"
            "+\n"
            "!!!\n"
            "@read_exact\n"
            "ACGTA\n"
            "+\n"
            "!!!!!\n"
            "@read_long\n"
            "ACGTACGT\n"
            "+\n"
            "!!!!!!!!\n"
        )

    kept = fastq_to_trimmed_fasta(fastq_fp, fasta_fp, trim_length)

    assert kept == 2
    assert fasta_fp.read_text() == ">read_exact\nACGTA\n>read_long\nACGTA\n"


def test_parse_deblur_clean_fasta_parses_sizes_and_counts_last_sequence(tmp_path: Path) -> None:
    clean_fp = tmp_path / "mock.clean"
    clean_fp.write_text(
        ">seqA;size=4;\n"
        "ACGT\n"
        ">seqB;size=2;\n"
        "TGCA\n"
        ">seqC;size=7;\n"
        "ACGT\n"
    )

    counts = parse_deblur_clean_fasta(clean_fp)

    assert counts["ACGT"] == 11
    assert counts["TGCA"] == 2


def read_timing_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_timing_recorder_writes_completed_and_skipped_steps(tmp_path: Path) -> None:
    timings_path = tmp_path / "timings.tsv"
    timing = TimingRecorder(timings_path, component="test")

    with timing.step("operation", item="sample-a"):
        running_rows = read_timing_rows(timings_path)
        assert running_rows[0]["status"] == "running"

    timing.skipped("cached_operation", message="output exists")

    rows = read_timing_rows(timings_path)
    assert [row["step"] for row in rows] == [
        "operation",
        "operation",
        "cached_operation",
    ]
    assert [row["status"] for row in rows] == ["running", "completed", "skipped"]
    assert rows[1]["exit_code"] == "0"
    assert float(rows[1]["elapsed_seconds"]) >= 0
    assert rows[2]["message"] == "output exists"


def test_timed_command_records_failure_before_reraising(tmp_path: Path) -> None:
    timings_path = tmp_path / "timings.tsv"
    timing = TimingRecorder(timings_path, component="test")
    command = [sys.executable, "-c", "raise SystemExit(7)"]

    with pytest.raises(subprocess.CalledProcessError):
        run_command(command, timing=timing, step="failing_command")

    rows = read_timing_rows(timings_path)
    assert rows[0]["status"] == "running"
    row = rows[1]
    assert row["status"] == "failed"
    assert row["exit_code"] == "7"
    assert row["command"]


def test_disabled_timing_recorder_does_not_write_a_file(tmp_path: Path) -> None:
    timing = TimingRecorder(None, component="test")
    with timing.step("operation"):
        pass
    timing.skipped("cached")

    assert list(tmp_path.iterdir()) == []


def test_timed_main_records_nonzero_return_as_failed(tmp_path: Path) -> None:
    timings_path = tmp_path / "timings.tsv"
    timing = TimingRecorder(timings_path, component="test")

    exit_code = run_timed_main(timing, lambda: 3)

    assert exit_code == 3
    rows = read_timing_rows(timings_path)
    assert rows[0]["status"] == "running"
    row = rows[1]
    assert row["step"] == "total"
    assert row["status"] == "failed"
    assert row["exit_code"] == "3"
