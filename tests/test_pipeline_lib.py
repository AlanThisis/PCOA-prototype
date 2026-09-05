from __future__ import annotations

import csv
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from pipeline_lib import (
    TimingRecorder,
    run_command,
    run_timed_main,
)


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


def test_timing_recorder_writes_valid_rows_from_concurrent_steps(tmp_path: Path) -> None:
    timings_path = tmp_path / "timings.tsv"
    timing = TimingRecorder(timings_path, component="concurrent-test")

    def record_step(index: int) -> None:
        with timing.step("operation", item=str(index)):
            pass

    with ThreadPoolExecutor(max_workers=8) as executor:
        list(executor.map(record_step, range(40)))

    rows = read_timing_rows(timings_path)
    assert len(rows) == 80
    assert {row["item"] for row in rows} == {str(index) for index in range(40)}
    assert all(row.get(None) is None for row in rows)
    assert [row["status"] for row in rows].count("running") == 40
    assert [row["status"] for row in rows].count("completed") == 40
