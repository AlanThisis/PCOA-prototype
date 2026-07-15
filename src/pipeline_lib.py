from __future__ import annotations

import argparse
import csv
import shlex
import shutil
import subprocess
import sys
import time
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterator


TIMING_COLUMNS = [
    "component",
    "step",
    "item",
    "start_utc",
    "end_utc",
    "elapsed_seconds",
    "status",
    "exit_code",
    "command",
    "message",
]


@dataclass
class TimingRecord:
    component: str
    step: str
    item: str
    start_utc: str
    end_utc: str = ""
    elapsed_seconds: str = ""
    status: str = "running"
    exit_code: str = ""
    command: str = ""
    message: str = ""


class TimingRecorder:
    """Persist optional operation-level wall-clock timing events to a TSV file."""

    def __init__(self, path: Path | None, component: str) -> None:
        self.path = path.resolve() if path is not None else None
        self.component = component
        if self.path is not None:
            self._initialize_file()

    @property
    def enabled(self) -> bool:
        return self.path is not None

    @contextmanager
    def step(
        self,
        name: str,
        *,
        item: str = "",
        command: str = "",
    ) -> Iterator[TimingRecord | None]:
        if not self.enabled:
            yield None
            return

        record = TimingRecord(
            component=self.component,
            step=name,
            item=item,
            start_utc=datetime.now(timezone.utc).isoformat(),
            command=command,
        )
        started = time.perf_counter()
        self._append(record)
        try:
            yield record
        except BaseException as exc:
            exception_exit_code = getattr(exc, "returncode", None)
            if exception_exit_code is None:
                exception_exit_code = getattr(exc, "code", 1)
            if exception_exit_code is None:
                exception_exit_code = 1
            self._finish_record(
                record,
                started,
                status="failed",
                exit_code=str(exception_exit_code),
                message=str(exc),
            )
            raise
        else:
            requested_exit_code = record.exit_code or "0"
            self._finish_record(
                record,
                started,
                status="failed" if requested_exit_code != "0" else "completed",
                exit_code=requested_exit_code,
                message=record.message,
            )

    def skipped(self, name: str, *, item: str = "", message: str = "") -> None:
        if not self.enabled:
            return
        now = datetime.now(timezone.utc).isoformat()
        record = TimingRecord(
            component=self.component,
            step=name,
            item=item,
            start_utc=now,
            end_utc=now,
            elapsed_seconds="0.000000",
            status="skipped",
            exit_code="0",
            message=message,
        )
        self._append(record)

    def _finish_record(
        self,
        record: TimingRecord,
        started: float,
        *,
        status: str,
        exit_code: str,
        message: str = "",
    ) -> None:
        record.end_utc = datetime.now(timezone.utc).isoformat()
        record.elapsed_seconds = f"{time.perf_counter() - started:.6f}"
        record.status = status
        record.exit_code = exit_code
        record.message = message
        self._append(record)

    def _initialize_file(self) -> None:
        if self.path is None:
            return
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=TIMING_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()

    def _append(self, record: TimingRecord) -> None:
        if self.path is None:
            return
        with self.path.open("a", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=TIMING_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writerow(asdict(record))


def add_timing_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--timings-tsv",
        type=Path,
        default=None,
        help="Optional TSV path for operation-level wall-clock timings.",
    )


def run_timed_main(timing: TimingRecorder, callback: Callable[[], int]) -> int:
    with timing.step("total") as record:
        exit_code = callback()
        if record is not None and exit_code != 0:
            record.exit_code = str(exit_code)
            record.message = f"process returned exit code {exit_code}"
        return exit_code


def resolve_executable(name: str) -> str:
    env_bin = Path(sys.executable).resolve().parent / name
    if env_bin.exists():
        return str(env_bin)

    resolved = shutil.which(name)
    if resolved is None:
        raise RuntimeError(f"Required executable not found in PATH or env bin: {name}")
    return resolved


def discover_inputs(data_dir: Path) -> list[Path]:
    patterns = ("*_1.fastq.gz", "*_R1_001.fastq.gz")
    discovered: set[Path] = set()
    for pattern in patterns:
        discovered.update(data_dir.rglob(pattern))
    run_files = sorted(discovered)
    if not run_files:
        raise FileNotFoundError(
            "No forward FASTQs matching '*_1.fastq.gz' or '*_R1_001.fastq.gz' "
            f"found under {data_dir}"
        )

    empty = [str(fp) for fp in run_files if fp.stat().st_size == 0]
    if empty:
        raise RuntimeError("Input FASTQ(s) are empty: " + ", ".join(empty))

    return run_files


def run_command(
    args: list[str],
    cwd: Path | None = None,
    *,
    timing: TimingRecorder | None = None,
    step: str | None = None,
    item: str = "",
) -> None:
    if timing is not None and step is not None:
        with timing.step(step, item=item, command=shlex.join(args)):
            subprocess.run(args, check=True, cwd=cwd)
        return
    subprocess.run(args, check=True, cwd=cwd)
