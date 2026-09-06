#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import shutil
import subprocess
import tempfile
from pathlib import Path

from pipeline_lib import (
    TimingRecorder,
    add_timing_argument,
    discover_inputs,
    resolve_executable,
    run_command,
    run_timed_main,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run forward-read Deblur processing for all discovered samples."
    )
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--work-dir", type=Path, default=Path("work/deblur"))
    parser.add_argument("--trim-length", type=int, default=150)
    parser.add_argument(
        "--error-dist",
        default="1,0.06,0.02,0.02,0.01,0.005,0.005,0.005,0.001,0.001,0.001,0.0005",
    )
    parser.add_argument(
        "--min-reads",
        type=int,
        default=0,
        help=(
            "Retain only features observed at least this many times across all "
            "samples. Use 0 to disable cross-sample filtering."
        ),
    )
    parser.add_argument("--jobs-to-start", type=int, default=1)
    parser.add_argument(
        "--keep-tmp-files",
        action="store_true",
        help="Retain Deblur workflow internals for debugging (large; disabled by default).",
    )
    add_timing_argument(parser)
    return parser.parse_args()


def run_deblur_workflow(
    seqs_dir: Path,
    work_dir: Path,
    trim_length: int,
    error_dist: str,
    min_reads: int,
    jobs_to_start: int,
    deblur_executable_path: str,
    timing: TimingRecorder | None = None,
    keep_tmp_files: bool = False,
) -> Path:
    timing = timing or TimingRecorder(None, component="run_deblur")
    workflow_output_dir = work_dir / "workflow"
    command = [
        deblur_executable_path,
        "workflow",
        "--seqs-fp",
        str(seqs_dir),
        "--output-dir",
        str(workflow_output_dir),
        "--trim-length",
        str(trim_length),
        "--error-dist",
        error_dist,
        "--min-reads",
        str(min_reads),
        "--jobs-to-start",
        str(jobs_to_start),
    ]
    if keep_tmp_files:
        command.append("--keep-tmp-files")
    command.append("--overwrite")
    run_command(
        command,
        timing=timing,
        step="deblur_workflow",
        item=str(workflow_output_dir),
    )
    if not workflow_output_dir.exists():
        raise FileNotFoundError(f"Deblur workflow output missing: {workflow_output_dir}")
    return workflow_output_dir


def filter_invalid_fastq_records(source: Path, destination: Path) -> tuple[int, int]:
    total = 0
    dropped = 0
    with gzip.open(source, "rb") as input_handle, gzip.open(destination, "wb") as output_handle:
        while True:
            record = [input_handle.readline() for _ in range(4)]
            if not record[0]:
                break
            total += 1
            header, sequence, separator, quality = record
            sequence = sequence.rstrip(b"\r\n")
            quality = quality.rstrip(b"\r\n")
            valid = (
                all(record)
                and header.startswith(b"@")
                and separator.startswith(b"+")
                and len(sequence) == len(quality)
                and all(33 <= value <= 95 for value in quality)
            )
            if valid:
                output_handle.writelines(record)
            else:
                dropped += 1
    return total, dropped


def stage_inputs_for_deblur(
    fastq_paths: list[Path],
    staging_dir: Path,
    *,
    sanitize_invalid_qualities: bool = False,
    sanitation_report: Path | None = None,
) -> Path:
    staging_dir.mkdir(parents=True, exist_ok=True)

    by_name: dict[str, Path] = {}
    for fastq_path in fastq_paths:
        existing = by_name.get(fastq_path.name)
        if existing is not None and existing != fastq_path:
            raise RuntimeError(
                "Duplicate FASTQ basenames discovered; cannot stage uniquely for Deblur: "
                f"{existing} and {fastq_path}"
            )
        by_name[fastq_path.name] = fastq_path

    sanitation_rows: list[tuple[str, int, int]] = []
    for basename, source_fastq_path in sorted(by_name.items()):
        staged_fastq_path = staging_dir / basename
        if sanitize_invalid_qualities:
            candidate = staging_dir / f".{basename}.filtered"
            total, dropped = filter_invalid_fastq_records(source_fastq_path, candidate)
            if dropped:
                candidate.replace(staged_fastq_path)
                sanitation_rows.append((basename, total, dropped))
                continue
            candidate.unlink()
        try:
            staged_fastq_path.symlink_to(source_fastq_path)
        except OSError:
            shutil.copy2(source_fastq_path, staged_fastq_path)

    if sanitation_report is not None:
        sanitation_report.parent.mkdir(parents=True, exist_ok=True)
        with sanitation_report.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(("fastq", "records_scanned", "records_dropped"))
            writer.writerows(sanitation_rows)

    return staging_dir


def run(args: argparse.Namespace, timing: TimingRecorder) -> int:
    args.data_dir = args.data_dir.resolve()
    args.work_dir = args.work_dir.resolve()
    args.work_dir.mkdir(parents=True, exist_ok=True)

    deblur_executable_path = resolve_executable("deblur")
    with timing.step("discover_forward_fastqs", item=str(args.data_dir)):
        fastq_paths = discover_inputs(args.data_dir)
    print(
        f"Discovered {len(fastq_paths)} forward-read FASTQs under {args.data_dir}.",
        flush=True,
    )
    print("Running deblur workflow directly on discovered FASTQs.", flush=True)

    with tempfile.TemporaryDirectory(
        prefix="deblur-inputs-", dir=str(args.work_dir)
    ) as tmp_dir:
        with timing.step("stage_deblur_inputs", item=f"{len(fastq_paths)} FASTQs"):
            staged_inputs_dir = stage_inputs_for_deblur(fastq_paths, Path(tmp_dir))
        try:
            workflow_output_dir = run_deblur_workflow(
                staged_inputs_dir,
                args.work_dir,
                args.trim_length,
                args.error_dist,
                args.min_reads,
                args.jobs_to_start,
                deblur_executable_path,
                timing,
                args.keep_tmp_files,
            )
        except subprocess.CalledProcessError:
            report = args.work_dir / "invalid_fastq_records.tsv"
            with tempfile.TemporaryDirectory(
                prefix="deblur-sanitized-inputs-", dir=str(args.work_dir)
            ) as sanitized_tmp_dir:
                sanitized_inputs = stage_inputs_for_deblur(
                    fastq_paths,
                    Path(sanitized_tmp_dir),
                    sanitize_invalid_qualities=True,
                    sanitation_report=report,
                )
                with report.open(encoding="utf-8") as handle:
                    dropped = sum(
                        int(row["records_dropped"])
                        for row in csv.DictReader(handle, delimiter="\t")
                    )
                if dropped == 0:
                    report.unlink()
                    raise
                print(
                    f"Deblur failed with {dropped} invalid FASTQ record(s); "
                    f"retrying with sanitized temporary inputs. Report: {report}",
                    flush=True,
                )
                workflow_output_dir = run_deblur_workflow(
                    sanitized_inputs,
                    args.work_dir,
                    args.trim_length,
                    args.error_dist,
                    args.min_reads,
                    args.jobs_to_start,
                    deblur_executable_path,
                    timing,
                    args.keep_tmp_files,
                )
    print(f"Finished. Deblur outputs written under: {workflow_output_dir}")
    return 0


def main() -> int:
    args = parse_args()
    timing = TimingRecorder(args.timings_tsv, component="run_deblur")
    return run_timed_main(timing, lambda: run(args, timing))


if __name__ == "__main__":
    raise SystemExit(main())
