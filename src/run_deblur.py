#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
import tempfile
from pathlib import Path

from pipeline_lib import discover_inputs, resolve_executable, run_command


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
    return parser.parse_args()


def run_deblur_workflow(
    seqs_dir: Path,
    work_dir: Path,
    trim_length: int,
    error_dist: str,
    min_reads: int,
    jobs_to_start: int,
    deblur_executable_path: str,
) -> Path:
    workflow_output_dir = work_dir / "workflow"
    run_command(
        [
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
            "--keep-tmp-files",
            "--overwrite",
        ]
    )
    if not workflow_output_dir.exists():
        raise FileNotFoundError(f"Deblur workflow output missing: {workflow_output_dir}")
    return workflow_output_dir


def stage_inputs_for_deblur(fastq_paths: list[Path], staging_dir: Path) -> Path:
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

    for basename, source_fastq_path in sorted(by_name.items()):
        staged_fastq_path = staging_dir / basename
        try:
            staged_fastq_path.symlink_to(source_fastq_path)
        except OSError:
            shutil.copy2(source_fastq_path, staged_fastq_path)

    return staging_dir


def main() -> int:
    args = parse_args()
    args.data_dir = args.data_dir.resolve()
    args.work_dir = args.work_dir.resolve()
    args.work_dir.mkdir(parents=True, exist_ok=True)

    deblur_executable_path = resolve_executable("deblur")
    fastq_paths = discover_inputs(args.data_dir)
    print(
        f"Discovered {len(fastq_paths)} forward-read FASTQs under {args.data_dir}.",
        flush=True,
    )
    print("Running deblur workflow directly on discovered FASTQs.", flush=True)

    with tempfile.TemporaryDirectory(
        prefix="deblur-inputs-", dir=str(args.work_dir)
    ) as tmp_dir:
        staged_inputs_dir = stage_inputs_for_deblur(fastq_paths, Path(tmp_dir))
        workflow_output_dir = run_deblur_workflow(
            staged_inputs_dir,
            args.work_dir,
            args.trim_length,
            args.error_dist,
            args.min_reads,
            args.jobs_to_start,
            deblur_executable_path,
        )
    print(f"Finished. Deblur outputs written under: {workflow_output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
