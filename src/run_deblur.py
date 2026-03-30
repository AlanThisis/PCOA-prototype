#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from pipeline_lib import discover_inputs, resolve_executable, run_command


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run forward-read Deblur processing for all discovered samples."
    )
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--work-dir", type=Path, default=Path("work/deblur"))
    parser.add_argument("--trim-length", type=int, default=250)
    parser.add_argument(
        "--error-dist",
        default="1,0.06,0.02,0.02,0.01,0.005,0.005,0.005,0.001,0.001,0.001,0.0005",
    )
    parser.add_argument("--jobs-to-start", type=int, default=1)
    return parser.parse_args()


def run_fastp_for_sample(
    fastq_path: Path,
    fastp_output_dir: Path,
    fastp_executable_path: str,
) -> Path:
    sample_id = fastq_path.name.replace("_1.fastq.gz", "")
    filtered_fastq_path = fastp_output_dir / f"{sample_id}_1.fastq.gz"
    run_command(
        [
            fastp_executable_path,
            "--in1",
            str(fastq_path),
            "--out1",
            str(filtered_fastq_path),
        ]
    )
    return filtered_fastq_path


def run_deblur_workflow(
    fastp_output_dir: Path,
    work_dir: Path,
    trim_length: int,
    error_dist: str,
    jobs_to_start: int,
    deblur_executable_path: str,
) -> Path:
    workflow_output_dir = work_dir / "workflow"
    run_command(
        [
            deblur_executable_path,
            "workflow",
            "--seqs-fp",
            str(fastp_output_dir),
            "--output-dir",
            str(workflow_output_dir),
            "--trim-length",
            str(trim_length),
            "--error-dist",
            error_dist,
            "--jobs-to-start",
            str(jobs_to_start),
            "--keep-tmp-files",
            "--overwrite",
        ]
    )
    if not workflow_output_dir.exists():
        raise FileNotFoundError(f"Deblur workflow output missing: {workflow_output_dir}")
    return workflow_output_dir


def main() -> int:
    args = parse_args()
    args.data_dir = args.data_dir.resolve()
    args.work_dir = args.work_dir.resolve()
    args.work_dir.mkdir(parents=True, exist_ok=True)

    fastp_executable_path = resolve_executable("fastp")
    deblur_executable_path = resolve_executable("deblur")
    fastq_paths = discover_inputs(args.data_dir)
    fastp_output_dir = args.work_dir / "fastp"
    fastp_output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Discovered {len(fastq_paths)} forward-read FASTQs under {args.data_dir}")

    for fastq_path in fastq_paths:
        print(f"Running fastp for: {fastq_path}", flush=True)
        filtered_fastq_path = run_fastp_for_sample(
            fastq_path,
            fastp_output_dir,
            fastp_executable_path,
        )
        print(f"  staged filtered reads at: {filtered_fastq_path}", flush=True)

    workflow_output_dir = run_deblur_workflow(
        fastp_output_dir,
        args.work_dir,
        args.trim_length,
        args.error_dist,
        args.jobs_to_start,
        deblur_executable_path,
    )
    print(f"Finished. Deblur outputs written under: {workflow_output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
