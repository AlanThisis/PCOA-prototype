#!/usr/bin/env python3
from __future__ import annotations

import argparse
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
        description="Subsample forward FASTQ reads with seqkit sample2."
    )
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--percent",
        type=float,
        required=True,
        help="Sampling percent in (0, 100]. Example: 50 for 50%%.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=11,
        help="Random seed passed to seqkit sample2.",
    )
    add_timing_argument(parser)
    return parser.parse_args()


def validate_percent(percent: float) -> float:
    if percent <= 0 or percent > 100:
        raise ValueError(f"--percent must be in (0, 100], got {percent}")
    return percent / 100.0


def subsample_fastqs(
    input_dir: Path,
    output_dir: Path,
    percent: float,
    seed: int,
    timing: TimingRecorder | None = None,
) -> list[Path]:
    timing = timing or TimingRecorder(None, component="subsample_fastq")
    seqkit_executable = resolve_executable("seqkit")
    input_dir = input_dir.resolve()
    output_dir = output_dir.resolve()
    fraction = validate_percent(percent)
    with timing.step("discover_forward_fastqs", item=str(input_dir)):
        fastq_paths = discover_inputs(input_dir)
        fastq_paths = [
            fastq_path
            for fastq_path in fastq_paths
            if not fastq_path.is_relative_to(output_dir)
        ]
    if not fastq_paths:
        raise FileNotFoundError(
            "No forward FASTQ inputs available after excluding output directory. "
            "Choose an output directory outside the input tree or clear old outputs."
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    sampled_paths: list[Path] = []
    for fastq_path in fastq_paths:
        relative_fastq_path = fastq_path.relative_to(input_dir)
        output_fastq = output_dir / relative_fastq_path
        output_fastq.parent.mkdir(parents=True, exist_ok=True)
        command = [
            seqkit_executable,
            "sample2",
            "-p",
            str(fraction),
            "-s",
            str(seed),
            "-2",
            str(fastq_path),
            "-o",
            str(output_fastq),
        ]
        run_command(
            command,
            timing=timing,
            step="seqkit_subsampling",
            item=relative_fastq_path.as_posix(),
        )
        sampled_paths.append(output_fastq)

    return sampled_paths


def run(args: argparse.Namespace, timing: TimingRecorder) -> int:
    sampled_paths = subsample_fastqs(
        args.input_dir, args.output_dir, args.percent, args.seed, timing
    )
    print(
        f"Subsampled {len(sampled_paths)} forward FASTQs from {args.input_dir.resolve()} "
        f"to {args.output_dir.resolve()} at {args.percent}% (seed={args.seed}, two-pass).",
        flush=True,
    )
    return 0


def main() -> int:
    args = parse_args()
    timing = TimingRecorder(args.timings_tsv, component="subsample_fastq")
    return run_timed_main(timing, lambda: run(args, timing))


if __name__ == "__main__":
    raise SystemExit(main())
