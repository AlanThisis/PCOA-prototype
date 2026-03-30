#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from pipeline_lib import discover_inputs, fastq_to_trimmed_fasta, resolve_executable, run_command


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
    return parser.parse_args()


def run_deblur_for_sample(
    fastq_fp: Path,
    work_dir: Path,
    trim_length: int,
    error_dist: str,
    deblur_exe: str,
) -> Path:
    sample_id = fastq_fp.name.replace("_1.fastq.gz", "")
    fasta_fp = work_dir / f"{sample_id}.trim{trim_length}.fasta"
    derep_fp = work_dir / f"{fasta_fp.name}.derep"
    clean_fp = Path(str(derep_fp) + ".clean")

    kept_reads = fastq_to_trimmed_fasta(fastq_fp, fasta_fp, trim_length)
    print(f"  {sample_id}: converted {kept_reads} reads to trimmed FASTA", flush=True)

    run_command(
        [
            deblur_exe,
            "dereplicate",
            fasta_fp.name,
            derep_fp.name,
            "--min-size",
            "2",
        ],
        cwd=work_dir,
    )
    run_command(
        [
            deblur_exe,
            "deblur-seqs",
            derep_fp.name,
            "--error-dist",
            error_dist,
        ],
        cwd=work_dir,
    )

    if not clean_fp.exists():
        raise FileNotFoundError(f"Deblur clean output missing: {clean_fp}")
    return clean_fp


def main() -> int:
    args = parse_args()
    args.data_dir = args.data_dir.resolve()
    args.work_dir = args.work_dir.resolve()
    args.work_dir.mkdir(parents=True, exist_ok=True)

    deblur_exe = resolve_executable("deblur")
    inputs = discover_inputs(args.data_dir)
    print(f"Discovered {len(inputs)} forward-read FASTQs under {args.data_dir}")

    for fp in inputs:
        print(f"Running Deblur for: {fp}", flush=True)
        run_deblur_for_sample(fp, args.work_dir, args.trim_length, args.error_dist, deblur_exe)

    print(f"Finished. Deblur outputs written under: {args.work_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
