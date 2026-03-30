from __future__ import annotations

import gzip
import re
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path


SIZE_PATTERN = re.compile(r";size=(\d+);?$")


def resolve_executable(name: str) -> str:
    env_bin = Path(sys.executable).resolve().parent / name
    if env_bin.exists():
        return str(env_bin)

    resolved = shutil.which(name)
    if resolved is None:
        raise RuntimeError(f"Required executable not found in PATH or env bin: {name}")
    return resolved


def discover_inputs(data_dir: Path) -> list[Path]:
    run_files = sorted(data_dir.rglob("*_1.fastq.gz"))
    if not run_files:
        raise FileNotFoundError(
            f"No forward FASTQs matching '*_1.fastq.gz' found under {data_dir}"
        )

    empty = [str(fp) for fp in run_files if fp.stat().st_size == 0]
    if empty:
        raise RuntimeError("Input FASTQ(s) are empty: " + ", ".join(empty))

    return run_files


def fastq_to_trimmed_fasta(fastq_fp: Path, fasta_fp: Path, trim_length: int) -> int:
    kept = 0
    with gzip.open(fastq_fp, "rt") as fin, fasta_fp.open("w") as fout:
        while True:
            header = fin.readline().strip()
            if not header:
                break
            seq = fin.readline().strip()
            fin.readline()
            fin.readline()
            if len(seq) < trim_length:
                continue
            kept += 1
            seq_id = header[1:].split()[0]
            fout.write(f">{seq_id}\n{seq[:trim_length]}\n")

    if kept == 0:
        raise RuntimeError(f"No reads >= {trim_length} bp in {fastq_fp}")
    return kept


def run_command(args: list[str], cwd: Path | None = None) -> None:
    subprocess.run(args, check=True, cwd=cwd)


def parse_deblur_clean_fasta(clean_fp: Path) -> Counter[str]:
    counts: Counter[str] = Counter()
    current_count = 0
    current_seq: str | None = None

    with clean_fp.open() as fin:
        for raw_line in fin:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_seq is not None:
                    counts[current_seq] += current_count
                match = SIZE_PATTERN.search(line)
                current_count = int(match.group(1)) if match else 1
                current_seq = None
            else:
                current_seq = line

    if current_seq is not None:
        counts[current_seq] += current_count

    return counts
