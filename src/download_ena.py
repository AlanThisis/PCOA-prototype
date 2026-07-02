#!/usr/bin/env python3
"""Download R1 FASTQs from ENA for a list of run accessions.

Usage
-----
# Single run accession:
  python src/download_ena.py ERR123456 --out-dir data/fastq_data/my_study

# Multiple accessions from a newline-separated file:
  python src/download_ena.py --accession-list accessions.txt --out-dir data/fastq_data/my_study

# Also download R2:
  python src/download_ena.py ERR123456 --out-dir data/fastq_data/my_study --r2

Accession file format: one run accession (ERR/SRR/DRR) per line; blank lines and
lines starting with '#' are ignored.

The script resolves each accession to its FASTQ FTP URL via the ENA search API,
then downloads with urllib. Files already present in --out-dir are skipped.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ENA_SEARCH_URL = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
    "?accession={acc}&result=read_run&fields=run_accession,fastq_ftp&format=json"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download R1 (and optionally R2) FASTQs from ENA by run accession."
    )
    acc_group = parser.add_mutually_exclusive_group(required=True)
    acc_group.add_argument(
        "accession",
        nargs="?",
        metavar="ACCESSION",
        help="Single run accession (ERR/SRR/DRR).",
    )
    acc_group.add_argument(
        "--accession-list",
        type=Path,
        metavar="FILE",
        help="Newline-separated file of run accessions.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="Directory to write downloaded FASTQs into.",
    )
    parser.add_argument(
        "--r2",
        action="store_true",
        default=False,
        help="Also download R2 (reverse) reads. Default: R1 only.",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=3,
        help="Number of download retries per file (default: 3).",
    )
    return parser.parse_args()


def load_accessions(args: argparse.Namespace) -> list[str]:
    if args.accession:
        return [args.accession.strip()]
    lines = args.accession_list.read_text().splitlines()
    return [ln.strip() for ln in lines if ln.strip() and not ln.startswith("#")]


def resolve_urls(acc: str, want_r2: bool, retries: int) -> list[tuple[str, str]]:
    """Return list of (url, filename) pairs for the given run accession."""
    url = ENA_SEARCH_URL.format(acc=acc)
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=20) as resp:
                data = json.loads(resp.read())
            break
        except Exception as e:
            if attempt == retries - 1:
                raise RuntimeError(f"Could not resolve URLs for {acc}: {e}") from e
            time.sleep(2)

    if not data:
        raise RuntimeError(f"ENA returned no records for accession: {acc}")

    fastq_ftp: str = data[0].get("fastq_ftp", "")
    if not fastq_ftp:
        raise RuntimeError(f"No fastq_ftp field for {acc} — may not be a run accession.")

    raw_urls = [u.strip() for u in fastq_ftp.split(";") if u.strip()]

    def make_https(u: str) -> str:
        return "https://" + u if not u.startswith("http") else u

    r1_urls = [
        u for u in raw_urls
        if u.endswith("_1.fastq.gz") or u.endswith("_R1_001.fastq.gz")
    ]
    r2_urls = [
        u for u in raw_urls
        if u.endswith("_2.fastq.gz") or u.endswith("_R2_001.fastq.gz")
    ]

    # Fall back to first URL if no explicit R1 found (single-end run)
    if not r1_urls:
        r1_urls = raw_urls[:1]

    selected = r1_urls
    if want_r2 and r2_urls:
        selected = r1_urls + r2_urls

    return [(make_https(u), Path(u).name) for u in selected]


def download_file(url: str, dest: Path, retries: int) -> None:
    for attempt in range(retries):
        try:
            urllib.request.urlretrieve(url, dest)
            return
        except Exception as e:
            if attempt == retries - 1:
                dest.unlink(missing_ok=True)
                raise RuntimeError(f"Failed to download {url}: {e}") from e
            time.sleep(3)


def main() -> int:
    args = parse_args()
    accessions = load_accessions(args)
    if not accessions:
        print("No accessions provided.", file=sys.stderr)
        return 1

    args.out_dir.mkdir(parents=True, exist_ok=True)

    failed: list[str] = []
    for acc in accessions:
        print(f"{acc}: resolving URLs ...", flush=True)
        try:
            file_pairs = resolve_urls(acc, args.r2, args.retries)
        except RuntimeError as e:
            print(f"  ERROR: {e}", flush=True)
            failed.append(acc)
            continue

        for url, filename in file_pairs:
            dest = args.out_dir / filename
            if dest.exists():
                print(f"  {filename}: already exists, skipping", flush=True)
                continue
            print(f"  {filename}: downloading ...", flush=True)
            try:
                download_file(url, dest, args.retries)
                print(f"  {filename}: done", flush=True)
            except RuntimeError as e:
                print(f"  {filename}: ERROR — {e}", flush=True)
                failed.append(acc)

    if failed:
        print(f"\nFailed accessions ({len(failed)}): {', '.join(failed)}", file=sys.stderr)
        return 1

    print(f"\nDone. Files written to {args.out_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
