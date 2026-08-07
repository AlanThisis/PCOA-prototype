#!/usr/bin/env python3
"""Download forward-read FASTQs from ENA AMPLICON runs in a study."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlparse

import requests


ENA_FILE_REPORT_BASE = "https://www.ebi.ac.uk/ena/portal/api/filereport"
FILE_REPORT_FIELDS = [
    "run_accession",
    "library_layout",
    "library_strategy",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
]
MANIFEST_ATTEMPTS = 5
MANIFEST_COLUMNS = [
    "project_accession",
    "run_accession",
    "library_layout",
    "library_strategy",
    "url",
    "md5",
    "bytes",
    "output_filename",
]
STATUS_COLUMNS = MANIFEST_COLUMNS + [
    "status",
    "elapsed_seconds",
    "downloaded_bytes",
    "message",
]


@dataclass(frozen=True)
class DownloadSpec:
    project_accession: str
    run_accession: str
    library_layout: str
    library_strategy: str
    url: str
    md5: str
    bytes: int
    output_filename: str


@dataclass(frozen=True)
class DownloadResult:
    spec: DownloadSpec
    status: str
    elapsed_seconds: float
    downloaded_bytes: int
    message: str = ""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def format_elapsed(seconds: float) -> str:
    whole_seconds = round(seconds)
    hours, remainder = divmod(whole_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{seconds:02d}"


def file_report_params(accession: str) -> dict[str, str]:
    return {
        "accession": accession,
        "result": "read_run",
        "fields": ",".join(FILE_REPORT_FIELDS),
        "format": "tsv",
    }


def https_url(raw_url: str) -> str:
    raw_url = raw_url.strip().rstrip("/")
    if raw_url.startswith("ftp://"):
        return "https://" + raw_url.removeprefix("ftp://")
    if "://" not in raw_url:
        return "https://" + raw_url
    return raw_url


def ftp_url(raw_url: str) -> str:
    parsed = urlparse(https_url(raw_url))
    return parsed._replace(scheme="ftp").geturl()


def split_field(record: dict[str, object], name: str) -> list[str]:
    value = str(record.get(name, "") or "")
    return value.split(";") if value else []


def select_forward_download(
    accession: str, record: dict[str, object]
) -> DownloadSpec | None:
    run_accession = str(record.get("run_accession", "") or "").strip()
    layout = str(record.get("library_layout", "") or "").strip().upper()
    strategy = str(record.get("library_strategy", "") or "").strip().upper()
    urls = split_field(record, "fastq_ftp")
    checksums = split_field(record, "fastq_md5")
    sizes = split_field(record, "fastq_bytes")

    if not run_accession or strategy != "AMPLICON":
        return None
    if not (len(urls) == len(checksums) == len(sizes)):
        raise ValueError(f"Mismatched FASTQ metadata fields for {run_accession}")

    candidates = list(zip(urls, checksums, sizes))
    selected: tuple[str, str, str] | None = None
    if layout == "PAIRED":
        expected_name = f"{run_accession}_1.fastq.gz"
        selected = next(
            (item for item in candidates if Path(urlparse(item[0]).path).name == expected_name),
            None,
        )
    elif layout == "SINGLE":
        accepted_names = {
            f"{run_accession}.fastq.gz",
            f"{run_accession}_1.fastq.gz",
        }
        selected = next(
            (
                item
                for item in candidates
                if Path(urlparse(item[0]).path).name in accepted_names
            ),
            None,
        )

    if selected is None:
        return None

    raw_url, expected_md5, expected_bytes = selected
    if not expected_md5:
        raise ValueError(f"Missing FASTQ MD5 for {run_accession}")
    try:
        byte_count = int(expected_bytes)
    except ValueError as exc:
        raise ValueError(f"Invalid FASTQ byte count for {run_accession}") from exc
    if byte_count <= 0:
        raise ValueError(f"Invalid FASTQ byte count for {run_accession}: {byte_count}")

    return DownloadSpec(
        project_accession=accession,
        run_accession=run_accession,
        library_layout=layout,
        library_strategy=strategy,
        url=https_url(raw_url),
        md5=expected_md5.lower(),
        bytes=byte_count,
        output_filename=f"{run_accession}_1.fastq.gz",
    )


def parse_file_report(report_text: str) -> list[dict[str, str]]:
    reader = csv.DictReader(report_text.splitlines(), delimiter="\t")
    fields = set(reader.fieldnames or [])
    missing = sorted(set(FILE_REPORT_FIELDS) - fields)
    if missing:
        preview = report_text[:500].replace("\n", " ")
        raise RuntimeError(
            f"ENA file report missing fields {missing}; response starts with {preview!r}"
        )
    records = list(reader)
    if any(None in record for record in records):
        raise RuntimeError("ENA file report contains a malformed TSV row")
    return records


def fetch_manifest(accession: str, timeout: float) -> tuple[list[DownloadSpec], list[str]]:
    records: list[dict[str, str]] | None = None
    last_error: Exception | None = None
    for attempt in range(1, MANIFEST_ATTEMPTS + 1):
        try:
            response = requests.get(
                ENA_FILE_REPORT_BASE,
                params=file_report_params(accession),
                headers={"User-Agent": "ena-amplicon-forward-downloader"},
                timeout=timeout,
            )
            response.raise_for_status()
            records = parse_file_report(response.text)
            break
        except (requests.RequestException, RuntimeError, csv.Error) as exc:
            last_error = exc
            if attempt == MANIFEST_ATTEMPTS:
                break
            delay = min(5 * (2 ** (attempt - 1)), 60)
            print(
                f"ENA manifest attempt {attempt}/{MANIFEST_ATTEMPTS} failed for "
                f"{accession}: {exc}; retrying in {delay}s",
                file=sys.stderr,
                flush=True,
            )
            time.sleep(delay)

    if records is None:
        raise RuntimeError(
            f"ENA manifest retrieval failed for {accession} after "
            f"{MANIFEST_ATTEMPTS} attempts: {last_error}"
        )

    specs: list[DownloadSpec] = []
    skipped: list[str] = []
    for record in records:
        if str(record.get("library_strategy", "")).strip().upper() != "AMPLICON":
            continue
        spec = select_forward_download(accession, record)
        if spec is None:
            run_accession = str(record.get("run_accession", "unknown"))
            skipped.append(
                f"{run_accession}: no archive-generated forward FASTQ for its layout"
            )
        else:
            specs.append(spec)
    return specs, skipped


def file_md5(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def completed_file_is_valid(path: Path, spec: DownloadSpec) -> bool:
    return (
        path.is_file()
        and path.stat().st_size == spec.bytes
        and file_md5(path) == spec.md5
    )


def verify_download(
    spec: DownloadSpec, part_path: Path, destination: Path
) -> None:
    if part_path.stat().st_size != spec.bytes:
        raise RuntimeError(
            f"size mismatch: expected {spec.bytes}, got {part_path.stat().st_size}"
        )
    actual_md5 = file_md5(part_path)
    if actual_md5 != spec.md5:
        part_path.unlink(missing_ok=True)
        raise RuntimeError(f"MD5 mismatch: expected {spec.md5}, got {actual_md5}")
    part_path.replace(destination)


def transfer_via_curl(spec: DownloadSpec, destination: Path, timeout: float) -> int:
    part_path = destination.with_name(destination.name + ".part")
    curl = shutil.which("curl")
    if curl is None:
        raise RuntimeError("curl is not installed")

    https = https_url(spec.url)
    urls = [https, https.replace("https://", "http://", 1), ftp_url(https)]
    errors: list[str] = []
    for url in urls:
        part_path.unlink(missing_ok=True)
        result = subprocess.run(
            [
                curl,
                "--fail",
                "--location",
                "--silent",
                "--show-error",
                "--connect-timeout",
                str(timeout),
                "--speed-limit",
                "1",
                "--speed-time",
                str(timeout),
                "--retry",
                "2",
                "--retry-delay",
                "2",
                "--output",
                str(part_path),
                url,
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            errors.append(f"{url}: curl exit {result.returncode}: {result.stderr.strip()}")
            continue
        try:
            verify_download(spec, part_path, destination)
            return spec.bytes
        except (OSError, RuntimeError) as exc:
            errors.append(f"{url}: {exc}")

    part_path.unlink(missing_ok=True)
    raise RuntimeError("; ".join(errors))


def transfer_via_https(spec: DownloadSpec, destination: Path, timeout: float) -> int:
    part_path = destination.with_name(destination.name + ".part")
    download_url = spec.url.strip().rstrip("/")
    existing_bytes = part_path.stat().st_size if part_path.exists() else 0
    headers = {"User-Agent": "ena-amplicon-forward-downloader"}
    if 0 < existing_bytes < spec.bytes:
        headers["Range"] = f"bytes={existing_bytes}-"
    elif existing_bytes >= spec.bytes:
        part_path.unlink()
        existing_bytes = 0

    response = requests.get(
        download_url,
        headers=headers,
        stream=True,
        timeout=(timeout, timeout),
        allow_redirects=False,
    )
    content_type = response.headers.get("Content-Type", "").lower()
    rejected_resume = response.status_code == 416 or "text/html" in content_type
    if rejected_resume and existing_bytes > 0:
        response.close()
        part_path.unlink()
        existing_bytes = 0
        headers.pop("Range", None)
        response = requests.get(
            download_url,
            headers=headers,
            stream=True,
            timeout=(timeout, timeout),
            allow_redirects=False,
        )

    with response:
        response.raise_for_status()
        resume_accepted = existing_bytes > 0 and response.status_code == 206
        content_type = response.headers.get("Content-Type", "").lower()
        if "text/html" in content_type:
            part_path.unlink(missing_ok=True)
            raise RuntimeError(
                f"ENA returned HTML instead of FASTQ for {spec.run_accession}"
            )
        content_length = response.headers.get("Content-Length")
        expected_response_bytes = spec.bytes - existing_bytes if resume_accepted else spec.bytes
        if content_length is not None and int(content_length) != expected_response_bytes:
            part_path.unlink(missing_ok=True)
            raise RuntimeError(
                f"HTTP content length mismatch: expected {expected_response_bytes}, "
                f"got {content_length}"
            )
        mode = "ab" if resume_accepted else "wb"
        downloaded_bytes = 0
        with part_path.open(mode) as handle:
            for chunk in response.iter_content(chunk_size=8 * 1024 * 1024):
                if chunk:
                    handle.write(chunk)
                    downloaded_bytes += len(chunk)

    verify_download(spec, part_path, destination)
    return downloaded_bytes


def transfer_once(spec: DownloadSpec, destination: Path, timeout: float) -> int:
    try:
        return transfer_via_https(spec, destination, timeout)
    except (requests.HTTPError, RuntimeError) as https_error:
        try:
            return transfer_via_curl(spec, destination, timeout)
        except (OSError, RuntimeError) as curl_error:
            raise RuntimeError(
                f"Python HTTPS failed ({https_error}); curl fallback failed ({curl_error})"
            ) from curl_error


def download_one(
    spec: DownloadSpec,
    output_dir: Path,
    retries: int,
    timeout: float,
) -> DownloadResult:
    started = time.perf_counter()
    destination = output_dir / spec.output_filename
    if completed_file_is_valid(destination, spec):
        destination.with_name(destination.name + ".part").unlink(missing_ok=True)
        return DownloadResult(
            spec=spec,
            status="skipped",
            elapsed_seconds=time.perf_counter() - started,
            downloaded_bytes=0,
            message="existing file passed size and MD5 checks",
        )

    last_error = ""
    downloaded_bytes = 0
    for attempt in range(1, retries + 1):
        try:
            downloaded_bytes += transfer_once(spec, destination, timeout)
            return DownloadResult(
                spec=spec,
                status="downloaded",
                elapsed_seconds=time.perf_counter() - started,
                downloaded_bytes=downloaded_bytes,
                message=f"completed in {attempt} attempt(s)",
            )
        except (OSError, requests.RequestException, RuntimeError) as exc:
            last_error = str(exc)
            if attempt < retries:
                time.sleep(min(2 ** (attempt - 1), 30))

    return DownloadResult(
        spec=spec,
        status="failed",
        elapsed_seconds=time.perf_counter() - started,
        downloaded_bytes=downloaded_bytes,
        message=last_error,
    )


def write_manifest(path: Path, specs: list[DownloadSpec]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=MANIFEST_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(asdict(spec) for spec in specs)


def write_status(path: Path, results: list[DownloadResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=STATUS_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for result in results:
            row = asdict(result.spec)
            row.update(
                status=result.status,
                elapsed_seconds=f"{result.elapsed_seconds:.6f}",
                downloaded_bytes=result.downloaded_bytes,
                message=result.message,
            )
            writer.writerow(row)


def write_summary(
    path: Path,
    accession: str,
    started_at: str,
    ended_at: str,
    elapsed_seconds: float,
    specs: list[DownloadSpec],
    results: list[DownloadResult],
    status: str,
) -> None:
    counts = {
        name: sum(result.status == name for result in results)
        for name in ("downloaded", "skipped", "failed")
    }
    columns = [
        "project_accession",
        "started_utc",
        "ended_utc",
        "elapsed_seconds",
        "elapsed_hms",
        "selected_runs",
        "downloaded_runs",
        "skipped_runs",
        "failed_runs",
        "status",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerow(
            {
                "project_accession": accession,
                "started_utc": started_at,
                "ended_utc": ended_at,
                "elapsed_seconds": f"{elapsed_seconds:.6f}",
                "elapsed_hms": format_elapsed(elapsed_seconds),
                "selected_runs": len(specs),
                "downloaded_runs": counts["downloaded"],
                "skipped_runs": counts["skipped"],
                "failed_runs": counts["failed"],
                "status": status,
            }
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Download only forward-read archive-generated FASTQs from AMPLICON "
            "runs in an ENA study/project."
        )
    )
    parser.add_argument("accession", help="ENA study/project accession")
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Destination (default: data/fastq_data/<accession>/full)",
    )
    parser.add_argument(
        "--workers", type=int, default=8, help="Concurrent downloads (default: 8)"
    )
    parser.add_argument(
        "--retries", type=int, default=3, help="Attempts per FASTQ (default: 3)"
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=60.0,
        help="Connection/read timeout in seconds (default: 60)",
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="Write the manifest without downloading"
    )
    return parser


def run(args: argparse.Namespace) -> int:
    accession = args.accession.strip().upper()
    if not accession or not accession.isalnum():
        raise ValueError(f"Invalid ENA accession: {args.accession!r}")
    if args.workers <= 0 or args.retries <= 0 or args.timeout <= 0:
        raise ValueError("--workers, --retries, and --timeout must be positive")

    output_dir = args.output_dir or Path("data/fastq_data") / accession / "full"
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "ena_amplicon_forward_manifest.tsv"
    status_path = output_dir / "ena_amplicon_forward_status.tsv"
    summary_path = output_dir / "ena_amplicon_forward_summary.tsv"

    started_at = utc_now()
    started = time.perf_counter()
    specs: list[DownloadSpec] = []
    results: list[DownloadResult] = []
    final_status = "failed"
    try:
        print(f"Querying ENA AMPLICON runs for {accession}...", file=sys.stderr)
        specs, skipped = fetch_manifest(accession, args.timeout)
        write_manifest(manifest_path, specs)
        for message in skipped:
            print(f"Skipping {message}", file=sys.stderr)
        if not specs:
            raise RuntimeError(f"No downloadable forward AMPLICON FASTQs found for {accession}")

        total_bytes = sum(spec.bytes for spec in specs)
        print(
            f"Selected {len(specs)} forward FASTQ(s), {total_bytes / 10**9:.2f} GB total.",
            file=sys.stderr,
        )
        print(f"Manifest: {manifest_path}", file=sys.stderr)

        if args.dry_run:
            final_status = "dry-run"
            return 0

        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = {
                executor.submit(
                    download_one, spec, output_dir, args.retries, args.timeout
                ): spec
                for spec in specs
            }
            result_by_run: dict[str, DownloadResult] = {}
            for future in as_completed(futures):
                result = future.result()
                result_by_run[result.spec.run_accession] = result
                print(
                    f"[{result.status}] {result.spec.run_accession} "
                    f"{format_elapsed(result.elapsed_seconds)} {result.message}",
                    file=sys.stderr,
                    flush=True,
                )

        results = [result_by_run[spec.run_accession] for spec in specs]
        write_status(status_path, results)
        failed = [result for result in results if result.status == "failed"]
        final_status = "failed" if failed else "completed"
        return 1 if failed else 0
    finally:
        elapsed = time.perf_counter() - started
        ended_at = utc_now()
        write_summary(
            summary_path,
            accession,
            started_at,
            ended_at,
            elapsed,
            specs,
            results,
            final_status,
        )
        print(
            f"Accession {accession} finished in {format_elapsed(elapsed)} "
            f"({elapsed:.2f} seconds).",
            file=sys.stderr,
        )
        print(f"Summary: {summary_path}", file=sys.stderr)


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    try:
        return run(args)
    except (OSError, ValueError, RuntimeError, requests.RequestException) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
