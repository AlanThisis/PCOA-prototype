#!/usr/bin/env python3
"""Download forward-read FASTQs from ENA AMPLICON runs in a study."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
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
DEFAULT_RECOVERY_ROUNDS = 4
DEFAULT_MAX_RUNTIME_SECONDS = 2 * 60 * 60
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


@dataclass(frozen=True)
class ValidationReport:
    valid_runs: frozenset[str]
    problems_by_run: dict[str, str]
    general_problems: tuple[str, ...] = ()

    @property
    def complete(self) -> bool:
        return not self.problems_by_run and not self.general_problems


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
        if selected is None:
            # ENA documents an unsuffixed archive FASTQ for paired-layout runs
            # as containing unpaired reads. Retain it only when no mate FASTQs
            # are present; when _1 exists it remains the forward-read choice.
            names = {Path(urlparse(item[0]).path).name for item in candidates}
            unsuffixed_name = f"{run_accession}.fastq.gz"
            mate_names = {
                f"{run_accession}_1.fastq.gz",
                f"{run_accession}_2.fastq.gz",
            }
            if unsuffixed_name in names and names.isdisjoint(mate_names):
                selected = next(
                    item
                    for item in candidates
                    if Path(urlparse(item[0]).path).name == unsuffixed_name
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


def fetch_manifest(
    accession: str, timeout: float
) -> tuple[list[DownloadSpec], list[tuple[str, str]]]:
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
            if not records:
                raise RuntimeError("ENA returned an empty file report")
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
    skipped: list[tuple[str, str]] = []
    for record in records:
        if str(record.get("library_strategy", "")).strip().upper() != "AMPLICON":
            continue
        spec = select_forward_download(accession, record)
        if spec is None:
            run_accession = str(record.get("run_accession", "unknown"))
            skipped.append(
                (
                    run_accession,
                    "no archive-generated forward FASTQ for its layout",
                )
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


def validate_downloads(
    specs: list[DownloadSpec], output_dir: Path, *, clean_stale_partials: bool = True
) -> ValidationReport:
    """Validate every selected FASTQ against the ENA manifest."""
    valid_runs: set[str] = set()
    problems_by_run: dict[str, str] = {}
    expected_parts: set[Path] = set()

    for spec in specs:
        destination = output_dir / spec.output_filename
        part_path = destination.with_name(destination.name + ".part")
        expected_parts.add(part_path)
        if not destination.is_file():
            problems_by_run[spec.run_accession] = f"missing {destination.name}"
            continue
        actual_size = destination.stat().st_size
        if actual_size != spec.bytes:
            problems_by_run[spec.run_accession] = (
                f"size mismatch: expected {spec.bytes}, got {actual_size}"
            )
            continue
        actual_md5 = file_md5(destination)
        if actual_md5 != spec.md5:
            problems_by_run[spec.run_accession] = (
                f"MD5 mismatch: expected {spec.md5}, got {actual_md5}"
            )
            continue
        valid_runs.add(spec.run_accession)
        if clean_stale_partials:
            part_path.unlink(missing_ok=True)

    unexpected_parts = [
        path for path in sorted(output_dir.glob("*.part")) if path not in expected_parts
    ]
    if clean_stale_partials:
        quarantine_dir = output_dir / ".ena_orphaned_partials"
        for path in unexpected_parts:
            quarantine_dir.mkdir(exist_ok=True)
            target = quarantine_dir / path.name
            suffix = 1
            while target.exists():
                target = quarantine_dir / f"{path.name}.{suffix}"
                suffix += 1
            path.replace(target)
        general_problems: tuple[str, ...] = ()
    else:
        general_problems = tuple(
            f"unexpected partial file: {path.name}" for path in unexpected_parts
        )
    return ValidationReport(
        valid_runs=frozenset(valid_runs),
        problems_by_run=problems_by_run,
        general_problems=general_problems,
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
    urls = [https, ftp_url(https)]
    errors: list[str] = []
    for url in urls:
        attempts = [True, False] if part_path.exists() else [False]
        for resume in attempts:
            if not resume:
                part_path.unlink(missing_ok=True)
            existing_bytes = part_path.stat().st_size if resume else 0
            command = [
                curl,
                "--fail",
                "--location",
                "--silent",
                "--show-error",
                "--connect-timeout",
                str(math.ceil(timeout)),
                "--speed-limit",
                "1",
                "--speed-time",
                str(math.ceil(timeout)),
                "--retry",
                "2",
                "--retry-delay",
                "2",
                "--output",
                str(part_path),
            ]
            if resume:
                command.extend(["--continue-at", "-"])
            command.append(url)
            result = subprocess.run(
                command,
                check=False,
                capture_output=True,
                text=True,
            )
            if result.returncode != 0:
                errors.append(
                    f"{url}: curl exit {result.returncode}: {result.stderr.strip()}"
                )
                continue
            try:
                verify_download(spec, part_path, destination)
                return spec.bytes - existing_bytes
            except (OSError, RuntimeError) as exc:
                errors.append(f"{url}: {exc}")
                if part_path.exists() and part_path.stat().st_size < spec.bytes:
                    continue
                part_path.unlink(missing_ok=True)
                break

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
        allow_redirects=True,
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
    except (requests.RequestException, RuntimeError) as https_error:
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
    part_path = destination.with_name(destination.name + ".part")
    if completed_file_is_valid(destination, spec):
        part_path.unlink(missing_ok=True)
        return DownloadResult(
            spec=spec,
            status="skipped",
            elapsed_seconds=time.perf_counter() - started,
            downloaded_bytes=0,
            message="existing file passed size and MD5 checks",
        )
    if part_path.is_file() and part_path.stat().st_size == spec.bytes:
        try:
            verify_download(spec, part_path, destination)
            return DownloadResult(
                spec=spec,
                status="downloaded",
                elapsed_seconds=time.perf_counter() - started,
                downloaded_bytes=0,
                message="promoted a complete partial after size and MD5 checks",
            )
        except (OSError, RuntimeError):
            # verify_download removes a full-length partial with the wrong MD5.
            pass

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


def write_unavailable(
    path: Path, accession: str, unavailable: list[tuple[str, str]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["project_accession", "run_accession", "reason"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for run_accession, reason in unavailable:
            writer.writerow(
                {
                    "project_accession": accession,
                    "run_accession": run_accession,
                    "reason": reason,
                }
            )


def write_summary(
    path: Path,
    accession: str,
    started_at: str,
    ended_at: str,
    elapsed_seconds: float,
    specs: list[DownloadSpec],
    results: list[DownloadResult],
    status: str,
    recovery_rounds: int,
    validation_problems: int,
    unavailable_runs: int,
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
        "unavailable_runs",
        "recovery_rounds",
        "validation_problems",
        "status",
    ]
    row: dict[str, object] = {
        "project_accession": accession,
        "started_utc": started_at,
        "ended_utc": ended_at,
        "elapsed_seconds": f"{elapsed_seconds:.6f}",
        "elapsed_hms": format_elapsed(elapsed_seconds),
        "selected_runs": len(specs),
        "downloaded_runs": counts["downloaded"],
        "skipped_runs": counts["skipped"],
        "failed_runs": counts["failed"],
        "unavailable_runs": unavailable_runs,
        "recovery_rounds": recovery_rounds,
        "validation_problems": validation_problems,
        "status": status,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    previous_rows: list[dict[str, str]] = []
    if path.is_file():
        with path.open(newline="", encoding="utf-8") as handle:
            previous_rows = list(csv.DictReader(handle, delimiter="\t"))

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerow(row)

    history_path = path.with_name("ena_amplicon_forward_history.tsv")
    history_rows: list[dict[str, str]] = []
    if history_path.is_file():
        with history_path.open(newline="", encoding="utf-8") as handle:
            history_rows = list(csv.DictReader(handle, delimiter="\t"))
    known_invocations = {
        (existing.get("started_utc", ""), existing.get("ended_utc", ""))
        for existing in history_rows
    }
    for previous in previous_rows:
        key = (previous.get("started_utc", ""), previous.get("ended_utc", ""))
        if key not in known_invocations:
            history_rows.append(previous)
            known_invocations.add(key)
    history_rows.append({name: str(row.get(name, "")) for name in columns})
    with history_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for history_row in history_rows:
            writer.writerow({name: history_row.get(name, "") for name in columns})


def write_completion_marker(
    path: Path,
    accession: str,
    specs: list[DownloadSpec],
    elapsed_seconds: float,
    unavailable_runs: int,
) -> None:
    payload = {
        "project_accession": accession,
        "completed_utc": utc_now(),
        "selected_fastqs": len(specs),
        "selected_bytes": sum(spec.bytes for spec in specs),
        "unavailable_amplicon_runs": unavailable_runs,
        "elapsed_seconds_this_invocation": round(elapsed_seconds, 6),
        "validation": "all files passed ENA byte-count and MD5 checks",
    }
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


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
    parser.add_argument(
        "--recovery-rounds",
        type=int,
        default=DEFAULT_RECOVERY_ROUNDS,
        help=(
            "Full passes over unresolved FASTQs before exiting incomplete "
            f"(default: {DEFAULT_RECOVERY_ROUNDS})"
        ),
    )
    parser.add_argument(
        "--retry-until-complete",
        action="store_true",
        help="Keep making recovery passes until complete or --max-runtime is reached",
    )
    parser.add_argument(
        "--max-runtime",
        type=float,
        default=DEFAULT_MAX_RUNTIME_SECONDS,
        help="Soft runtime limit in seconds (default: 7200; 0 disables the limit)",
    )
    parser.add_argument(
        "--allow-missing-forward-fastq",
        action="store_true",
        help=(
            "Permit completion when ENA lists AMPLICON runs without a downloadable "
            "forward/unpaired archive FASTQ"
        ),
    )
    return parser


def run(args: argparse.Namespace) -> int:
    accession = args.accession.strip().upper()
    if not accession or not accession.isalnum():
        raise ValueError(f"Invalid ENA accession: {args.accession!r}")
    if args.workers <= 0 or args.retries <= 0 or args.timeout <= 0:
        raise ValueError("--workers, --retries, and --timeout must be positive")
    if args.recovery_rounds <= 0:
        raise ValueError("--recovery-rounds must be positive")
    if args.max_runtime < 0:
        raise ValueError("--max-runtime cannot be negative")

    output_dir = args.output_dir or Path("data/fastq_data") / accession / "full"
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "ena_amplicon_forward_manifest.tsv"
    status_path = output_dir / "ena_amplicon_forward_status.tsv"
    summary_path = output_dir / "ena_amplicon_forward_summary.tsv"
    unavailable_path = output_dir / "ena_amplicon_forward_unavailable.tsv"
    completion_path = output_dir / ".ena_download_complete.json"

    started_at = utc_now()
    started = time.perf_counter()
    specs: list[DownloadSpec] = []
    results: list[DownloadResult] = []
    result_by_run: dict[str, DownloadResult] = {}
    validation = ValidationReport(frozenset(), {})
    recovery_rounds = 0
    unavailable: list[tuple[str, str]] = []
    final_status = "failed"
    try:
        print(f"Querying ENA AMPLICON runs for {accession}...", file=sys.stderr)
        specs, unavailable = fetch_manifest(accession, args.timeout)
        write_manifest(manifest_path, specs)
        write_unavailable(unavailable_path, accession, unavailable)
        for run_accession, reason in unavailable:
            print(f"Unavailable {run_accession}: {reason}", file=sys.stderr)
        if not args.dry_run:
            completion_path.unlink(missing_ok=True)
        if not specs:
            raise RuntimeError(
                f"No downloadable forward AMPLICON FASTQs found for {accession}"
            )

        total_bytes = sum(spec.bytes for spec in specs)
        print(
            f"Selected {len(specs)} forward FASTQ(s), {total_bytes / 10**9:.2f} GB total.",
            file=sys.stderr,
        )
        print(f"Manifest: {manifest_path}", file=sys.stderr)

        if args.dry_run:
            final_status = "dry-run"
            return 0

        validation = validate_downloads(specs, output_dir)
        initially_valid = validation.valid_runs
        selection_complete = not unavailable or args.allow_missing_forward_fastq
        if validation.complete:
            results = [
                DownloadResult(
                    spec=spec,
                    status="skipped",
                    elapsed_seconds=0.0,
                    downloaded_bytes=0,
                    message="existing file passed size and MD5 checks",
                )
                for spec in specs
            ]
            write_status(status_path, results)
        else:
            print(
                f"Initial validation: {len(initially_valid)}/{len(specs)} FASTQ(s) "
                "complete.",
                file=sys.stderr,
            )
            results = []
            for spec in specs:
                if spec.run_accession in validation.valid_runs:
                    result = DownloadResult(
                        spec=spec,
                        status="skipped",
                        elapsed_seconds=0.0,
                        downloaded_bytes=0,
                        message="existing file passed size and MD5 checks",
                    )
                else:
                    result = DownloadResult(
                        spec=spec,
                        status="failed",
                        elapsed_seconds=0.0,
                        downloaded_bytes=0,
                        message=validation.problems_by_run.get(
                            spec.run_accession, "initial validation failed"
                        ),
                    )
                result_by_run[spec.run_accession] = result
                results.append(result)
            write_status(status_path, results)

        while not validation.complete:
            elapsed = time.perf_counter() - started
            if args.max_runtime and elapsed >= args.max_runtime:
                print(
                    f"Stopping recovery after reaching the {format_elapsed(args.max_runtime)} "
                    "soft runtime limit.",
                    file=sys.stderr,
                )
                break
            if not args.retry_until_complete and recovery_rounds >= args.recovery_rounds:
                break

            recovery_rounds += 1
            unresolved = [
                spec
                for spec in specs
                if spec.run_accession not in validation.valid_runs
            ]
            if recovery_rounds > 1:
                delay = min(15 * (2 ** (recovery_rounds - 2)), 300)
                if args.max_runtime:
                    remaining = args.max_runtime - (time.perf_counter() - started)
                    if remaining <= 0:
                        break
                    delay = min(delay, max(0, remaining))
                print(
                    f"Recovery round {recovery_rounds}: retrying {len(unresolved)} "
                    f"FASTQ(s) after {delay:.0f}s...",
                    file=sys.stderr,
                    flush=True,
                )
                time.sleep(delay)
                if (
                    args.max_runtime
                    and time.perf_counter() - started >= args.max_runtime
                ):
                    break
            else:
                print(
                    f"Download round {recovery_rounds}: {len(unresolved)} FASTQ(s) "
                    "need transfer or repair.",
                    file=sys.stderr,
                )

            with ThreadPoolExecutor(max_workers=args.workers) as executor:
                futures = {
                    executor.submit(
                        download_one, spec, output_dir, args.retries, args.timeout
                    ): spec
                    for spec in unresolved
                }
                for future in as_completed(futures):
                    spec = futures[future]
                    try:
                        result = future.result()
                    except Exception as exc:  # preserve the rest of the batch
                        result = DownloadResult(
                            spec=spec,
                            status="failed",
                            elapsed_seconds=0.0,
                            downloaded_bytes=0,
                            message=f"unexpected worker error: {exc}",
                        )
                    result_by_run[result.spec.run_accession] = result
                    print(
                        f"[{result.status}] {result.spec.run_accession} "
                        f"{format_elapsed(result.elapsed_seconds)} {result.message}",
                        file=sys.stderr,
                        flush=True,
                    )

            validation = validate_downloads(specs, output_dir)
            for spec in specs:
                if spec.run_accession in validation.valid_runs:
                    result_by_run.setdefault(
                        spec.run_accession,
                        DownloadResult(
                            spec=spec,
                            status="skipped",
                            elapsed_seconds=0.0,
                            downloaded_bytes=0,
                            message="existing file passed size and MD5 checks",
                        ),
                    )
                else:
                    previous = result_by_run.get(spec.run_accession)
                    validation_message = validation.problems_by_run.get(
                        spec.run_accession, "final validation failed"
                    )
                    if previous and previous.message:
                        message = (
                            f"{previous.message}; validation: {validation_message}"
                        )
                    else:
                        message = validation_message
                    result_by_run[spec.run_accession] = DownloadResult(
                        spec=spec,
                        status="failed",
                        elapsed_seconds=previous.elapsed_seconds if previous else 0.0,
                        downloaded_bytes=previous.downloaded_bytes if previous else 0,
                        message=message,
                    )
            results = [result_by_run[spec.run_accession] for spec in specs]
            write_status(status_path, results)
            print(
                f"Validation after round {recovery_rounds}: "
                f"{len(validation.valid_runs)}/{len(specs)} FASTQ(s) complete.",
                file=sys.stderr,
                flush=True,
            )

        if validation.complete and selection_complete:
            final_status = "completed"
            write_completion_marker(
                completion_path,
                accession,
                specs,
                time.perf_counter() - started,
                len(unavailable),
            )
            print(
                f"VALIDATION PASSED: all {len(specs)} selected FASTQ(s) passed "
                "size and MD5 checks.",
                file=sys.stderr,
            )
            print(f"Completion marker: {completion_path}", file=sys.stderr)
            return 0

        final_status = "incomplete"
        problem_count = len(validation.problems_by_run) + len(
            validation.general_problems
        ) + (0 if args.allow_missing_forward_fastq else len(unavailable))
        print(
            f"INCOMPLETE: {len(validation.valid_runs)}/{len(specs)} FASTQ(s) "
            f"validated; {problem_count} problem(s) remain.",
            file=sys.stderr,
        )
        for run_accession, message in validation.problems_by_run.items():
            print(f"  - {run_accession}: {message}", file=sys.stderr)
        for message in validation.general_problems:
            print(f"  - {message}", file=sys.stderr)
        if unavailable and not args.allow_missing_forward_fastq:
            print(
                f"  - {len(unavailable)} AMPLICON run(s) lack a downloadable "
                f"forward/unpaired archive FASTQ; see {unavailable_path}",
                file=sys.stderr,
            )
        print(
            "Safe to rerun the same command; verified files will be skipped.",
            file=sys.stderr,
        )
        return 1
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
            recovery_rounds,
            len(validation.problems_by_run)
            + len(validation.general_problems)
            + (0 if args.allow_missing_forward_fastq else len(unavailable)),
            len(unavailable),
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
