#!/usr/bin/env python3
"""Fetch ENA sample metadata for a study/project accession and write CSV.

Examples
--------
Print CSV to stdout:
    python src/get_ENA_metadata.py PRJEB46665

Write CSV to a file:
    python src/get_ENA_metadata.py PRJEB46665 -o samples.csv
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import time
import xml.etree.ElementTree as ET
from collections import OrderedDict

import requests


ENA_PORTAL_SEARCH_BASE = "https://www.ebi.ac.uk/ena/portal/api/search"
ENA_BROWSER_XML_BASE = "https://www.ebi.ac.uk/ena/browser/api/xml"
TOOL_NAME = "ena_metadata_harmonization"

# I did not find a clear public ENA RPS cap in the current docs, so use a
# conservative default for repeated XML fetches.
DEFAULT_REQUESTS_PER_SECOND = 2.0
MIN_REQUEST_INTERVAL_SECONDS = 1.0 / DEFAULT_REQUESTS_PER_SECOND
PROJECT_RE = re.compile(r"^[A-Z]{4,}\d+$", re.IGNORECASE)
XML_BATCH_SIZE = 200

CORE_COLUMNS = [
    "project_accession",
    "sample_accession",
    "sample_alias",
    "sample_title",
    "center_name",
    "primary_id",
    "secondary_id",
    "taxon_id",
    "scientific_name",
    "description",
    "error",
]


class ENAClient:
    def __init__(self, tool: str = TOOL_NAME) -> None:
        self.tool = tool
        self._last_request_time = 0.0

    def _throttle(self) -> None:
        elapsed = time.monotonic() - self._last_request_time
        if elapsed < MIN_REQUEST_INTERVAL_SECONDS:
            time.sleep(MIN_REQUEST_INTERVAL_SECONDS - elapsed)

    def _get(self, url: str, params: dict[str, str] | None = None) -> requests.Response:
        self._throttle()
        try:
            response = requests.get(
                url,
                params=params,
                headers={"User-Agent": self.tool},
                timeout=30,
            )
            response.raise_for_status()
        except requests.HTTPError as exc:
            detail = exc.response.text if exc.response is not None else ""
            raise RuntimeError(f"HTTP error for {response.url}: {detail}") from exc
        except requests.RequestException as exc:
            raise RuntimeError(f"Request failed for {url}: {exc}") from exc

        self._last_request_time = time.monotonic()
        return response

    def fetch_sample_accessions(self, project_accession: str) -> list[str]:
        response = self._get(
            ENA_PORTAL_SEARCH_BASE,
            {
                "result": "sample",
                "query": f'study_accession="{project_accession}"',
                "fields": "sample_accession",
                "format": "json",
                "limit": "0",
            },
        )
        records = response.json()
        sample_accessions = []
        for record in records:
            accession = str(record.get("sample_accession", "")).strip()
            if accession:
                sample_accessions.append(accession)
        return sample_accessions

    def fetch_samples_xml_batch(self, sample_accessions: list[str]) -> str:
        response = self._get(f"{ENA_BROWSER_XML_BASE}/{','.join(sample_accessions)}")
        return response.text


def normalize_column_name(raw_name: str) -> str:
    normalized = raw_name.strip().lower()
    normalized = re.sub(r"[^a-z0-9]+", "_", normalized)
    normalized = normalized.strip("_")
    return normalized or "unnamed_field"


def add_value(row: dict[str, str], column_name: str, value: str) -> None:
    clean_value = " ".join(value.split())
    if not clean_value:
        return
    existing = row.get(column_name, "")
    if not existing:
        row[column_name] = clean_value
        return
    existing_parts = existing.split("; ")
    if clean_value not in existing_parts:
        row[column_name] = f"{existing}; {clean_value}"


def unique_preserving_order(values: list[str]) -> list[str]:
    seen: set[str] = set()
    ordered = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return ordered


def assign_dynamic_column(
    raw_tag: str,
    value: str,
    row: dict[str, str],
    tag_to_column: OrderedDict[str, str],
    column_to_tag: dict[str, str],
) -> None:
    if not raw_tag.strip():
        return
    if raw_tag not in tag_to_column:
        base = normalize_column_name(raw_tag)
        candidate = base
        suffix = 2
        while candidate in column_to_tag and column_to_tag[candidate] != raw_tag:
            candidate = f"{base}_{suffix}"
            suffix += 1
        tag_to_column[raw_tag] = candidate
        column_to_tag[candidate] = raw_tag
    add_value(row, tag_to_column[raw_tag], value)


def parse_sample_element(
    project_accession: str,
    sample: ET.Element,
    fallback_accession: str,
    tag_to_column: OrderedDict[str, str],
    column_to_tag: dict[str, str],
) -> dict[str, str]:
    row: dict[str, str] = {
        "project_accession": project_accession,
        "sample_accession": sample.get("accession", fallback_accession),
        "sample_alias": sample.get("alias", ""),
        "center_name": sample.get("center_name", ""),
        "sample_title": sample.findtext("./TITLE", default="").strip(),
        "primary_id": sample.findtext("./IDENTIFIERS/PRIMARY_ID", default="").strip(),
        "secondary_id": sample.findtext("./IDENTIFIERS/SECONDARY_ID", default="").strip(),
        "taxon_id": sample.findtext("./SAMPLE_NAME/TAXON_ID", default="").strip(),
        "scientific_name": sample.findtext("./SAMPLE_NAME/SCIENTIFIC_NAME", default="").strip(),
        "description": sample.findtext("./DESCRIPTION", default="").strip(),
        "error": "",
    }
    for sample_attribute in sample.findall("./SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"):
        raw_tag = sample_attribute.findtext("./TAG", default="")
        value = sample_attribute.findtext("./VALUE", default="")
        assign_dynamic_column(raw_tag, value, row, tag_to_column, column_to_tag)
    return row


def chunked(values: list[str], size: int) -> list[list[str]]:
    return [values[i:i + size] for i in range(0, len(values), size)]


def fetch_rows(project_accession: str, max_samples: int | None = None) -> tuple[list[dict[str, str]], list[str]]:
    import sys
    client = ENAClient()
    print(f"Fetching samples for {project_accession}...", file=sys.stderr, flush=True)
    sample_accessions = unique_preserving_order(client.fetch_sample_accessions(project_accession))
    if not sample_accessions:
        raise RuntimeError(f"No sample accessions found for project {project_accession}")
    if max_samples is not None:
        sample_accessions = sample_accessions[:max_samples]

    total = len(sample_accessions)
    batches = chunked(sample_accessions, XML_BATCH_SIZE)
    print(f"  Found {total} sample(s), fetching in {len(batches)} batch(es).", file=sys.stderr, flush=True)

    rows: list[dict[str, str]] = []
    tag_to_column: OrderedDict[str, str] = OrderedDict()
    column_to_tag: dict[str, str] = {}
    fetched = 0

    for batch_i, batch in enumerate(batches, 1):
        start = (batch_i - 1) * XML_BATCH_SIZE + 1
        end = start + len(batch) - 1
        print(f"  Batch [{batch_i}/{len(batches)}] samples {start}–{end}...", file=sys.stderr, end="", flush=True)
        try:
            batch_xml = client.fetch_samples_xml_batch(batch)
            root = ET.fromstring(batch_xml)
            sample_elements = root.findall("./SAMPLE")
            by_accession = {el.get("accession", ""): el for el in sample_elements}
            print(f" got {len(sample_elements)}", file=sys.stderr, flush=True)
        except Exception as exc:
            print(f" error: {exc}", file=sys.stderr, flush=True)
            for acc in batch:
                rows.append({"project_accession": project_accession, "sample_accession": acc, "error": str(exc)})
            continue

        for acc in batch:
            sample_el = by_accession.get(acc)
            if sample_el is None:
                rows.append({"project_accession": project_accession, "sample_accession": acc, "error": "not returned in batch XML"})
            else:
                rows.append(parse_sample_element(project_accession, sample_el, acc, tag_to_column, column_to_tag))
        fetched += len(sample_elements)

    errors = sum(1 for r in rows if r.get("error"))
    suffix = f", {errors} error(s)" if errors else ""
    print(f"Done: {total - errors}/{total} samples fetched{suffix}", file=sys.stderr, flush=True)

    columns = CORE_COLUMNS + [column for column in tag_to_column.values() if column not in CORE_COLUMNS]
    return rows, columns


def write_csv(rows: list[dict[str, str]], columns: list[str], output_path: str | None) -> None:
    if output_path:
        with open(output_path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
            writer.writeheader()
            for row in rows:
                writer.writerow({column: row.get(column, "") for column in columns})
        return

    writer = csv.DictWriter(sys.stdout, fieldnames=columns, extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow({column: row.get(column, "") for column in columns})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Fetch ENA sample metadata for a project/study accession and flatten it into CSV.",
    )
    parser.add_argument(
        "project_accession",
        help="ENA project/study accession, for example PRJEB46665.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Write CSV to this path. If omitted, CSV is printed to stdout.",
    )
    parser.add_argument(
        "--max-samples",
        type=int,
        help="Only fetch the first N sample accessions from the project. Useful for quick tests on large studies.",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    project_accession = args.project_accession.strip().upper()

    if not PROJECT_RE.fullmatch(project_accession):
        print(f"Error: invalid project accession: {args.project_accession!r}", file=sys.stderr)
        return 1
    if args.max_samples is not None and args.max_samples <= 0:
        print("Error: --max-samples must be a positive integer", file=sys.stderr)
        return 1

    try:
        rows, columns = fetch_rows(project_accession, max_samples=args.max_samples)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    write_csv(rows, columns, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
