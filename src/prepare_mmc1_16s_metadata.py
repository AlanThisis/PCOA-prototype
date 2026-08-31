#!/usr/bin/env python3
"""Build human-only, run-resolved MMC1 metadata and an ENA run allowlist."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


PROJECT_EXCLUSIONS = {
    "PRJNA339931": "study cohort collapsed into one umbrella BioSample",
    "PRJNA564314": (
        "multiple respiratory sites collapsed into one BioSample and Oxford "
        "Nanopore reads are incompatible with the current Deblur workflow"
    ),
    "PRJNA985881": (
        "infant feces, maternal feces, and breast-milk material collapsed into "
        "one BioSample with an incorrect uniform body-site label"
    ),
    "PRJNA1212986": (
        "PacBio subreads are incompatible with the current fixed-length Deblur "
        "workflow and are not ENA archive FASTQs with supported names"
    ),
    "PRJNA381268": "ENA FASTQs are protected by dbGaP and cannot be downloaded",
    "PRJEB75547": (
        "Oxford Nanopore GridION full-length amplicon reads use quality scores "
        "that are incompatible with the current short-read Deblur workflow"
    ),
}
RUN_EXCLUSIONS = {
    "SRR10377565": (
        "ENA advertises an archive FASTQ but repeatedly serves an HTML error "
        "document instead of the file"
    ),
}
LAB_CONTROL_PATTERN = re.compile(
    r"(?:pcr\s*neg|pcrneg|negative\s*control|extraction\s*negative|"
    r"extr(?:action)?nc|nc[._-]?pcr|(?:^|[_-])blank(?:[_-]|$))",
    re.IGNORECASE,
)
NON_16S_PATTERN = re.compile(r"(?:\bITS\b|\b18S\b|fungal)", re.IGNORECASE)
RUN_FIELDS = (
    "experiment_accession",
    "experiment_alias",
    "experiment_title",
    "run_alias",
    "library_name",
    "read_count",
    "base_count",
)


@dataclass(frozen=True)
class Exclusion:
    project_accession: str
    sample_accession: str
    run_accession: str
    reason: str


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv", type=Path)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/mmc1_16s_cleaning"),
    )
    return parser.parse_args(argv)


def make_unique_headers(headers: list[str]) -> list[str]:
    seen: Counter[str] = Counter()
    unique: list[str] = []
    for header in headers:
        seen[header] += 1
        if seen[header] == 1:
            unique.append(header)
        elif header == "study_title" and seen[header] == 2:
            unique.append("ena_study_title")
        else:
            unique.append(f"{header}_{seen[header]}")
    return unique


def split_values(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def value_for_run(row: dict[str, str], field: str, index: int, count: int) -> str:
    values = split_values(row.get(field, ""))
    if len(values) == count:
        return values[index]
    if len(values) == 1:
        return values[0]
    return row.get(field, "").strip()


def run_text(row: dict[str, str], index: int, count: int) -> str:
    return " ".join(value_for_run(row, field, index, count) for field in RUN_FIELDS)


def exclusion_reason(
    row: dict[str, str],
    run_index: int,
    run_count: int,
    project_exclusions: dict[str, str] | None = None,
) -> str | None:
    project = row["source_study"].strip()
    project_exclusions = project_exclusions or PROJECT_EXCLUSIONS
    if project in project_exclusions:
        return project_exclusions[project]
    run_accession = split_values(row.get("run_accession", ""))[run_index]
    if run_accession in RUN_EXCLUSIONS:
        return RUN_EXCLUSIONS[run_accession]
    if row.get("study_hmr", "").strip() != "Human":
        return "study is not exclusively human"
    if row.get("host_species_resolved", "").strip() != "Human":
        return "specimen is not human"

    text = run_text(row, run_index, run_count)
    if LAB_CONTROL_PATTERN.search(text):
        return "explicit laboratory negative/blank control"
    if NON_16S_PATTERN.search(text):
        return "experiment is explicitly ITS/18S/fungal rather than bacterial 16S"
    if project == "PRJNA597168" and not re.search(
        r"\bforward\b|(?:^|[_-])\d+F(?:\.|[_-]|$)", text, re.IGNORECASE
    ):
        return "separately deposited reverse-direction run"
    return None


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        raw_headers = next(reader)
        headers = make_unique_headers(raw_headers)
        rows = []
        for line_number, values in enumerate(reader, start=2):
            if len(values) != len(headers):
                raise ValueError(
                    f"Malformed CSV row {line_number}: expected {len(headers)} "
                    f"columns, found {len(values)}"
                )
            rows.append(dict(zip(headers, values)))
    return headers, rows


def write_csv(path: Path, headers: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)


def prepare_metadata(input_csv: Path, output_dir: Path) -> dict[str, int]:
    input_csv = input_csv.resolve()
    if not input_csv.is_file() or input_csv.stat().st_size == 0:
        raise FileNotFoundError(f"Input metadata not found or empty: {input_csv}")
    output_dir.mkdir(parents=True, exist_ok=True)

    headers, source_rows = read_rows(input_csv)
    if "run_accession" not in headers or "source_study" not in headers:
        raise ValueError("Input must contain run_accession and source_study columns")

    output_headers = ["sample-id", *headers]
    cleaned: list[dict[str, str]] = []
    exclusions: list[Exclusion] = []
    project_study_labels: dict[str, set[str]] = {}
    for row in source_rows:
        project_study_labels.setdefault(row["source_study"].strip(), set()).add(
            row.get("study_hmr", "").strip()
        )
    project_level_exclusions = dict(PROJECT_EXCLUSIONS)
    for project, labels in project_study_labels.items():
        if labels != {"Human"}:
            project_level_exclusions.setdefault(
                project, "study is not exclusively human"
            )

    for row in source_rows:
        project = row["source_study"].strip()
        sample = row.get("sample_accession", "").strip()
        runs = split_values(row["run_accession"])
        if not runs:
            exclusions.append(Exclusion(project, sample, "", "missing run accession"))
            continue
        for index, run_accession in enumerate(runs):
            reason = exclusion_reason(
                row, index, len(runs), project_level_exclusions
            )
            if reason:
                exclusions.append(Exclusion(project, sample, run_accession, reason))
                continue

            resolved = dict(row)
            resolved["run_accession"] = run_accession
            for field in RUN_FIELDS:
                if field in resolved:
                    resolved[field] = value_for_run(row, field, index, len(runs))
            resolved["sample-id"] = run_accession
            cleaned.append(resolved)

    sample_ids = [row["sample-id"] for row in cleaned]
    if len(sample_ids) != len(set(sample_ids)):
        duplicates = sorted(key for key, n in Counter(sample_ids).items() if n > 1)
        raise ValueError(f"Duplicate selected run/sample IDs: {duplicates[:10]}")
    if any(row["study_hmr"] != "Human" for row in cleaned):
        raise ValueError("Non-human study survived cleanup")
    if any(row["host_species_resolved"] != "Human" for row in cleaned):
        raise ValueError("Non-human specimen survived cleanup")

    project_counts = Counter(row["source_study"] for row in cleaned)
    if not project_counts:
        raise ValueError("Cleanup selected no projects")
    source_projects = set(project_study_labels)
    excluded_project_reasons = {
        project: project_level_exclusions.get(project, "all runs failed run-level filters")
        for project in source_projects - set(project_counts)
    }

    write_csv(output_dir / "cleaned_metadata.csv", output_headers, cleaned)
    with (output_dir / "approved_runs.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["project_accession", "run_accession", "sample_id"])
        for row in cleaned:
            writer.writerow(
                [row["source_study"], row["run_accession"], row["sample-id"]]
            )
    (output_dir / "projects.txt").write_text(
        "".join(f"{project}\n" for project in sorted(project_counts)),
        encoding="utf-8",
    )

    with (output_dir / "excluded_runs.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["project_accession", "sample_accession", "run_accession", "reason"]
        )
        for item in exclusions:
            writer.writerow(
                [
                    item.project_accession,
                    item.sample_accession,
                    item.run_accession,
                    item.reason,
                ]
            )
    with (output_dir / "excluded_projects.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["project_accession", "reason"])
        for project, reason in sorted(excluded_project_reasons.items()):
            writer.writerow([project, reason])

    summary = {
        "source_rows": len(source_rows),
        "selected_projects": len(project_counts),
        "selected_runs": len(cleaned),
        "excluded_run_records": len(exclusions),
        "excluded_projects": len(excluded_project_reasons),
        "explicitly_excluded_projects": sum(
            project in excluded_project_reasons for project in PROJECT_EXCLUSIONS
        ),
    }
    (output_dir / "cleanup_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    summary = prepare_metadata(args.input_csv, args.output_dir)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
