#!/usr/bin/env python3
"""Prepare deduplicated MMC sample-label runs for ENA download and PCoA."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import TextIO


PROJECT_PREFIXES = ("PRJNA", "PRJEB", "PRJDB")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("samples_tsv", type=Path)
    parser.add_argument("studies_tsv", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--tier",
        action="append",
        default=[],
        help="Tier to retain; repeat as needed (default: 1)",
    )
    return parser.parse_args(argv)


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open(encoding="utf-8", newline="")


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Input not found or empty: {path}")
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        return fields, list(reader)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def canonical_key(row: dict[str, str]) -> tuple[int, str, str]:
    accession = row["study_accession"].strip()
    return (
        0 if accession.startswith(PROJECT_PREFIXES) else 1,
        row["record_id"].strip(),
        accession,
    )


def prepare_sample_labels(
    samples_tsv: Path,
    studies_tsv: Path,
    output_dir: Path,
    tiers: set[str],
) -> dict[str, object]:
    sample_fields, sample_rows = read_tsv(samples_tsv)
    study_fields, study_rows = read_tsv(studies_tsv)
    required_samples = {
        "record_id",
        "study_accession",
        "sample_accession",
        "tier",
        "run_accession",
        "library_strategy",
        "instrument_platform",
    }
    required_studies = {"record_id", "body_site", "disease_group", "title"}
    if missing := sorted(required_samples - set(sample_fields)):
        raise ValueError(f"Sample table missing columns: {', '.join(missing)}")
    if missing := sorted(required_studies - set(study_fields)):
        raise ValueError(f"Study table missing columns: {', '.join(missing)}")

    studies = {row["record_id"].strip(): row for row in study_rows}
    candidates = [
        row
        for row in sample_rows
        if row["tier"].strip() in tiers
        and row["library_strategy"].strip().upper() == "AMPLICON"
        and row["run_accession"].strip()
    ]
    missing_studies = sorted(
        {row["record_id"] for row in candidates if row["record_id"] not in studies}
    )
    if missing_studies:
        raise ValueError(f"Sample rows reference missing studies: {missing_studies[:10]}")

    by_run: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in candidates:
        by_run[row["run_accession"].strip()].append(row)

    selected: list[dict[str, str]] = []
    duplicate_rows: list[dict[str, str]] = []
    conflict_counts: Counter[str] = Counter()
    conflict_fields = (
        "study_accession",
        "record_id",
        "case_control",
        "mondo_label",
    )
    for run_accession in sorted(by_run):
        alternatives = sorted(by_run[run_accession], key=canonical_key)
        chosen = alternatives[0]
        study = studies[chosen["record_id"]]
        combined = dict(chosen)
        combined.update(
            {
                "sample-id": run_accession,
                "source_study": chosen["study_accession"].strip(),
                "source_record_id": chosen["record_id"].strip(),
                "environment_harmonized": study["body_site"].strip(),
                "study_level_disease_group": study["disease_group"].strip(),
                "study_title": study["title"].strip(),
                "citation_count": str(len(alternatives)),
            }
        )
        selected.append(combined)
        for field in conflict_fields:
            if len({row.get(field, "").strip() for row in alternatives}) > 1:
                conflict_counts[field] += 1
        if len(alternatives) > 1:
            for alternative in alternatives:
                duplicate_rows.append(
                    {
                        "run_accession": run_accession,
                        "selected": "yes" if alternative is chosen else "no",
                        "record_id": alternative["record_id"],
                        "study_accession": alternative["study_accession"],
                        "sample_accession": alternative["sample_accession"],
                        "case_control": alternative.get("case_control", ""),
                        "mondo_label": alternative.get("mondo_label", ""),
                    }
                )

    output_dir.mkdir(parents=True, exist_ok=True)
    metadata_fields = [
        "sample-id",
        "run_accession",
        "sample_accession",
        "source_study",
        "source_record_id",
        "tier",
        "instrument_platform",
        "environment_harmonized",
        "case_control",
        "mondo_id",
        "mondo_label",
        "study_level_disease_group",
        "host_age",
        "host_sex",
        "mapping_method",
        "study_title",
        "citation_count",
    ]
    write_tsv(output_dir / "all_amplicon_metadata.tsv", metadata_fields, selected)
    illumina = [
        row for row in selected if row["instrument_platform"].strip() == "ILLUMINA"
    ]
    write_tsv(output_dir / "illumina_pipeline_metadata.tsv", metadata_fields, illumina)

    allowlist_fields = ["project_accession", "run_accession", "sample_id"]
    allowlist = [
        {
            "project_accession": row["source_study"],
            "run_accession": row["run_accession"],
            "sample_id": row["sample-id"],
        }
        for row in selected
    ]
    write_tsv(output_dir / "approved_runs.tsv", allowlist_fields, allowlist)
    illumina_allowlist = [
        row
        for row, metadata in zip(allowlist, selected)
        if metadata["instrument_platform"].strip() == "ILLUMINA"
    ]
    write_tsv(
        output_dir / "illumina_approved_runs.tsv",
        allowlist_fields,
        illumina_allowlist,
    )
    write_tsv(
        output_dir / "duplicate_citations.tsv",
        [
            "run_accession",
            "selected",
            "record_id",
            "study_accession",
            "sample_accession",
            "case_control",
            "mondo_label",
        ],
        duplicate_rows,
    )
    non_illumina = [
        {
            "run_accession": row["run_accession"],
            "sample_accession": row["sample_accession"],
            "project_accession": row["source_study"],
            "instrument_platform": row["instrument_platform"],
        }
        for row in selected
        if row["instrument_platform"].strip() != "ILLUMINA"
    ]
    write_tsv(
        output_dir / "non_illumina_runs.tsv",
        [
            "run_accession",
            "sample_accession",
            "project_accession",
            "instrument_platform",
        ],
        non_illumina,
    )
    projects = sorted({row["source_study"] for row in selected})
    (output_dir / "projects.txt").write_text(
        "".join(f"{project}\n" for project in projects), encoding="utf-8"
    )

    summary: dict[str, object] = {
        "tiers": sorted(tiers),
        "source_sample_rows": len(sample_rows),
        "candidate_amplicon_rows": len(candidates),
        "unique_amplicon_runs": len(selected),
        "illumina_runs": len(illumina),
        "non_illumina_runs": len(non_illumina),
        "project_accessions": len(projects),
        "runs_with_multiple_citations": sum(len(rows) > 1 for rows in by_run.values()),
        "duplicate_citation_rows": len(candidates) - len(selected),
        "conflicting_duplicate_groups": dict(sorted(conflict_counts.items())),
    }
    (output_dir / "preparation_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    tiers = set(args.tier or ["1"])
    summary = prepare_sample_labels(
        args.samples_tsv, args.studies_tsv, args.output_dir, tiers
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
