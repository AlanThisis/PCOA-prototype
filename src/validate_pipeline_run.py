#!/usr/bin/env python3
"""Validate that a pipeline run completed with internally consistent outputs."""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    return parser.parse_args(argv)


def read_json(path: Path) -> dict:
    if not path.is_file():
        raise ValueError(f"Missing required file: {path}")
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def validate_run(run_dir: Path) -> list[str]:
    run_dir = run_dir.resolve()
    problems: list[str] = []
    state = read_json(run_dir / "run_state.json")
    manifest = read_json(run_dir / "run_manifest.json")
    attempts = state.get("attempts", [])
    if not attempts or attempts[-1].get("status") not in {
        "completed",
        "completed_with_exclusions",
    }:
        problems.append("latest pipeline attempt is not complete")
    for name, stage in state.get("stages", {}).items():
        if stage.get("status") not in {"completed", "excluded_empty"}:
            problems.append(f"stage {name} has status {stage.get('status')}")

    results = run_dir / "results"
    for filename in (
        "pcoa_coordinates_unweighted_unifrac.txt",
        "pcoa_plot_unweighted_unifrac.png",
        "analysis_summary.json",
        "pipeline_summary.json",
    ):
        path = results / filename
        if not path.is_file() or path.stat().st_size == 0:
            problems.append(f"missing or empty result: {filename}")
    for column in manifest.get("color_by", []):
        safe_column = re.sub(r"[^A-Za-z0-9._-]+", "_", column).strip("._")
        plot = results / f"pcoa_{safe_column}.png"
        if not plot.is_file() or plot.stat().st_size == 0:
            problems.append(f"missing or empty colored plot: {plot.name}")

    sample_status = results / "sample_processing_status.tsv"
    deblur_summary_path = results / "deblur_processing_summary.json"
    if sample_status.is_file() or deblur_summary_path.is_file():
        if not sample_status.is_file() or not deblur_summary_path.is_file():
            problems.append("balanced Deblur status outputs are incomplete")
        else:
            with sample_status.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            summary = read_json(deblur_summary_path)
            counts = summary.get("counts", {})
            if len(rows) != counts.get("expected"):
                problems.append("sample status row count does not match expected count")
            observed_counts: dict[str, int] = {}
            for row in rows:
                observed_counts[row["status"]] = observed_counts.get(row["status"], 0) + 1
            for status in (
                "completed",
                "completed_after_sanitation",
                "zero_features",
                "failed",
                "excluded_assay",
            ):
                if observed_counts.get(status, 0) != counts.get(status, 0):
                    problems.append(f"sample count mismatch for status {status}")
            if counts.get("failed", 0) > summary.get("failed_sample_limit", 0):
                problems.append("failed sample count exceeds configured tolerance")

    analysis = read_json(results / "analysis_summary.json")
    if analysis.get("mapped_samples", 0) <= 0:
        problems.append("analysis reports no GG2-mapped samples")
    if analysis.get("rarefied_samples", 0) <= 0:
        problems.append("analysis reports no samples retained after rarefaction")
    distance = analysis.get("distance_tsv", {})
    # Older QIIME summaries used qza_path; current QIIME and DART summaries
    # record the native distance artifact under backend_path.
    backend_path = distance.get("backend_path") or distance.get("qza_path")
    distance_artifact = Path(backend_path) if backend_path else None
    if (
        distance_artifact is None
        or not distance_artifact.is_file()
        or distance_artifact.stat().st_size == 0
    ):
        problems.append("UniFrac distance artifact is missing or empty")
    if distance.get("exported"):
        distance_tsv = results / "distance_matrix_unweighted_unifrac.tsv"
        if not distance_tsv.is_file() or distance_tsv.stat().st_size == 0:
            problems.append("distance TSV was declared exported but is missing")
    return problems


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        problems = validate_run(args.run_dir)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"VALIDATION FAILED: {exc}")
        return 1
    if problems:
        print(f"VALIDATION FAILED: {len(problems)} problem(s)")
        for problem in problems:
            print(f"  - {problem}")
        return 1
    print(f"VALIDATION OK: {args.run_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
