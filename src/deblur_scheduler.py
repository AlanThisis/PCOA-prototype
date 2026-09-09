#!/usr/bin/env python3
"""Run workload-balanced, failure-isolating Deblur shards.

The input TSV is produced by ``run_pipeline.py`` and contains one FASTQ per
row. Initial shards run concurrently. A persistently failing shard is split by
estimated workload until the failure is isolated to a single FASTQ, allowing
all unaffected samples to finish without a manual pipeline restart.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import biom
import biom.util

from merge_biom import load_seqs
from pipeline_lib import (
    TimingRecorder,
    add_timing_argument,
    run_command,
    run_timed_main,
)


@dataclass(frozen=True)
class FastqEntry:
    sample_id: str
    study: str
    fastq_path: Path
    size_bytes: int
    weight_bytes: int
    shard_id: int
    excluded_assay: bool = False


@dataclass(frozen=True)
class LeafResult:
    node_id: str
    entries: tuple[FastqEntry, ...]
    workflow_dir: Path | None
    error: str | None = None


class InfrastructureFailure(RuntimeError):
    """A shared execution failure that sample-level bisection cannot repair."""


INFRASTRUCTURE_ERROR_MARKERS = (
    "cannot allocate memory",
    "disk quota exceeded",
    "input/output error",
    "modulenotfounderror",
    "no space left on device",
    "permission denied",
    "required executable not found",
)


def is_infrastructure_failure(exc: BaseException, log_text: str = "") -> bool:
    """Return whether retrying progressively smaller data shards is futile."""
    if isinstance(exc, OSError):
        return True
    if isinstance(exc, RuntimeError) and not isinstance(
        exc, subprocess.CalledProcessError
    ):
        return "without a complete workflow" in str(exc)
    return_code = getattr(exc, "returncode", None)
    if return_code in {126, 127, -9, -15, 137, 143}:
        return True
    combined = f"{exc}\n{log_text}".lower()
    return any(marker in combined for marker in INFRASTRUCTURE_ERROR_MARKERS)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shard-manifest", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--trim-length", type=int, default=150)
    parser.add_argument(
        "--error-dist",
        default="1,0.06,0.02,0.02,0.01,0.005,0.005,0.005,0.001,0.001,0.001,0.0005",
    )
    parser.add_argument("--min-reads", type=int, default=0)
    parser.add_argument("--shard-workers", type=int, default=2)
    parser.add_argument("--concurrent-shards", type=int, default=1)
    parser.add_argument("--max-failed-fraction", type=float, default=0.02)
    parser.add_argument("--max-failed-samples", type=int)
    parser.add_argument("--keep-tmp-files", action="store_true")
    add_timing_argument(parser)
    return parser.parse_args(argv)


def read_shard_manifest(path: Path) -> list[FastqEntry]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "sample_id",
        "study",
        "fastq_path",
        "size_bytes",
        "weight_bytes",
        "shard_id",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Invalid or empty shard manifest: {path}")
    entries = [
        FastqEntry(
            sample_id=row["sample_id"],
            study=row["study"],
            fastq_path=Path(row["fastq_path"]),
            size_bytes=int(row["size_bytes"]),
            weight_bytes=int(row["weight_bytes"]),
            shard_id=int(row["shard_id"]),
            excluded_assay=row.get("excluded_assay", "0") in {"1", "true", "True"},
        )
        for row in rows
    ]
    if len({entry.sample_id for entry in entries}) != len(entries):
        raise ValueError("Shard manifest contains duplicate sample IDs")
    return entries


def split_entries(
    entries: tuple[FastqEntry, ...],
) -> tuple[tuple[FastqEntry, ...], tuple[FastqEntry, ...]]:
    """Split a failed shard into workload-balanced, deterministic children."""
    left: list[FastqEntry] = []
    right: list[FastqEntry] = []
    totals = [0, 0]
    for entry in sorted(entries, key=lambda item: (-item.weight_bytes, item.sample_id)):
        target = 0 if totals[0] <= totals[1] else 1
        (left if target == 0 else right).append(entry)
        totals[target] += entry.weight_bytes
    return tuple(sorted(left, key=lambda item: item.sample_id)), tuple(
        sorted(right, key=lambda item: item.sample_id)
    )


def write_json_atomic(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def workflow_is_complete(workflow_dir: Path) -> bool:
    biom_path = workflow_dir / "all.biom"
    fasta_path = workflow_dir / "all.seqs.fa"
    return (
        biom_path.is_file()
        and biom_path.stat().st_size > 0
        and fasta_path.is_file()
    )


def write_fastq_manifest(path: Path, entries: tuple[FastqEntry, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{entry.fastq_path}\n" for entry in entries), encoding="utf-8")


def run_leaf(
    node_id: str,
    entries: tuple[FastqEntry, ...],
    args: argparse.Namespace,
    timing: TimingRecorder,
) -> list[LeafResult]:
    node_dir = args.work_dir / "nodes" / node_id
    workflow_dir = node_dir / "workflow"
    if workflow_is_complete(workflow_dir):
        timing.skipped("deblur_shard", item=node_id, message="valid cached workflow")
        return [LeafResult(node_id, entries, workflow_dir)]
    split_marker = node_dir / "split.json"
    if split_marker.is_file():
        split = json.loads(split_marker.read_text(encoding="utf-8"))
        current_ids = {entry.sample_id for entry in entries}
        if set(split["left"] + split["right"]) == current_ids:
            by_id = {entry.sample_id: entry for entry in entries}
            results: list[LeafResult] = []
            for suffix, sample_ids in (("L", split["left"]), ("R", split["right"])):
                child = tuple(by_id[sample_id] for sample_id in sample_ids)
                if child:
                    results.extend(run_node(f"{node_id}-{suffix}", child, args, timing))
            return results

    input_manifest = node_dir / "fastqs.txt"
    write_fastq_manifest(input_manifest, entries)
    command = [
        sys.executable,
        str(Path(__file__).with_name("run_deblur.py")),
        "--input-manifest",
        str(input_manifest),
        "--work-dir",
        str(node_dir),
        "--trim-length",
        str(args.trim_length),
        "--error-dist",
        args.error_dist,
        "--min-reads",
        str(args.min_reads),
        "--jobs-to-start",
        str(args.shard_workers),
        "--timings-tsv",
        str(node_dir / "timings.tsv"),
    ]
    if args.keep_tmp_files:
        command.append("--keep-tmp-files")
    log_path = node_dir / "deblur.log"
    try:
        with log_path.open("w", encoding="utf-8") as log_handle:
            run_command(
                command,
                timing=timing,
                step="deblur_shard",
                item=node_id,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
            )
        if not workflow_is_complete(workflow_dir):
            raise RuntimeError("Deblur returned without a complete workflow")
        return [LeafResult(node_id, entries, workflow_dir)]
    except (subprocess.CalledProcessError, RuntimeError, OSError) as exc:
        log_text = (
            log_path.read_text(encoding="utf-8", errors="replace")
            if log_path.is_file()
            else ""
        )
        if is_infrastructure_failure(exc, log_text):
            raise InfrastructureFailure(
                f"Deblur infrastructure failure in {node_id}; refusing sample-level "
                f"bisection. See {log_path}: {type(exc).__name__}: {exc}"
            ) from exc
        if len(entries) == 1:
            return [
                LeafResult(node_id, entries, None, f"{type(exc).__name__}: {exc}")
            ]
        left, right = split_entries(entries)
        write_json_atomic(
            split_marker,
            {
                "left": [entry.sample_id for entry in left],
                "right": [entry.sample_id for entry in right],
            },
        )
        results: list[LeafResult] = []
        for suffix, child in (("L", left), ("R", right)):
            if child:
                results.extend(run_node(f"{node_id}-{suffix}", child, args, timing))
        return results


def asdict_leaf(result: LeafResult) -> dict[str, Any]:
    return {
        "node_id": result.node_id,
        "sample_ids": [entry.sample_id for entry in result.entries],
        "workflow_dir": str(result.workflow_dir) if result.workflow_dir else None,
        "error": result.error,
    }


def run_node(
    node_id: str,
    entries: tuple[FastqEntry, ...],
    args: argparse.Namespace,
    timing: TimingRecorder,
) -> list[LeafResult]:
    return run_leaf(node_id, entries, args, timing)


def deblur_ids_for_entry(entry: FastqEntry) -> tuple[str, str]:
    name = entry.fastq_path.name
    return entry.sample_id, name[: -len(".fastq.gz")]


def sanitation_samples(workflow_dir: Path, entries: tuple[FastqEntry, ...]) -> set[str]:
    report = workflow_dir.parent / "invalid_fastq_records.tsv"
    if not report.is_file():
        return set()
    basename_to_sample = {entry.fastq_path.name: entry.sample_id for entry in entries}
    with report.open(newline="", encoding="utf-8") as handle:
        return {
            basename_to_sample[row["fastq"]]
            for row in csv.DictReader(handle, delimiter="\t")
            if row["fastq"] in basename_to_sample and int(row["records_dropped"]) > 0
        }


def merge_workflows(leaves: list[LeafResult], output_dir: Path) -> set[str]:
    tables = []
    sequences: dict[str, str] = {}
    observed_samples: set[str] = set()
    for leaf in leaves:
        if leaf.workflow_dir is None:
            continue
        table = biom.load_table(str(leaf.workflow_dir / "all.biom"))
        if table.shape == (0, 0):
            continue
        tables.append(table)
        observed_samples.update(str(value) for value in table.ids(axis="sample"))
        sequences.update(load_seqs(leaf.workflow_dir / "all.seqs.fa"))
    if not tables:
        raise ValueError("No non-empty Deblur shard outputs are available to merge")
    merged = tables[0]
    for table in tables[1:]:
        merged = merged.merge(table)
    output_dir.mkdir(parents=True, exist_ok=True)
    with biom.util.biom_open(str(output_dir / "all.biom"), "w") as handle:
        merged.to_hdf5(handle, "deblur_scheduler.py")
    feature_ids = set(str(value) for value in merged.ids(axis="observation"))
    with (output_dir / "all.seqs.fa").open("w", encoding="utf-8") as handle:
        for sequence_id in sorted(feature_ids):
            if sequence_id in sequences:
                handle.write(f">{sequence_id}\n{sequences[sequence_id]}\n")
    return observed_samples


def allowed_failures(total: int, fraction: float, absolute: int | None) -> int:
    fraction_limit = int(total * fraction)
    return absolute if absolute is not None else fraction_limit


def write_sample_status(
    path: Path, entries: list[FastqEntry], leaves: list[LeafResult], observed: set[str]
) -> dict[str, int]:
    leaf_by_sample = {
        entry.sample_id: leaf for leaf in leaves for entry in leaf.entries
    }
    sanitized = set()
    for leaf in leaves:
        if leaf.workflow_dir:
            sanitized.update(sanitation_samples(leaf.workflow_dir, leaf.entries))
    counts = {
        "expected": len(entries),
        "completed": 0,
        "completed_after_sanitation": 0,
        "zero_features": 0,
        "failed": 0,
        "excluded_assay": 0,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ("sample_id", "study", "fastq_path", "status", "node_id", "message")
        )
        for entry in sorted(entries, key=lambda item: (item.study, item.sample_id)):
            if entry.excluded_assay:
                status = "excluded_assay"
                message = "Excluded by exact sample-ID assay manifest"
                node_id = ""
            else:
                leaf = leaf_by_sample[entry.sample_id]
                node_id = leaf.node_id
            if not entry.excluded_assay and leaf.workflow_dir is None:
                status = "failed"
                message = leaf.error or "persistent Deblur failure"
            elif not entry.excluded_assay and any(
                identifier in observed for identifier in deblur_ids_for_entry(entry)
            ):
                status = (
                    "completed_after_sanitation"
                    if entry.sample_id in sanitized
                    else "completed"
                )
                message = ""
            elif not entry.excluded_assay:
                status = "zero_features"
                message = "Deblur emitted no features for this sample"
            counts[status] += 1
            writer.writerow(
                (
                    entry.sample_id,
                    entry.study,
                    entry.fastq_path,
                    status,
                    node_id,
                    message,
                )
            )
    return counts


def run(args: argparse.Namespace, timing: TimingRecorder) -> int:
    if args.min_reads != 0:
        raise ValueError("Balanced shard scheduling requires --min-reads 0")
    if args.shard_workers <= 0 or args.concurrent_shards <= 0:
        raise ValueError("Shard worker counts must be positive")
    if not 0 <= args.max_failed_fraction <= 1:
        raise ValueError("--max-failed-fraction must be between 0 and 1")
    if args.max_failed_samples is not None and args.max_failed_samples < 0:
        raise ValueError("--max-failed-samples cannot be negative")

    args.work_dir = args.work_dir.resolve()
    args.results_dir = args.results_dir.resolve()
    # Each Deblur process is an explicit worker. Prevent numerical libraries
    # loaded by those processes from silently spawning additional CPU pools.
    for env_name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[env_name] = "1"
    entries = read_shard_manifest(args.shard_manifest.resolve())
    grouped: dict[int, list[FastqEntry]] = {}
    for entry in entries:
        if not entry.excluded_assay:
            grouped.setdefault(entry.shard_id, []).append(entry)

    leaves: list[LeafResult] = []
    with ThreadPoolExecutor(
        max_workers=min(args.concurrent_shards, len(grouped))
    ) as executor:
        futures = {
            executor.submit(
                run_node,
                f"shard-{shard_id:04d}",
                tuple(sorted(shard_entries, key=lambda item: item.sample_id)),
                args,
                timing,
            ): shard_id
            for shard_id, shard_entries in sorted(grouped.items())
        }
        for future in as_completed(futures):
            leaves.extend(future.result())

    leaves.sort(key=lambda item: item.node_id)
    merge_error: ValueError | None = None
    try:
        with timing.step("merge_deblur_shards", item=f"{len(leaves)} leaf workflows"):
            observed = merge_workflows(leaves, args.work_dir / "workflow")
    except ValueError as exc:
        if "No non-empty Deblur shard outputs" not in str(exc):
            raise
        observed = set()
        merge_error = exc
    counts = write_sample_status(
        args.results_dir / "sample_processing_status.tsv", entries, leaves, observed
    )
    failed_limit = allowed_failures(
        counts["expected"] - counts["excluded_assay"],
        args.max_failed_fraction,
        args.max_failed_samples,
    )
    if merge_error is not None:
        status = "failed_no_features"
    elif counts["failed"] > failed_limit:
        status = "failed_tolerance_exceeded"
    elif counts["failed"] or counts["zero_features"]:
        status = "completed_with_exclusions"
    else:
        status = "completed"
    summary = {
        "status": status,
        "counts": counts,
        "observed_deblur_sample_ids": len(observed),
        "failed_sample_limit": failed_limit,
        "failure_fraction": counts["failed"]
        / (counts["expected"] - counts["excluded_assay"]),
        "leaf_nodes": [asdict_leaf(leaf) for leaf in leaves],
    }
    write_json_atomic(args.results_dir / "deblur_processing_summary.json", summary)
    if merge_error is not None:
        raise merge_error
    if counts["failed"] > failed_limit:
        raise RuntimeError(
            f"{counts['failed']} samples failed Deblur, exceeding allowed limit {failed_limit}"
        )
    print(
        f"Deblur complete: {counts['completed']} completed, "
        f"{counts['completed_after_sanitation']} sanitized, "
        f"{counts['zero_features']} zero-feature, {counts['failed']} failed.",
        flush=True,
    )
    return 0


def main() -> int:
    args = parse_args()
    timing = TimingRecorder(args.timings_tsv, component="deblur_scheduler")
    return run_timed_main(timing, lambda: run(args, timing))


if __name__ == "__main__":
    raise SystemExit(main())
