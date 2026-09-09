#!/usr/bin/env python3
"""Run the prepared-FASTQ-to-PCoA workflow in a self-contained directory."""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import socket
import statistics
import subprocess
import sys
import threading
from concurrent.futures import CancelledError, Future, ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from pipeline_lib import TimingRecorder, discover_inputs, resolve_executable, run_command
SCHEMA_VERSION = 1
STUDY_NAME_PATTERN = re.compile(r"^[A-Za-z0-9._-]+$")
GG2_BACKBONE_FILENAME = "2024.09.backbone.full-length.fna.qza"
GG2_ID_TREE_FILENAME = "2024.09.phylogeny.id.nwk.qza"
UNIFRAC_OUTPUT_NAMES = (
    "pcoa_coordinates_unweighted_unifrac.txt",
    "pcoa_plot_unweighted_unifrac.png",
    "analysis_summary.json",
)


@dataclass(frozen=True)
class Study:
    name: str
    fastq_dir: Path
    fastq_paths: tuple[Path, ...]


@dataclass(frozen=True)
class Stage:
    name: str
    command: tuple[str, ...]
    expected_outputs: tuple[Path, ...]
    dependencies: tuple[str, ...] = ()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run Deblur, optional cross-study BIOM merge, GG2/UniFrac, and "
            "metadata-colored PCoA plots. Inputs must be prepared forward-read FASTQs."
        )
    )
    parser.add_argument(
        "--study",
        action="append",
        required=True,
        metavar="NAME=FASTQ_DIR",
        help="Study name and prepared FASTQ directory; repeat for cross-study runs.",
    )
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument(
        "--exclude-samples",
        type=Path,
        help=(
            "Optional TXT/TSV/CSV exact sample-ID manifest to exclude known "
            "non-target assays before Deblur."
        ),
    )
    parser.add_argument(
        "--color-by",
        action="append",
        default=[],
        help="Metadata column used for a colored PCoA plot; repeat as needed.",
    )
    parser.add_argument("--trim-length", type=int, default=150)
    parser.add_argument(
        "--error-dist",
        default="1,0.06,0.02,0.02,0.01,0.005,0.005,0.005,0.001,0.001,0.001,0.0005",
    )
    parser.add_argument("--min-reads", type=int, default=0)
    parser.add_argument(
        "--sampling-depth",
        type=int,
        default=1000,
        help="Rarefaction depth after Greengenes2 mapping (default: 1000).",
    )
    parser.add_argument("--gg2-dir", type=Path, default=Path("data/gg2"))
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--deblur-scheduler",
        choices=("balanced-shards", "study"),
        default="balanced-shards",
        help="Deblur scheduling policy (default: balanced-shards).",
    )
    parser.add_argument(
        "--deblur-shard-workers",
        type=int,
        default=2,
        help="Deblur processes within each balanced shard (default: 2).",
    )
    parser.add_argument(
        "--deblur-shard-count",
        type=int,
        help="Initial balanced shard count (default: 2 x --threads).",
    )
    parser.add_argument(
        "--deblur-study-workers",
        type=int,
        default=1,
        help=(
            "Number of independent study-level Deblur workflows to run concurrently "
            "(default: 1)."
        ),
    )
    parser.add_argument(
        "--deblur-jobs-per-study",
        type=int,
        default=None,
        help=(
            "Deblur sample jobs allowed within each concurrent study. By default, "
            "divide --threads across the active study workers."
        ),
    )
    parser.add_argument(
        "--max-failed-studies",
        type=int,
        default=0,
        help=(
            "Maximum number of persistently failed Deblur studies that may be "
            "excluded before merge (default: 0). Exclusions are recorded."
        ),
    )
    parser.add_argument(
        "--max-failed-fraction",
        type=float,
        default=0.02,
        help="Maximum persistently failed sample fraction in balanced mode (default: 0.02).",
    )
    parser.add_argument(
        "--max-failed-samples",
        type=int,
        help="Absolute failed-sample limit; overrides --max-failed-fraction.",
    )
    parser.add_argument(
        "--pcoa-method",
        choices=("auto", "eigh", "fsvd"),
        default="auto",
    )
    parser.add_argument("--pcoa-dimensions", type=int, default=10)
    parser.add_argument("--pcoa-memory-budget-gb", type=float)
    parser.add_argument(
        "--export-distance-tsv",
        choices=("auto", "always", "never"),
        default="auto",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume a compatible existing run; completed stages are validated before skipping.",
    )
    parser.add_argument(
        "--keep-deblur-tmp-files",
        action="store_true",
        help="Retain Deblur workflow internals for debugging (large; disabled by default).",
    )
    return parser.parse_args(argv)


def parse_studies(values: list[str]) -> list[Study]:
    studies: list[Study] = []
    seen_names: set[str] = set()
    seen_sample_ids: dict[str, Path] = {}

    for value in values:
        if "=" not in value:
            raise ValueError(f"Invalid --study {value!r}; expected NAME=FASTQ_DIR")
        name, raw_dir = value.split("=", 1)
        if not name or not raw_dir:
            raise ValueError(f"Invalid --study {value!r}; expected NAME=FASTQ_DIR")
        if name in {".", ".."} or not STUDY_NAME_PATTERN.fullmatch(name):
            raise ValueError(
                f"Invalid study name {name!r}; use letters, numbers, '.', '_', or '-', "
                "but not '.' or '..'"
            )
        if name in seen_names:
            raise ValueError(f"Duplicate study name: {name}")

        fastq_dir = Path(raw_dir).expanduser().resolve()
        if not fastq_dir.is_dir():
            raise FileNotFoundError(f"Study FASTQ directory not found: {fastq_dir}")
        fastq_paths = tuple(path.resolve() for path in discover_inputs(fastq_dir))
        for fastq_path in fastq_paths:
            sample_id = sample_id_from_fastq(fastq_path)
            previous = seen_sample_ids.get(sample_id)
            if previous is not None:
                raise ValueError(
                    f"Duplicate sample identifier {sample_id!r}: {previous} and {fastq_path}"
                )
            seen_sample_ids[sample_id] = fastq_path

        studies.append(Study(name, fastq_dir, fastq_paths))
        seen_names.add(name)

    return studies


def sample_id_from_fastq(path: Path) -> str:
    for suffix in ("_R1_001.fastq.gz", "_1.fastq.gz"):
        if path.name.endswith(suffix):
            return path.name[: -len(suffix)]
    raise ValueError(f"Unsupported forward FASTQ name: {path.name}")


def safe_output_name(value: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._")
    if not safe:
        raise ValueError(f"Metadata column cannot form a safe output filename: {value!r}")
    return safe


def file_fingerprint(path: Path) -> dict[str, str | int]:
    stat = path.stat()
    return {
        "path": str(path),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def validate_metadata(metadata: Path | None, color_by: list[str]) -> Path | None:
    if color_by and metadata is None:
        raise ValueError("--metadata is required when --color-by is supplied")
    if metadata is None:
        return None

    metadata = metadata.expanduser().resolve()
    if not metadata.is_file() or metadata.stat().st_size == 0:
        raise FileNotFoundError(f"Metadata file not found or empty: {metadata}")

    delimiter = "\t" if metadata.suffix.lower() == ".tsv" else ","
    with metadata.open(newline="", encoding="utf-8-sig") as handle:
        fields = csv.DictReader(handle, delimiter=delimiter).fieldnames or []
    if "run_accessions" not in fields and "sample-id" not in fields:
        raise ValueError(
            "Metadata must contain either a 'run_accessions' or 'sample-id' column"
        )
    missing = [column for column in color_by if column not in fields]
    if missing:
        raise ValueError(
            f"Metadata column(s) not found: {', '.join(missing)}. "
            f"Available: {', '.join(fields)}"
        )

    output_names = [safe_output_name(column) for column in color_by]
    if len(set(output_names)) != len(output_names):
        raise ValueError("--color-by columns produce duplicate output filenames")
    return metadata


def load_sample_id_manifest(path: Path) -> set[str]:
    path = path.expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Sample exclusion manifest not found or empty: {path}")
    if path.suffix.lower() in {".csv", ".tsv"}:
        delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
        with path.open(newline="", encoding="utf-8-sig") as handle:
            rows = list(csv.DictReader(handle, delimiter=delimiter))
        if not rows:
            raise ValueError(f"Sample exclusion manifest has no data rows: {path}")
        fields = rows[0].keys()
        id_column = next(
            (
                column
                for column in ("sample_id", "sample-id", "run_accessions", "run_accession")
                if column in fields
            ),
            None,
        )
        if id_column is None:
            raise ValueError(
                "Sample exclusion table needs sample_id, sample-id, run_accessions, "
                "or run_accession"
            )
        values = {row[id_column].strip() for row in rows if row[id_column].strip()}
    else:
        values = {
            line.strip()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        }
    if not values:
        raise ValueError(f"Sample exclusion manifest contains no IDs: {path}")
    return values


def validate_args(
    args: argparse.Namespace,
) -> tuple[list[Study], Path | None, dict[str, str]]:
    if args.trim_length <= 0:
        raise ValueError("--trim-length must be greater than zero")
    if args.min_reads < 0:
        raise ValueError("--min-reads cannot be negative")
    if args.sampling_depth is not None and args.sampling_depth <= 0:
        raise ValueError("--sampling-depth must be greater than zero")
    if args.threads <= 0:
        raise ValueError("--threads must be greater than zero")
    if args.deblur_study_workers <= 0:
        raise ValueError("--deblur-study-workers must be greater than zero")
    if args.deblur_jobs_per_study is not None and args.deblur_jobs_per_study <= 0:
        raise ValueError("--deblur-jobs-per-study must be greater than zero")
    if args.max_failed_studies < 0:
        raise ValueError("--max-failed-studies cannot be negative")
    if args.deblur_shard_workers <= 0:
        raise ValueError("--deblur-shard-workers must be greater than zero")
    if args.deblur_shard_count is not None and args.deblur_shard_count <= 0:
        raise ValueError("--deblur-shard-count must be greater than zero")
    if not 0 <= args.max_failed_fraction <= 1:
        raise ValueError("--max-failed-fraction must be between zero and one")
    if args.max_failed_samples is not None and args.max_failed_samples < 0:
        raise ValueError("--max-failed-samples cannot be negative")
    if args.pcoa_dimensions < 2:
        raise ValueError("--pcoa-dimensions must be at least two")
    if args.pcoa_memory_budget_gb is not None and args.pcoa_memory_budget_gb <= 0:
        raise ValueError("--pcoa-memory-budget-gb must be greater than zero")
    if args.deblur_scheduler == "balanced-shards" and args.min_reads != 0:
        raise ValueError("Balanced Deblur scheduling requires --min-reads 0")
    if not args.error_dist.strip():
        raise ValueError("--error-dist cannot be empty")

    studies = parse_studies(args.study)
    excluded_samples: set[str] = set()
    if args.exclude_samples is not None:
        if args.deblur_scheduler != "balanced-shards":
            raise ValueError(
                "--exclude-samples requires --deblur-scheduler balanced-shards"
            )
        args.exclude_samples = args.exclude_samples.expanduser().resolve()
        excluded_samples = load_sample_id_manifest(args.exclude_samples)
        available = {
            sample_id_from_fastq(path)
            for study in studies
            for path in study.fastq_paths
        }
        unknown = sorted(excluded_samples - available)
        if unknown:
            raise ValueError(
                "Sample exclusion manifest contains IDs absent from FASTQ inputs: "
                + ", ".join(unknown[:10])
            )
        if excluded_samples == available:
            raise ValueError("Sample exclusion manifest excludes every FASTQ input")
    if args.deblur_scheduler == "study":
        resolve_deblur_parallelism(args, len(studies))
    else:
        resolve_balanced_parallelism(
            args,
            sum(len(study.fastq_paths) for study in studies) - len(excluded_samples),
        )
    metadata = validate_metadata(args.metadata, args.color_by)
    gg2_dir = args.gg2_dir.expanduser().resolve()
    missing_gg2 = [
        gg2_dir / filename
        for filename in (GG2_BACKBONE_FILENAME, GG2_ID_TREE_FILENAME)
        if not (gg2_dir / filename).is_file()
        or (gg2_dir / filename).stat().st_size == 0
    ]
    if missing_gg2:
        raise FileNotFoundError(
            "Missing or empty Greengenes2 artifact(s): "
            + ", ".join(str(path) for path in missing_gg2)
        )

    executables = {
        "python": str(Path(sys.executable).resolve()),
        "deblur": resolve_executable("deblur"),
        "qiime": resolve_executable("qiime"),
    }
    validate_qiime_gg2(executables["qiime"])
    args.run_dir = args.run_dir.expanduser().resolve()
    args.gg2_dir = gg2_dir
    return studies, metadata, executables


def resolve_deblur_parallelism(
    args: argparse.Namespace, study_count: int
) -> tuple[int, int]:
    if study_count <= 0:
        raise ValueError("At least one study is required")
    study_workers = min(args.deblur_study_workers, study_count)
    jobs_per_study = args.deblur_jobs_per_study
    if jobs_per_study is None:
        jobs_per_study = max(1, args.threads // study_workers)
    requested_jobs = study_workers * jobs_per_study
    if requested_jobs > args.threads:
        raise ValueError(
            "Deblur parallelism oversubscribes --threads: "
            f"{study_workers} active studies x {jobs_per_study} jobs per study = "
            f"{requested_jobs}, but --threads={args.threads}"
        )
    return study_workers, jobs_per_study


def resolve_balanced_parallelism(
    args: argparse.Namespace, sample_count: int
) -> tuple[int, int, int]:
    if sample_count <= 0:
        raise ValueError("At least one FASTQ is required")
    shard_workers = args.deblur_shard_workers
    concurrent_shards = max(1, args.threads // shard_workers)
    if concurrent_shards * shard_workers > args.threads:
        raise ValueError("Balanced Deblur parallelism oversubscribes --threads")
    requested_shards = args.deblur_shard_count or (2 * args.threads)
    return min(requested_shards, sample_count), concurrent_shards, shard_workers


def eligible_sample_count(args: argparse.Namespace, studies: list[Study]) -> int:
    total = sum(len(study.fastq_paths) for study in studies)
    excluded = (
        len(load_sample_id_manifest(args.exclude_samples))
        if args.exclude_samples
        else 0
    )
    return total - excluded


def balanced_shard_rows(
    studies: list[Study], shard_count: int, excluded_samples: set[str] | None = None
) -> list[dict[str, str | int]]:
    excluded_samples = excluded_samples or set()
    paths = [path for study in studies for path in study.fastq_paths]
    median_size = int(statistics.median(path.stat().st_size for path in paths))
    weighted = [
        (
            path.stat().st_size + median_size,
            sample_id_from_fastq(path),
            study.name,
            path,
            path.stat().st_size,
        )
        for study in studies
        for path in study.fastq_paths
    ]
    totals = [0] * shard_count
    rows: list[dict[str, str | int]] = []
    for weight, sample_id, study_name, path, size in sorted(
        weighted, key=lambda item: (-item[0], item[1], item[2])
    ):
        excluded = sample_id in excluded_samples
        shard_id = -1 if excluded else min(
            range(shard_count), key=lambda index: (totals[index], index)
        )
        if not excluded:
            totals[shard_id] += weight
        rows.append(
            {
                "sample_id": sample_id,
                "study": study_name,
                "fastq_path": str(path),
                "size_bytes": size,
                "weight_bytes": weight,
                "shard_id": shard_id,
                "excluded_assay": int(excluded),
            }
        )
    return sorted(rows, key=lambda row: (int(row["shard_id"]), str(row["sample_id"])))


def write_shard_manifest(path: Path, rows: list[dict[str, str | int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def validate_qiime_gg2(qiime: str) -> None:
    try:
        subprocess.run(
            [qiime, "greengenes2", "--help"],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        detail = getattr(exc, "stderr", "") or str(exc)
        raise RuntimeError(
            "QIIME2 could not load the greengenes2 plugin: " + detail.strip()
        ) from exc


def git_info(repo_dir: Path) -> tuple[str | None, bool | None]:
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_dir,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=no"],
            cwd=repo_dir,
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None, None
    return commit, bool(status.strip())


def build_manifest(
    args: argparse.Namespace,
    studies: list[Study],
    metadata: Path | None,
    executables: dict[str, str],
    repo_dir: Path,
) -> dict[str, Any]:
    commit, tracked_dirty = git_info(repo_dir)
    sample_count = sum(len(study.fastq_paths) for study in studies)
    excluded_samples = (
        load_sample_id_manifest(args.exclude_samples) if args.exclude_samples else set()
    )
    balanced = None
    if args.deblur_scheduler == "balanced-shards":
        shard_count, concurrent_shards, shard_workers = resolve_balanced_parallelism(
            args, sample_count - len(excluded_samples)
        )
        balanced = {
            "shard_count": shard_count,
            "concurrent_shards": concurrent_shards,
            "shard_workers": shard_workers,
            "assignments": balanced_shard_rows(
                studies, shard_count, excluded_samples
            ),
            "max_failed_fraction": args.max_failed_fraction,
            "max_failed_samples": args.max_failed_samples,
        }
    return {
        "schema_version": SCHEMA_VERSION,
        "created_utc": utc_now(),
        "git_commit": commit,
        "git_tracked_dirty": tracked_dirty,
        "hostname": socket.gethostname(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "conda_environment": os.environ.get("CONDA_DEFAULT_ENV"),
        "conda_prefix": os.environ.get("CONDA_PREFIX"),
        "python_version": sys.version,
        "executables": executables,
        "studies": [
            {
                "name": study.name,
                "fastq_dir": str(study.fastq_dir),
                "fastqs": [file_fingerprint(path) for path in study.fastq_paths],
            }
            for study in studies
        ],
        "metadata": file_fingerprint(metadata) if metadata else None,
        "excluded_samples": (
            {
                "file": file_fingerprint(args.exclude_samples),
                "sample_ids": sorted(excluded_samples),
            }
            if args.exclude_samples
            else None
        ),
        "color_by": list(args.color_by),
        "gg2_dir": str(args.gg2_dir),
        "gg2_artifacts": [
            file_fingerprint(args.gg2_dir / filename)
            for filename in (GG2_BACKBONE_FILENAME, GG2_ID_TREE_FILENAME)
        ],
        "scientific_parameters": {
            "trim_length": args.trim_length,
            "error_dist": args.error_dist,
            "min_reads": args.min_reads,
            "sampling_depth": args.sampling_depth,
            "pcoa_method": args.pcoa_method,
            "pcoa_dimensions": args.pcoa_dimensions,
        },
        "deblur_scheduler": args.deblur_scheduler,
        "balanced_deblur": balanced,
        "pcoa_memory_budget_gb": args.pcoa_memory_budget_gb,
        "export_distance_tsv": args.export_distance_tsv,
        "keep_deblur_tmp_files": args.keep_deblur_tmp_files,
    }


def compatibility_view(manifest: dict[str, Any]) -> dict[str, Any]:
    keys = (
        "studies",
        "metadata",
        "excluded_samples",
        "color_by",
        "gg2_dir",
        "gg2_artifacts",
        "scientific_parameters",
        "deblur_scheduler",
        "balanced_deblur",
        "pcoa_memory_budget_gb",
        "export_distance_tsv",
        "keep_deblur_tmp_files",
    )
    return {key: manifest.get(key) for key in keys}


def write_json_atomic(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def read_json(path: Path) -> dict[str, Any]:
    try:
        with path.open(encoding="utf-8") as handle:
            value = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Cannot read run metadata {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise RuntimeError(f"Run metadata is not a JSON object: {path}")
    return value


def stage_output_paths(
    run_dir: Path,
    studies: list[Study],
    color_by: list[str],
    deblur_scheduler: str = "study",
) -> dict[str, tuple[Path, ...]]:
    outputs: dict[str, tuple[Path, ...]] = {}
    if deblur_scheduler == "balanced-shards":
        workflow_dir = run_dir / "work" / "deblur-balanced" / "workflow"
        outputs["deblur:balanced"] = (
            workflow_dir / "all.biom",
            workflow_dir / "all.seqs.fa",
            run_dir / "results" / "sample_processing_status.tsv",
            run_dir / "results" / "deblur_processing_summary.json",
        )
    else:
        for study in studies:
            workflow_dir = run_dir / "work" / "deblur" / study.name / "workflow"
            outputs[f"deblur:{study.name}"] = (
                workflow_dir / "all.biom",
                workflow_dir / "all.seqs.fa",
            )
    if deblur_scheduler == "study" and len(studies) > 1:
        merged_dir = run_dir / "work" / "merged"
        outputs["merge"] = (merged_dir / "all.biom", merged_dir / "all.seqs.fa")
    results_dir = run_dir / "results"
    outputs["unifrac"] = tuple(results_dir / name for name in UNIFRAC_OUTPUT_NAMES)
    for column in color_by:
        outputs[f"plot:{column}"] = (results_dir / f"pcoa_{safe_output_name(column)}.png",)
    return outputs


def build_stages(
    args: argparse.Namespace,
    studies: list[Study],
    metadata: Path | None,
    attempt_dir: Path,
    repo_dir: Path,
) -> list[Stage]:
    run_dir = args.run_dir
    output_paths = stage_output_paths(
        run_dir, studies, args.color_by, args.deblur_scheduler
    )
    stages: list[Stage] = []
    deblur_workflows: list[Path] = []
    if args.deblur_scheduler == "balanced-shards":
        shard_count, concurrent_shards, shard_workers = resolve_balanced_parallelism(
            args, eligible_sample_count(args, studies)
        )
        balanced_dir = run_dir / "work" / "deblur-balanced"
        command = [
            sys.executable,
            str(repo_dir / "src" / "deblur_scheduler.py"),
            "--shard-manifest",
            str(balanced_dir / "shard_manifest.tsv"),
            "--work-dir",
            str(balanced_dir),
            "--results-dir",
            str(run_dir / "results"),
            "--trim-length",
            str(args.trim_length),
            "--error-dist",
            args.error_dist,
            "--min-reads",
            str(args.min_reads),
            "--shard-workers",
            str(shard_workers),
            "--concurrent-shards",
            str(concurrent_shards),
            "--max-failed-fraction",
            str(args.max_failed_fraction),
            "--timings-tsv",
            str(attempt_dir / "deblur-balanced.tsv"),
        ]
        if args.max_failed_samples is not None:
            command.extend(("--max-failed-samples", str(args.max_failed_samples)))
        if args.keep_deblur_tmp_files:
            command.append("--keep-tmp-files")
        stages.append(
            Stage("deblur:balanced", tuple(command), output_paths["deblur:balanced"])
        )
        unifrac_input_dir = balanced_dir / "workflow"
        unifrac_dependencies = ("deblur:balanced",)
    else:
        _, deblur_jobs_per_study = resolve_deblur_parallelism(args, len(studies))

        for study in studies:
            deblur_dir = run_dir / "work" / "deblur" / study.name
            workflow_dir = deblur_dir / "workflow"
            deblur_workflows.append(workflow_dir)
            command = [
                sys.executable,
                str(repo_dir / "src" / "run_deblur.py"),
                "--data-dir",
                str(study.fastq_dir),
                "--work-dir",
                str(deblur_dir),
                "--trim-length",
                str(args.trim_length),
                "--error-dist",
                args.error_dist,
                "--min-reads",
                str(args.min_reads),
                "--jobs-to-start",
                str(deblur_jobs_per_study),
                "--timings-tsv",
                str(attempt_dir / f"deblur-{study.name}.tsv"),
            ]
            if args.keep_deblur_tmp_files:
                command.append("--keep-tmp-files")
            stages.append(
                Stage(
                    f"deblur:{study.name}",
                    tuple(command),
                    output_paths[f"deblur:{study.name}"],
                )
            )

    if args.deblur_scheduler == "study" and len(studies) > 1:
        merged_dir = run_dir / "work" / "merged"
        command = [
            sys.executable,
            str(repo_dir / "src" / "merge_biom.py"),
            "--deblur-dirs",
            *(str(path) for path in deblur_workflows),
            "--out-dir",
            str(merged_dir),
            "--skip-empty",
            "--run-state",
            str(run_dir / "run_state.json"),
            "--timings-tsv",
            str(attempt_dir / "merge.tsv"),
        ]
        stages.append(
            Stage(
                "merge",
                tuple(command),
                output_paths["merge"],
                tuple(f"deblur:{study.name}" for study in studies),
            )
        )
        unifrac_input_dir = merged_dir
        unifrac_dependencies = ("merge",)
    elif args.deblur_scheduler == "study":
        unifrac_input_dir = deblur_workflows[0]
        unifrac_dependencies = (f"deblur:{studies[0].name}",)

    results_dir = run_dir / "results"
    command = [
        sys.executable,
        str(repo_dir / "src" / "unifrac.py"),
        "--deblur-dir",
        str(unifrac_input_dir),
        "--results-dir",
        str(results_dir),
        "--gg2-dir",
        str(args.gg2_dir),
        "--threads",
        str(args.threads),
        "--work-dir",
        str(run_dir / "work" / "qiime2"),
        "--timings-tsv",
        str(attempt_dir / "unifrac.tsv"),
        "--pcoa-method",
        args.pcoa_method,
        "--pcoa-dimensions",
        str(args.pcoa_dimensions),
        "--export-distance-tsv",
        args.export_distance_tsv,
    ]
    if args.pcoa_memory_budget_gb is not None:
        command.extend(("--pcoa-memory-budget-gb", str(args.pcoa_memory_budget_gb)))
    if args.sampling_depth is not None:
        command.extend(("--sampling-depth", str(args.sampling_depth)))
    stages.append(
        Stage(
            "unifrac",
            tuple(command),
            output_paths["unifrac"],
            unifrac_dependencies,
        )
    )

    for column in args.color_by:
        assert metadata is not None
        safe_column = safe_output_name(column)
        command = [
            sys.executable,
            str(repo_dir / "src" / "plot_pcoa.py"),
            "--pcoa",
            str(results_dir / "pcoa_coordinates_unweighted_unifrac.txt"),
            "--metadata",
            str(metadata),
            "--color-by",
            column,
            "--out",
            str(results_dir / f"pcoa_{safe_column}.png"),
            "--title",
            f"{run_dir.name} - {column}",
            "--timings-tsv",
            str(attempt_dir / f"plot-{safe_column}.tsv"),
        ]
        stages.append(
            Stage(
                f"plot:{column}",
                tuple(command),
                output_paths[f"plot:{column}"],
                ("unifrac",),
            )
        )
    return stages


def outputs_exist(paths: tuple[Path, ...]) -> bool:
    return all(path.is_file() and path.stat().st_size > 0 for path in paths)


def initialize_run(
    args: argparse.Namespace,
    manifest: dict[str, Any],
    stage_names: list[str],
) -> tuple[dict[str, Any], int]:
    run_dir = args.run_dir
    manifest_path = run_dir / "run_manifest.json"
    state_path = run_dir / "run_state.json"

    if run_dir.exists() and not run_dir.is_dir():
        raise RuntimeError(f"Run path exists but is not a directory: {run_dir}")

    if args.resume:
        if not run_dir.is_dir():
            raise FileNotFoundError(f"Cannot resume missing run directory: {run_dir}")
        if not manifest_path.is_file() or not state_path.is_file():
            raise RuntimeError("Resume requires run_manifest.json and run_state.json")
        saved_manifest = read_json(manifest_path)
        if saved_manifest.get("schema_version") != SCHEMA_VERSION:
            raise RuntimeError("Unsupported run_manifest.json schema version")
        if compatibility_view(saved_manifest) != compatibility_view(manifest):
            raise RuntimeError(
                "Resume configuration does not match run_manifest.json; use a new --run-dir"
            )
        state = read_json(state_path)
        if state.get("schema_version") != SCHEMA_VERSION:
            raise RuntimeError("Unsupported run_state.json schema version")
        saved_stages = state.get("stages")
        if not isinstance(saved_stages, dict) or set(saved_stages) != set(stage_names):
            raise RuntimeError("Saved stage layout is incompatible with this invocation")
        attempts = state.get("attempts")
        if not isinstance(attempts, list):
            raise RuntimeError("Invalid attempts list in run_state.json")
        attempt_number = max(
            (int(attempt.get("number", 0)) for attempt in attempts if isinstance(attempt, dict)),
            default=0,
        ) + 1
        return state, attempt_number

    if run_dir.exists() and any(run_dir.iterdir()):
        raise RuntimeError(
            f"Run directory is not empty: {run_dir}. Use --resume for a compatible run."
        )
    run_dir.mkdir(parents=True, exist_ok=True)
    state = {
        "schema_version": SCHEMA_VERSION,
        "stages": {
            name: {
                "status": "pending",
                "attempt": None,
                "started_utc": None,
                "ended_utc": None,
                "error": None,
            }
            for name in stage_names
        },
        "attempts": [],
    }
    write_json_atomic(manifest_path, manifest)
    write_json_atomic(state_path, state)
    return state, 1


def run_stage(
    stage: Stage,
    state: dict[str, Any],
    state_path: Path,
    attempt_number: int,
    timing: TimingRecorder,
    state_lock: threading.Lock | None = None,
) -> None:
    state_lock = state_lock or threading.Lock()
    with state_lock:
        stage_state = state["stages"][stage.name]
        stage_state.update(
            {
                "status": "running",
                "attempt": attempt_number,
                "started_utc": utc_now(),
                "ended_utc": None,
                "error": None,
            }
        )
        write_json_atomic(state_path, state)
    try:
        run_command(
            list(stage.command),
            timing=timing,
            step=stage.name,
            item=", ".join(str(path) for path in stage.expected_outputs),
        )
        missing = [path for path in stage.expected_outputs if not outputs_exist((path,))]
        if missing:
            if is_empty_deblur_output(stage):
                message = "Deblur produced a zero-feature table; excluding study downstream"
                with state_lock:
                    stage_state = state["stages"][stage.name]
                    stage_state.update(
                        {
                            "status": "excluded_empty",
                            "ended_utc": utc_now(),
                            "error": message,
                        }
                    )
                    write_json_atomic(state_path, state)
                print(f"Excluded empty stage: {stage.name}", flush=True)
                timing.skipped(
                    stage.name,
                    item=str(stage.expected_outputs[0]),
                    message=message,
                )
                return
            raise RuntimeError(
                f"Stage {stage.name} completed without expected non-empty output(s): "
                + ", ".join(str(path) for path in missing)
            )
    except BaseException as exc:
        with state_lock:
            stage_state = state["stages"][stage.name]
            stage_state.update(
                {
                    "status": "failed",
                    "ended_utc": utc_now(),
                    "error": f"{type(exc).__name__}: {exc}",
                }
            )
            write_json_atomic(state_path, state)
        raise
    with state_lock:
        stage_state = state["stages"][stage.name]
        stage_state.update(
            {"status": "completed", "ended_utc": utc_now(), "error": None}
        )
        write_json_atomic(state_path, state)


def should_skip_stage(
    stage: Stage, state: dict[str, Any], rerun_stages: set[str]
) -> bool:
    stage_state = state["stages"][stage.name]
    if stage_state.get("status") == "excluded_empty":
        return is_empty_deblur_output(stage)
    return (
        stage_state.get("status") == "completed"
        and outputs_exist(stage.expected_outputs)
        and not rerun_stages.intersection(stage.dependencies)
    )


def is_empty_deblur_output(stage: Stage) -> bool:
    if not stage.name.startswith("deblur:") or len(stage.expected_outputs) != 2:
        return False
    biom_path, fasta_path = stage.expected_outputs
    return (
        biom_path.is_file()
        and biom_path.stat().st_size > 0
        and fasta_path.is_file()
        and fasta_path.stat().st_size == 0
    )


def record_skipped_stage(stage: Stage, timing: TimingRecorder) -> None:
    print(f"Skipping completed stage: {stage.name}", flush=True)
    timing.skipped(
        stage.name,
        item=", ".join(str(path) for path in stage.expected_outputs),
        message="state completed and expected outputs exist",
    )


def run_deblur_stages(
    stages: list[Stage],
    state: dict[str, Any],
    state_path: Path,
    attempt_number: int,
    timing: TimingRecorder,
    study_workers: int,
) -> tuple[set[str], set[str]]:
    rerun_stages: set[str] = set()
    failed_stages: set[str] = set()
    stages_to_run: list[Stage] = []
    for stage in stages:
        if should_skip_stage(stage, state, rerun_stages):
            record_skipped_stage(stage, timing)
        else:
            stages_to_run.append(stage)

    if not stages_to_run:
        return rerun_stages, failed_stages

    if study_workers == 1:
        for stage in stages_to_run:
            print(f"Running stage: {stage.name}", flush=True)
            try:
                run_stage(stage, state, state_path, attempt_number, timing)
            except BaseException as exc:
                failed_stages.add(stage.name)
                print(f"Failed independent stage: {stage.name}: {exc}", flush=True)
            else:
                rerun_stages.add(stage.name)
        return rerun_stages, failed_stages

    state_lock = threading.Lock()
    futures: dict[Future[None], Stage] = {}
    with ThreadPoolExecutor(max_workers=min(study_workers, len(stages_to_run))) as executor:
        for stage in stages_to_run:
            print(f"Submitting concurrent stage: {stage.name}", flush=True)
            future = executor.submit(
                run_stage,
                stage,
                state,
                state_path,
                attempt_number,
                timing,
                state_lock,
            )
            futures[future] = stage
        for future in as_completed(futures):
            stage = futures[future]
            try:
                future.result()
            except CancelledError:
                continue
            except BaseException as exc:
                failed_stages.add(stage.name)
                print(f"Failed independent stage: {stage.name}: {exc}", flush=True)
            else:
                rerun_stages.add(stage.name)
                print(f"Completed concurrent stage: {stage.name}", flush=True)
    return rerun_stages, failed_stages


def write_study_processing_status(
    run_dir: Path, studies: list[Study], state: dict[str, Any]
) -> Path:
    path = run_dir / "results" / "study_processing_status.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("study", "status", "attempt", "message"))
        for study in studies:
            stage = state["stages"][f"deblur:{study.name}"]
            writer.writerow(
                (
                    study.name,
                    stage.get("status", "pending"),
                    stage.get("attempt") or "",
                    stage.get("error") or "",
                )
            )
    return path


def write_pipeline_summary(run_dir: Path, attempt_status: str) -> Path:
    results_dir = run_dir / "results"
    summary: dict[str, Any] = {"status": attempt_status}
    deblur_summary = results_dir / "deblur_processing_summary.json"
    analysis_summary = results_dir / "analysis_summary.json"
    if deblur_summary.is_file():
        summary["deblur"] = read_json(deblur_summary)
    if analysis_summary.is_file():
        summary["analysis"] = read_json(analysis_summary)
    deblur_counts = summary.get("deblur", {}).get("counts", {})
    analysis_counts = summary.get("analysis", {})
    if deblur_counts:
        summary["counts"] = {
            "expected": deblur_counts.get("expected", 0),
            "processed": sum(
                deblur_counts.get(key, 0)
                for key in ("completed", "completed_after_sanitation", "zero_features")
            ),
            "sanitized": deblur_counts.get("completed_after_sanitation", 0),
            "zero_features": deblur_counts.get("zero_features", 0),
            "failed": deblur_counts.get("failed", 0),
            "excluded_assay": deblur_counts.get("excluded_assay", 0),
            "mapped": analysis_counts.get("mapped_samples"),
            "rarefied": analysis_counts.get("rarefied_samples"),
        }
    path = results_dir / "pipeline_summary.json"
    write_json_atomic(path, summary)
    return path


def execute_pipeline(args: argparse.Namespace) -> Path:
    repo_dir = Path(__file__).resolve().parent.parent
    studies, metadata, executables = validate_args(args)
    manifest = build_manifest(args, studies, metadata, executables, repo_dir)
    output_paths = stage_output_paths(
        args.run_dir, studies, args.color_by, args.deblur_scheduler
    )
    state, attempt_number = initialize_run(args, manifest, list(output_paths))

    if args.deblur_scheduler == "balanced-shards":
        balanced = manifest["balanced_deblur"]
        assert isinstance(balanced, dict)
        write_shard_manifest(
            args.run_dir / "work" / "deblur-balanced" / "shard_manifest.tsv",
            balanced["assignments"],
        )

    attempt_dir = args.run_dir / "timings" / f"attempt-{attempt_number:03d}"
    while attempt_dir.exists():
        attempt_number += 1
        attempt_dir = args.run_dir / "timings" / f"attempt-{attempt_number:03d}"
    attempt_dir.mkdir(parents=True)
    stages = build_stages(args, studies, metadata, attempt_dir, repo_dir)
    state_path = args.run_dir / "run_state.json"
    current_commit, tracked_dirty = git_info(repo_dir)
    if args.deblur_scheduler == "study":
        deblur_stage_workers, deblur_jobs_per_stage = resolve_deblur_parallelism(
            args, len(studies)
        )
    else:
        shard_count, concurrent_shards, shard_workers = resolve_balanced_parallelism(
            args, eligible_sample_count(args, studies)
        )
        deblur_stage_workers = 1
        deblur_jobs_per_stage = shard_workers
    attempt = {
        "number": attempt_number,
        "status": "running",
        "started_utc": utc_now(),
        "ended_utc": None,
        "threads": args.threads,
        "deblur_scheduler": args.deblur_scheduler,
        "deblur_study_workers": deblur_stage_workers if args.deblur_scheduler == "study" else None,
        "deblur_jobs_per_study": deblur_jobs_per_stage if args.deblur_scheduler == "study" else None,
        "deblur_shard_count": shard_count if args.deblur_scheduler == "balanced-shards" else None,
        "deblur_concurrent_shards": concurrent_shards if args.deblur_scheduler == "balanced-shards" else None,
        "deblur_shard_workers": deblur_jobs_per_stage if args.deblur_scheduler == "balanced-shards" else None,
        "max_failed_studies": args.max_failed_studies,
        "max_failed_fraction": args.max_failed_fraction,
        "max_failed_samples": args.max_failed_samples,
        "git_commit": current_commit,
        "git_tracked_dirty": tracked_dirty,
        "hostname": socket.gethostname(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "python_executable": str(Path(sys.executable).resolve()),
        "timings_dir": str(attempt_dir),
    }
    state["attempts"].append(attempt)
    write_json_atomic(state_path, state)

    timing = TimingRecorder(attempt_dir / "pipeline.tsv", component="run_pipeline")
    rerun_stages: set[str] = set()
    try:
        with timing.step("total", item=str(args.run_dir)):
            deblur_stages = [
                stage for stage in stages if stage.name.startswith("deblur:")
            ]
            downstream_stages = [
                stage for stage in stages if not stage.name.startswith("deblur:")
            ]
            rerun_deblur, failed_deblur = run_deblur_stages(
                deblur_stages,
                state,
                state_path,
                attempt_number,
                timing,
                deblur_stage_workers,
            )
            rerun_stages.update(rerun_deblur)
            if args.deblur_scheduler == "study":
                status_path = write_study_processing_status(args.run_dir, studies, state)
                print(f"Study processing status: {status_path}", flush=True)
            if args.deblur_scheduler == "study" and len(failed_deblur) > args.max_failed_studies:
                names = ", ".join(sorted(failed_deblur))
                raise RuntimeError(
                    f"{len(failed_deblur)} Deblur study stage(s) failed, exceeding "
                    f"--max-failed-studies={args.max_failed_studies}: {names}"
                )
            if args.deblur_scheduler == "balanced-shards" and failed_deblur:
                raise RuntimeError("Balanced Deblur scheduler failed; see its processing summary")
            if args.deblur_scheduler == "study" and failed_deblur:
                print(
                    f"Excluding {len(failed_deblur)} failed Deblur study stage(s) "
                    f"within configured limit: {', '.join(sorted(failed_deblur))}",
                    flush=True,
                )
            for stage in downstream_stages:
                if should_skip_stage(stage, state, rerun_stages):
                    record_skipped_stage(stage, timing)
                    continue

                print(f"Running stage: {stage.name}", flush=True)
                stage_to_run = stage
                if stage.name == "unifrac" and rerun_stages.intersection(
                    stage.dependencies
                ):
                    stage_to_run = Stage(
                        stage.name,
                        stage.command + ("--refresh-input-artifacts",),
                        stage.expected_outputs,
                        stage.dependencies,
                    )
                run_stage(stage_to_run, state, state_path, attempt_number, timing)
                rerun_stages.add(stage.name)
    except BaseException:
        attempt.update({"status": "failed", "ended_utc": utc_now()})
        write_json_atomic(state_path, state)
        raise

    final_status = "completed"
    deblur_summary_path = args.run_dir / "results" / "deblur_processing_summary.json"
    if deblur_summary_path.is_file():
        deblur_status = read_json(deblur_summary_path).get("status")
        if deblur_status == "completed_with_exclusions":
            final_status = "completed_with_exclusions"
    attempt.update({"status": final_status, "ended_utc": utc_now()})
    write_json_atomic(state_path, state)
    summary_path = write_pipeline_summary(args.run_dir, final_status)
    print(f"Pipeline summary: {summary_path}", flush=True)
    print(f"Pipeline completed. Results: {args.run_dir / 'results'}", flush=True)
    return args.run_dir


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        execute_pipeline(args)
    except (FileNotFoundError, RuntimeError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
