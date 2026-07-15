#!/usr/bin/env python3
"""Run the prepared-FASTQ-to-PCoA workflow in a self-contained directory."""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import socket
import subprocess
import sys
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
    "distance_matrix_unweighted_unifrac.tsv",
    "pcoa_coordinates_unweighted_unifrac.txt",
    "pcoa_plot_unweighted_unifrac.png",
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
    parser.add_argument("--sampling-depth", type=int)
    parser.add_argument("--gg2-dir", type=Path, default=Path("data/gg2"))
    parser.add_argument("--threads", type=int, default=4)
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


def validate_args(args: argparse.Namespace) -> tuple[list[Study], Path | None, dict[str, str]]:
    if args.trim_length <= 0:
        raise ValueError("--trim-length must be greater than zero")
    if args.min_reads < 0:
        raise ValueError("--min-reads cannot be negative")
    if args.sampling_depth is not None and args.sampling_depth <= 0:
        raise ValueError("--sampling-depth must be greater than zero")
    if args.threads <= 0:
        raise ValueError("--threads must be greater than zero")
    if not args.error_dist.strip():
        raise ValueError("--error-dist cannot be empty")

    studies = parse_studies(args.study)
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
        },
        "keep_deblur_tmp_files": args.keep_deblur_tmp_files,
    }


def compatibility_view(manifest: dict[str, Any]) -> dict[str, Any]:
    keys = (
        "git_commit",
        "studies",
        "metadata",
        "color_by",
        "gg2_dir",
        "gg2_artifacts",
        "scientific_parameters",
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


def stage_output_paths(run_dir: Path, studies: list[Study], color_by: list[str]) -> dict[str, tuple[Path, ...]]:
    outputs: dict[str, tuple[Path, ...]] = {}
    for study in studies:
        workflow_dir = run_dir / "work" / "deblur" / study.name / "workflow"
        outputs[f"deblur:{study.name}"] = (
            workflow_dir / "all.biom",
            workflow_dir / "all.seqs.fa",
        )
    if len(studies) > 1:
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
    output_paths = stage_output_paths(run_dir, studies, args.color_by)
    stages: list[Stage] = []
    deblur_workflows: list[Path] = []

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
            str(args.threads),
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

    if len(studies) > 1:
        merged_dir = run_dir / "work" / "merged"
        command = [
            sys.executable,
            str(repo_dir / "src" / "merge_biom.py"),
            "--deblur-dirs",
            *(str(path) for path in deblur_workflows),
            "--out-dir",
            str(merged_dir),
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
    else:
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
    ]
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
) -> None:
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
            raise RuntimeError(
                f"Stage {stage.name} completed without expected non-empty output(s): "
                + ", ".join(str(path) for path in missing)
            )
    except BaseException as exc:
        stage_state.update(
            {
                "status": "failed",
                "ended_utc": utc_now(),
                "error": f"{type(exc).__name__}: {exc}",
            }
        )
        write_json_atomic(state_path, state)
        raise
    stage_state.update({"status": "completed", "ended_utc": utc_now(), "error": None})
    write_json_atomic(state_path, state)


def execute_pipeline(args: argparse.Namespace) -> Path:
    repo_dir = Path(__file__).resolve().parent.parent
    studies, metadata, executables = validate_args(args)
    manifest = build_manifest(args, studies, metadata, executables, repo_dir)
    output_paths = stage_output_paths(args.run_dir, studies, args.color_by)
    state, attempt_number = initialize_run(args, manifest, list(output_paths))

    attempt_dir = args.run_dir / "timings" / f"attempt-{attempt_number:03d}"
    while attempt_dir.exists():
        attempt_number += 1
        attempt_dir = args.run_dir / "timings" / f"attempt-{attempt_number:03d}"
    attempt_dir.mkdir(parents=True)
    stages = build_stages(args, studies, metadata, attempt_dir, repo_dir)
    state_path = args.run_dir / "run_state.json"
    current_commit, tracked_dirty = git_info(repo_dir)
    attempt = {
        "number": attempt_number,
        "status": "running",
        "started_utc": utc_now(),
        "ended_utc": None,
        "threads": args.threads,
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
            for stage in stages:
                stage_state = state["stages"][stage.name]
                can_skip = (
                    stage_state.get("status") == "completed"
                    and outputs_exist(stage.expected_outputs)
                    and not rerun_stages.intersection(stage.dependencies)
                )
                if can_skip:
                    print(f"Skipping completed stage: {stage.name}", flush=True)
                    timing.skipped(
                        stage.name,
                        item=", ".join(str(path) for path in stage.expected_outputs),
                        message="state completed and expected outputs exist",
                    )
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

    attempt.update({"status": "completed", "ended_utc": utc_now()})
    write_json_atomic(state_path, state)
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
