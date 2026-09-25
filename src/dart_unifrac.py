#!/usr/bin/env python3
"""Helpers for running DartUniFrac on QIIME 2 rarefied tables and trees."""
from __future__ import annotations

import shutil
import subprocess
import zipfile
from pathlib import Path

import biom
import skbio

from pipeline_lib import TimingRecorder, run_command


DART_METHOD = "dmh"
DART_PCOA_DIMENSIONS = 10
DART_REQUIRED_VERSION = "0.3.0"


def dart_version(executable: str, required: str = DART_REQUIRED_VERSION) -> str:
    """Return the DART version string and enforce the validated release."""
    completed = subprocess.run(
        [executable, "--version"],
        check=True,
        capture_output=True,
        text=True,
    )
    version = " ".join(
        part.strip() for part in (completed.stdout, completed.stderr) if part.strip()
    )
    if required not in version:
        raise RuntimeError(
            f"DartUniFrac {required} is required; executable reported {version!r}"
        )
    return version


def extract_qza_member(qza: Path, member_suffix: str, destination: Path) -> Path:
    """Extract exactly one QZA data member matching *member_suffix*."""
    with zipfile.ZipFile(qza) as archive:
        matches = [
            name
            for name in archive.namelist()
            if name.endswith(f"/data/{member_suffix}") and not name.endswith("/")
        ]
        if len(matches) != 1:
            raise ValueError(
                f"Expected exactly one data/{member_suffix} in {qza}; "
                f"found {len(matches)}"
            )
        destination.parent.mkdir(parents=True, exist_ok=True)
        with archive.open(matches[0]) as source, destination.open("wb") as target:
            shutil.copyfileobj(source, target)
    if destination.stat().st_size == 0:
        raise ValueError(f"Extracted an empty data/{member_suffix} from {qza}")
    return destination


def prepare_dart_inputs(
    rarefied_table_qza: Path,
    tree_qza: Path,
    work_dir: Path,
) -> tuple[Path, Path, dict[str, int]]:
    """Extract and validate the BIOM table and Newick tree required by DART."""
    input_dir = work_dir / "inputs"
    biom_fp = extract_qza_member(
        rarefied_table_qza, "feature-table.biom", input_dir / "feature-table.biom"
    )
    tree_fp = extract_qza_member(tree_qza, "tree.nwk", input_dir / "tree.nwk")
    table = biom.load_table(str(biom_fp))
    if table.shape[1] < 2:
        raise ValueError(
            "DartUniFrac PCoA requires at least two samples after rarefaction; "
            f"found {table.shape[1]}"
        )
    return biom_fp, tree_fp, {
        "sample_count": int(table.shape[1]),
        "feature_count": int(table.shape[0]),
    }


def build_dart_command(
    executable: str,
    tree_fp: Path,
    biom_fp: Path,
    output_fp: Path,
    *,
    threads: int,
    sketch_size: int,
    seed: int,
    bbits: int,
    compress: bool,
    weighted: bool = False,
) -> list[str]:
    command = [
        executable,
        "-t",
        str(tree_fp),
        "-b",
        str(biom_fp),
        "-m",
        DART_METHOD,
        "-s",
        str(sketch_size),
        "-o",
        str(output_fp),
        "-T",
        str(threads),
        "--seed",
        str(seed),
        "--bbits",
        str(bbits),
        "--pcoa",
    ]
    if weighted:
        command.append("--weighted")
    if compress:
        command.append("--compress")
    return command


def validate_dart_parameters(
    *, sketch_size: int, seed: int, bbits: int, pcoa_dimensions: int
) -> None:
    if sketch_size <= 0:
        raise ValueError("--dart-sketch-size must be greater than zero")
    if seed < 0:
        raise ValueError("--dart-seed cannot be negative")
    if bbits not in {16, 32, 64}:
        raise ValueError("--dart-bbits must be one of 16, 32, or 64")
    if pcoa_dimensions != DART_PCOA_DIMENSIONS:
        raise ValueError(
            "DartUniFrac v0.3.0 emits 10 fPCoA dimensions; "
            "set --pcoa-dimensions 10"
        )


def _find_distance_output(requested: Path, compress: bool) -> Path:
    candidates = [requested]
    if compress:
        candidates = [
            requested.with_suffix(requested.suffix + ".zst"),
            Path(str(requested) + ".zst"),
            requested,
        ]
    for candidate in candidates:
        if candidate.is_file() and candidate.stat().st_size > 0:
            return candidate
    raise RuntimeError(
        "DartUniFrac completed without producing its requested distance matrix: "
        + ", ".join(str(path) for path in candidates)
    )


def run_dart_unifrac(
    *,
    executable: str,
    biom_fp: Path,
    tree_fp: Path,
    work_dir: Path,
    threads: int,
    sketch_size: int,
    seed: int,
    bbits: int,
    compress: bool,
    timing: TimingRecorder,
    weighted: bool = False,
) -> dict[str, object]:
    """Run DART DMH UniFrac plus fPCoA and validate its stable outputs."""
    work_dir.mkdir(parents=True, exist_ok=True)
    metric = "weighted" if weighted else "unweighted"
    distance_requested = work_dir / f"{metric}_unifrac_distance_matrix.tsv"
    for candidate in (
        distance_requested,
        Path(str(distance_requested) + ".zst"),
        work_dir / "pcoa.txt",
        work_dir / "ordination.txt",
    ):
        candidate.unlink(missing_ok=True)
    command = build_dart_command(
        executable,
        tree_fp,
        biom_fp,
        distance_requested,
        threads=threads,
        sketch_size=sketch_size,
        seed=seed,
        bbits=bbits,
        compress=compress,
        weighted=weighted,
    )
    run_command(
        command,
        cwd=work_dir,
        timing=timing,
        step=f"dart_{metric}_unifrac_fpcoa",
        item=str(distance_requested),
    )
    distance_fp = _find_distance_output(distance_requested, compress)
    ordination_fp = work_dir / "ordination.txt"
    pcoa_fp = work_dir / "pcoa.txt"
    for output in (ordination_fp, pcoa_fp):
        if not output.is_file() or output.stat().st_size == 0:
            raise RuntimeError(f"DartUniFrac did not produce a valid {output.name}")
    ordination = skbio.io.read(
        str(ordination_fp),
        format="ordination",
        into=skbio.stats.ordination.OrdinationResults,
    )
    if ordination.samples.shape[1] < 2:
        raise RuntimeError("DartUniFrac ordination contains fewer than two axes")
    return {
        "distance_path": distance_fp,
        "ordination_path": ordination_fp,
        "raw_pcoa_path": pcoa_fp,
        "sample_count": int(ordination.samples.shape[0]),
        "dimensions": int(ordination.samples.shape[1]),
        "compressed": distance_fp.suffix == ".zst",
    }
