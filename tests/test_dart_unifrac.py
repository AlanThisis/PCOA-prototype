from __future__ import annotations

import zipfile
import subprocess
from pathlib import Path

import biom
import numpy as np
import pytest

import dart_unifrac
from pipeline_lib import TimingRecorder


def write_qza(path: Path, member: str, contents: bytes) -> None:
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr(f"uuid/data/{member}", contents)


def test_extract_qza_member_requires_exactly_one_match(tmp_path: Path) -> None:
    qza = tmp_path / "tree.qza"
    write_qza(qza, "tree.nwk", b"(a,b);\n")

    output = dart_unifrac.extract_qza_member(qza, "tree.nwk", tmp_path / "tree.nwk")

    assert output.read_text() == "(a,b);\n"
    with pytest.raises(ValueError, match="found 0"):
        dart_unifrac.extract_qza_member(qza, "missing", tmp_path / "missing")


def test_prepare_dart_inputs_extracts_and_validates_biom(tmp_path: Path) -> None:
    table = biom.Table(
        np.array([[1, 0], [0, 2]]),
        observation_ids=["a", "b"],
        sample_ids=["s1", "s2"],
    )
    biom_source = tmp_path / "source.biom"
    with biom.util.biom_open(biom_source, "w") as handle:
        table.to_hdf5(handle, "test")
    table_qza = tmp_path / "table.qza"
    write_qza(table_qza, "feature-table.biom", biom_source.read_bytes())
    tree_qza = tmp_path / "tree.qza"
    write_qza(tree_qza, "tree.nwk", b"(a:1,b:1);\n")

    biom_fp, tree_fp, stats = dart_unifrac.prepare_dart_inputs(
        table_qza, tree_qza, tmp_path / "work"
    )

    assert biom_fp.is_file()
    assert tree_fp.read_text() == "(a:1,b:1);\n"
    assert stats == {"sample_count": 2, "feature_count": 2}


def test_dart_command_uses_reproducible_dmh_defaults() -> None:
    command = dart_unifrac.build_dart_command(
        "/env/bin/dartunifrac",
        Path("tree.nwk"),
        Path("table.biom"),
        Path("distance.tsv"),
        threads=8,
        sketch_size=2048,
        seed=1337,
        bbits=16,
        compress=True,
    )

    assert command[command.index("-m") + 1] == "dmh"
    assert command[command.index("-s") + 1] == "2048"
    assert command[command.index("-T") + 1] == "8"
    assert command[command.index("--seed") + 1] == "1337"
    assert command[command.index("--bbits") + 1] == "16"
    assert "--pcoa" in command
    assert "--compress" in command


def test_validate_dart_parameters_requires_the_v030_dimension_count() -> None:
    dart_unifrac.validate_dart_parameters(
        sketch_size=2048, seed=1337, bbits=16, pcoa_dimensions=10
    )
    with pytest.raises(ValueError, match="10 fPCoA dimensions"):
        dart_unifrac.validate_dart_parameters(
            sketch_size=2048, seed=1337, bbits=16, pcoa_dimensions=5
        )


def test_dart_version_requires_pinned_release(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        dart_unifrac.subprocess,
        "run",
        lambda *args, **kwargs: subprocess.CompletedProcess(
            args[0], 0, stdout="dartunifrac 0.3.0\n", stderr=""
        ),
    )
    assert dart_unifrac.dart_version("dartunifrac") == "dartunifrac 0.3.0"

    monkeypatch.setattr(
        dart_unifrac.subprocess,
        "run",
        lambda *args, **kwargs: subprocess.CompletedProcess(
            args[0], 0, stdout="dartunifrac 0.4.0\n", stderr=""
        ),
    )
    with pytest.raises(RuntimeError, match="0.3.0 is required"):
        dart_unifrac.dart_version("dartunifrac")


def test_run_dart_requires_distance_and_both_pcoa_outputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def fake_run(command: list[str], *, cwd: Path, **_: object) -> None:
        requested = Path(command[command.index("-o") + 1])
        Path(str(requested) + ".zst").write_bytes(b"compressed")
        (cwd / "pcoa.txt").write_text("raw\n")
        # Deliberately omit ordination.txt to exercise the output contract.

    monkeypatch.setattr(dart_unifrac, "run_command", fake_run)
    with pytest.raises(RuntimeError, match="ordination.txt"):
        dart_unifrac.run_dart_unifrac(
            executable="dartunifrac",
            biom_fp=tmp_path / "table.biom",
            tree_fp=tmp_path / "tree.nwk",
            work_dir=tmp_path / "work",
            threads=4,
            sketch_size=2048,
            seed=1337,
            bbits=16,
            compress=True,
            timing=TimingRecorder(None, component="test"),
        )
