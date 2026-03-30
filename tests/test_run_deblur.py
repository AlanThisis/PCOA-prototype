from __future__ import annotations

from pathlib import Path

import pytest

from run_deblur import discover_inputs


def test_discover_inputs_finds_all_matching_forward_fastqs(tmp_path: Path) -> None:
    for name in ("SRR100_1.fastq.gz", "SRR101_1.fastq.gz", "SRR999_2.fastq.gz"):
        (tmp_path / name).write_bytes(b"non-empty")

    found = discover_inputs(tmp_path)

    assert [fp.name for fp in found] == ["SRR100_1.fastq.gz", "SRR101_1.fastq.gz"]


def test_discover_inputs_raises_for_empty_matching_file(tmp_path: Path) -> None:
    (tmp_path / "SRR100_1.fastq.gz").write_bytes(b"")

    with pytest.raises(RuntimeError, match="empty"):
        discover_inputs(tmp_path)


def test_discover_inputs_raises_when_no_matching_files(tmp_path: Path) -> None:
    (tmp_path / "SRR100_2.fastq.gz").write_bytes(b"non-empty")

    with pytest.raises(FileNotFoundError, match=r"\*_1.fastq.gz"):
        discover_inputs(tmp_path)
