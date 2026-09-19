from __future__ import annotations

import argparse
from pathlib import Path

import biom
import biom.util
import numpy as np

import merge_biom
from pipeline_lib import TimingRecorder


def write_deblur_output(path: Path) -> None:
    path.mkdir(parents=True)
    table = biom.Table(
        np.array([[5, 3], [0, 0]]),
        observation_ids=["observed-feature", "empty-feature"],
        sample_ids=["sample-a", "sample-b"],
        input_is_dense=True,
    )
    with biom.util.biom_open(str(path / "all.biom"), "w") as handle:
        table.to_hdf5(handle, "test")
    (path / "all.seqs.fa").write_text(
        ">observed-feature\nACGT\n>empty-feature\nTGCA\n",
        encoding="utf-8",
    )


def test_merge_prunes_empty_features_without_changing_samples(tmp_path: Path) -> None:
    source = tmp_path / "source"
    output = tmp_path / "merged"
    write_deblur_output(source)
    args = argparse.Namespace(
        deblur_dirs=[source],
        out_dir=output,
        skip_empty=False,
        run_state=None,
    )

    assert merge_biom.run(args, TimingRecorder(None, "test")) == 0

    merged = biom.load_table(str(output / "all.biom"))
    assert list(merged.ids(axis="sample")) == ["sample-a", "sample-b"]
    assert list(merged.ids(axis="observation")) == ["observed-feature"]
    assert merged.sum(axis="sample").tolist() == [5, 3]
    assert (output / "all.seqs.fa").read_text() == (
        ">observed-feature\nACGT\n"
    )
