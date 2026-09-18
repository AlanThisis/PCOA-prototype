from __future__ import annotations

import argparse
import csv
import json
import zipfile
from pathlib import Path

import biom
import numpy as np
import pandas as pd
import pytest

import alpha_rarefaction
from pipeline_lib import TimingRecorder


def test_choose_max_depth_uses_percentile_and_includes_reference() -> None:
    depths = [100, 500, 800, 1200, 5000]
    assert alpha_rarefaction.choose_max_depth(
        depths, requested=None, percentile=60, reference_depth=1000
    ) == (1000, "nearest-rank-p60")


def test_explicit_max_depth_cannot_exceed_observed_maximum() -> None:
    with pytest.raises(ValueError, match="exceeds the maximum mapped"):
        alpha_rarefaction.choose_max_depth(
            [10, 20, 30], requested=1000, percentile=90, reference_depth=1000
        )


def test_evenly_spaced_depths_include_endpoints() -> None:
    assert alpha_rarefaction.evenly_spaced_depths(1, 1000, 5) == [
        1,
        250,
        500,
        750,
        1000,
    ]


def test_evenly_spaced_depths_matches_qiime_step_validation() -> None:
    with pytest.raises(ValueError, match="exceeds the possible steps"):
        alpha_rarefaction.evenly_spaced_depths(1, 5, 10)


def test_prepare_qiime_metadata_matches_forward_read_suffixes(tmp_path: Path) -> None:
    source = tmp_path / "metadata.csv"
    source.write_text(
        "sample-id,group\nSRR1,case\nSRR2,control\n", encoding="utf-8"
    )
    output = tmp_path / "qiime_metadata.tsv"

    alpha_rarefaction.prepare_qiime_metadata(
        source, ["SRR1_1", "SRR2_R1_001"], output
    )

    with output.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows == [
        {"sample-id": "SRR1_1", "group": "case"},
        {"sample-id": "SRR2_R1_001", "group": "control"},
    ]


def test_prepare_qiime_metadata_reports_missing_ids(tmp_path: Path) -> None:
    source = tmp_path / "metadata.tsv"
    source.write_text("sample-id\nSRR1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="missing 1 of 2"):
        alpha_rarefaction.prepare_qiime_metadata(
            source, ["SRR1_1", "SRR2_1"], tmp_path / "out.tsv"
        )


def test_run_invokes_qiime_faith_pd_and_writes_audit_outputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    mapped_table = tmp_path / "mapped.qza"
    phylogeny = tmp_path / "tree.qza"
    metadata = tmp_path / "metadata.tsv"
    for path in (mapped_table, phylogeny):
        path.write_text("fixture", encoding="utf-8")
    metadata.write_text("sample-id\tgroup\ns1\tA\ns2\tA\ns3\tB\n", encoding="utf-8")
    output_dir = tmp_path / "output"
    commands: list[tuple[list[str], str]] = []

    def fake_run_command(args: list[str], **kwargs: object) -> None:
        commands.append((args, str(kwargs["step"])))
        if kwargs["step"] == "export_mapped_table":
            export_dir = Path(args[args.index("--output-path") + 1])
            export_dir.mkdir(parents=True)
            table = biom.Table(
                np.array([[10, 500, 2000], [0, 500, 1000]]),
                observation_ids=["f1", "f2"],
                sample_ids=["s1", "s2", "s3"],
                input_is_dense=True,
            )
            with biom.util.biom_open(export_dir / "feature-table.biom", "w") as handle:
                table.to_hdf5(handle, "test")
        elif kwargs["step"] == "faith_pd_alpha_rarefaction":
            output = Path(args[args.index("--o-visualization") + 1])
            output.write_text("qzv", encoding="utf-8")
        else:
            export_dir = Path(args[args.index("--output-path") + 1])
            export_dir.mkdir(parents=True)
            pd.DataFrame(
                {
                    "sample-id": ["s1", "s2", "s3"],
                    "depth-1_iter-1": [1.0, 2.0, 3.0],
                    "depth-1_iter-2": [1.2, 2.2, 3.2],
                    "depth-3000_iter-1": [np.nan, np.nan, 8.0],
                    "depth-3000_iter-2": [np.nan, np.nan, 8.2],
                }
            ).to_csv(export_dir / "faith_pd.csv", index=False)

    monkeypatch.setattr(alpha_rarefaction, "run_command", fake_run_command)
    monkeypatch.setattr(alpha_rarefaction, "resolve_executable", lambda _: "/env/bin/qiime")
    args = argparse.Namespace(
        mapped_table=mapped_table,
        phylogeny=phylogeny,
        output_dir=output_dir,
        metadata=metadata,
        max_depth=None,
        max_depth_percentile=90.0,
        reference_depth=1000,
        min_depth=1,
        steps=10,
        iterations=4,
        overwrite=False,
    )

    assert alpha_rarefaction.run(args, TimingRecorder(None, "test")) == 0

    alpha_command = next(command for command, step in commands if step == "faith_pd_alpha_rarefaction")
    assert alpha_command[1:3] == ["diversity", "alpha-rarefaction"]
    assert alpha_command[alpha_command.index("--p-metrics") + 1] == "faith_pd"
    assert alpha_command[alpha_command.index("--p-max-depth") + 1] == "3000"
    assert alpha_command[alpha_command.index("--p-iterations") + 1] == "4"
    qiime_metadata = Path(
        alpha_command[alpha_command.index("--m-metadata-file") + 1]
    )
    assert qiime_metadata == output_dir / "qiime_metadata.tsv"
    with qiime_metadata.open(newline="", encoding="utf-8") as handle:
        assert [row["sample-id"] for row in csv.DictReader(handle, delimiter="\t")] == [
            "s1",
            "s2",
            "s3",
        ]

    with (output_dir / "sample_retention_by_depth.tsv").open(
        newline="", encoding="utf-8"
    ) as handle:
        retention = list(csv.DictReader(handle, delimiter="\t"))
    assert retention[-1]["depth"] == "3000"
    assert retention[-1]["retained_samples"] == "1"

    summary = json.loads(
        (output_dir / "alpha_rarefaction_summary.json").read_text(encoding="utf-8")
    )
    assert summary["metric"] == "faith_pd"
    assert summary["samples_retained_at_reference_depth"] == 2
    assert summary["positive_sample_count"] == 3
    assert (output_dir / "faith_pd_alpha_rarefaction.png").is_file()
    curve = pd.read_csv(output_dir / "faith_pd_curve_summary.tsv", sep="\t")
    assert curve["retained_samples"].tolist() == [3, 1]
    assert curve["faith_pd_median"].tolist() == pytest.approx([2.1, 8.1])


def test_plot_curve_comparison_requires_matching_depths(tmp_path: Path) -> None:
    first = pd.DataFrame(
        {
            "depth": [1, 1000],
            "faith_pd_median": [1.0, 3.0],
            "retained_percent": [100.0, 60.0],
        }
    )
    second = pd.DataFrame(
        {
            "depth": [1, 500],
            "faith_pd_median": [0.8, 2.0],
            "retained_percent": [100.0, 70.0],
        }
    )
    with pytest.raises(ValueError, match="different depth grid"):
        alpha_rarefaction.plot_curve_comparison(
            [("Full", first), ("10%", second)], tmp_path / "combined.png", 1000
        )


def test_plot_curve_comparison_writes_png(tmp_path: Path) -> None:
    curve = pd.DataFrame(
        {
            "depth": [1, 500, 1000],
            "faith_pd_median": [1.0, 2.0, 2.4],
            "retained_percent": [100.0, 80.0, 60.0],
        }
    )
    output = tmp_path / "combined.png"
    alpha_rarefaction.plot_curve_comparison(
        [("Full", curve), ("10%", curve.copy())], output, 1000
    )
    assert output.is_file()


def test_load_and_plot_group_summary_from_qzv(tmp_path: Path) -> None:
    qzv = tmp_path / "alpha.qzv"
    payload = {
        "columns": [
            "environment_harmonized",
            "_alpha_rarefaction_depth_column_",
            "25%",
            "50%",
            "75%",
            "count",
        ],
        "index": [0, 1, 2, 3],
        "data": [
            ["Gut", 1, 1.0, 2.0, 3.0, 12],
            ["Gut", 1000, 3.0, 4.0, 5.0, 10],
            ["Oral", 1, 0.5, 1.5, 2.5, 8],
            ["Oral", 1000, 2.0, 3.0, 4.0, 6],
        ],
    }
    with zipfile.ZipFile(qzv, "w") as archive:
        archive.writestr(
            "uuid/data/faith_pd-environment_harmonized.jsonp",
            f"load_data('faith_pd', 'environment_harmonized',{json.dumps(payload)})",
        )

    frame = alpha_rarefaction.load_group_summary_from_qzv(
        qzv, "environment_harmonized"
    )
    assert frame["group"].tolist() == ["Gut", "Gut", "Oral", "Oral"]
    assert frame["faith_pd_median"].tolist() == [2.0, 4.0, 1.5, 3.0]

    output = tmp_path / "grouped.png"
    alpha_rarefaction.plot_group_summary(
        frame, output, "environment_harmonized", 1000
    )
    assert output.is_file()


def test_existing_qzv_requires_overwrite(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    mapped_table = tmp_path / "mapped.qza"
    phylogeny = tmp_path / "tree.qza"
    mapped_table.write_text("fixture", encoding="utf-8")
    phylogeny.write_text("fixture", encoding="utf-8")
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    (output_dir / "faith_pd_alpha_rarefaction.qzv").write_text("old", encoding="utf-8")

    monkeypatch.setattr(
        alpha_rarefaction,
        "run_command",
        lambda *args, **kwargs: pytest.fail("existing output should fail before export"),
    )
    monkeypatch.setattr(alpha_rarefaction, "resolve_executable", lambda _: "/env/bin/qiime")
    args = argparse.Namespace(
        mapped_table=mapped_table,
        phylogeny=phylogeny,
        output_dir=output_dir,
        metadata=None,
        max_depth=10,
        max_depth_percentile=90.0,
        reference_depth=5,
        min_depth=1,
        steps=2,
        iterations=1,
        overwrite=False,
    )

    with pytest.raises(FileExistsError, match="--overwrite"):
        alpha_rarefaction.run(args, TimingRecorder(None, "test"))
