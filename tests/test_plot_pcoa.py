import argparse
from pathlib import Path

import pytest

from plot_pcoa import (
    UNKNOWN_LABEL,
    deterministic_color,
    load_coordinates,
    load_id_to_label,
    resolve_style,
    run,
)
from pipeline_lib import TimingRecorder


def style_args(**overrides):
    values = {
        "legend": "auto",
        "legend_columns": None,
        "point_size": None,
        "alpha": None,
        "fig_width": None,
        "fig_height": None,
    }
    values.update(overrides)
    return argparse.Namespace(**values)


def test_adaptive_style_scales_points_and_legend():
    small = resolve_style(200, 5, style_args())
    large = resolve_style(20_000, 100, style_args())

    assert small.point_size > large.point_size
    assert small.legend == "show"
    assert large.legend == "separate"


def test_deterministic_colors_do_not_depend_on_category_order():
    assert deterministic_color("Gut") == deterministic_color("Gut")
    assert deterministic_color("Gut") != deterministic_color("Oral")
    assert deterministic_color(UNKNOWN_LABEL) == "#C9CDD3"


def test_load_compact_coordinates_and_variance(tmp_path: Path):
    coordinates = tmp_path / "pcoa_coordinates_pc1_pc3.tsv"
    coordinates.write_text("sample-id\tPC1\tPC2\tPC3\nS1\t1\t2\t3\n", encoding="utf-8")
    variance = tmp_path / "pcoa_variance_pc1_pc3.tsv"
    variance.write_text(
        "axis\tproportion_explained\nPC1\t0.4\nPC2\t0.3\nPC3\t0.2\n",
        encoding="utf-8",
    )

    frame, proportions = load_coordinates(coordinates, (1, 2))

    assert frame.to_dict("records") == [{"sample-id": "S1", "x": 1, "y": 2}]
    assert proportions == pytest.approx((0.4, 0.3))


def test_blank_metadata_labels_become_unknown(tmp_path: Path):
    metadata = tmp_path / "metadata.csv"
    metadata.write_text(
        "sample-id,environment_harmonized\nS1,Gut\nS2,\n", encoding="utf-8"
    )

    labels = load_id_to_label(metadata, "environment_harmonized")

    assert labels == {"S1": "Gut", "S2": UNKNOWN_LABEL}


def test_metadata_forward_read_suffixes_are_normalized(tmp_path: Path):
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "sample-id\tenvironment_harmonized\nDRR100552_1\tGut\nSRR20_R1_001\tOral\n",
        encoding="utf-8",
    )

    labels = load_id_to_label(metadata, "environment_harmonized")

    assert labels == {"DRR100552": "Gut", "SRR20": "Oral"}


def test_conflicting_normalized_metadata_ids_are_rejected(tmp_path: Path):
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "sample-id\tenvironment_harmonized\nS1\tGut\nS1_1\tOral\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="Conflicting metadata labels"):
        load_id_to_label(metadata, "environment_harmonized")


def test_plot_rejects_zero_metadata_id_matches(tmp_path: Path):
    coordinates = tmp_path / "coordinates.tsv"
    coordinates.write_text(
        "sample-id\tPC1\tPC2\nS1\t0\t1\nS2\t1\t0\n", encoding="utf-8"
    )
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "sample-id\tenvironment_harmonized\nOTHER_1\tGut\n", encoding="utf-8"
    )
    args = argparse.Namespace(
        pc=[1, 2], pcoa=coordinates, variance=None, metadata=metadata,
        color_by="environment_harmonized", palette_file=None,
        out=tmp_path / "plot.png", color_key=None, report=None,
        title=None, legend_title=None, dpi=72, **style_args().__dict__,
    )

    with pytest.raises(ValueError, match="all-Unknown"):
        run(args, TimingRecorder(tmp_path / "timings.tsv", component="test"))
