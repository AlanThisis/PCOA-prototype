import argparse
from pathlib import Path

import pytest

from plot_pcoa import (
    UNKNOWN_LABEL,
    deterministic_color,
    load_coordinates,
    load_id_to_label,
    resolve_style,
)


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
