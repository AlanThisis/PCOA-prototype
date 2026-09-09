from __future__ import annotations

import json
from pathlib import Path

import validate_pipeline_run


def write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


def test_validate_complete_balanced_run(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    results = run_dir / "results"
    results.mkdir(parents=True)
    qza = run_dir / "work" / "qiime2" / "unweighted_unifrac_distance_matrix.qza"
    qza.parent.mkdir(parents=True)
    qza.write_bytes(b"artifact")
    write_json(
        run_dir / "run_state.json",
        {
            "attempts": [{"status": "completed_with_exclusions"}],
            "stages": {"deblur:balanced": {"status": "completed"}, "unifrac": {"status": "completed"}},
        },
    )
    write_json(run_dir / "run_manifest.json", {"color_by": []})
    for filename in (
        "pcoa_coordinates_unweighted_unifrac.txt",
        "pcoa_plot_unweighted_unifrac.png",
        "pipeline_summary.json",
    ):
        (results / filename).write_bytes(b"result")
    (results / "sample_processing_status.tsv").write_text(
        "sample_id\tstudy\tfastq_path\tstatus\tnode_id\tmessage\n"
        "A\tstudy\t/a\tcompleted\tshard-0\t\n"
    )
    write_json(
        results / "deblur_processing_summary.json",
        {
            "counts": {
                "expected": 1,
                "completed": 1,
                "completed_after_sanitation": 0,
                "zero_features": 0,
                "failed": 0,
                "excluded_assay": 0,
            },
            "failed_sample_limit": 0,
        },
    )
    write_json(
        results / "analysis_summary.json",
        {
            "mapped_samples": 1,
            "rarefied_samples": 1,
            "distance_tsv": {"exported": False, "qza_path": str(qza)},
        },
    )

    assert validate_pipeline_run.validate_run(run_dir) == []


def test_validator_reports_incomplete_stage_and_missing_qza(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    results = run_dir / "results"
    results.mkdir(parents=True)
    write_json(
        run_dir / "run_state.json",
        {"attempts": [{"status": "failed"}], "stages": {"unifrac": {"status": "failed"}}},
    )
    write_json(run_dir / "run_manifest.json", {"color_by": []})
    for filename in (
        "pcoa_coordinates_unweighted_unifrac.txt",
        "pcoa_plot_unweighted_unifrac.png",
        "pipeline_summary.json",
    ):
        (results / filename).write_bytes(b"result")
    write_json(
        results / "analysis_summary.json",
        {
            "mapped_samples": 1,
            "rarefied_samples": 1,
            "distance_tsv": {
                "exported": False,
                "qza_path": str(run_dir / "missing.qza"),
            },
        },
    )

    problems = validate_pipeline_run.validate_run(run_dir)

    assert "latest pipeline attempt is not complete" in problems
    assert "stage unifrac has status failed" in problems
    assert "UniFrac distance QZA is missing or empty" in problems
