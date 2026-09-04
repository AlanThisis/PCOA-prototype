from __future__ import annotations

import csv
import json
from pathlib import Path

from prepare_mmc_sample_labels import prepare_sample_labels


SAMPLE_FIELDS = [
    "record_id",
    "study_accession",
    "sample_accession",
    "tier",
    "raw_disease_value",
    "disease_source_field",
    "case_control",
    "mondo_id",
    "mondo_label",
    "doid",
    "mesh",
    "host_age",
    "host_sex",
    "mapping_method",
    "study_mesh_heading",
    "custom_attrs_available",
    "name_derived",
    "unexpected_annotation",
    "run_accession",
    "library_strategy",
    "instrument_platform",
]


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def sample(**updates: str) -> dict[str, str]:
    row = {field: "" for field in SAMPLE_FIELDS}
    row.update(
        {
            "record_id": "MMC2-00001",
            "study_accession": "PRJNA1",
            "sample_accession": "SAMN1",
            "tier": "1",
            "case_control": "case",
            "mondo_label": "disease one",
            "run_accession": "SRR1",
            "library_strategy": "AMPLICON",
            "instrument_platform": "ILLUMINA",
        }
    )
    row.update(updates)
    return row


def test_prepares_lossless_download_and_illumina_analysis_views(tmp_path: Path) -> None:
    samples = tmp_path / "samples.tsv"
    studies = tmp_path / "studies.tsv"
    output = tmp_path / "out"
    write_tsv(
        samples,
        SAMPLE_FIELDS,
        [
            sample(),
            sample(
                record_id="MMC2-00002",
                study_accession="SRP1",
                case_control="control",
                mondo_label="",
            ),
            sample(
                record_id="MMC2-00003",
                study_accession="PRJEB2",
                sample_accession="SAMEA2",
                run_accession="ERR2",
                instrument_platform="LS454",
            ),
            sample(
                record_id="MMC2-00004",
                study_accession="PRJNA3",
                sample_accession="SAMN3",
                run_accession="SRR3",
                library_strategy="WGS",
            ),
            sample(
                record_id="MMC2-00005",
                study_accession="PRJNA4",
                sample_accession="SAMN4",
                run_accession="SRR4",
                tier="2",
            ),
        ],
    )
    write_tsv(
        studies,
        ["record_id", "body_site", "disease_group", "title"],
        [
            {
                "record_id": f"MMC2-{number:05d}",
                "body_site": "Gut",
                "disease_group": "Cancer",
                "title": f"Study {number}",
            }
            for number in range(1, 6)
        ],
    )

    summary = prepare_sample_labels(samples, studies, output, {"1"})

    assert summary["unique_amplicon_runs"] == 2
    assert summary["illumina_runs"] == 1
    assert summary["runs_with_multiple_citations"] == 1
    assert summary["conflicting_duplicate_groups"] == {
        "case_control": 1,
        "mondo_label": 1,
        "record_id": 1,
        "study_accession": 1,
    }
    with (output / "all_amplicon_metadata.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["sample-id"] for row in rows] == ["ERR2", "SRR1"]
    assert rows[1]["source_study"] == "PRJNA1"
    with (output / "illumina_pipeline_metadata.tsv").open() as handle:
        illumina = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["sample-id"] for row in illumina] == ["SRR1"]
    assert json.loads((output / "preparation_summary.json").read_text()) == summary


def test_can_prepare_combined_tiers(tmp_path: Path) -> None:
    samples = tmp_path / "samples.tsv"
    studies = tmp_path / "studies.tsv"
    output = tmp_path / "out"
    write_tsv(
        samples,
        SAMPLE_FIELDS,
        [sample(), sample(run_accession="SRR2", sample_accession="SAMN2", tier="2")],
    )
    write_tsv(
        studies,
        ["record_id", "body_site", "disease_group", "title"],
        [
            {
                "record_id": "MMC2-00001",
                "body_site": "Gut",
                "disease_group": "Cancer",
                "title": "Study",
            }
        ],
    )

    summary = prepare_sample_labels(samples, studies, output, {"1", "2"})

    assert summary["unique_amplicon_runs"] == 2
    assert summary["tiers"] == ["1", "2"]
