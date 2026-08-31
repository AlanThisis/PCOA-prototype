from __future__ import annotations

import csv
import json
from pathlib import Path

from prepare_mmc1_16s_metadata import prepare_metadata


HEADERS = [
    "sample_accession",
    "source_study",
    "study_hmr",
    "host_species_resolved",
    "run_accession",
    "experiment_accession",
    "experiment_alias",
    "experiment_title",
    "run_alias",
    "library_name",
    "read_count",
    "base_count",
    "environment_harmonized",
    "disease_harmonized",
]


def write_source(path: Path, rows: list[list[str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(HEADERS)
        writer.writerows(rows)


def base_row(project: str, sample: str, runs: str) -> list[str]:
    return [
        sample,
        project,
        "Human",
        "Human",
        runs,
        "ERX1",
        "ordinary sample",
        "16S amplicon",
        "ordinary.fastq.gz",
        "ordinary",
        "100",
        "1000",
        "Gut",
        "healthy",
    ]


def read_csv(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def test_cleanup_filters_human_forward_16s_controls_and_umbrella_runs(
    tmp_path: Path,
) -> None:
    rows = [
        base_row("PRJ1", "SAM1", "SRR1"),
        base_row("PRJNA339931", "SAM_BAD", "SRR_BAD"),
        [
            "SAM_EYE",
            "PRJNA597168",
            "Human",
            "Human",
            "SRR_REVERSE; SRR_FORWARD",
            "ERX_R; ERX_F",
            "15 Reverse; 15 Forward",
            "Reverse; Forward",
            "15R.fastq.gz; 15F.fastq.gz",
            "15 Reverse; 15 Forward",
            "10; 20",
            "100; 200",
            "Ocular",
            "disease",
        ],
        [
            "SAM_ORAL",
            "PRJNA1158838",
            "Human",
            "Human",
            "SRR_16S; SRR_ITS",
            "ERX_16S; ERX_ITS",
            "oral; oral",
            "Cancer mucosities 16S Amplicons; Cancer mucosities ITS Amplicons",
            "a.fastq.gz; b.fastq.gz",
            "oral; oral",
            "10; 20",
            "100; 200",
            "Oral",
            "disease",
        ],
        [
            "SAM_UMBRELLA",
            "PRJ2",
            "Human",
            "Human",
            "SRR2; SRR3",
            "ERX2; ERX3",
            "specimen 2; specimen 3",
            "16S; 16S",
            "two.fastq.gz; three.fastq.gz",
            "two; three",
            "10; 20",
            "100; 200",
            "Pulmonary",
            "disease",
        ],
        [
            "SAM_CTRL",
            "PRJ2",
            "Human",
            "Human",
            "SRR_CTRL",
            "ERX_CTRL",
            "PCRNegCtrl",
            "16S",
            "PCRNegCtrl.fastq.gz",
            "control",
            "10",
            "100",
            "Pulmonary",
            "NA",
        ],
        [
            "SAM_MOUSE",
            "PRJ_MOUSE",
            "Mouse/Rat",
            "Mouse",
            "SRR_MOUSE",
            "ERX_MOUSE",
            "mouse",
            "16S",
            "mouse.fastq.gz",
            "mouse",
            "10",
            "100",
            "Gut",
            "disease",
        ],
    ]
    source = tmp_path / "source.csv"
    output = tmp_path / "cleaned"
    write_source(source, rows)

    summary = prepare_metadata(source, output)

    cleaned = read_csv(output / "cleaned_metadata.csv")
    assert [row["sample-id"] for row in cleaned] == [
        "SRR1",
        "SRR_FORWARD",
        "SRR_16S",
        "SRR2",
        "SRR3",
    ]
    assert all(row["study_hmr"] == "Human" for row in cleaned)
    assert all(row["host_species_resolved"] == "Human" for row in cleaned)
    assert cleaned[-2]["experiment_accession"] == "ERX2"
    assert cleaned[-1]["experiment_accession"] == "ERX3"
    assert summary == {
        "source_rows": 7,
        "selected_projects": 4,
        "selected_runs": 5,
        "excluded_run_records": 5,
        "excluded_projects": 2,
        "explicitly_excluded_projects": 1,
    }
    assert json.loads((output / "cleanup_summary.json").read_text()) == summary


def test_cleanup_rejects_duplicate_selected_runs(tmp_path: Path) -> None:
    source = tmp_path / "source.csv"
    write_source(
        source,
        [base_row("PRJ1", "SAM1", "SRR1"), base_row("PRJ2", "SAM2", "SRR1")],
    )

    try:
        prepare_metadata(source, tmp_path / "out")
    except ValueError as exc:
        assert "Duplicate selected run/sample IDs" in str(exc)
    else:
        raise AssertionError("Duplicate selected runs were accepted")


def test_cleanup_excludes_known_unavailable_ena_run(tmp_path: Path) -> None:
    source = tmp_path / "source.csv"
    write_source(
        source,
        [
            base_row("PRJNA580506", "SAM_BAD", "SRR10377565"),
            base_row("PRJNA580506", "SAM_GOOD", "SRR10377566"),
        ],
    )

    prepare_metadata(source, tmp_path / "out")

    cleaned = read_csv(tmp_path / "out" / "cleaned_metadata.csv")
    assert [row["sample-id"] for row in cleaned] == ["SRR10377566"]
    excluded = read_csv(tmp_path / "out" / "excluded_runs.tsv", delimiter="\t")
    assert excluded[0]["run_accession"] == "SRR10377565"
    assert "HTML error" in excluded[0]["reason"]


def test_cleanup_excludes_known_nanopore_project(tmp_path: Path) -> None:
    source = tmp_path / "source.csv"
    write_source(
        source,
        [
            base_row("PRJEB75547", "SAM_NANOPORE", "ERR13071189"),
            base_row("PRJ_ILLUMINA", "SAM_ILLUMINA", "SRR_ILLUMINA"),
        ],
    )

    summary = prepare_metadata(source, tmp_path / "out")

    assert summary["selected_runs"] == 1
    assert summary["explicitly_excluded_projects"] == 1
    excluded = read_csv(tmp_path / "out" / "excluded_projects.tsv", delimiter="\t")
    assert excluded == [
        {
            "project_accession": "PRJEB75547",
            "reason": (
                "Oxford Nanopore GridION full-length amplicon reads use quality "
                "scores that are incompatible with the current short-read Deblur workflow"
            ),
        }
    ]
