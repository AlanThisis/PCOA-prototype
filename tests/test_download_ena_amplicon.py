from __future__ import annotations

from download_ena_amplicon import (
    format_elapsed,
    portal_params,
    select_forward_download,
)


def test_portal_query_filters_primary_project_to_amplicon_runs() -> None:
    params = portal_params("PRJEB46665")

    assert params["query"] == (
        'study_accession="PRJEB46665" AND library_strategy="AMPLICON"'
    )
    assert params["limit"] == "0"


def test_portal_query_uses_secondary_study_accession_for_erp() -> None:
    params = portal_params("ERP005534")

    assert params["query"] == (
        'secondary_study_accession="ERP005534" AND library_strategy="AMPLICON"'
    )


def test_select_forward_download_chooses_only_first_mate() -> None:
    record = {
        "run_accession": "ERR100",
        "library_layout": "PAIRED",
        "library_strategy": "AMPLICON",
        "target_gene": "16S rRNA",
        "fastq_ftp": "ftp.sra.ebi.ac.uk/ERR100_1.fastq.gz;ftp.sra.ebi.ac.uk/ERR100_2.fastq.gz",
        "fastq_md5": "aaa;bbb",
        "fastq_bytes": "100;120",
    }

    spec = select_forward_download("ERP1", record)

    assert spec is not None
    assert spec.url == "https://ftp.sra.ebi.ac.uk/ERR100_1.fastq.gz"
    assert spec.md5 == "aaa"
    assert spec.bytes == 100
    assert spec.output_filename == "ERR100_1.fastq.gz"


def test_select_forward_download_normalizes_single_end_filename() -> None:
    record = {
        "run_accession": "ERR200",
        "library_layout": "SINGLE",
        "library_strategy": "AMPLICON",
        "target_gene": "",
        "fastq_ftp": "ftp.sra.ebi.ac.uk/ERR200.fastq.gz",
        "fastq_md5": "abc123",
        "fastq_bytes": "99",
    }

    spec = select_forward_download("ERP1", record)

    assert spec is not None
    assert spec.output_filename == "ERR200_1.fastq.gz"


def test_select_forward_download_rejects_wgs_and_reverse_only_records() -> None:
    wgs_record = {
        "run_accession": "ERR300",
        "library_layout": "PAIRED",
        "library_strategy": "WGS",
        "fastq_ftp": "ftp.sra.ebi.ac.uk/ERR300_1.fastq.gz",
        "fastq_md5": "abc",
        "fastq_bytes": "100",
    }
    reverse_only_record = {
        "run_accession": "ERR400",
        "library_layout": "PAIRED",
        "library_strategy": "AMPLICON",
        "fastq_ftp": "ftp.sra.ebi.ac.uk/ERR400_2.fastq.gz",
        "fastq_md5": "abc",
        "fastq_bytes": "100",
    }

    assert select_forward_download("ERP1", wgs_record) is None
    assert select_forward_download("ERP1", reverse_only_record) is None


def test_format_elapsed_supports_accession_runtime_over_24_hours() -> None:
    assert format_elapsed(90061) == "25:01:01"
