from __future__ import annotations

from download_ena_amplicon import (
    file_report_params,
    format_elapsed,
    parse_file_report,
    select_forward_download,
)


def test_file_report_requests_run_download_fields_for_project() -> None:
    params = file_report_params("PRJEB46665")

    assert params["accession"] == "PRJEB46665"
    assert params["result"] == "read_run"
    assert params["format"] == "tsv"
    assert "fastq_ftp" in params["fields"]


def test_parse_file_report_rejects_ena_error_response() -> None:
    malformed_response = '[\n{"message":"Query failed"}'

    try:
        parse_file_report(malformed_response)
    except RuntimeError as exc:
        assert "missing fields" in str(exc)
    else:
        raise AssertionError("Malformed ENA response was accepted")


def test_select_forward_download_chooses_only_first_mate() -> None:
    record = {
        "run_accession": "ERR100",
        "library_layout": "PAIRED",
        "library_strategy": "AMPLICON",
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
