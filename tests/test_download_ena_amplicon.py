from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

import download_ena_amplicon
from download_ena_amplicon import (
    DownloadSpec,
    download_one,
    file_report_params,
    format_elapsed,
    ftp_url,
    parse_file_report,
    select_forward_download,
    transfer_once,
    https_url,
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


def test_https_url_strips_whitespace_and_trailing_slash() -> None:
    assert https_url(" ftp.sra.ebi.ac.uk/run.fastq.gz/ \n") == (
        "https://ftp.sra.ebi.ac.uk/run.fastq.gz"
    )


def test_ftp_url_restores_ena_native_transport() -> None:
    assert ftp_url("https://ftp.sra.ebi.ac.uk/run.fastq.gz") == (
        "ftp://ftp.sra.ebi.ac.uk/run.fastq.gz"
    )


class FakeDownloadResponse:
    def __init__(
        self,
        status_code: int,
        content: bytes = b"",
        content_type: str = "application/x-gzip",
    ) -> None:
        self.status_code = status_code
        self.content = content
        self.closed = False
        self.headers = {
            "Content-Type": content_type,
            "Content-Length": str(len(content)),
        }

    def close(self) -> None:
        self.closed = True

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise download_ena_amplicon.requests.HTTPError(
                f"HTTP {self.status_code}"
            )

    def iter_content(self, chunk_size: int):
        del chunk_size
        yield self.content

    def __enter__(self):
        return self

    def __exit__(self, *_: object) -> None:
        self.close()


class FakeFtpResponse:
    def __init__(self, content: bytes) -> None:
        self.content = content
        self.offset = 0
        self.headers = {"Content-Length": str(len(content))}

    def read(self, _: int) -> bytes:
        if self.offset:
            return b""
        self.offset = len(self.content)
        return self.content

    def __enter__(self):
        return self

    def __exit__(self, *_: object) -> None:
        return None


def make_spec(content: bytes) -> DownloadSpec:
    return DownloadSpec(
        project_accession="PRJ1",
        run_accession="SRR1",
        library_layout="PAIRED",
        library_strategy="AMPLICON",
        url="https://ftp.sra.ebi.ac.uk/SRR1_1.fastq.gz",
        md5=hashlib.md5(content).hexdigest(),
        bytes=len(content),
        output_filename="SRR1_1.fastq.gz",
    )


def test_transfer_restarts_without_range_after_416(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    content = b"complete-fastq"
    spec = make_spec(content)
    destination = tmp_path / spec.output_filename
    part_path = tmp_path / f"{spec.output_filename}.part"
    part_path.write_bytes(b"partial")
    responses = [FakeDownloadResponse(416), FakeDownloadResponse(200, content)]
    request_headers: list[dict[str, str]] = []

    def fake_get(*_: object, **kwargs: object) -> FakeDownloadResponse:
        request_headers.append(dict(kwargs["headers"]))
        return responses.pop(0)

    monkeypatch.setattr(download_ena_amplicon.requests, "get", fake_get)

    downloaded = transfer_once(spec, destination, timeout=60)

    assert downloaded == len(content)
    assert destination.read_bytes() == content
    assert not part_path.exists()
    assert request_headers[0]["Range"] == "bytes=7-"
    assert "Range" not in request_headers[1]


def test_transfer_discards_stale_partial_and_restarts_after_html(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    content = b"complete-fastq"
    spec = make_spec(content)
    destination = tmp_path / spec.output_filename
    part_path = destination.with_name(destination.name + ".part")
    part_path.write_bytes(b"<html>")
    html = b"<html>not a fastq</html>"
    responses = [
        FakeDownloadResponse(200, html, "text/html"),
        FakeDownloadResponse(200, content),
    ]
    request_headers: list[dict[str, str]] = []

    def fake_get(*_: object, **kwargs: object) -> FakeDownloadResponse:
        request_headers.append(dict(kwargs["headers"]))
        return responses.pop(0)

    monkeypatch.setattr(download_ena_amplicon.requests, "get", fake_get)

    downloaded = transfer_once(spec, destination, timeout=60)

    assert downloaded == len(content)
    assert destination.read_bytes() == content
    assert not part_path.exists()
    assert request_headers[0]["Range"] == "bytes=6-"
    assert "Range" not in request_headers[1]


def test_transfer_falls_back_to_ftp_after_clean_https_returns_html(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    content = b"complete-fastq"
    spec = make_spec(content)
    destination = tmp_path / spec.output_filename
    response = FakeDownloadResponse(200, b"<html>error</html>", "text/html")
    monkeypatch.setattr(
        download_ena_amplicon.requests,
        "get",
        lambda *_args, **_kwargs: response,
    )
    ftp_requests: list[str] = []

    def fake_urlopen(url: str, **_: object) -> FakeFtpResponse:
        ftp_requests.append(url)
        return FakeFtpResponse(content)

    monkeypatch.setattr(download_ena_amplicon, "urlopen", fake_urlopen)

    downloaded = transfer_once(spec, destination, timeout=60)

    assert downloaded == len(content)
    assert destination.read_bytes() == content
    assert not destination.with_name(destination.name + ".part").exists()
    assert ftp_requests == ["ftp://ftp.sra.ebi.ac.uk/SRR1_1.fastq.gz"]


def test_valid_final_file_removes_stale_partial(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    content = b"complete-fastq"
    spec = make_spec(content)
    destination = tmp_path / spec.output_filename
    part_path = tmp_path / f"{spec.output_filename}.part"
    destination.write_bytes(content)
    part_path.write_bytes(b"stale")
    monkeypatch.setattr(
        download_ena_amplicon.requests,
        "get",
        lambda *_args, **_kwargs: pytest.fail("download should be skipped"),
    )

    result = download_one(spec, tmp_path, retries=3, timeout=60)

    assert result.status == "skipped"
    assert destination.read_bytes() == content
    assert not part_path.exists()
