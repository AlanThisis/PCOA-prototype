from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest

import download_ena_amplicon
from download_ena_amplicon import (
    DownloadSpec,
    DownloadResult,
    download_one,
    fetch_manifest,
    file_report_params,
    format_elapsed,
    ftp_url,
    load_approved_runs,
    parse_file_report,
    quarantine_unapproved_fastqs,
    select_forward_download,
    select_approved_manifest,
    transfer_once,
    validate_downloads,
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


def test_fetch_manifest_retries_a_transient_empty_report(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    header = "\t".join(download_ena_amplicon.FILE_REPORT_FIELDS)
    complete_report = "\n".join(
        [
            header,
            "\t".join(
                [
                    "SRR1",
                    "SINGLE",
                    "AMPLICON",
                    "ftp.sra.ebi.ac.uk/SRR1.fastq.gz",
                    "abc",
                    "100",
                ]
            ),
        ]
    )
    responses = [header + "\n", complete_report]

    class FakeManifestResponse:
        def __init__(self, text: str) -> None:
            self.text = text

        def raise_for_status(self) -> None:
            return None

    def fake_get(*_: object, **__: object) -> FakeManifestResponse:
        return FakeManifestResponse(responses.pop(0))

    monkeypatch.setattr(download_ena_amplicon.requests, "get", fake_get)
    monkeypatch.setattr(download_ena_amplicon.time, "sleep", lambda _: None)

    specs, unavailable = fetch_manifest("PRJ1", timeout=1)

    assert len(specs) == 1
    assert specs[0].run_accession == "SRR1"
    assert unavailable == []
    assert responses == []


def test_exact_run_manifest_selects_only_allowlisted_project_runs(
    tmp_path: Path,
) -> None:
    manifest = tmp_path / "approved.tsv"
    manifest.write_text(
        "project_accession\trun_accession\tsample_id\n"
        "PRJ1\tSRR1\tSRR1\n"
        "PRJ2\tSRR9\tSRR9\n",
        encoding="utf-8",
    )
    approved = load_approved_runs(manifest, "PRJ1")
    specs = [
        replace(make_spec(b"one"), run_accession="SRR1"),
        replace(make_spec(b"two"), run_accession="SRR2"),
    ]

    selected, unavailable = select_approved_manifest(specs, [], approved)

    assert [spec.run_accession for spec in selected] == ["SRR1"]
    assert unavailable == []


def test_exact_run_manifest_fails_when_approved_run_is_missing(
    tmp_path: Path,
) -> None:
    manifest = tmp_path / "approved.tsv"
    manifest.write_text(
        "project_accession\trun_accession\nPRJ1\tSRR_MISSING\n",
        encoding="utf-8",
    )
    approved = load_approved_runs(manifest, "PRJ1")

    with pytest.raises(RuntimeError, match="SRR_MISSING"):
        select_approved_manifest([], [], approved)


def test_exact_run_download_quarantines_preexisting_unapproved_fastq(
    tmp_path: Path,
) -> None:
    spec = make_spec(b"approved")
    extra = tmp_path / "SRR_EXTRA_1.fastq.gz"
    extra.write_bytes(b"not approved")

    moved = quarantine_unapproved_fastqs(tmp_path, [spec])

    assert moved == [tmp_path / "unapproved_fastqs" / extra.name]
    assert moved[0].read_bytes() == b"not approved"
    assert not extra.exists()


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


def test_select_forward_download_accepts_only_unpaired_file_from_paired_layout() -> None:
    record = {
        "run_accession": "ERR500",
        "library_layout": "PAIRED",
        "library_strategy": "AMPLICON",
        "fastq_ftp": "ftp.sra.ebi.ac.uk/ERR500.fastq.gz",
        "fastq_md5": "abc",
        "fastq_bytes": "100",
    }

    spec = select_forward_download("ERP1", record)

    assert spec is not None
    assert spec.url == "https://ftp.sra.ebi.ac.uk/ERR500.fastq.gz"
    assert spec.output_filename == "ERR500_1.fastq.gz"


def test_select_forward_download_prefers_forward_mate_over_unpaired_file() -> None:
    record = {
        "run_accession": "ERR600",
        "library_layout": "PAIRED",
        "library_strategy": "AMPLICON",
        "fastq_ftp": (
            "ftp.sra.ebi.ac.uk/ERR600.fastq.gz;"
            "ftp.sra.ebi.ac.uk/ERR600_1.fastq.gz;"
            "ftp.sra.ebi.ac.uk/ERR600_2.fastq.gz"
        ),
        "fastq_md5": "unpaired;forward;reverse",
        "fastq_bytes": "50;100;100",
    }

    spec = select_forward_download("ERP1", record)

    assert spec is not None
    assert spec.md5 == "forward"
    assert spec.output_filename == "ERR600_1.fastq.gz"


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


def test_transfer_falls_back_to_verified_curl_after_python_https_returns_html(
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
    curl_requests: list[str] = []

    def fake_run(command: list[str], **_: object) -> object:
        curl_requests.append(command[-1])
        output_path = Path(command[command.index("--output") + 1])
        output_path.write_bytes(content)
        return type("Result", (), {"returncode": 0, "stderr": ""})()

    monkeypatch.setattr(download_ena_amplicon.shutil, "which", lambda _: "/usr/bin/curl")
    monkeypatch.setattr(download_ena_amplicon.subprocess, "run", fake_run)

    downloaded = transfer_once(spec, destination, timeout=60)

    assert downloaded == len(content)
    assert destination.read_bytes() == content
    assert not destination.with_name(destination.name + ".part").exists()
    assert curl_requests == ["https://ftp.sra.ebi.ac.uk/SRR1_1.fastq.gz"]


def test_transfer_falls_back_to_curl_after_requests_timeout(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    content = b"complete-fastq"
    spec = make_spec(content)
    destination = tmp_path / spec.output_filename
    monkeypatch.setattr(
        download_ena_amplicon.requests,
        "get",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            download_ena_amplicon.requests.Timeout("timed out")
        ),
    )

    def fake_curl(
        received_spec: DownloadSpec, received_destination: Path, timeout: float
    ) -> int:
        assert received_spec == spec
        assert timeout == 60
        received_destination.write_bytes(content)
        return len(content)

    monkeypatch.setattr(download_ena_amplicon, "transfer_via_curl", fake_curl)

    assert transfer_once(spec, destination, timeout=60) == len(content)
    assert destination.read_bytes() == content


def test_curl_receives_integer_timeout_and_resumes_partial(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    content = b"complete-fastq"
    spec = make_spec(content)
    destination = tmp_path / spec.output_filename
    part_path = destination.with_name(destination.name + ".part")
    part_path.write_bytes(content[:4])
    commands: list[list[str]] = []

    def fake_run(command: list[str], **_: object) -> object:
        commands.append(command)
        part_path.write_bytes(content)
        return type("Result", (), {"returncode": 0, "stderr": ""})()

    monkeypatch.setattr(download_ena_amplicon.shutil, "which", lambda _: "/usr/bin/curl")
    monkeypatch.setattr(download_ena_amplicon.subprocess, "run", fake_run)

    downloaded = download_ena_amplicon.transfer_via_curl(
        spec, destination, timeout=60.0
    )

    assert downloaded == len(content) - 4
    assert "--continue-at" in commands[0]
    assert commands[0][commands[0].index("--speed-time") + 1] == "60"
    assert destination.read_bytes() == content


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


def test_complete_partial_is_verified_and_promoted_without_network(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    content = b"complete-fastq"
    spec = make_spec(content)
    destination = tmp_path / spec.output_filename
    part_path = destination.with_name(destination.name + ".part")
    part_path.write_bytes(content)
    monkeypatch.setattr(
        download_ena_amplicon.requests,
        "get",
        lambda *_args, **_kwargs: pytest.fail("download should not be attempted"),
    )

    result = download_one(spec, tmp_path, retries=3, timeout=60)

    assert result.status == "downloaded"
    assert "promoted" in result.message
    assert destination.read_bytes() == content
    assert not part_path.exists()


def test_validation_checks_all_files_and_cleans_stale_partials(tmp_path: Path) -> None:
    content = b"complete-fastq"
    spec = make_spec(content)
    destination = tmp_path / spec.output_filename
    destination.write_bytes(content)
    destination.with_name(destination.name + ".part").write_bytes(b"stale")
    (tmp_path / "obsolete.fastq.gz.part").write_bytes(b"stale")

    report = validate_downloads([spec], tmp_path)

    assert report.complete
    assert report.valid_runs == frozenset({spec.run_accession})
    assert not list(tmp_path.glob("*.part"))
    assert (tmp_path / ".ena_orphaned_partials" / "obsolete.fastq.gz.part").is_file()


def make_run_args(output_dir: Path, **overrides: object) -> SimpleNamespace:
    values: dict[str, object] = {
        "accession": "PRJ1",
        "output_dir": output_dir,
        "workers": 1,
        "retries": 1,
        "timeout": 1.0,
        "dry_run": False,
        "recovery_rounds": 2,
        "retry_until_complete": False,
        "max_runtime": 60.0,
        "allow_missing_forward_fastq": False,
    }
    values.update(overrides)
    return SimpleNamespace(**values)


def test_run_exits_success_only_after_final_validation_and_writes_marker(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    content = b"complete-fastq"
    spec = make_spec(content)
    monkeypatch.setattr(
        download_ena_amplicon, "fetch_manifest", lambda *_: ([spec], [])
    )

    def fake_download(
        received_spec: DownloadSpec,
        output_dir: Path,
        retries: int,
        timeout: float,
    ) -> DownloadResult:
        del retries, timeout
        (output_dir / received_spec.output_filename).write_bytes(content)
        return DownloadResult(received_spec, "downloaded", 0.1, len(content), "ok")

    monkeypatch.setattr(download_ena_amplicon, "download_one", fake_download)

    assert download_ena_amplicon.run(make_run_args(tmp_path)) == 0
    marker_path = tmp_path / ".ena_download_complete.json"
    marker = json.loads(marker_path.read_text())
    assert marker["project_accession"] == "PRJ1"
    assert marker["selected_fastqs"] == 1
    assert marker["validation"].startswith("all files passed")
    assert "\tcompleted\n" in (
        tmp_path / "ena_amplicon_forward_summary.tsv"
    ).read_text()


def test_run_retries_revalidates_and_exits_incomplete(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    content = b"complete-fastq"
    spec = make_spec(content)
    attempts = 0
    monkeypatch.setattr(
        download_ena_amplicon, "fetch_manifest", lambda *_: ([spec], [])
    )
    monkeypatch.setattr(download_ena_amplicon.time, "sleep", lambda _: None)

    def fake_failure(*_: object, **__: object) -> DownloadResult:
        nonlocal attempts
        attempts += 1
        return DownloadResult(spec, "failed", 0.1, 0, "ENA unavailable")

    monkeypatch.setattr(download_ena_amplicon, "download_one", fake_failure)

    assert download_ena_amplicon.run(make_run_args(tmp_path)) == 1
    assert attempts == 2
    assert not (tmp_path / ".ena_download_complete.json").exists()
    status = (tmp_path / "ena_amplicon_forward_status.tsv").read_text()
    assert "failed" in status
    assert "ENA unavailable" in status
    assert "\tincomplete\n" in (
        tmp_path / "ena_amplicon_forward_summary.tsv"
    ).read_text()


def test_runtime_limit_before_first_recovery_writes_current_status(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    content = b"complete-fastq"
    spec = make_spec(content)
    monkeypatch.setattr(
        download_ena_amplicon, "fetch_manifest", lambda *_: ([spec], [])
    )
    monkeypatch.setattr(
        download_ena_amplicon,
        "download_one",
        lambda *_args, **_kwargs: pytest.fail("download should not start"),
    )

    args = make_run_args(tmp_path, max_runtime=1e-12)
    assert download_ena_amplicon.run(args) == 1

    status_path = tmp_path / "ena_amplicon_forward_status.tsv"
    with status_path.open(newline="", encoding="utf-8") as handle:
        statuses = list(csv.DictReader(handle, delimiter="\t"))
    assert len(statuses) == 1
    assert statuses[0]["run_accession"] == spec.run_accession
    assert statuses[0]["status"] == "failed"
    assert "missing" in statuses[0]["message"]


def test_run_recovers_on_a_later_full_pass(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    content = b"complete-fastq"
    spec = make_spec(content)
    attempts = 0
    monkeypatch.setattr(
        download_ena_amplicon, "fetch_manifest", lambda *_: ([spec], [])
    )
    monkeypatch.setattr(download_ena_amplicon.time, "sleep", lambda _: None)

    def fail_then_succeed(
        received_spec: DownloadSpec,
        output_dir: Path,
        *_: object,
    ) -> DownloadResult:
        nonlocal attempts
        attempts += 1
        if attempts == 1:
            return DownloadResult(
                received_spec, "failed", 0.1, 0, "temporary ENA failure"
            )
        (output_dir / received_spec.output_filename).write_bytes(content)
        return DownloadResult(
            received_spec, "downloaded", 0.1, len(content), "recovered"
        )

    monkeypatch.setattr(download_ena_amplicon, "download_one", fail_then_succeed)

    assert download_ena_amplicon.run(make_run_args(tmp_path)) == 0
    assert attempts == 2
    assert (tmp_path / ".ena_download_complete.json").exists()
    summary = (tmp_path / "ena_amplicon_forward_summary.tsv").read_text()
    assert "\t2\t0\tcompleted\n" in summary


def test_run_refuses_completion_when_amplicon_run_has_no_forward_fastq(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    content = b"complete-fastq"
    spec = make_spec(content)
    (tmp_path / spec.output_filename).write_bytes(content)
    monkeypatch.setattr(
        download_ena_amplicon,
        "fetch_manifest",
        lambda *_: ([spec], [("SRR2", "no forward archive FASTQ")]),
    )

    assert download_ena_amplicon.run(make_run_args(tmp_path)) == 1
    assert not (tmp_path / ".ena_download_complete.json").exists()
    unavailable = (tmp_path / "ena_amplicon_forward_unavailable.tsv").read_text()
    assert "SRR2" in unavailable
    assert "no forward archive FASTQ" in unavailable


def test_run_can_explicitly_allow_amplicon_run_without_forward_fastq(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    content = b"complete-fastq"
    spec = make_spec(content)
    (tmp_path / spec.output_filename).write_bytes(content)
    monkeypatch.setattr(
        download_ena_amplicon,
        "fetch_manifest",
        lambda *_: ([spec], [("SRR2", "no forward archive FASTQ")]),
    )

    args = make_run_args(tmp_path, allow_missing_forward_fastq=True)
    assert download_ena_amplicon.run(args) == 0
    marker = json.loads((tmp_path / ".ena_download_complete.json").read_text())
    assert marker["unavailable_amplicon_runs"] == 1


def test_summary_history_preserves_multiple_invocations(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    content = b"complete-fastq"
    spec = make_spec(content)
    (tmp_path / spec.output_filename).write_bytes(content)
    monkeypatch.setattr(
        download_ena_amplicon, "fetch_manifest", lambda *_: ([spec], [])
    )

    assert download_ena_amplicon.run(make_run_args(tmp_path)) == 0
    assert download_ena_amplicon.run(make_run_args(tmp_path)) == 0

    history_path = tmp_path / "ena_amplicon_forward_history.tsv"
    with history_path.open(newline="", encoding="utf-8") as handle:
        history = list(csv.DictReader(handle, delimiter="\t"))
    assert len(history) == 2
    assert all(row["status"] == "completed" for row in history)


def test_non_dry_run_with_no_downloadable_fastqs_invalidates_old_marker(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    marker = tmp_path / ".ena_download_complete.json"
    marker.write_text("stale")
    monkeypatch.setattr(
        download_ena_amplicon, "fetch_manifest", lambda *_: ([], [])
    )

    with pytest.raises(RuntimeError, match="No downloadable"):
        download_ena_amplicon.run(make_run_args(tmp_path))

    assert not marker.exists()
