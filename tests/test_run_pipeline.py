from __future__ import annotations

import json
import threading
import time
from pathlib import Path

import pytest

import run_pipeline


def write_fastq(directory: Path, sample_id: str) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{sample_id}_1.fastq.gz"
    path.write_bytes(b"prepared-fastq")
    return path


def write_gg2(directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for filename in (
        run_pipeline.GG2_BACKBONE_FILENAME,
        run_pipeline.GG2_ID_TREE_FILENAME,
    ):
        artifact = directory / filename
        if not artifact.exists():
            artifact.write_bytes(b"artifact")


def command_value(command: list[str], option: str) -> Path:
    return Path(command[command.index(option) + 1])


def install_fake_environment(
    monkeypatch: pytest.MonkeyPatch,
    commands: list[list[str]],
    *,
    fail_script: str | None = None,
    empty_study: str | None = None,
) -> None:
    monkeypatch.setattr(
        run_pipeline,
        "resolve_executable",
        lambda name: f"/environment/bin/{name}",
    )
    monkeypatch.setattr(run_pipeline, "git_info", lambda _: ("abc123", False))
    monkeypatch.setattr(run_pipeline, "validate_qiime_gg2", lambda _: None)

    def fake_run_command(command: list[str], **_: object) -> None:
        commands.append(command)
        script = Path(command[1]).name
        if script == fail_script:
            raise RuntimeError(f"injected failure in {script}")
        if script == "run_deblur.py":
            output_dir = command_value(command, "--work-dir") / "workflow"
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / "all.biom").write_bytes(b"biom")
            is_empty = command_value(command, "--data-dir").name == empty_study
            (output_dir / "all.seqs.fa").write_bytes(b"" if is_empty else b"fasta")
        elif script == "merge_biom.py":
            output_dir = command_value(command, "--out-dir")
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / "all.biom").write_bytes(b"biom")
            (output_dir / "all.seqs.fa").write_bytes(b"fasta")
        elif script == "unifrac.py":
            output_dir = command_value(command, "--results-dir")
            output_dir.mkdir(parents=True, exist_ok=True)
            for filename in run_pipeline.UNIFRAC_OUTPUT_NAMES:
                (output_dir / filename).write_bytes(b"result")
        elif script == "plot_pcoa.py":
            output_path = command_value(command, "--out")
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_bytes(b"plot")
        else:
            raise AssertionError(f"Unexpected component command: {command}")

    monkeypatch.setattr(run_pipeline, "run_command", fake_run_command)


def make_args(
    tmp_path: Path,
    studies: list[tuple[str, Path]],
    *,
    metadata: Path | None = None,
    color_by: list[str] | None = None,
    run_dir: Path | None = None,
    extra: list[str] | None = None,
) -> object:
    gg2_dir = tmp_path / "gg2"
    write_gg2(gg2_dir)
    argv: list[str] = []
    for name, fastq_dir in studies:
        argv.extend(("--study", f"{name}={fastq_dir}"))
    argv.extend(
        (
            "--run-dir",
            str(run_dir or tmp_path / "run"),
            "--gg2-dir",
            str(gg2_dir),
            "--trim-length",
            "120",
            "--min-reads",
            "0",
            "--threads",
            "16",
        )
    )
    if metadata:
        argv.extend(("--metadata", str(metadata)))
    for column in color_by or []:
        argv.extend(("--color-by", column))
    argv.extend(extra or [])
    return run_pipeline.parse_args(argv)


def test_one_study_run_skips_merge_and_writes_manifest_and_state(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    study_dir = tmp_path / "study"
    write_fastq(study_dir, "ERR1")
    commands: list[list[str]] = []
    install_fake_environment(monkeypatch, commands)
    args = make_args(tmp_path, [("ERP", study_dir)])

    run_dir = run_pipeline.execute_pipeline(args)

    assert [Path(command[1]).name for command in commands] == [
        "run_deblur.py",
        "unifrac.py",
    ]
    unifrac_command = commands[1]
    assert command_value(unifrac_command, "--deblur-dir") == (
        run_dir / "work" / "deblur" / "ERP" / "workflow"
    )
    assert unifrac_command[unifrac_command.index("--sampling-depth") + 1] == "1000"
    manifest = json.loads((run_dir / "run_manifest.json").read_text())
    assert manifest["git_commit"] == "abc123"
    assert manifest["scientific_parameters"]["trim_length"] == 120
    assert manifest["scientific_parameters"]["sampling_depth"] == 1000
    assert manifest["studies"][0]["fastqs"][0]["path"] == str(
        (study_dir / "ERR1_1.fastq.gz").resolve()
    )
    assert manifest["studies"][0]["fastqs"][0]["size"] > 0
    state = json.loads((run_dir / "run_state.json").read_text())
    assert set(state["stages"]) == {"deblur:ERP", "unifrac"}
    assert all(stage["status"] == "completed" for stage in state["stages"].values())
    assert state["attempts"][0]["threads"] == 16
    assert (run_dir / "timings" / "attempt-001" / "pipeline.tsv").is_file()


def test_cross_study_run_merges_and_generates_metadata_plots(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    write_fastq(first, "ERR1")
    write_fastq(second, "ERR2")
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "sample-id\tdisease status\tstudy\nERR1\tcancer\tERP\nERR2\thealthy\tPRJ\n"
    )
    commands: list[list[str]] = []
    install_fake_environment(monkeypatch, commands)
    args = make_args(
        tmp_path,
        [("ERP", first), ("PRJ", second)],
        metadata=metadata,
        color_by=["disease status", "study"],
        extra=["--keep-deblur-tmp-files", "--sampling-depth", "100"],
    )

    run_dir = run_pipeline.execute_pipeline(args)

    assert [Path(command[1]).name for command in commands] == [
        "run_deblur.py",
        "run_deblur.py",
        "merge_biom.py",
        "unifrac.py",
        "plot_pcoa.py",
        "plot_pcoa.py",
    ]
    assert all("--keep-tmp-files" in command for command in commands[:2])
    merge_command = commands[2]
    merge_inputs = merge_command[
        merge_command.index("--deblur-dirs") + 1 : merge_command.index("--out-dir")
    ]
    assert merge_inputs == [
        str(run_dir / "work" / "deblur" / "ERP" / "workflow"),
        str(run_dir / "work" / "deblur" / "PRJ" / "workflow"),
    ]
    assert command_value(commands[3], "--deblur-dir") == run_dir / "work" / "merged"
    assert commands[3][commands[3].index("--sampling-depth") + 1] == "100"
    assert "--refresh-input-artifacts" in commands[3]
    assert (run_dir / "results" / "pcoa_disease_status.png").is_file()
    assert (run_dir / "results" / "pcoa_study.png").is_file()


def test_cross_study_run_excludes_zero_feature_deblur_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    empty = tmp_path / "empty-study"
    retained = tmp_path / "retained-study"
    write_fastq(empty, "ERR0")
    write_fastq(retained, "ERR1")
    commands: list[list[str]] = []
    install_fake_environment(monkeypatch, commands, empty_study=empty.name)

    run_dir = run_pipeline.execute_pipeline(
        make_args(tmp_path, [("empty", empty), ("retained", retained)])
    )

    state = json.loads((run_dir / "run_state.json").read_text())
    assert state["stages"]["deblur:empty"]["status"] == "excluded_empty"
    assert state["stages"]["deblur:retained"]["status"] == "completed"
    merge_command = next(
        command for command in commands if Path(command[1]).name == "merge_biom.py"
    )
    assert "--skip-empty" in merge_command
    assert state["stages"]["merge"]["status"] == "completed"
    assert state["stages"]["unifrac"]["status"] == "completed"


def test_cross_study_deblur_runs_concurrently_with_bounded_inner_jobs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    studies = []
    for index in range(4):
        study_dir = tmp_path / f"study-{index}"
        write_fastq(study_dir, f"ERR{index}")
        studies.append((f"study-{index}", study_dir))

    commands: list[list[str]] = []
    install_fake_environment(monkeypatch, commands)
    active = 0
    maximum_active = 0
    lock = threading.Lock()
    original_run_command = run_pipeline.run_command

    def delayed_run_command(command: list[str], **kwargs: object) -> None:
        nonlocal active, maximum_active
        if Path(command[1]).name == "run_deblur.py":
            with lock:
                active += 1
                maximum_active = max(maximum_active, active)
            time.sleep(0.03)
            try:
                original_run_command(command, **kwargs)
            finally:
                with lock:
                    active -= 1
        else:
            original_run_command(command, **kwargs)

    monkeypatch.setattr(run_pipeline, "run_command", delayed_run_command)
    args = make_args(
        tmp_path,
        studies,
        extra=[
            "--deblur-study-workers",
            "2",
            "--deblur-jobs-per-study",
            "3",
        ],
    )

    run_dir = run_pipeline.execute_pipeline(args)

    assert maximum_active == 2
    deblur_commands = [
        command for command in commands if Path(command[1]).name == "run_deblur.py"
    ]
    assert len(deblur_commands) == 4
    assert all(
        command[command.index("--jobs-to-start") + 1] == "3"
        for command in deblur_commands
    )
    state = json.loads((run_dir / "run_state.json").read_text())
    assert all(stage["status"] == "completed" for stage in state["stages"].values())
    assert state["attempts"][0]["deblur_study_workers"] == 2
    assert state["attempts"][0]["deblur_jobs_per_study"] == 3


def test_deblur_parallelism_auto_divides_threads_across_effective_workers() -> None:
    class Args:
        threads = 16
        deblur_study_workers = 8
        deblur_jobs_per_study = None

    assert run_pipeline.resolve_deblur_parallelism(Args(), 2) == (2, 8)


def test_deblur_parallelism_rejects_oversubscription() -> None:
    class Args:
        threads = 16
        deblur_study_workers = 8
        deblur_jobs_per_study = 4

    with pytest.raises(ValueError, match="oversubscribes"):
        run_pipeline.resolve_deblur_parallelism(Args(), 8)


def test_concurrent_deblur_failure_leaves_resumable_state(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    studies = []
    for index in range(3):
        study_dir = tmp_path / f"study-{index}"
        write_fastq(study_dir, f"ERR{index}")
        studies.append((f"study-{index}", study_dir))
    run_dir = tmp_path / "run"
    commands: list[list[str]] = []
    install_fake_environment(monkeypatch, commands)
    original_run_command = run_pipeline.run_command

    def fail_one_study(command: list[str], **kwargs: object) -> None:
        if Path(command[1]).name == "run_deblur.py":
            data_dir = command_value(command, "--data-dir")
            if data_dir.name == "study-0":
                time.sleep(0.01)
                raise RuntimeError("injected concurrent Deblur failure")
            time.sleep(0.05)
        original_run_command(command, **kwargs)

    monkeypatch.setattr(run_pipeline, "run_command", fail_one_study)
    args = make_args(
        tmp_path,
        studies,
        run_dir=run_dir,
        extra=["--deblur-study-workers", "2"],
    )

    with pytest.raises(RuntimeError, match="concurrent Deblur failure"):
        run_pipeline.execute_pipeline(args)

    state = json.loads((run_dir / "run_state.json").read_text())
    deblur_statuses = {
        name: stage["status"]
        for name, stage in state["stages"].items()
        if name.startswith("deblur:")
    }
    assert deblur_statuses["deblur:study-0"] == "failed"
    assert "running" not in deblur_statuses.values()
    assert state["attempts"][0]["status"] == "failed"

    resumed_commands: list[list[str]] = []
    install_fake_environment(monkeypatch, resumed_commands)
    resume_args = make_args(
        tmp_path,
        studies,
        run_dir=run_dir,
        extra=["--resume", "--deblur-study-workers", "2"],
    )
    run_pipeline.execute_pipeline(resume_args)

    resumed_state = json.loads((run_dir / "run_state.json").read_text())
    assert all(
        stage["status"] == "completed"
        for stage in resumed_state["stages"].values()
    )
    assert resumed_state["attempts"][1]["status"] == "completed"


def test_failed_stage_is_persisted_and_can_be_resumed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    study_dir = tmp_path / "study"
    write_fastq(study_dir, "ERR1")
    run_dir = tmp_path / "run"
    failed_commands: list[list[str]] = []
    install_fake_environment(monkeypatch, failed_commands, fail_script="unifrac.py")
    args = make_args(tmp_path, [("ERP", study_dir)], run_dir=run_dir)

    with pytest.raises(RuntimeError, match="injected failure"):
        run_pipeline.execute_pipeline(args)

    state = json.loads((run_dir / "run_state.json").read_text())
    assert state["stages"]["deblur:ERP"]["status"] == "completed"
    assert state["stages"]["unifrac"]["status"] == "failed"
    assert state["attempts"][0]["status"] == "failed"

    resumed_commands: list[list[str]] = []
    install_fake_environment(monkeypatch, resumed_commands)
    resume_args = make_args(
        tmp_path,
        [("ERP", study_dir)],
        run_dir=run_dir,
        extra=["--resume", "--threads", "2"],
    )
    run_pipeline.execute_pipeline(resume_args)

    assert [Path(command[1]).name for command in resumed_commands] == ["unifrac.py"]
    assert "--refresh-input-artifacts" not in resumed_commands[0]
    state = json.loads((run_dir / "run_state.json").read_text())
    assert state["attempts"][1]["status"] == "completed"
    assert state["attempts"][1]["threads"] == 2
    assert (run_dir / "timings" / "attempt-002" / "pipeline.tsv").is_file()


def test_missing_completed_output_reruns_stage_and_downstream(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    write_fastq(first, "ERR1")
    write_fastq(second, "ERR2")
    run_dir = tmp_path / "run"
    commands: list[list[str]] = []
    install_fake_environment(monkeypatch, commands)
    args = make_args(tmp_path, [("ERP", first), ("PRJ", second)], run_dir=run_dir)
    run_pipeline.execute_pipeline(args)
    (run_dir / "work" / "merged" / "all.biom").unlink()

    resumed_commands: list[list[str]] = []
    install_fake_environment(monkeypatch, resumed_commands)
    resume_args = make_args(
        tmp_path,
        [("ERP", first), ("PRJ", second)],
        run_dir=run_dir,
        extra=["--resume"],
    )
    run_pipeline.execute_pipeline(resume_args)

    assert [Path(command[1]).name for command in resumed_commands] == [
        "merge_biom.py",
        "unifrac.py",
    ]
    assert "--refresh-input-artifacts" in resumed_commands[1]


def test_missing_deblur_output_does_not_rerun_independent_study(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    write_fastq(first, "ERR1")
    write_fastq(second, "ERR2")
    run_dir = tmp_path / "run"
    commands: list[list[str]] = []
    install_fake_environment(monkeypatch, commands)
    run_pipeline.execute_pipeline(
        make_args(tmp_path, [("ERP", first), ("PRJ", second)], run_dir=run_dir)
    )
    (run_dir / "work" / "deblur" / "ERP" / "workflow" / "all.biom").unlink()

    resumed_commands: list[list[str]] = []
    install_fake_environment(monkeypatch, resumed_commands)
    run_pipeline.execute_pipeline(
        make_args(
            tmp_path,
            [("ERP", first), ("PRJ", second)],
            run_dir=run_dir,
            extra=["--resume"],
        )
    )

    assert [Path(command[1]).name for command in resumed_commands] == [
        "run_deblur.py",
        "merge_biom.py",
        "unifrac.py",
    ]
    assert command_value(resumed_commands[0], "--data-dir") == first.resolve()
    assert "--refresh-input-artifacts" in resumed_commands[2]


def test_resume_rejects_manifest_mismatch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    study_dir = tmp_path / "study"
    write_fastq(study_dir, "ERR1")
    run_dir = tmp_path / "run"
    commands: list[list[str]] = []
    install_fake_environment(monkeypatch, commands)
    run_pipeline.execute_pipeline(
        make_args(tmp_path, [("ERP", study_dir)], run_dir=run_dir)
    )

    mismatch = make_args(
        tmp_path,
        [("ERP", study_dir)],
        run_dir=run_dir,
        extra=["--resume", "--trim-length", "121"],
    )
    with pytest.raises(RuntimeError, match="does not match"):
        run_pipeline.execute_pipeline(mismatch)


def test_nonempty_run_requires_explicit_resume(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    study_dir = tmp_path / "study"
    write_fastq(study_dir, "ERR1")
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "unrelated.txt").write_text("occupied")
    commands: list[list[str]] = []
    install_fake_environment(monkeypatch, commands)

    with pytest.raises(RuntimeError, match="not empty"):
        run_pipeline.execute_pipeline(
            make_args(tmp_path, [("ERP", study_dir)], run_dir=run_dir)
        )


def test_validation_rejects_duplicate_sample_ids_across_studies(tmp_path: Path) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    write_fastq(first, "ERR1")
    write_fastq(second, "ERR1")

    with pytest.raises(ValueError, match="Duplicate sample identifier"):
        run_pipeline.parse_studies([f"ERP={first}", f"PRJ={second}"])


@pytest.mark.parametrize("name", [".", ".."])
def test_validation_rejects_dot_path_study_names(tmp_path: Path, name: str) -> None:
    study_dir = tmp_path / "study"
    write_fastq(study_dir, "ERR1")

    with pytest.raises(ValueError, match="Invalid study name"):
        run_pipeline.parse_studies([f"{name}={study_dir}"])


@pytest.mark.parametrize("changed_input", ["fastq", "metadata", "gg2"])
def test_resume_rejects_modified_input_files(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    changed_input: str,
) -> None:
    study_dir = tmp_path / "study"
    fastq = write_fastq(study_dir, "ERR1")
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("sample-id\tgroup\nERR1\tcase\n")
    run_dir = tmp_path / "run"
    commands: list[list[str]] = []
    install_fake_environment(monkeypatch, commands)
    run_pipeline.execute_pipeline(
        make_args(
            tmp_path,
            [("ERP", study_dir)],
            metadata=metadata,
            color_by=["group"],
            run_dir=run_dir,
        )
    )

    if changed_input == "fastq":
        fastq.write_bytes(b"modified-prepared-fastq")
    elif changed_input == "metadata":
        metadata.write_text("sample-id\tgroup\nERR1\tcontrol\n")
    else:
        (tmp_path / "gg2" / run_pipeline.GG2_BACKBONE_FILENAME).write_bytes(
            b"modified-artifact"
        )

    with pytest.raises(RuntimeError, match="does not match"):
        run_pipeline.execute_pipeline(
            make_args(
                tmp_path,
                [("ERP", study_dir)],
                metadata=metadata,
                color_by=["group"],
                run_dir=run_dir,
                extra=["--resume"],
            )
        )


def test_qiime_preflight_checks_greengenes2_plugin(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    commands: list[list[str]] = []

    def fake_run(command: list[str], **_: object) -> object:
        commands.append(command)
        return object()

    monkeypatch.setattr(run_pipeline.subprocess, "run", fake_run)

    run_pipeline.validate_qiime_gg2("/environment/bin/qiime")

    assert commands == [["/environment/bin/qiime", "greengenes2", "--help"]]


def test_qiime_preflight_reports_broken_greengenes2_plugin(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail_run(command: list[str], **_: object) -> object:
        raise run_pipeline.subprocess.CalledProcessError(
            1,
            command,
            stderr="No such command: greengenes2",
        )

    monkeypatch.setattr(run_pipeline.subprocess, "run", fail_run)

    with pytest.raises(RuntimeError, match="No such command: greengenes2"):
        run_pipeline.validate_qiime_gg2("/environment/bin/qiime")
