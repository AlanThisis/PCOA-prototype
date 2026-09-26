from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

import faith_pd


def write_metadata(path: Path, groups: dict[str, str]) -> Path:
    rows = ["sample-id\tbody_site"] + [f"{sample}\t{group}" for sample, group in groups.items()]
    path.write_text("\n".join(rows) + "\n")
    return path


def test_run_exports_values_and_plots_each_grouping(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    values = {"A_1": 10.0, "B_1": 12.0, "C_1": 30.0, "D_1": 33.0, "E_1": 5.0}
    commands: list[list[str]] = []

    def fake_run_command(command: list[str], **_: object) -> None:
        commands.append(command)
        if command[1:3] == ["tools", "export"]:
            out = Path(command[command.index("--output-path") + 1])
            out.mkdir(parents=True)
            lines = ["\tfaith_pd"] + [f"{s}\t{v}" for s, v in values.items()]
            (out / "alpha-diversity.tsv").write_text("\n".join(lines) + "\n")

    monkeypatch.setattr(faith_pd, "run_command", fake_run_command)
    monkeypatch.setattr(faith_pd, "resolve_executable", lambda name: f"/env/bin/{name}")
    table = tmp_path / "rarefied.qza"
    tree = tmp_path / "tree.qza"
    table.write_bytes(b"qza")
    tree.write_bytes(b"qza")
    # E has no label; read suffixes are stripped to match metadata IDs.
    metadata = write_metadata(tmp_path / "meta.tsv", {"A": "gut", "B": "gut", "C": "oral", "D": "oral"})
    out = tmp_path / "results"

    assert faith_pd.main([
        "run", "--rarefied-table", str(table), "--phylogeny", str(tree),
        "--output-dir", str(out), "--metadata", str(metadata), "--group-by", "body_site",
    ]) == 0

    beta = commands[0]
    assert beta[1:3] == ["diversity", "alpha-phylogenetic"]
    assert beta[beta.index("--p-metric") + 1] == "faith_pd"
    assert beta[beta.index("--i-table") + 1] == str(table.resolve())
    written = pd.read_csv(out / "faith_pd.tsv", sep="\t")
    assert list(written.columns) == ["sample-id", "faith_pd"]
    assert len(written) == 5
    assert (out / "faith_pd_body_site.png").stat().st_size > 0
    summary = json.loads((out / "faith_pd_summary.json").read_text())
    grouping = summary["groupings"]["body_site"]
    assert grouping["n"] == 4
    assert grouping["unlabelled_samples_excluded"] == 1
    assert grouping["group_medians"] == {"oral": 31.5, "gut": 11.0}
    assert grouping["epsilon2"] == pytest.approx(grouping["kruskal_h"] / 3)


def test_run_requires_metadata_for_grouping(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="--group-by requires --metadata"):
        faith_pd.run(
            faith_pd.parse_args([
                "run", "--rarefied-table", "t.qza", "--phylogeny", "p.qza",
                "--output-dir", str(tmp_path), "--group-by", "body_site",
            ]),
            timing=None,
        )


def write_values(path: Path, values: dict[str, float]) -> Path:
    pd.Series(values, name="faith_pd").rename_axis("sample-id").to_frame().to_csv(path, sep="\t")
    return path


def test_compare_uses_shared_samples_and_reports_drift(tmp_path: Path) -> None:
    full = write_values(tmp_path / "full.tsv", {"A": 10, "B": 20, "C": 30, "D": 40, "X": 50})
    sub = write_values(tmp_path / "sub.tsv", {"A": 9, "B": 18, "C": 27, "D": 36})
    metadata = write_metadata(tmp_path / "meta.tsv", {"A": "gut", "B": "gut", "C": "oral", "D": "oral", "X": "gut"})
    out = tmp_path / "compare"

    assert faith_pd.main([
        "compare", "--run", f"full={full}", "--run", f"sub={sub}",
        "--metadata", str(metadata), "--group-by", "body_site", "--output-dir", str(out),
    ]) == 0

    summary = pd.read_csv(out / "faith_pd_levels_summary.tsv", sep="\t").set_index("run")
    assert summary.loc["full", "n"] == 4  # X is missing from sub, so it is excluded
    assert summary.loc["sub", "spearman_vs_reference"] == pytest.approx(1.0)
    assert summary.loc["sub", "median_pct_change_vs_reference"] == pytest.approx(-10.0)
    assert (out / "faith_pd_levels_body_site.png").stat().st_size > 0
    assert (out / "faith_pd_levels_agreement.png").stat().st_size > 0


def test_compare_rejects_cohort_ids_missing_from_a_run(tmp_path: Path) -> None:
    full = write_values(tmp_path / "full.tsv", {"A": 10, "B": 20})
    sub = write_values(tmp_path / "sub.tsv", {"A": 9})
    cohort = tmp_path / "cohort.txt"
    cohort.write_text("A\nB\n")
    with pytest.raises(ValueError, match="1 cohort IDs lack"):
        faith_pd.compare_frames({"full": full, "sub": sub}, cohort)
