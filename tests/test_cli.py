"""CLI-boundary tests. These exercise the Typer layer, where flag defaults
turn into CohortBuildOptions — the service layer is tested in test_services.py.
"""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from tcga_pull import cli, services
from tcga_pull.config import CohortSpec

runner = CliRunner()


def _capture_spec(monkeypatch) -> list[CohortSpec]:
    """Patch services.run_cohort so `pull` builds a spec but does no I/O."""
    captured: list[CohortSpec] = []
    monkeypatch.setattr(services, "run_cohort", lambda spec, **_: captured.append(spec))
    return captured


def test_pull_without_flags_preserves_yaml_incremental(tmp_path: Path, monkeypatch):
    """Regression: omitting --incremental/--delete-raw must not force standard mode.

    The flags are tri-state (None = unset); a bare `tcga-pull pull <yaml>` has to
    leave the YAML's processing block untouched.
    """
    yaml_path = tmp_path / "cohort.yaml"
    yaml_path.write_text(
        "name: c\n"
        "filters: {project: TCGA-BRCA}\n"
        "processing:\n"
        "  mode: incremental\n"
        "  batch_size: 200\n"
        "  delete_raw_after_processing: true\n"
    )
    captured = _capture_spec(monkeypatch)

    result = runner.invoke(cli.app, ["pull", str(yaml_path)])

    assert result.exit_code == 0, result.output
    assert len(captured) == 1
    assert captured[0].processing.mode == "incremental"
    assert captured[0].processing.batch_size == 200
    assert captured[0].processing.delete_raw_after_processing is True


def test_pull_incremental_flag_overrides_standard_yaml(tmp_path: Path, monkeypatch):
    yaml_path = tmp_path / "cohort.yaml"
    yaml_path.write_text("name: c\nfilters: {project: TCGA-BRCA}\n")
    captured = _capture_spec(monkeypatch)

    result = runner.invoke(
        cli.app,
        ["pull", str(yaml_path), "--incremental", "--processing-batch-size", "50"],
    )

    assert result.exit_code == 0, result.output
    assert captured[0].processing.mode == "incremental"
    assert captured[0].processing.batch_size == 50


def test_pull_no_incremental_flag_forces_standard_over_yaml(tmp_path: Path, monkeypatch):
    yaml_path = tmp_path / "cohort.yaml"
    yaml_path.write_text(
        "name: c\nfilters: {project: TCGA-BRCA}\nprocessing: {mode: incremental}\n"
    )
    captured = _capture_spec(monkeypatch)

    result = runner.invoke(cli.app, ["pull", str(yaml_path), "--no-incremental"])

    assert result.exit_code == 0, result.output
    assert captured[0].processing.mode == "standard"
