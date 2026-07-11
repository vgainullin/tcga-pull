"""CLI-boundary tests. These exercise the Typer layer, where flag defaults
turn into CohortBuildOptions — the service layer is tested in test_services.py.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from typer.testing import CliRunner

from tcga_pull import cli, services
from tcga_pull.config import CohortSpec
from tcga_pull.coverage import CoverageMatrix, CoverageOutputs
from tcga_pull.overlap import (
    IntersectionSummary,
    OverlapOutputs,
    OverlapReport,
    SelectionSummary,
)

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


def test_dataset_command_passes_options_to_service(tmp_path: Path, monkeypatch):
    cohort = tmp_path / "cohort"
    cohort.mkdir()
    yaml_path = tmp_path / "cohort.yaml"
    yaml_path.write_text(
        "name: c\n"
        "filters: {project: TCGA-BRCA}\n"
        "recipe_options:\n"
        "  model_dataset:\n"
        "    label_column: lineage\n"
    )
    captured: dict[str, Any] = {}

    def fake_write_model_dataset_recipe(*args: Any, **kwargs: Any):
        captured["args"] = args
        captured["kwargs"] = kwargs
        out = tmp_path / "model_dataset"
        return services.ModelDatasetOutputs(
            path=out,
            manifest=out / "manifest.json",
            samples=out / "samples.parquet",
            feature_index=out / "feature_index.parquet",
            matrices={"snv": out / "snv.parquet"},
        )

    monkeypatch.setattr(services, "write_model_dataset_recipe", fake_write_model_dataset_recipe)

    result = runner.invoke(
        cli.app,
        [
            "dataset",
            str(cohort),
            "--config",
            str(yaml_path),
            "--modality",
            "snv",
            "--min-class-count",
            "1",
            "--seed",
            "42",
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured["args"] == (cohort,)
    assert captured["kwargs"]["recipe_options"]["model_dataset"]["label_column"] == "lineage"
    assert captured["kwargs"]["modalities"] == ["snv"]
    assert captured["kwargs"]["min_class_count"] == 1
    assert captured["kwargs"]["seed"] == 42


def test_coverage_command_writes_matrix_outputs(tmp_path: Path, monkeypatch):
    captured: dict[str, object] = {}
    matrix = CoverageMatrix(
        program="TCGA",
        project_ids=("TCGA-BRCA",),
        rows=(
            {"support_status": "supported"},
            {"support_status": "raw_only"},
        ),
    )

    def fake_build_coverage_matrix(*, program: str, projects: list[str] | None = None):
        captured["program"] = program
        captured["projects"] = projects
        return matrix

    def fake_write_coverage_outputs(matrix_arg: CoverageMatrix, out_dir: Path):
        captured["matrix"] = matrix_arg
        captured["out_dir"] = out_dir
        return CoverageOutputs(
            parquet_path=tmp_path / "tcga_open_access_coverage_matrix.parquet",
            markdown_path=tmp_path / "tcga_open_access_coverage_matrix.md",
            n_rows=2,
            n_projects=1,
        )

    monkeypatch.setattr(services, "build_coverage_matrix", fake_build_coverage_matrix)
    monkeypatch.setattr(services, "write_coverage_outputs", fake_write_coverage_outputs)

    result = runner.invoke(
        cli.app,
        [
            "coverage",
            "--program",
            "TCGA",
            "--project",
            "TCGA-BRCA",
            "--out-dir",
            str(tmp_path),
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured["program"] == "TCGA"
    assert captured["projects"] == ["TCGA-BRCA"]
    assert captured["matrix"] == matrix
    assert captured["out_dir"] == tmp_path
    assert "supported=1" in result.output
    assert "raw_only=1" in result.output


def test_overlap_command_passes_omics_and_writes_outputs(tmp_path: Path, monkeypatch):
    config = tmp_path / "cohort.yaml"
    config.write_text(
        "name: c\nfilters: {project: TCGA-CHOL, data_format: MAF}\n"
        "optional_omics:\n  - name: rna\n    filters: {data_type: RNA}\n"
    )
    selection = SelectionSummary("primary", {}, 2, 3, 100, {}, {})
    intersection = IntersectionSummary(("primary", "rna"), 2, {})
    report = OverlapReport(
        "c", "2026-07-11T00:00:00+00:00", (selection,), (intersection,), intersection
    )
    captured: dict[str, object] = {}

    def fake_build(spec, *, omics):
        captured["omics"] = omics
        return report

    def fake_write(report_arg, *, json_path=None, parquet_path=None):
        captured["report"] = report_arg
        captured["json_path"] = json_path
        return OverlapOutputs(json_path=json_path)

    monkeypatch.setattr(services, "build_overlap_report", fake_build)
    monkeypatch.setattr(services, "write_overlap_outputs", fake_write)

    result = runner.invoke(
        cli.app, ["overlap", str(config), "--omics", "rna", "--json", str(tmp_path / "o.json")]
    )

    assert result.exit_code == 0, result.output
    assert captured["omics"] == ["rna"]
    assert captured["report"] == report
    assert "all selected" in result.output
    assert "queried_at=2026-07-11" in result.output
