"""Offline tests for metadata-only multi-modality overlap."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import polars as pl

from tcga_pull.config import CohortSpec, OptionalOmicsSpec
from tcga_pull.overlap import build_overlap_report, write_overlap_outputs


class FakeClient:
    def __init__(self, responses: list[list[dict]]) -> None:
        self.responses = iter(responses)
        self.filters: list[dict] = []

    def fetch_files(self, filters: dict, fields: list[str] | None = None) -> list[dict]:
        self.filters.append(filters)
        return next(self.responses)


def _file(file_id: str, size: int, *cases: tuple[str, str, str]) -> dict:
    return {
        "file_id": file_id,
        "file_size": size,
        "cases": [
            {
                "case_id": case_id,
                "project": {"project_id": project},
                "samples": [{"sample_type": sample_type}],
            }
            for case_id, project, sample_type in cases
        ],
    }


def _spec(tmp_path: Path) -> CohortSpec:
    return CohortSpec(
        name="multi",
        out_dir=tmp_path,
        filters={"project": ["TCGA-A", "TCGA-B"], "data_format": ["MAF"]},
        optional_omics=[
            OptionalOmicsSpec(name="rna", filters={"data_type": ["RNA"]}),
            OptionalOmicsSpec(name="protein", filters={"data_type": ["Protein"]}),
        ],
    )


def test_build_overlap_report_counts_selections_intersections_and_breakdowns(tmp_path: Path):
    client = FakeClient(
        [
            [
                _file("snv-1", 10, ("A", "TCGA-A", "Primary Tumor")),
                _file(
                    "snv-2",
                    20,
                    ("B", "TCGA-A", "Primary Tumor"),
                    ("C", "TCGA-B", "Metastatic"),
                ),
            ],
            [
                _file("rna-1", 30, ("B", "TCGA-A", "Primary Tumor")),
                _file("rna-2", 40, ("C", "TCGA-B", "Metastatic")),
            ],
            [_file("protein-1", 50, ("C", "TCGA-B", "Metastatic"))],
        ]
    )

    report = build_overlap_report(
        _spec(tmp_path),
        omics=["rna", "protein", "rna"],
        client=client,
        queried_at=datetime(2026, 7, 11, tzinfo=timezone.utc),
    )

    assert [item.name for item in report.selections] == ["primary", "rna", "protein"]
    assert report.selections[0].n_files == 2
    assert report.selections[0].n_cases == 3
    assert report.selections[0].total_size == 30
    assert report.selections[0].project_breakdown == {
        "TCGA-A": {"n_files": 2, "n_cases": 2, "total_size": 30},
        "TCGA-B": {"n_files": 1, "n_cases": 1, "total_size": 20},
    }
    assert report.selections[0].cases_by_sample_type == {"Metastatic": 1, "Primary Tumor": 2}
    assert [item.n_cases for item in report.pairwise] == [2, 1, 1]
    assert report.all_selected.n_cases == 1
    assert report.all_selected.cases_by_project == {"TCGA-B": 1}
    assert report.queried_at == "2026-07-11T00:00:00+00:00"
    assert len(client.filters) == 3


def test_build_overlap_report_requires_an_optional_selection(tmp_path: Path):
    try:
        build_overlap_report(_spec(tmp_path), omics=[], client=FakeClient([]))
    except ValueError as exc:
        assert "select at least one" in str(exc)
    else:
        raise AssertionError("expected ValueError")


def test_write_overlap_outputs_preserves_provenance(tmp_path: Path):
    report = build_overlap_report(
        _spec(tmp_path),
        omics=["rna"],
        client=FakeClient(
            [
                [_file("snv", 10, ("A", "TCGA-A", "Primary Tumor"))],
                [_file("rna", 20, ("A", "TCGA-A", "Primary Tumor"))],
            ]
        ),
        queried_at=datetime(2026, 7, 11, tzinfo=timezone.utc),
    )
    json_path = tmp_path / "report.json"
    parquet_path = tmp_path / "report.parquet"

    outputs = write_overlap_outputs(report, json_path=json_path, parquet_path=parquet_path)

    assert outputs.json_path == json_path.resolve()
    payload = json.loads(json_path.read_text())
    assert payload["selections"][0]["resolved_filter"]
    assert payload["all_selected"]["n_cases"] == 1
    df = pl.read_parquet(parquet_path)
    assert set(df["row_type"]) == {"selection", "pairwise", "all_selected"}
    assert df.filter(pl.col("row_type") == "selection").height == 2
