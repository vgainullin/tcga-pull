"""Offline tests for shared case selection across optional omics."""

from __future__ import annotations

import io
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import cast

import pytest
from rich.console import Console

from tcga_pull.case_set import (
    CASE_SET_ORDERING,
    build_shared_case_set,
    load_shared_case_set,
    write_shared_case_set,
)
from tcga_pull.config import CohortSpec, LimitSpec, OptionalOmicsSpec
from tcga_pull.gdc import GDCClient
from tcga_pull.pipeline import fetch_preview
from tcga_pull.pipeline import run as pipeline_run
from tcga_pull.services import SpecBuildError, apply_shared_case_set


def _leaf_values(node: dict, field: str) -> list[str] | None:
    content = node.get("content")
    if isinstance(content, list):
        for child in content:
            values = _leaf_values(child, field)
            if values is not None:
                return values
    if isinstance(content, dict) and content.get("field") == field:
        return list(content["value"])
    return None


def _file(modality: str, case_id: str, submitter_id: str, project_id: str) -> dict:
    return {
        "file_id": f"{modality}-{case_id}",
        "file_name": f"{modality}-{case_id}.tsv",
        "file_size": 100,
        "cases": [
            {
                "case_id": case_id,
                "submitter_id": submitter_id,
                "project": {"project_id": project_id},
            }
        ],
    }


class RoutingClient:
    def __init__(self) -> None:
        self.files = {
            "SNV": [
                _file("snv", "A1", "PA-001", "TCGA-A"),
                _file("snv", "A2", "PA-002", "TCGA-A"),
                _file("snv", "A3", "PA-003", "TCGA-A"),
                _file("snv", "B1", "PB-001", "TCGA-B"),
                _file("snv", "B2", "PB-002", "TCGA-B"),
            ],
            "RNA": [
                _file("rna", "A2", "PA-002", "TCGA-A"),
                _file("rna", "A3", "PA-003", "TCGA-A"),
                _file("rna", "B1", "PB-001", "TCGA-B"),
                _file("rna", "B2", "PB-002", "TCGA-B"),
                _file("rna", "B3", "PB-003", "TCGA-B"),
            ],
            "Protein": [
                _file("protein", "A1", "PA-001", "TCGA-A"),
                _file("protein", "A2", "PA-002", "TCGA-A"),
                _file("protein", "A3", "PA-003", "TCGA-A"),
                _file("protein", "B2", "PB-002", "TCGA-B"),
            ],
        }

    def fetch_files(self, filters: dict, fields: list[str] | None = None) -> list[dict]:
        modalities = _leaf_values(filters, "data_type")
        assert modalities
        hits = list(self.files[modalities[0]])
        selected = _leaf_values(filters, "cases.case_id")
        if selected is not None:
            selected_ids = set(selected)
            hits = [hit for hit in hits if hit["cases"][0]["case_id"] in selected_ids]
        return hits

    def fetch_clinical(self, filters: dict) -> list[dict]:
        selected = set(_leaf_values(filters, "cases.case_id") or [])
        cases = [hit["cases"][0] for hit in self.files["SNV"]]
        return [case for case in cases if case["case_id"] in selected]


def _spec(tmp_path: Path, *, limit: int | None = 2) -> CohortSpec:
    return CohortSpec(
        name="shared",
        out_dir=tmp_path,
        filters={"project": ["TCGA-A", "TCGA-B"], "data_type": "SNV"},
        limit=LimitSpec(per_project=limit),
        optional_omics=[
            OptionalOmicsSpec(name="rna", filters={"data_type": "RNA"}),
            OptionalOmicsSpec(name="protein", filters={"data_type": "Protein"}),
        ],
    )


def test_build_shared_case_set_selects_deterministic_n_way_intersection(tmp_path: Path):
    case_set = build_shared_case_set(
        _spec(tmp_path),
        omics=["rna", "protein", "rna"],
        client=RoutingClient(),
        created_at=datetime(2026, 7, 19, tzinfo=timezone.utc),
    )

    assert case_set.selections == ("primary", "rna", "protein")
    assert case_set.candidate_count == 3
    assert case_set.selected_count == 3
    assert case_set.candidate_counts_by_project == {"TCGA-A": 2, "TCGA-B": 1}
    assert case_set.selected_counts_by_project == {"TCGA-A": 2, "TCGA-B": 1}
    assert case_set.shortfalls_by_project == {"TCGA-B": {"requested": 2, "available": 1}}
    assert case_set.ordering == CASE_SET_ORDERING
    assert case_set.case_ids == ("A2", "A3", "B2")
    assert set(case_set.resolved_filters) == {"primary", "rna", "protein"}
    assert case_set.created_at == "2026-07-19T00:00:00+00:00"


def test_shared_case_set_caps_each_project_by_submitter_order(tmp_path: Path):
    case_set = build_shared_case_set(
        _spec(tmp_path, limit=1),
        omics=["rna", "protein"],
        client=RoutingClient(),
    )

    assert case_set.candidate_count == 3
    assert case_set.selected_count == 2
    assert case_set.case_ids == ("A2", "B2")
    assert case_set.selected_counts_by_project == {"TCGA-A": 1, "TCGA-B": 1}


def test_case_set_roundtrip_and_exact_preview_alignment(tmp_path: Path):
    client = RoutingClient()
    spec = _spec(tmp_path)
    case_set = build_shared_case_set(spec, omics=["rna", "protein"], client=client)
    artifact = write_shared_case_set(case_set, tmp_path / "shared-cases.json")

    assert load_shared_case_set(artifact) == case_set
    expected = {"A2", "A3", "B2"}
    for selection in case_set.selections:
        selected_spec = apply_shared_case_set(spec, selection=selection, path=artifact)
        preview = fetch_preview(selected_spec, client=cast(GDCClient, client))
        actual = {case["case_id"] for hit in preview.file_hits for case in hit.get("cases") or []}
        assert actual == expected
        assert selected_spec.limit.per_project is None
        assert selected_spec.case_set_provenance is not None
        assert selected_spec.case_set_provenance["selection"] == selection
        assert selected_spec.case_set_provenance["selected_count"] == 3
        assert len(selected_spec.case_set_provenance["sha256"]) == 64


def test_pipeline_records_case_set_provenance(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    client = RoutingClient()
    spec = _spec(tmp_path)
    case_set = build_shared_case_set(spec, omics=["rna", "protein"], client=client)
    artifact = write_shared_case_set(case_set, tmp_path / "shared-cases.json")
    selected_spec = apply_shared_case_set(spec, selection="primary", path=artifact)

    def fake_bulk(file_hits: list[dict], download_dir: Path, **kwargs) -> None:
        for hit in file_hits:
            path = download_dir / hit["file_id"] / hit["file_name"]
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("fixture")

    monkeypatch.setattr("tcga_pull.pipeline.bulk_download_via_api", fake_bulk)
    cohort_dir = pipeline_run(
        selected_spec,
        client=cast(GDCClient, client),
        console=Console(file=io.StringIO()),
    )

    provenance = json.loads((cohort_dir / "cohort.json").read_text())
    assert provenance["n_cases"] == 3
    assert provenance["case_set"]["selection"] == "primary"
    assert provenance["case_set"]["selected_count"] == 3
    assert provenance["case_set"]["source"] == str(artifact)


def test_case_set_preview_rejects_missing_selected_cases(tmp_path: Path):
    client = RoutingClient()
    spec = _spec(tmp_path)
    case_set = build_shared_case_set(spec, omics=["rna", "protein"], client=client)
    artifact = write_shared_case_set(case_set, tmp_path / "shared-cases.json")
    selected_spec = apply_shared_case_set(spec, selection="primary", path=artifact)
    client.files["SNV"] = [hit for hit in client.files["SNV"] if hit["cases"][0]["case_id"] != "B2"]

    with pytest.raises(ValueError, match=r"missing=\['B2'\]"):
        fetch_preview(selected_spec, client=cast(GDCClient, client))


def test_shared_case_set_requires_limit(tmp_path: Path):
    with pytest.raises(ValueError, match=r"limit\.per_project is required"):
        build_shared_case_set(_spec(tmp_path, limit=None), omics=["rna"], client=RoutingClient())


def test_shared_case_set_reports_projects_with_zero_intersection(tmp_path: Path):
    client = RoutingClient()
    client.files["Protein"] = [
        hit
        for hit in client.files["Protein"]
        if hit["cases"][0]["project"]["project_id"] == "TCGA-A"
    ]

    case_set = build_shared_case_set(_spec(tmp_path), omics=["rna", "protein"], client=client)

    assert case_set.candidate_counts_by_project == {"TCGA-A": 2, "TCGA-B": 0}
    assert case_set.selected_counts_by_project == {"TCGA-A": 2, "TCGA-B": 0}
    assert case_set.shortfalls_by_project["TCGA-B"] == {"requested": 2, "available": 0}


def test_apply_shared_case_set_rejects_wrong_cohort_or_selection(tmp_path: Path):
    spec = _spec(tmp_path)
    case_set = build_shared_case_set(spec, omics=["rna"], client=RoutingClient())
    artifact = write_shared_case_set(case_set, tmp_path / "shared-cases.json")

    with pytest.raises(SpecBuildError, match="not in case set"):
        apply_shared_case_set(spec, selection="protein", path=artifact)

    payload = json.loads(artifact.read_text())
    payload["cohort"] = "other"
    artifact.write_text(json.dumps(payload))
    with pytest.raises(SpecBuildError, match="belongs to cohort"):
        apply_shared_case_set(spec, selection="primary", path=artifact)


def test_load_shared_case_set_rejects_unknown_schema_version(tmp_path: Path):
    artifact = tmp_path / "future.json"
    artifact.write_text(json.dumps({"schema_version": 999}))

    with pytest.raises(ValueError, match="unsupported case-set schema_version"):
        load_shared_case_set(artifact)


def test_load_shared_case_set_rejects_inconsistent_selected_count(tmp_path: Path):
    artifact = write_shared_case_set(_shared_case_set_fixture(), tmp_path / "invalid.json")
    payload = json.loads(artifact.read_text())
    payload["selected_count"] = 99
    artifact.write_text(json.dumps(payload))

    with pytest.raises(ValueError, match="selected_count does not match"):
        load_shared_case_set(artifact)


def _shared_case_set_fixture():
    return build_shared_case_set(
        _spec(Path("/tmp"), limit=1),
        omics=["rna"],
        client=RoutingClient(),
    )
