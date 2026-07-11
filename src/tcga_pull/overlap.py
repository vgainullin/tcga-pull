"""Metadata-only case overlap reporting across cohort selections."""

from __future__ import annotations

import json
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Any, Protocol

import polars as pl

from .config import CohortSpec
from .gdc import FILE_FIELDS, GDCClient


class OverlapClient(Protocol):
    def fetch_files(self, filters: dict, fields: list[str] | None = None) -> list[dict]: ...


@dataclass(frozen=True)
class SelectionSummary:
    name: str
    resolved_filter: dict
    n_files: int
    n_cases: int
    total_size: int
    project_breakdown: dict[str, dict[str, int]]
    cases_by_sample_type: dict[str, int]


@dataclass(frozen=True)
class IntersectionSummary:
    selections: tuple[str, ...]
    n_cases: int
    cases_by_project: dict[str, int]


@dataclass(frozen=True)
class OverlapReport:
    cohort: str
    queried_at: str
    selections: tuple[SelectionSummary, ...]
    pairwise: tuple[IntersectionSummary, ...]
    all_selected: IntersectionSummary


@dataclass(frozen=True)
class OverlapOutputs:
    json_path: Path | None = None
    parquet_path: Path | None = None


@dataclass
class _SelectionData:
    summary: SelectionSummary
    case_ids: set[str]
    projects_by_case: dict[str, set[str]]


def build_overlap_report(
    spec: CohortSpec,
    *,
    omics: Iterable[str],
    client: OverlapClient | None = None,
    queried_at: datetime | None = None,
) -> OverlapReport:
    """Resolve primary plus named optional omics and compare unique case IDs."""
    names = list(dict.fromkeys(omics))
    if not names:
        raise ValueError("select at least one optional omics entry")

    client = client or GDCClient()
    resolved: list[tuple[str, CohortSpec]] = [("primary", spec)]
    resolved.extend((name, spec.optional_omics_cohort(name)) for name in names)
    data = [_fetch_selection(name, selected, client) for name, selected in resolved]

    pairwise = tuple(_intersection(combo, data) for combo in combinations(range(len(data)), 2))
    all_selected = _intersection(tuple(range(len(data))), data)
    timestamp = (queried_at or datetime.now(timezone.utc)).astimezone(timezone.utc).isoformat()
    return OverlapReport(
        cohort=spec.name,
        queried_at=timestamp,
        selections=tuple(item.summary for item in data),
        pairwise=pairwise,
        all_selected=all_selected,
    )


def write_overlap_outputs(
    report: OverlapReport,
    *,
    json_path: Path | None = None,
    parquet_path: Path | None = None,
) -> OverlapOutputs:
    if json_path is None and parquet_path is None:
        raise ValueError("provide json_path and/or parquet_path")
    resolved_json = _resolved_output_path(json_path) if json_path else None
    resolved_parquet = _resolved_output_path(parquet_path) if parquet_path else None
    if resolved_json:
        resolved_json.parent.mkdir(parents=True, exist_ok=True)
        resolved_json.write_text(json.dumps(asdict(report), indent=2, sort_keys=True) + "\n")
    if resolved_parquet:
        resolved_parquet.parent.mkdir(parents=True, exist_ok=True)
        pl.DataFrame(_output_rows(report)).write_parquet(resolved_parquet)
    return OverlapOutputs(json_path=resolved_json, parquet_path=resolved_parquet)


def _fetch_selection(name: str, spec: CohortSpec, client: OverlapClient) -> _SelectionData:
    resolved_filter = spec.resolve_filter()
    hits = client.fetch_files(resolved_filter, fields=FILE_FIELDS)
    files: dict[str, dict] = {}
    case_ids: set[str] = set()
    projects_by_case: dict[str, set[str]] = defaultdict(set)
    sample_types_by_case: dict[str, set[str]] = defaultdict(set)
    file_ids_by_project: dict[str, set[str]] = defaultdict(set)
    bytes_by_project: dict[str, int] = defaultdict(int)

    for hit in hits:
        file_id = str(hit.get("file_id") or hit.get("id") or "")
        if file_id:
            files[file_id] = hit
        hit_projects: set[str] = set()
        for case in hit.get("cases") or []:
            case_id = case.get("case_id")
            if not case_id:
                continue
            case_ids.add(case_id)
            project_id = (case.get("project") or {}).get("project_id")
            if project_id:
                projects_by_case[case_id].add(project_id)
                hit_projects.add(project_id)
            for sample in case.get("samples") or []:
                sample_type = sample.get("sample_type")
                if sample_type:
                    sample_types_by_case[case_id].add(sample_type)
        for project_id in hit_projects:
            if file_id:
                file_ids_by_project[project_id].add(file_id)
            bytes_by_project[project_id] += int(hit.get("file_size") or 0)

    summary = SelectionSummary(
        name=name,
        resolved_filter=resolved_filter,
        n_files=len(files),
        n_cases=len(case_ids),
        total_size=sum(int(hit.get("file_size") or 0) for hit in files.values()),
        project_breakdown={
            project_id: {
                "n_files": len(file_ids_by_project[project_id]),
                "n_cases": n_cases,
                "total_size": bytes_by_project[project_id],
            }
            for project_id, n_cases in _group_counts(case_ids, projects_by_case).items()
        },
        cases_by_sample_type=_group_counts(case_ids, sample_types_by_case),
    )
    return _SelectionData(summary=summary, case_ids=case_ids, projects_by_case=projects_by_case)


def _intersection(indices: tuple[int, ...], data: list[_SelectionData]) -> IntersectionSummary:
    case_ids = set.intersection(*(data[index].case_ids for index in indices))
    projects: dict[str, set[str]] = defaultdict(set)
    for case_id in case_ids:
        for index in indices:
            projects[case_id].update(data[index].projects_by_case.get(case_id, set()))
    return IntersectionSummary(
        selections=tuple(data[index].summary.name for index in indices),
        n_cases=len(case_ids),
        cases_by_project=_group_counts(case_ids, projects),
    )


def _group_counts(case_ids: set[str], values_by_case: dict[str, set[str]]) -> dict[str, int]:
    counts: dict[str, int] = defaultdict(int)
    for case_id in case_ids:
        for value in values_by_case.get(case_id, set()):
            counts[value] += 1
    return dict(sorted(counts.items()))


def _output_rows(report: OverlapReport) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for item in report.selections:
        rows.append(
            {
                "row_type": "selection",
                "selections": item.name,
                "n_files": item.n_files,
                "n_cases": item.n_cases,
                "total_size": item.total_size,
                "queried_at": report.queried_at,
                "resolved_filters_json": json.dumps(item.resolved_filter, sort_keys=True),
                "project_breakdown_json": json.dumps(item.project_breakdown, sort_keys=True),
                "cases_by_sample_type_json": json.dumps(item.cases_by_sample_type, sort_keys=True),
            }
        )
    for row_type, items in (
        ("pairwise", report.pairwise),
        ("all_selected", (report.all_selected,)),
    ):
        for intersection in items:
            rows.append(
                {
                    "row_type": row_type,
                    "selections": ",".join(intersection.selections),
                    "n_files": None,
                    "n_cases": intersection.n_cases,
                    "total_size": None,
                    "queried_at": report.queried_at,
                    "resolved_filters_json": None,
                    "project_breakdown_json": json.dumps(
                        intersection.cases_by_project, sort_keys=True
                    ),
                    "cases_by_sample_type_json": None,
                }
            )
    return rows


def _resolved_output_path(path: Path) -> Path:
    return Path(path).expanduser().resolve()
