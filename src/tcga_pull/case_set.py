"""Deterministic shared case selection across cohort modalities."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Protocol

from .config import CohortSpec
from .gdc import FILE_FIELDS, GDCClient

CASE_SET_SCHEMA_VERSION = 1
CASE_SET_ORDERING = "project_id, submitter_id, case_id ascending"


class CaseSetClient(Protocol):
    def fetch_files(self, filters: dict, fields: list[str] | None = None) -> list[dict]: ...


@dataclass(frozen=True)
class SharedCase:
    case_id: str
    submitter_id: str
    project_id: str


@dataclass(frozen=True)
class SharedCaseSet:
    schema_version: int
    cohort: str
    created_at: str
    selections: tuple[str, ...]
    resolved_filters: dict[str, dict]
    requested_per_project: int
    candidate_count: int
    selected_count: int
    candidate_counts_by_project: dict[str, int]
    selected_counts_by_project: dict[str, int]
    shortfalls_by_project: dict[str, dict[str, int]]
    ordering: str
    cases: tuple[SharedCase, ...]

    @property
    def case_ids(self) -> tuple[str, ...]:
        return tuple(case.case_id for case in self.cases)


def build_shared_case_set(
    spec: CohortSpec,
    *,
    omics: Iterable[str],
    client: CaseSetClient | None = None,
    created_at: datetime | None = None,
) -> SharedCaseSet:
    """Select the first N cases per project from an N-way modality intersection."""
    names = list(dict.fromkeys(omics))
    if not names:
        raise ValueError("select at least one optional omics entry")
    requested = spec.limit.per_project
    if requested is None:
        raise ValueError("limit.per_project is required for shared case selection")

    client = client or GDCClient()
    resolved: list[tuple[str, CohortSpec]] = [("primary", spec)]
    resolved.extend((name, spec.optional_omics_cohort(name)) for name in names)

    resolved_filters: dict[str, dict] = {}
    cases_by_selection: list[dict[str, SharedCase]] = []
    for name, selected_spec in resolved:
        selected_filter = selected_spec.resolve_filter()
        resolved_filters[name] = selected_filter
        hits = client.fetch_files(selected_filter, fields=FILE_FIELDS)
        cases_by_selection.append(_cases_from_hits(hits))

    candidate_ids = set.intersection(*(set(cases) for cases in cases_by_selection))
    candidates = [_merge_case(case_id, cases_by_selection) for case_id in candidate_ids]
    candidates.sort(key=_case_sort_key)

    project_ids = sorted({case.project_id for case in cases_by_selection[0].values()})
    candidate_counts_raw = _counts_by_project(candidates)
    candidate_counts = {
        project_id: candidate_counts_raw.get(project_id, 0)
        for project_id in sorted(set(project_ids) | set(candidate_counts_raw))
    }
    selected: list[SharedCase] = []
    selected_per_project: Counter[str] = Counter()
    for case in candidates:
        if selected_per_project[case.project_id] >= requested:
            continue
        selected.append(case)
        selected_per_project[case.project_id] += 1

    selected_counts_raw = _counts_by_project(selected)
    selected_counts = {
        project_id: selected_counts_raw.get(project_id, 0) for project_id in candidate_counts
    }
    shortfalls = {
        project_id: {"requested": requested, "available": available}
        for project_id, available in candidate_counts.items()
        if available < requested
    }
    timestamp = (created_at or datetime.now(timezone.utc)).astimezone(timezone.utc).isoformat()
    return SharedCaseSet(
        schema_version=CASE_SET_SCHEMA_VERSION,
        cohort=spec.name,
        created_at=timestamp,
        selections=tuple(name for name, _ in resolved),
        resolved_filters=resolved_filters,
        requested_per_project=requested,
        candidate_count=len(candidates),
        selected_count=len(selected),
        candidate_counts_by_project=candidate_counts,
        selected_counts_by_project=selected_counts,
        shortfalls_by_project=shortfalls,
        ordering=CASE_SET_ORDERING,
        cases=tuple(selected),
    )


def write_shared_case_set(case_set: SharedCaseSet, path: Path) -> Path:
    _validate_shared_case_set(case_set)
    resolved = Path(path).expanduser().resolve()
    resolved.parent.mkdir(parents=True, exist_ok=True)
    resolved.write_text(json.dumps(asdict(case_set), indent=2, sort_keys=True) + "\n")
    return resolved


def load_shared_case_set(path: Path) -> SharedCaseSet:
    resolved = Path(path).expanduser().resolve()
    payload = json.loads(resolved.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"invalid case-set artifact: {resolved}")
    version = payload.get("schema_version")
    if version != CASE_SET_SCHEMA_VERSION:
        raise ValueError(
            f"unsupported case-set schema_version {version!r}; expected {CASE_SET_SCHEMA_VERSION}"
        )
    try:
        cases = tuple(SharedCase(**case) for case in payload["cases"])
        case_set = SharedCaseSet(
            schema_version=version,
            cohort=payload["cohort"],
            created_at=payload["created_at"],
            selections=tuple(payload["selections"]),
            resolved_filters=payload["resolved_filters"],
            requested_per_project=int(payload["requested_per_project"]),
            candidate_count=int(payload["candidate_count"]),
            selected_count=int(payload["selected_count"]),
            candidate_counts_by_project=payload["candidate_counts_by_project"],
            selected_counts_by_project=payload["selected_counts_by_project"],
            shortfalls_by_project=payload["shortfalls_by_project"],
            ordering=payload["ordering"],
            cases=cases,
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"invalid case-set artifact: {resolved}") from exc
    _validate_shared_case_set(case_set)
    return case_set


def case_set_sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).expanduser().resolve().read_bytes()).hexdigest()


def _validate_shared_case_set(case_set: SharedCaseSet) -> None:
    if not case_set.cohort or not case_set.selections or case_set.selections[0] != "primary":
        raise ValueError("case set must identify a cohort and start with the primary selection")
    if len(set(case_set.selections)) != len(case_set.selections):
        raise ValueError("case set contains duplicate selections")
    if case_set.requested_per_project <= 0:
        raise ValueError("case set requested_per_project must be positive")
    if case_set.selected_count != len(case_set.cases):
        raise ValueError("case set selected_count does not match cases")
    if case_set.candidate_count < case_set.selected_count:
        raise ValueError("case set candidate_count is smaller than selected_count")
    if sum(case_set.candidate_counts_by_project.values()) != case_set.candidate_count:
        raise ValueError("case set candidate project counts do not match candidate_count")
    if sum(case_set.selected_counts_by_project.values()) != case_set.selected_count:
        raise ValueError("case set selected project counts do not match selected_count")
    if len(set(case_set.case_ids)) != len(case_set.case_ids):
        raise ValueError("case set contains duplicate case IDs")
    if any(
        not case.case_id or not case.submitter_id or not case.project_id for case in case_set.cases
    ):
        raise ValueError("case set contains incomplete case metadata")
    if set(case_set.selections) != set(case_set.resolved_filters):
        raise ValueError("case set selections do not match resolved_filters")
    if case_set.ordering != CASE_SET_ORDERING:
        raise ValueError(f"unsupported case set ordering: {case_set.ordering!r}")


def _cases_from_hits(hits: list[dict]) -> dict[str, SharedCase]:
    cases: dict[str, SharedCase] = {}
    for hit in hits:
        for case in hit.get("cases") or []:
            case_id = case.get("case_id")
            if not case_id:
                continue
            submitter_id = case.get("submitter_id") or case_id
            project_id = (case.get("project") or {}).get("project_id") or "unknown"
            cases[case_id] = SharedCase(
                case_id=str(case_id),
                submitter_id=str(submitter_id),
                project_id=str(project_id),
            )
    return cases


def _merge_case(case_id: str, cases_by_selection: list[dict[str, SharedCase]]) -> SharedCase:
    records = [cases[case_id] for cases in cases_by_selection]
    submitters = [record.submitter_id for record in records if record.submitter_id != case_id]
    projects = [record.project_id for record in records if record.project_id != "unknown"]
    submitter_id = min(submitters) if submitters else case_id
    project_id = min(projects) if projects else "unknown"
    return SharedCase(case_id=case_id, submitter_id=submitter_id, project_id=project_id)


def _case_sort_key(case: SharedCase) -> tuple[str, str, str]:
    return (case.project_id, case.submitter_id, case.case_id)


def _counts_by_project(cases: Iterable[SharedCase]) -> dict[str, int]:
    counts = Counter(case.project_id for case in cases)
    return dict(sorted(counts.items()))
