"""Cohort spec: YAML / dict-based config with a small sugar layer over GDC filters."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from .gdc import f_and, f_in, open_access

# sugar key -> GDC field name on the /files endpoint.
# Convention: bare names for file-rooted fields, `cases.X` for case joins,
# `analysis.X` for analysis joins. for_cases_endpoint() does the inverse
# translation when the filter is sent to /cases. The bare form for file
# fields matters: GDC's /files endpoint 500s on faceted queries when
# file-rooted clauses carry the `files.` prefix.
SUGAR_FIELDS: dict[str, str] = {
    "project": "cases.project.project_id",
    "data_type": "data_type",
    "data_category": "data_category",
    "data_format": "data_format",
    "experimental_strategy": "experimental_strategy",
    "workflow": "analysis.workflow_type",
    "sample_type": "cases.samples.sample_type",
    "tissue_type": "cases.samples.tissue_type",
    "primary_site": "cases.primary_site",
    "disease_type": "cases.disease_type",
    "case_id": "cases.case_id",
    "submitter_id": "cases.submitter_id",
}


@dataclass
class CohortSpec:
    name: str
    out_dir: Path
    filters: dict[str, Any] = field(default_factory=dict)
    gdc_filter: dict | None = None
    n_processes: int = 4

    @property
    def cohort_dir(self) -> Path:
        return self.out_dir / self.name

    def resolve_filter(self) -> dict:
        """Merge sugar `filters` with `gdc_filter` and force open access."""
        if self.gdc_filter:
            return open_access(self.gdc_filter)
        clauses = []
        for key, value in self.filters.items():
            if key not in SUGAR_FIELDS:
                raise ValueError(
                    f"unknown filter key: {key!r}. "
                    f"Known: {sorted(SUGAR_FIELDS)}. Use `gdc_filter` for raw fields."
                )
            clauses.append(f_in(SUGAR_FIELDS[key], value))
        return open_access(f_and(*clauses) if clauses else None)


def load_yaml(path: Path) -> CohortSpec:
    data = yaml.safe_load(path.read_text())
    return CohortSpec(
        name=data["name"],
        out_dir=Path(data.get("out_dir", "./cohorts")).expanduser().resolve(),
        filters=data.get("filters", {}) or {},
        gdc_filter=data.get("gdc_filter"),
        n_processes=int((data.get("download") or {}).get("n_processes", 4)),
    )


def from_flags(
    name: str,
    out_dir: Path,
    *,
    project: list[str] | None = None,
    data_type: list[str] | None = None,
    data_category: list[str] | None = None,
    data_format: list[str] | None = None,
    experimental_strategy: list[str] | None = None,
    workflow: list[str] | None = None,
    sample_type: list[str] | None = None,
    n_processes: int = 4,
) -> CohortSpec:
    filters: dict[str, Any] = {}
    if project:
        filters["project"] = project
    if data_type:
        filters["data_type"] = data_type
    if data_category:
        filters["data_category"] = data_category
    if data_format:
        filters["data_format"] = data_format
    if experimental_strategy:
        filters["experimental_strategy"] = experimental_strategy
    if workflow:
        filters["workflow"] = workflow
    if sample_type:
        filters["sample_type"] = sample_type
    return CohortSpec(name=name, out_dir=out_dir, filters=filters, n_processes=n_processes)


def read_projects_file(path: Path) -> list[str]:
    """One project_id per line. Lines starting with '#' and blank lines are ignored."""
    out: list[str] = []
    for raw in Path(path).read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        out.append(line)
    return out
