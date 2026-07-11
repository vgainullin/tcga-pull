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


# Recipes that can run after a successful pull. Names map to functions in
# pipeline.RECIPE_REGISTRY (registered there to avoid circular imports).
KNOWN_RECIPES: tuple[str, ...] = (
    "variants",
    "samples",
    "frequency",
    "rna_expression",
    "mirna_expression",
    "methylation",
    "copy_number",
    "protein_expression",
    "multiomics",
    "model_dataset",
)


@dataclass
class LimitSpec:
    """Optional per-cohort sampling limits applied at file-list stage,
    before any bytes are downloaded.

    `per_project`: keep at most N unique cases per GDC project_id, chosen
    deterministically by submitter_id sort. Useful for prototyping or for
    building balanced training subsets across heterogeneous lineages.
    """

    per_project: int | None = None

    def __post_init__(self) -> None:
        if self.per_project is not None and self.per_project <= 0:
            raise ValueError(f"limit.per_project must be > 0, got {self.per_project!r}")


@dataclass
class ProcessingSpec:
    """Post-download processing behavior.

    Standard mode downloads/restructures the whole cohort, then runs recipes.
    Incremental mode downloads in batches, lets batch-aware recipes process each
    batch, and can delete handled raw files before downloading the next batch.
    """

    mode: str = "standard"
    batch_size: int = 200
    delete_raw_after_processing: bool = False

    def __post_init__(self) -> None:
        if self.mode not in {"standard", "incremental"}:
            raise ValueError("processing.mode must be 'standard' or 'incremental'")
        if self.batch_size <= 0:
            raise ValueError(f"processing.batch_size must be > 0, got {self.batch_size!r}")


@dataclass
class OptionalOmicsSpec:
    """Named add-on dataset that can be pulled alongside a primary cohort.

    Optional omics inherit the parent cohort's project filter when they do not
    declare their own `project` filter. They are otherwise ordinary cohort
    specs, so they preview/download into a separate cohort directory and can be
    merged downstream by case_id or submitter_id.
    """

    name: str
    filters: dict[str, Any] = field(default_factory=dict)
    gdc_filter: dict | None = None
    notes: str | None = None
    recipes: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        if not self.filters and not self.gdc_filter:
            raise ValueError(f"optional omics {self.name!r} must define filters or gdc_filter")
        bad = [r for r in self.recipes if r not in KNOWN_RECIPES]
        if bad:
            raise ValueError(f"unknown recipe(s): {bad}. Known: {list(KNOWN_RECIPES)}")


@dataclass
class CohortSpec:
    name: str
    out_dir: Path
    filters: dict[str, Any] = field(default_factory=dict)
    gdc_filter: dict | None = None
    n_processes: int = 4
    recipes: list[str] = field(default_factory=list)
    limit: LimitSpec = field(default_factory=LimitSpec)
    processing: ProcessingSpec = field(default_factory=ProcessingSpec)
    recipe_options: dict[str, Any] = field(default_factory=dict)
    optional_omics: list[OptionalOmicsSpec] = field(default_factory=list)

    def __post_init__(self) -> None:
        bad = [r for r in self.recipes if r not in KNOWN_RECIPES]
        if bad:
            raise ValueError(f"unknown recipe(s): {bad}. Known: {list(KNOWN_RECIPES)}")
        names = [omics.name for omics in self.optional_omics]
        dupes = sorted({name for name in names if names.count(name) > 1})
        if dupes:
            raise ValueError(f"duplicate optional omics name(s): {dupes}")

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

    def optional_omics_cohort(self, name: str) -> CohortSpec:
        """Return a concrete CohortSpec for one optional omics add-on."""
        matches = [omics for omics in self.optional_omics if omics.name == name]
        if not matches:
            known = ", ".join(sorted(omics.name for omics in self.optional_omics)) or "none"
            raise ValueError(f"unknown optional omics: {name!r}. Known: {known}")
        omics = matches[0]
        filters = dict(omics.filters)
        if omics.gdc_filter is None and "project" not in filters and "project" in self.filters:
            filters = {"project": self.filters["project"], **filters}
        return CohortSpec(
            name=f"{self.name}__{omics.name}",
            out_dir=self.out_dir,
            filters=filters,
            gdc_filter=omics.gdc_filter,
            n_processes=self.n_processes,
            recipes=list(omics.recipes),
            limit=self.limit,
            processing=self.processing,
            recipe_options=self.recipe_options,
        )


def load_yaml(path: Path, *, out_dir_override: Path | None = None) -> CohortSpec:
    """Load a cohort YAML. If `out_dir_override` is given (e.g. from a CLI
    `--out` flag), it wins over the YAML's `out_dir` field."""
    data = yaml.safe_load(path.read_text())
    if out_dir_override is not None:
        out_dir = out_dir_override.expanduser().resolve()
    else:
        out_dir = Path(data.get("out_dir", "./cohorts")).expanduser().resolve()
    limit_block = data.get("limit") or {}
    limit = LimitSpec(per_project=limit_block.get("per_project"))
    processing_block = data.get("processing") or {}
    processing = ProcessingSpec(
        mode=processing_block.get("mode", "standard"),
        batch_size=int(processing_block.get("batch_size", 200)),
        delete_raw_after_processing=bool(
            processing_block.get("delete_raw_after_processing", False)
        ),
    )
    optional_omics = [
        OptionalOmicsSpec(
            name=item["name"],
            filters=item.get("filters", {}) or {},
            gdc_filter=item.get("gdc_filter"),
            notes=item.get("notes") or item.get("description"),
            recipes=list(item.get("recipes") or []),
        )
        for item in data.get("optional_omics") or []
    ]
    return CohortSpec(
        name=data["name"],
        out_dir=out_dir,
        filters=data.get("filters", {}) or {},
        gdc_filter=data.get("gdc_filter"),
        n_processes=int((data.get("download") or {}).get("n_processes", 4)),
        recipes=list(data.get("recipes") or []),
        limit=limit,
        processing=processing,
        recipe_options=data.get("recipe_options", {}) or {},
        optional_omics=optional_omics,
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
    processing: ProcessingSpec | None = None,
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
    return CohortSpec(
        name=name,
        out_dir=out_dir,
        filters=filters,
        n_processes=n_processes,
        processing=processing or ProcessingSpec(),
    )


def read_projects_file(path: Path) -> list[str]:
    """One project_id per line. Lines starting with '#' and blank lines are ignored."""
    out: list[str] = []
    for raw in Path(path).read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        out.append(line)
    return out
