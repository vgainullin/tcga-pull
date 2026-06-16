"""Open-access GDC coverage matrix for TCGA projects.

The matrix is an inventory, not a downloader. It queries GDC file metadata,
groups files by project and data-type fields, and classifies each group against
the recipes tcga-pull can normalize today.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol

import polars as pl

from .gdc import GDCClient, f_in, open_access

COVERAGE_FILE_FIELDS: list[str] = [
    "file_id",
    "file_name",
    "data_category",
    "data_type",
    "data_format",
    "experimental_strategy",
    "analysis.workflow_type",
    "access",
    "file_size",
    "cases.case_id",
    "cases.submitter_id",
    "cases.project.project_id",
]

COVERAGE_SCHEMA: dict[str, type[pl.DataType]] = {
    "program": pl.Utf8,
    "project_id": pl.Utf8,
    "data_category": pl.Utf8,
    "data_type": pl.Utf8,
    "data_format": pl.Utf8,
    "experimental_strategy": pl.Utf8,
    "workflow_type": pl.Utf8,
    "access": pl.Utf8,
    "n_files": pl.Int64,
    "n_cases": pl.Int64,
    "total_size": pl.Int64,
    "support_status": pl.Utf8,
    "recipe": pl.Utf8,
    "outputs": pl.Utf8,
    "notes": pl.Utf8,
}


class CoverageClient(Protocol):
    def list_projects(self, program: str = "TCGA") -> list[dict]: ...

    def fetch_files(self, filters: dict, fields: list[str] | None = None) -> list[dict]: ...


@dataclass(frozen=True)
class SupportClassification:
    status: str
    recipe: str | None
    outputs: tuple[str, ...]
    notes: str


@dataclass(frozen=True)
class CoverageMatrix:
    program: str
    project_ids: tuple[str, ...]
    rows: tuple[dict[str, Any], ...]


@dataclass(frozen=True)
class CoverageOutputs:
    parquet_path: Path
    markdown_path: Path
    n_rows: int
    n_projects: int


def classify_file_group(
    *,
    data_category: str | None,
    data_type: str | None,
    data_format: str | None,
    experimental_strategy: str | None,
    workflow_type: str | None,
) -> SupportClassification:
    """Classify a GDC file group against currently implemented recipes."""
    category = data_category or ""
    dtype = data_type or ""
    fmt = data_format or ""
    strategy = experimental_strategy or ""
    workflow = workflow_type or ""

    if category == "Simple Nucleotide Variation" and fmt == "MAF":
        return SupportClassification(
            "supported",
            "variants",
            ("variants.parquet",),
            "Masked somatic mutation MAFs are normalized into variant rows.",
        )
    if (
        category == "Transcriptome Profiling"
        and dtype == "Gene Expression Quantification"
        and strategy == "RNA-Seq"
        and workflow == "STAR - Counts"
    ):
        return SupportClassification(
            "supported",
            "rna_expression",
            ("rna_expression.parquet",),
            "STAR count files are normalized to long-form gene expression.",
        )
    if (
        category == "Transcriptome Profiling"
        and dtype == "miRNA Expression Quantification"
        and strategy == "miRNA-Seq"
    ):
        return SupportClassification(
            "supported",
            "mirna_expression",
            ("mirna_expression.parquet",),
            "miRNA quantification files are normalized to long-form miRNA expression.",
        )
    if category == "DNA Methylation" and dtype == "Methylation Beta Value":
        return SupportClassification(
            "supported",
            "methylation",
            ("methylation_beta.parquet",),
            "Methylation beta-value tables are normalized by probe.",
        )
    if category == "Copy Number Variation" and dtype in {
        "Copy Number Segment",
        "Masked Copy Number Segment",
    }:
        return SupportClassification(
            "supported",
            "copy_number",
            ("copy_number_segments.parquet",),
            "Segment-level CNV files are normalized by genomic segment.",
        )
    if category == "Copy Number Variation" and dtype == "Gene Level Copy Number":
        return SupportClassification(
            "supported",
            "copy_number",
            ("gene_copy_number.parquet",),
            "Gene-level CNV files are normalized by gene.",
        )
    if category == "Proteome Profiling" and dtype == "Protein Expression Quantification":
        return SupportClassification(
            "supported",
            "protein_expression",
            ("protein_expression.parquet",),
            "RPPA protein expression files are normalized by target.",
        )
    if category == "Clinical":
        return SupportClassification(
            "raw_only",
            None,
            (),
            "Core case clinical metadata is fetched from /cases; clinical supplement files are retained raw.",
        )
    return SupportClassification(
        "raw_only",
        None,
        (),
        "Open-access files can be downloaded and tracked, but no normalizing recipe exists yet.",
    )


def build_coverage_matrix(
    *,
    program: str = "TCGA",
    projects: Iterable[str] | None = None,
    client: CoverageClient | None = None,
) -> CoverageMatrix:
    client = client or GDCClient()
    project_ids = _project_ids(program=program, projects=projects, client=client)
    if not project_ids:
        return CoverageMatrix(program=program, project_ids=(), rows=())

    filters = open_access(f_in("cases.project.project_id", project_ids))
    file_hits = client.fetch_files(filters, fields=COVERAGE_FILE_FIELDS)
    rows = _matrix_rows(program, project_ids, file_hits)
    return CoverageMatrix(program=program, project_ids=tuple(project_ids), rows=tuple(rows))


def coverage_dataframe(matrix: CoverageMatrix) -> pl.DataFrame:
    return pl.DataFrame(matrix.rows, schema=COVERAGE_SCHEMA)


def write_coverage_outputs(matrix: CoverageMatrix, out_dir: Path) -> CoverageOutputs:
    out_dir = Path(out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{matrix.program.lower()}_open_access_coverage_matrix"
    parquet_path = out_dir / f"{stem}.parquet"
    markdown_path = out_dir / f"{stem}.md"

    df = coverage_dataframe(matrix)
    df.write_parquet(parquet_path)
    markdown_path.write_text(render_coverage_markdown(matrix))
    return CoverageOutputs(
        parquet_path=parquet_path,
        markdown_path=markdown_path,
        n_rows=len(matrix.rows),
        n_projects=len(matrix.project_ids),
    )


def render_coverage_markdown(matrix: CoverageMatrix) -> str:
    rows = list(matrix.rows)
    supported = sum(1 for row in rows if row["support_status"] == "supported")
    raw_only = sum(1 for row in rows if row["support_status"] == "raw_only")

    lines = [
        f"# {matrix.program} Open-Access Coverage Matrix",
        "",
        f"- Projects: {len(matrix.project_ids)}",
        f"- File groups: {len(rows)}",
        f"- Supported groups: {supported}",
        f"- Raw-only groups: {raw_only}",
        "",
        "| project_id | data_category | data_type | data_format | strategy | workflow | files | cases | status | recipe |",
        "|---|---|---|---|---|---|---:|---:|---|---|",
    ]
    for row in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    _md(row["project_id"]),
                    _md(row["data_category"]),
                    _md(row["data_type"]),
                    _md(row["data_format"]),
                    _md(row["experimental_strategy"]),
                    _md(row["workflow_type"]),
                    str(row["n_files"]),
                    str(row["n_cases"]),
                    _md(row["support_status"]),
                    _md(row["recipe"]),
                ]
            )
            + " |"
        )
    lines.append("")
    return "\n".join(lines)


def _project_ids(
    *,
    program: str,
    projects: Iterable[str] | None,
    client: CoverageClient,
) -> list[str]:
    if projects is not None:
        return sorted(dict.fromkeys(p for p in projects if p))
    rows = client.list_projects(program=program)
    return sorted(
        project_id
        for project_id in (row.get("project_id") for row in rows)
        if isinstance(project_id, str) and project_id
    )


def _matrix_rows(program: str, project_ids: list[str], file_hits: list[dict]) -> list[dict]:
    project_filter = set(project_ids)
    grouped: dict[tuple, dict[str, Any]] = {}
    for hit in file_hits:
        workflow = (hit.get("analysis") or {}).get("workflow_type")
        base = (
            _clean(hit.get("data_category")),
            _clean(hit.get("data_type")),
            _clean(hit.get("data_format")),
            _clean(hit.get("experimental_strategy")),
            _clean(workflow),
            _clean(hit.get("access")) or "open",
        )
        file_size = int(hit.get("file_size") or 0)
        file_id = str(hit.get("file_id") or "")
        for project_id, case_ids in _project_case_ids(hit, project_filter):
            key = (project_id, *base)
            bucket = grouped.setdefault(
                key,
                {
                    "program": program,
                    "project_id": project_id,
                    "data_category": base[0],
                    "data_type": base[1],
                    "data_format": base[2],
                    "experimental_strategy": base[3],
                    "workflow_type": base[4],
                    "access": base[5],
                    "_file_ids": set(),
                    "_case_ids": set(),
                    "total_size": 0,
                },
            )
            if file_id not in bucket["_file_ids"]:
                bucket["_file_ids"].add(file_id)
                bucket["total_size"] += file_size
            bucket["_case_ids"].update(case_ids)

    rows: list[dict] = []
    for bucket in grouped.values():
        support = classify_file_group(
            data_category=bucket["data_category"],
            data_type=bucket["data_type"],
            data_format=bucket["data_format"],
            experimental_strategy=bucket["experimental_strategy"],
            workflow_type=bucket["workflow_type"],
        )
        rows.append(
            {
                "program": bucket["program"],
                "project_id": bucket["project_id"],
                "data_category": bucket["data_category"],
                "data_type": bucket["data_type"],
                "data_format": bucket["data_format"],
                "experimental_strategy": bucket["experimental_strategy"],
                "workflow_type": bucket["workflow_type"],
                "access": bucket["access"],
                "n_files": len(bucket["_file_ids"]),
                "n_cases": len(bucket["_case_ids"]),
                "total_size": bucket["total_size"],
                "support_status": support.status,
                "recipe": support.recipe,
                "outputs": ", ".join(support.outputs),
                "notes": support.notes,
            }
        )

    return sorted(
        rows,
        key=lambda row: (
            row["project_id"],
            row["data_category"],
            row["data_type"],
            row["data_format"],
            row["experimental_strategy"],
            row["workflow_type"],
        ),
    )


def _project_case_ids(hit: dict, project_filter: set[str]) -> list[tuple[str, set[str]]]:
    by_project: dict[str, set[str]] = {}
    for case in hit.get("cases") or []:
        project_id = ((case.get("project") or {}).get("project_id")) or ""
        if not project_id or project_id not in project_filter:
            continue
        case_id = case.get("case_id")
        by_project.setdefault(project_id, set())
        if case_id:
            by_project[project_id].add(str(case_id))
    return sorted(by_project.items())


def _clean(value: object) -> str | None:
    if value is None:
        return None
    text = str(value)
    return text if text else None


def _md(value: object) -> str:
    if value is None:
        return ""
    return str(value).replace("|", "\\|")
