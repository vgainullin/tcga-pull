"""Find open-access omics that overlap a pancancer SNV cohort by patient.

This is metadata-only: it queries GDC /cases for case IDs, but does not
download any data files.

Example:
    uv run python examples/pancancer_omics_overlap.py examples/pancancer_snv.yaml
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from tcga_pull.config import load_yaml
from tcga_pull.gdc import GDCClient, f_and, f_in, open_access


@dataclass(frozen=True)
class OmicsCandidate:
    name: str
    filters: dict[str, Any]
    notes: str = ""


DEFAULT_CANDIDATES: tuple[OmicsCandidate, ...] = (
    OmicsCandidate(
        "rna_gene_expression_star_counts",
        {
            "data_category": "Transcriptome Profiling",
            "data_type": "Gene Expression Quantification",
            "experimental_strategy": "RNA-Seq",
            "analysis.workflow_type": "STAR - Counts",
        },
        "GDC harmonized gene-level RNA-seq counts.",
    ),
    OmicsCandidate(
        "mirna_expression",
        {
            "data_category": "Transcriptome Profiling",
            "data_type": "miRNA Expression Quantification",
            "experimental_strategy": "miRNA-Seq",
        },
        "Patient-level miRNA abundance; workflow names vary by project.",
    ),
    OmicsCandidate(
        "dna_methylation_beta",
        {
            "data_category": "DNA Methylation",
            "data_type": "Methylation Beta Value",
        },
        "Beta values; TCGA includes 27K/450K-era array data and newer harmonized releases.",
    ),
    OmicsCandidate(
        "masked_copy_number_segments",
        {
            "data_category": "Copy Number Variation",
            "data_type": "Masked Copy Number Segment",
        },
        "Segment-level CNV suitable for broad/focal event derivation.",
    ),
    OmicsCandidate(
        "copy_number_segments",
        {
            "data_category": "Copy Number Variation",
            "data_type": "Copy Number Segment",
        },
        "Unmasked/legacy segment-level CNV where available.",
    ),
    OmicsCandidate(
        "gene_level_copy_number",
        {
            "data_category": "Copy Number Variation",
            "data_type": "Gene Level Copy Number",
        },
        "Gene-level CNV scores where GDC exposes them as open access.",
    ),
    OmicsCandidate(
        "protein_expression_rppa",
        {
            "data_category": "Proteome Profiling",
            "data_type": "Protein Expression Quantification",
            "experimental_strategy": "Reverse Phase Protein Array",
        },
        "TCGA RPPA protein expression.",
    ),
)


def _filter_from_dict(filters: dict[str, Any]) -> dict:
    return f_and(*(f_in(field, value) for field, value in filters.items()))


def _project_filter(projects: list[str]) -> dict:
    return f_in("cases.project.project_id", projects)


def _case_rows(client: GDCClient, filters: dict) -> list[dict[str, Any]]:
    fields = ["case_id", "submitter_id", "project.project_id"]
    return list(client.iter_cases(filters, fields=fields, expand=[], page_size=1000))


def _case_ids(rows: list[dict[str, Any]]) -> set[str]:
    return {row["case_id"] for row in rows if row.get("case_id")}


def _submitter_ids(rows: list[dict[str, Any]]) -> set[str]:
    return {row["submitter_id"] for row in rows if row.get("submitter_id")}


def _summarize_candidate(
    *,
    client: GDCClient,
    projects: list[str],
    baseline_case_ids: set[str],
    baseline_submitter_ids: set[str],
    candidate: OmicsCandidate,
) -> dict[str, Any]:
    candidate_filter = open_access(
        f_and(_project_filter(projects), _filter_from_dict(candidate.filters))
    )
    return _summarize_filter(
        client=client,
        baseline_case_ids=baseline_case_ids,
        baseline_submitter_ids=baseline_submitter_ids,
        name=candidate.name,
        filters=candidate_filter,
        notes=candidate.notes,
        candidate_filters=candidate.filters,
    )


def _summarize_filter(
    *,
    client: GDCClient,
    baseline_case_ids: set[str],
    baseline_submitter_ids: set[str],
    name: str,
    filters: dict,
    notes: str,
    candidate_filters: dict[str, Any],
) -> dict[str, Any]:
    candidate_filter = filters
    rows = _case_rows(client, candidate_filter)
    candidate_case_ids = _case_ids(rows)
    candidate_submitter_ids = _submitter_ids(rows)
    case_overlap = candidate_case_ids & baseline_case_ids
    submitter_overlap = candidate_submitter_ids & baseline_submitter_ids
    return {
        "name": name,
        "overlap_cases": len(case_overlap),
        "overlap_submitter_ids": len(submitter_overlap),
        "available_cases_in_projects": len(candidate_case_ids),
        "baseline_cases": len(baseline_case_ids),
        "overlap_pct": round(100 * len(case_overlap) / len(baseline_case_ids), 2)
        if baseline_case_ids
        else 0.0,
        "filters": candidate_filters,
        "notes": notes,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cohort_yaml", type=Path)
    parser.add_argument("--json", action="store_true", help="Emit JSON instead of TSV.")
    args = parser.parse_args()

    spec = load_yaml(args.cohort_yaml)
    raw_projects = spec.filters.get("project") or []
    projects = [raw_projects] if isinstance(raw_projects, str) else list(raw_projects)
    if not projects:
        raise SystemExit("cohort YAML must include filters.project")

    client = GDCClient(timeout=120)
    baseline_rows = _case_rows(client, spec.resolve_filter())
    baseline_case_ids = _case_ids(baseline_rows)
    baseline_submitter_ids = _submitter_ids(baseline_rows)

    if spec.optional_omics:
        rows = []
        for item in spec.optional_omics:
            omics = spec.optional_omics_cohort(item.name)
            rows.append(
                _summarize_filter(
                    client=client,
                    baseline_case_ids=baseline_case_ids,
                    baseline_submitter_ids=baseline_submitter_ids,
                    name=item.name,
                    filters=omics.resolve_filter(),
                    notes=item.notes or "",
                    candidate_filters=item.filters,
                )
            )
    else:
        rows = [
            _summarize_candidate(
                client=client,
                projects=projects,
                baseline_case_ids=baseline_case_ids,
                baseline_submitter_ids=baseline_submitter_ids,
                candidate=candidate,
            )
            for candidate in DEFAULT_CANDIDATES
        ]
    rows.sort(key=lambda row: row["overlap_cases"], reverse=True)

    if args.json:
        print(json.dumps(rows, indent=2))
        return

    columns = [
        "name",
        "overlap_cases",
        "overlap_submitter_ids",
        "available_cases_in_projects",
        "baseline_cases",
        "overlap_pct",
        "notes",
    ]
    print("\t".join(columns))
    for row in rows:
        print("\t".join(str(row[col]) for col in columns))


if __name__ == "__main__":
    main()
