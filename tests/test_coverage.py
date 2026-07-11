"""Offline tests for the open-access coverage matrix."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import polars as pl

from tcga_pull.coverage import (
    COVERAGE_FILE_FIELDS,
    build_coverage_matrix,
    classify_file_group,
    render_coverage_markdown,
    write_coverage_outputs,
)


def _hit(
    file_id: str,
    *,
    project_id: str,
    case_id: str,
    data_category: str,
    data_type: str,
    data_format: str,
    experimental_strategy: str | None = None,
    workflow_type: str | None = None,
    file_size: int = 10,
) -> dict[str, Any]:
    hit: dict[str, Any] = {
        "file_id": file_id,
        "file_name": f"{file_id}.tsv",
        "data_category": data_category,
        "data_type": data_type,
        "data_format": data_format,
        "experimental_strategy": experimental_strategy,
        "access": "open",
        "file_size": file_size,
        "cases": [
            {
                "case_id": case_id,
                "submitter_id": f"S-{case_id}",
                "project": {"project_id": project_id},
            }
        ],
    }
    if workflow_type is not None:
        hit["analysis"] = {"workflow_type": workflow_type}
    return hit


def test_classify_file_group_matches_current_recipe_surface():
    assert (
        classify_file_group(
            data_category="Simple Nucleotide Variation",
            data_type="Masked Somatic Mutation",
            data_format="MAF",
            experimental_strategy="WXS",
            workflow_type=None,
        ).recipe
        == "variants"
    )
    assert classify_file_group(
        data_category="Transcriptome Profiling",
        data_type="Gene Expression Quantification",
        data_format="TSV",
        experimental_strategy="RNA-Seq",
        workflow_type="STAR - Counts",
    ).outputs == ("rna_expression.parquet",)
    assert (
        classify_file_group(
            data_category="Biospecimen",
            data_type="Biospecimen Supplement",
            data_format="BCR XML",
            experimental_strategy=None,
            workflow_type=None,
        ).status
        == "raw_only"
    )


def test_build_coverage_matrix_aggregates_file_hits_and_classifies():
    class FakeClient:
        def __init__(self) -> None:
            self.filters: dict | None = None
            self.fields: list[str] | None = None

        def list_projects(self, program: str = "TCGA") -> list[dict]:
            assert program == "TCGA"
            return [{"project_id": "TCGA-LUAD"}, {"project_id": "TCGA-BRCA"}]

        def fetch_files(self, filters: dict, fields: list[str] | None = None) -> list[dict]:
            self.filters = filters
            self.fields = fields
            return [
                _hit(
                    "maf-1",
                    project_id="TCGA-BRCA",
                    case_id="case-1",
                    data_category="Simple Nucleotide Variation",
                    data_type="Masked Somatic Mutation",
                    data_format="MAF",
                    experimental_strategy="WXS",
                    file_size=10,
                ),
                _hit(
                    "maf-2",
                    project_id="TCGA-BRCA",
                    case_id="case-2",
                    data_category="Simple Nucleotide Variation",
                    data_type="Masked Somatic Mutation",
                    data_format="MAF",
                    experimental_strategy="WXS",
                    file_size=20,
                ),
                _hit(
                    "rna-1",
                    project_id="TCGA-LUAD",
                    case_id="case-3",
                    data_category="Transcriptome Profiling",
                    data_type="Gene Expression Quantification",
                    data_format="TSV",
                    experimental_strategy="RNA-Seq",
                    workflow_type="STAR - Counts",
                    file_size=30,
                ),
                _hit(
                    "clinical-1",
                    project_id="TCGA-LUAD",
                    case_id="case-3",
                    data_category="Clinical",
                    data_type="Clinical Supplement",
                    data_format="BCR XML",
                    file_size=40,
                ),
            ]

    client = FakeClient()
    matrix = build_coverage_matrix(client=client)

    assert client.fields == COVERAGE_FILE_FIELDS
    assert matrix.project_ids == ("TCGA-BRCA", "TCGA-LUAD")
    assert [row["project_id"] for row in matrix.rows] == [
        "TCGA-BRCA",
        "TCGA-LUAD",
        "TCGA-LUAD",
    ]

    brca = matrix.rows[0]
    assert brca["data_category"] == "Simple Nucleotide Variation"
    assert brca["n_files"] == 2
    assert brca["n_cases"] == 2
    assert brca["total_size"] == 30
    assert brca["support_status"] == "supported"
    assert brca["recipe"] == "variants"

    clinical = matrix.rows[1]
    assert clinical["data_category"] == "Clinical"
    assert clinical["support_status"] == "raw_only"

    rna = matrix.rows[2]
    assert rna["recipe"] == "rna_expression"


def test_build_coverage_matrix_project_override_skips_project_listing():
    class FakeClient:
        def list_projects(self, program: str = "TCGA") -> list[dict]:
            raise AssertionError("project override should not list projects")

        def fetch_files(self, filters: dict, fields: list[str] | None = None) -> list[dict]:
            return [
                _hit(
                    "cnv-1",
                    project_id="TCGA-LUAD",
                    case_id="case-1",
                    data_category="Copy Number Variation",
                    data_type="Copy Number Segment",
                    data_format="TSV",
                    file_size=50,
                )
            ]

    matrix = build_coverage_matrix(projects=["TCGA-LUAD", "TCGA-LUAD"], client=FakeClient())

    assert matrix.project_ids == ("TCGA-LUAD",)
    assert matrix.rows[0]["recipe"] == "copy_number"
    assert matrix.rows[0]["outputs"] == "copy_number_segments.parquet"


def test_write_coverage_outputs_writes_parquet_and_markdown(tmp_path: Path):
    matrix = build_coverage_matrix(
        projects=["TCGA-BRCA"],
        client=type(
            "FakeClient",
            (),
            {
                "fetch_files": lambda self, filters, fields=None: [
                    _hit(
                        "protein-1",
                        project_id="TCGA-BRCA",
                        case_id="case-1",
                        data_category="Proteome Profiling",
                        data_type="Protein Expression Quantification",
                        data_format="TSV",
                        experimental_strategy="Reverse Phase Protein Array",
                    )
                ]
            },
        )(),
    )

    outputs = write_coverage_outputs(matrix, tmp_path)

    assert outputs.n_projects == 1
    assert outputs.n_rows == 1
    assert outputs.parquet_path.name == "tcga_open_access_coverage_matrix.parquet"
    assert outputs.markdown_path.name == "tcga_open_access_coverage_matrix.md"
    df = pl.read_parquet(outputs.parquet_path)
    assert df.select("recipe").item() == "protein_expression"
    assert "TCGA Open-Access Coverage Matrix" in outputs.markdown_path.read_text()
    assert "protein_expression" in render_coverage_markdown(matrix)
