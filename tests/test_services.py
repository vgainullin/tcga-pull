"""Tests for UI-neutral application services."""

from __future__ import annotations

from pathlib import Path

import pytest

from tcga_pull.services import (
    CohortBuildOptions,
    SpecBuildError,
    build_cohort_spec,
    default_cohort_name,
    list_projects,
    resolve_samples_writer,
    resolve_variants_writer,
)


def test_build_cohort_spec_from_direct_inputs_dedupes_projects(tmp_path: Path):
    projects_file = tmp_path / "projects.txt"
    projects_file.write_text("TCGA-LUAD\nTCGA-BRCA\n")

    spec = build_cohort_spec(
        CohortBuildOptions(
            name="mixed",
            out=tmp_path,
            project=["TCGA-BRCA"],
            projects_file=projects_file,
            data_category=["Simple Nucleotide Variation"],
            data_format=["MAF"],
            n_processes=8,
            limit_per_project=3,
        )
    )

    assert spec.name == "mixed"
    assert spec.out_dir == tmp_path.resolve()
    assert spec.filters["project"] == ["TCGA-BRCA", "TCGA-LUAD"]
    assert spec.filters["data_category"] == ["Simple Nucleotide Variation"]
    assert spec.filters["data_format"] == ["MAF"]
    assert spec.n_processes == 8
    assert spec.limit.per_project == 3


def test_build_cohort_spec_loads_yaml_and_applies_overrides(tmp_path: Path):
    yaml_path = tmp_path / "cohort.yaml"
    yaml_path.write_text(
        "name: yaml_cohort\n"
        "out_dir: /tmp/ignored\n"
        "filters: {project: TCGA-BRCA}\n"
        "limit: {per_project: 10}\n"
    )

    spec = build_cohort_spec(
        CohortBuildOptions(
            config=yaml_path,
            out=tmp_path / "out",
            limit_per_project=2,
        )
    )

    assert spec.name == "yaml_cohort"
    assert spec.out_dir == (tmp_path / "out").resolve()
    assert spec.filters["project"] == "TCGA-BRCA"
    assert spec.limit.per_project == 2


def test_build_cohort_spec_requires_project_when_not_loading_yaml():
    with pytest.raises(SpecBuildError, match="Need either"):
        build_cohort_spec(CohortBuildOptions())


def test_build_cohort_spec_rejects_non_positive_limit_override(tmp_path: Path):
    with pytest.raises(SpecBuildError, match="must be > 0"):
        build_cohort_spec(
            CohortBuildOptions(
                project=["TCGA-BRCA"],
                out=tmp_path,
                limit_per_project=0,
            )
        )


def test_default_cohort_name_matches_existing_cli_behavior():
    assert default_cohort_name(["TCGA-BRCA"], None) == "brca"
    assert default_cohort_name(["TCGA-BRCA"], ["Gene Expression Quantification"]) == (
        "brca_gene_expression_quantification"
    )
    assert default_cohort_name(["TCGA-BRCA", "TCGA-LUAD"], None) == "multi_2_projects"


def test_list_projects_sorts_by_project_id():
    class FakeClient:
        def list_projects(self, program: str = "TCGA") -> list[dict]:
            assert program == "TCGA"
            return [
                {"project_id": "TCGA-LUAD"},
                {"project_id": "TCGA-BRCA"},
                {"project_id": None},
            ]

    assert [row.get("project_id") for row in list_projects(client=FakeClient())] == [
        None,
        "TCGA-BRCA",
        "TCGA-LUAD",
    ]


def test_engine_resolvers_reject_unknown_engine():
    with pytest.raises(ValueError, match="unknown engine"):
        resolve_variants_writer("nope")
    with pytest.raises(ValueError, match="unknown engine"):
        resolve_samples_writer("nope")
