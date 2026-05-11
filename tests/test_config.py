"""Tests for config: YAML loading, multi-project filters, projects-file parsing."""

from __future__ import annotations

from pathlib import Path

import pytest

from tcga_pull.config import CohortSpec, from_flags, load_yaml, read_projects_file


def _find_leaf(flt: dict, field: str) -> dict | None:
    """Walk an AND/OR tree and return the leaf clause for the given field (or None)."""
    if not isinstance(flt, dict):
        return None
    op = flt.get("op")
    content = flt.get("content")
    if op in ("and", "or") and isinstance(content, list):
        for child in content:
            hit = _find_leaf(child, field)
            if hit is not None:
                return hit
        return None
    if op in ("in", "=") and isinstance(content, dict) and content.get("field") == field:
        return flt
    return None


def test_read_projects_file_ignores_comments_and_blanks(tmp_path: Path):
    p = tmp_path / "projects.txt"
    p.write_text(
        "TCGA-BRCA\n"
        "# a comment\n"
        "\n"
        "TCGA-LUAD  \n"
        "  TARGET-AML\n"  # leading whitespace stripped
        "# trailing comment\n"
    )
    assert read_projects_file(p) == ["TCGA-BRCA", "TCGA-LUAD", "TARGET-AML"]


def test_from_flags_multi_project():
    spec = from_flags(
        name="multi",
        out_dir=Path("/tmp"),
        project=["TCGA-BRCA", "TCGA-LUAD", "TARGET-AML"],
        data_category=["Simple Nucleotide Variation"],
        data_format=["MAF"],
    )
    flt = spec.resolve_filter()
    project_clause = _find_leaf(flt, "cases.project.project_id")
    assert project_clause is not None
    assert project_clause["content"]["value"] == ["TCGA-BRCA", "TCGA-LUAD", "TARGET-AML"]
    # File-rooted fields live as bare names on /files (no `files.` prefix)
    format_clause = _find_leaf(flt, "data_format")
    assert format_clause is not None
    assert format_clause["content"]["value"] == ["MAF"]


def test_yaml_multi_project(tmp_path: Path):
    yaml_path = tmp_path / "cohort.yaml"
    yaml_path.write_text(
        "name: pancancer_snv\n"
        "out_dir: /tmp/out\n"
        "filters:\n"
        "  project:\n"
        "    - TCGA-BRCA\n"
        "    - TCGA-LUAD\n"
        "  data_category: Simple Nucleotide Variation\n"
        "  data_format: MAF\n"
    )
    spec = load_yaml(yaml_path)
    assert spec.name == "pancancer_snv"
    assert spec.filters["project"] == ["TCGA-BRCA", "TCGA-LUAD"]
    # Single-string sugar still works
    assert spec.filters["data_category"] == "Simple Nucleotide Variation"

    flt = spec.resolve_filter()
    project_clause = _find_leaf(flt, "cases.project.project_id")
    assert project_clause is not None
    assert project_clause["content"]["value"] == ["TCGA-BRCA", "TCGA-LUAD"]


def test_yaml_unknown_filter_raises(tmp_path: Path):
    yaml_path = tmp_path / "cohort.yaml"
    yaml_path.write_text("name: bad\nfilters:\n  not_a_real_field: yes\n")
    spec = load_yaml(yaml_path)
    with pytest.raises(ValueError, match="unknown filter key"):
        spec.resolve_filter()


def test_cohort_spec_resolve_filter_wraps_open_access():
    spec = CohortSpec(name="x", out_dir=Path("/tmp"), filters={"project": "TCGA-CHOL"})
    flt = spec.resolve_filter()
    # Bare `access` — file-rooted fields are unprefixed on /files
    access_clause = _find_leaf(flt, "access")
    assert access_clause is not None
    assert access_clause["content"]["value"] == ["open"]


def test_yaml_with_recipes(tmp_path: Path):
    yaml_path = tmp_path / "cohort.yaml"
    yaml_path.write_text(
        "name: x\n"
        "filters: {project: TCGA-CHOL}\n"
        "recipes:\n"
        "  - variants\n"
        "  - samples\n"
        "  - frequency\n"
    )
    spec = load_yaml(yaml_path)
    assert spec.recipes == ["variants", "samples", "frequency"]


def test_cohort_spec_rejects_unknown_recipe():
    with pytest.raises(ValueError, match="unknown recipe"):
        CohortSpec(name="x", out_dir=Path("/tmp"), recipes=["variants", "no_such_recipe"])


def test_cohort_spec_recipes_default_empty():
    spec = CohortSpec(name="x", out_dir=Path("/tmp"))
    assert spec.recipes == []
