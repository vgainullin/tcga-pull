"""End-to-end tests against the live GDC.

Two tiers, gated by markers:

  @pytest.mark.network   — read-only GDC API queries. Fast (< 5 s), small payloads.
  @pytest.mark.download  — actually downloads + processes a tiny cohort.
                            ~30 s wall, ~3 MB on the wire, hits the bulk /data
                            endpoint (no gdc-client needed).

Run locally:
    uv run pytest -m network         # API queries only
    uv run pytest -m download        # downloads
    uv run pytest -m "not network and not download"  # default offline subset
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest
from rich.console import Console

from tcga_pull import load_cohort
from tcga_pull.config import load_yaml
from tcga_pull.gdc import GDCClient
from tcga_pull.pipeline import fetch_preview
from tcga_pull.pipeline import run as pipeline_run

# ---------------------------------------------------------------- @network


@pytest.mark.network
def test_list_projects_includes_tcga_brca():
    rows = GDCClient().list_projects(program="TCGA")
    project_ids = {r.get("project_id") for r in rows}
    assert "TCGA-BRCA" in project_ids
    assert "TCGA-LUAD" in project_ids
    # 33 TCGA projects expected; tolerate small drift
    assert len(project_ids) >= 30, f"expected >= 30 TCGA projects, got {len(project_ids)}"


@pytest.mark.network
def test_preview_tiny_filter_returns_nonzero_counts(tmp_path: Path):
    yaml = tmp_path / "preview.yaml"
    yaml.write_text(
        "name: preview\n"
        "filters:\n"
        "  project: TCGA-CHOL\n"
        "  data_category: Simple Nucleotide Variation\n"
        "  data_format: MAF\n"
    )
    spec = load_yaml(yaml)
    preview = fetch_preview(spec)
    assert preview.n_files > 0, "TCGA-CHOL should have SNV MAFs"
    assert preview.n_cases > 0
    assert preview.total_size > 0


# --------------------------------------------------------------- @download


@pytest.mark.download
def test_pull_with_recipes_end_to_end(tmp_path: Path):
    """Full pipeline against live GDC: pull -> variants -> samples -> frequency.

    Uses --limit-per-project=1 across 3 small projects → ~3-5 MAFs, ~30 s wall.
    """
    yaml_path = tmp_path / "ci_tiny.yaml"
    yaml_path.write_text(
        "name: ci_tiny\n"
        "filters:\n"
        "  project: [TCGA-CHOL, TCGA-DLBC, TCGA-ACC]\n"
        "  data_category: Simple Nucleotide Variation\n"
        "  data_format: MAF\n"
        "limit:\n"
        "  per_project: 1\n"
        "recipes:\n"
        "  - variants\n"
        "  - samples\n"
        "  - frequency\n"
    )
    spec = load_yaml(yaml_path, out_dir_override=tmp_path)
    # Silence the pull's Rich output during the test
    import io

    pipeline_run(spec, yes=True, console=Console(file=io.StringIO()))

    cohort = load_cohort(tmp_path / "ci_tiny")

    # All four required parquets land
    assert len(cohort.clinical) >= 3, "at least 3 cases (one per project)"
    assert len(cohort.manifest) >= 3, "at least 3 files"
    assert len(cohort.variants) > 0, "at least one somatic mutation across the cases"
    assert len(cohort.samples) >= 3

    # Recipes wrote their optional parquets
    assert cohort.gene_frequency is not None
    assert cohort.variant_frequency is not None
    assert len(cohort.gene_frequency) > 0
    assert len(cohort.variant_frequency) > 0

    # Required columns present (sanity-check schema didn't drift)
    for col in (
        "chrom",
        "pos",
        "ref",
        "alt",
        "hugo_symbol",
        "is_coding",
        "is_high_impact",
        "is_rare",
        "primary_aliquot",
    ):
        assert col in cohort.variants.columns, f"variants.parquet missing {col!r}"
    for col in (
        "case_id",
        "submitter_id",
        "project_id",
        "lineage",
        "n_variants_total",
        "n_variants_coding",
    ):
        assert col in cohort.samples.columns, f"samples.parquet missing {col!r}"

    # Lineage was derived (not just project_id)
    lineages = set(cohort.samples["lineage"].to_list())
    # TCGA-CHOL → bile_duct, TCGA-DLBC → lymph_node, TCGA-ACC → adrenal
    assert lineages & {"bile_duct", "lymph_node", "adrenal"}, (
        f"expected curated tissue labels, got {lineages}"
    )

    # primary_aliquot is exactly one per patient
    primaries_per_patient = (
        cohort.variants.filter(pl.col("primary_aliquot"))
        .group_by("submitter_id")
        .agg(pl.col("tumor_barcode").n_unique().alias("n"))
    )
    assert (primaries_per_patient["n"] == 1).all()

    # Provenance sidecar is populated
    prov = cohort.provenance
    assert "filter" in prov
    assert prov["n_files"] >= 3
