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
from typer.testing import CliRunner

from tcga_pull import load_cohort
from tcga_pull.cli import app
from tcga_pull.config import load_yaml
from tcga_pull.coverage import build_coverage_matrix
from tcga_pull.gdc import GDCClient
from tcga_pull.pipeline import fetch_preview
from tcga_pull.pipeline import run as pipeline_run

runner = CliRunner()

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


@pytest.mark.network
def test_coverage_matrix_live_gdc_chol_classifies_supported_and_raw_only():
    matrix = build_coverage_matrix(projects=["TCGA-CHOL"])

    assert matrix.project_ids == ("TCGA-CHOL",)
    assert len(matrix.rows) >= 10
    assert {row["access"] for row in matrix.rows} == {"open"}

    snv_maf = [
        row
        for row in matrix.rows
        if row["data_category"] == "Simple Nucleotide Variation" and row["data_format"] == "MAF"
    ]
    assert snv_maf, "TCGA-CHOL should expose open-access SNV MAF files"
    assert {row["support_status"] for row in snv_maf} == {"supported"}
    assert {row["recipe"] for row in snv_maf} == {"variants"}

    raw_only = [row for row in matrix.rows if row["support_status"] == "raw_only"]
    assert raw_only, "coverage matrix should expose gaps as raw_only rows"
    assert any(row["data_category"] in {"Biospecimen", "Clinical"} for row in raw_only)


@pytest.mark.network
def test_coverage_command_writes_live_gdc_outputs(tmp_path: Path):
    result = runner.invoke(
        app,
        ["coverage", "--program", "TCGA", "--project", "TCGA-CHOL", "--out-dir", str(tmp_path)],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    assert "projects=1" in result.output
    assert "supported=" in result.output
    assert "raw_only=" in result.output

    parquet_path = tmp_path / "tcga_open_access_coverage_matrix.parquet"
    markdown_path = tmp_path / "tcga_open_access_coverage_matrix.md"
    assert parquet_path.exists()
    assert markdown_path.exists()

    df = pl.read_parquet(parquet_path)
    assert df.height >= 10
    assert set(df["project_id"].to_list()) == {"TCGA-CHOL"}
    assert "supported" in set(df["support_status"].to_list())
    assert "raw_only" in set(df["support_status"].to_list())
    assert "TCGA-CHOL" in markdown_path.read_text()


# --------------------------------------------------------------- @download


@pytest.mark.download
def test_pull_with_recipes_end_to_end(tmp_path: Path):
    """CLI pull against live GDC: pull -> variants -> samples -> frequency.

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
        "processing:\n"
        "  mode: incremental\n"
        "  batch_size: 1\n"
        "  delete_raw_after_processing: true\n"
        "recipes:\n"
        "  - variants\n"
        "  - samples\n"
        "  - frequency\n"
    )

    result = runner.invoke(
        app,
        ["pull", str(yaml_path), "--out", str(tmp_path)],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output

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
    expected_lineages = {"bile_duct", "lymph_node", "adrenal"}
    assert expected_lineages <= lineages, f"expected curated tissue labels, got {lineages}"

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
    assert prov["processing"]["mode"] == "incremental"
    assert prov["processing"]["delete_raw_after_processing"] is True


@pytest.mark.download
def test_pull_multiomics_smoke_end_to_end(tmp_path: Path):
    """Tiny live-GDC smoke test for incremental multiomics recipes.

    The config is generated in tmp_path, not checked into examples. CHOL has
    small SNV/RNA/CNV/methylation coverage, and per_project=1 keeps the
    download bounded for CI.
    """
    yaml_path = tmp_path / "ci_multiomics_smoke.yaml"
    yaml_path.write_text(
        "name: ci_multiomics_smoke\n"
        "filters:\n"
        "  project: TCGA-CHOL\n"
        "limit:\n"
        "  per_project: 1\n"
        "processing:\n"
        "  mode: incremental\n"
        "  batch_size: 1\n"
        "  delete_raw_after_processing: true\n"
        "recipe_options:\n"
        "  rna_expression:\n"
        "    columns: [gene_id, gene_name, unstranded]\n"
        "  copy_number:\n"
        "    outputs: [segments]\n"
        "  model_dataset:\n"
        "    label_column: lineage\n"
        "    modalities: [snv, rna_expression, methylation_beta, gene_copy_number]\n"
        "    min_class_count: 1\n"
        "recipes:\n"
        "  - variants\n"
        "  - samples\n"
        "  - multiomics\n"
        "gdc_filter:\n"
        "  op: and\n"
        "  content:\n"
        "    - op: in\n"
        "      content:\n"
        "        field: cases.project.project_id\n"
        "        value: [TCGA-CHOL]\n"
        "    - op: or\n"
        "      content:\n"
        "        - op: and\n"
        "          content:\n"
        "            - op: in\n"
        "              content: {field: data_category, value: [Simple Nucleotide Variation]}\n"
        "            - op: in\n"
        "              content: {field: data_format, value: [MAF]}\n"
        "        - op: and\n"
        "          content:\n"
        "            - op: in\n"
        "              content: {field: data_category, value: [Transcriptome Profiling]}\n"
        "            - op: in\n"
        "              content: {field: data_type, value: [Gene Expression Quantification]}\n"
        "            - op: in\n"
        "              content: {field: experimental_strategy, value: [RNA-Seq]}\n"
        "            - op: in\n"
        "              content: {field: analysis.workflow_type, value: [STAR - Counts]}\n"
        "        - op: and\n"
        "          content:\n"
        "            - op: in\n"
        "              content: {field: data_category, value: [Copy Number Variation]}\n"
        "            - op: in\n"
        "              content: {field: data_type, value: [Copy Number Segment]}\n"
        "        - op: and\n"
        "          content:\n"
        "            - op: in\n"
        "              content: {field: data_category, value: [DNA Methylation]}\n"
        "            - op: in\n"
        "              content: {field: data_type, value: [Methylation Beta Value]}\n"
    )
    spec = load_yaml(yaml_path, out_dir_override=tmp_path)

    import io

    pipeline_run(spec, console=Console(file=io.StringIO()))
    cohort = load_cohort(tmp_path / "ci_multiomics_smoke")

    assert len(cohort.variants) > 0
    assert len(cohort.samples) >= 1
    assert cohort.rna_expression is not None
    assert cohort.copy_number_segments is not None
    assert cohort.methylation_beta is not None
    assert cohort.gene_copy_number is None

    assert len(cohort.rna_expression) > 0
    assert len(cohort.copy_number_segments) > 0
    assert len(cohort.methylation_beta) > 0

    assert cohort.rna_expression.columns == [
        "case_id",
        "submitter_id",
        "file_id",
        "file_name",
        "data_type",
        "experimental_strategy",
        "workflow_type",
        "gene_id",
        "gene_name",
        "unstranded",
    ]
    for col in ("case_id", "submitter_id", "chrom", "segment_mean"):
        assert col in cohort.copy_number_segments.columns
    for col in ("case_id", "submitter_id", "probe_id", "beta_value"):
        assert col in cohort.methylation_beta.columns

    manifest = cohort.manifest
    handled_omics = manifest.filter(
        pl.col("data_category").is_in(
            ["Transcriptome Profiling", "Copy Number Variation", "DNA Methylation"]
        )
    )
    assert handled_omics.height > 0
    assert set(handled_omics["status"].to_list()) == {"processed_deleted"}
    assert handled_omics["local_path"].null_count() == handled_omics.height

    prov = cohort.provenance
    assert prov["processing"]["mode"] == "incremental"
    assert prov["processing"]["delete_raw_after_processing"] is True

    result = runner.invoke(
        app,
        ["dataset", str(tmp_path / "ci_multiomics_smoke"), "--config", str(yaml_path)],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output

    model_dataset = load_cohort(tmp_path / "ci_multiomics_smoke").model_dataset
    assert model_dataset is not None
    assert (model_dataset.path / "manifest.json").exists()
    assert (model_dataset.path / "samples.parquet").exists()
    assert (model_dataset.path / "feature_index.parquet").exists()

    assert len(model_dataset.samples) == len(cohort.samples)
    assert model_dataset.manifest["n_samples"] == len(cohort.samples)
    assert model_dataset.manifest["label_column"] == "lineage"

    assert model_dataset.snv is not None
    assert model_dataset.rna_expression is not None
    assert model_dataset.gene_copy_number is None
    assert len(model_dataset.feature_index) > 0
    assert set(model_dataset.feature_index["modality"].to_list()) >= {
        "snv",
        "rna_expression",
    }
    assert "gene_copy_number" in model_dataset.manifest["skipped_modalities"]
    if "methylation_beta" in model_dataset.manifest["modalities"]:
        assert model_dataset.methylation_beta is not None
    else:
        assert "methylation_beta" in model_dataset.manifest["skipped_modalities"]
