"""Offline tests for model dataset matrix export."""

from __future__ import annotations

import math
from pathlib import Path

import polars as pl
import pytest

from tcga_pull.model_dataset import write_model_dataset


def _fixture_cohort(tmp_path: Path) -> Path:
    cohort = tmp_path / "cohort"
    cohort.mkdir()

    samples = pl.DataFrame(
        [
            {
                "case_id": "c1",
                "submitter_id": "P1",
                "lineage": "breast",
                "oncotree_tissue": "Breast",
            },
            {
                "case_id": "c2",
                "submitter_id": "P2",
                "lineage": "breast",
                "oncotree_tissue": "Breast",
            },
            {
                "case_id": "c3",
                "submitter_id": "P3",
                "lineage": "breast",
                "oncotree_tissue": "Breast",
            },
            {
                "case_id": "c4",
                "submitter_id": "P4",
                "lineage": "lung",
                "oncotree_tissue": "Lung",
            },
            {
                "case_id": "c5",
                "submitter_id": "P5",
                "lineage": "lung",
                "oncotree_tissue": "Lung",
            },
            {
                "case_id": "c6",
                "submitter_id": "P6",
                "lineage": "lung",
                "oncotree_tissue": "Lung",
            },
        ]
    )
    samples.write_parquet(cohort / "samples.parquet")
    pl.DataFrame(
        {"case_id": samples["case_id"], "submitter_id": samples["submitter_id"]}
    ).write_parquet(cohort / "clinical.parquet")

    manifest = pl.DataFrame(
        [
            {
                "file_id": f"maf-{submitter}",
                "submitter_id": submitter,
                "data_category": "Simple Nucleotide Variation",
                "data_format": "MAF",
            }
            for submitter in ("P1", "P2", "P3", "P4", "P5", "P6")
        ]
    )
    manifest.write_parquet(cohort / "manifest.parquet")

    pl.DataFrame(
        [
            {
                "submitter_id": "P1",
                "hugo_symbol": "TP53",
                "primary_aliquot": True,
                "is_coding": True,
                "is_rare": True,
            },
            {
                "submitter_id": "P2",
                "hugo_symbol": "PIK3CA",
                "primary_aliquot": True,
                "is_coding": True,
                "is_rare": True,
            },
            {
                "submitter_id": "P4",
                "hugo_symbol": "TP53",
                "primary_aliquot": False,
                "is_coding": True,
                "is_rare": True,
            },
            {
                "submitter_id": "P5",
                "hugo_symbol": "KRAS",
                "primary_aliquot": True,
                "is_coding": True,
                "is_rare": True,
            },
        ]
    ).write_parquet(cohort / "variants.parquet")

    pl.DataFrame(
        [
            {"submitter_id": "P1", "gene_id": "ENSG1", "unstranded": 9},
            {"submitter_id": "P1", "gene_id": "ENSG2", "unstranded": 3},
            {"submitter_id": "P4", "gene_id": "ENSG1", "unstranded": 0},
        ]
    ).write_parquet(cohort / "rna_expression.parquet")
    pl.DataFrame(
        [
            {"submitter_id": "P1", "probe_id": "cg0001", "beta_value": 0.25},
            {"submitter_id": "P4", "probe_id": "cg0001", "beta_value": 0.75},
        ]
    ).write_parquet(cohort / "methylation_beta.parquet")
    pl.DataFrame(
        [
            {"submitter_id": "P1", "gene_name": "TP53", "copy_number": 1.0},
            {"submitter_id": "P4", "gene_name": "TP53", "copy_number": 2.5},
        ]
    ).write_parquet(cohort / "gene_copy_number.parquet")
    pl.DataFrame(
        [
            {
                "submitter_id": "P1",
                "mirna_id": "hsa-let-7a",
                "reads_per_million_mirna_mapped": 10.0,
            },
            {
                "submitter_id": "P4",
                "mirna_id": "hsa-let-7a",
                "reads_per_million_mirna_mapped": 20.0,
            },
        ]
    ).write_parquet(cohort / "mirna_expression.parquet")
    pl.DataFrame(
        [
            {"submitter_id": "P1", "protein_id": "TP53|p53", "expression_value": -0.2},
            {"submitter_id": "P4", "protein_id": "TP53|p53", "expression_value": 0.4},
        ]
    ).write_parquet(cohort / "protein_expression.parquet")
    return cohort


def _row(df: pl.DataFrame, submitter_id: str) -> dict:
    rows = df.filter(pl.col("submitter_id") == submitter_id).to_dicts()
    assert len(rows) == 1
    return rows[0]


def test_write_model_dataset_exports_case_aligned_matrices(tmp_path: Path):
    cohort = _fixture_cohort(tmp_path)

    outputs = write_model_dataset(cohort, {"model_dataset": {"min_class_count": 2, "seed": 7}})

    assert outputs.path == cohort / "model_dataset"
    assert outputs.manifest.exists()
    assert set(outputs.matrices) == {
        "gene_copy_number",
        "methylation_beta",
        "mirna_expression",
        "protein_expression",
        "rna_expression",
        "snv",
    }

    sample_index = pl.read_parquet(outputs.samples)
    assert sample_index.height == 6
    assert set(sample_index["split"].to_list()) == {"train", "val", "test"}
    for col in (
        "has_snv",
        "has_rna_expression",
        "has_methylation_beta",
        "has_gene_copy_number",
        "has_mirna_expression",
        "has_protein_expression",
    ):
        assert col in sample_index.columns

    feature_index = pl.read_parquet(outputs.feature_index)
    assert set(feature_index["modality"].to_list()) == set(outputs.matrices)

    snv = pl.read_parquet(outputs.matrices["snv"])
    assert _row(snv, "P1")["snv__tp53"] == 1
    assert _row(snv, "P4")["snv__tp53"] == 0
    assert "snv__kras" in snv.columns

    rna = pl.read_parquet(outputs.matrices["rna_expression"])
    assert math.isclose(_row(rna, "P1")["rna_expression__ensg1"], math.log1p(9))
    assert math.isclose(_row(rna, "P4")["rna_expression__ensg1"], 0.0)
    assert _row(rna, "P2")["rna_expression__ensg1"] is None

    methylation = pl.read_parquet(outputs.matrices["methylation_beta"])
    assert math.isclose(_row(methylation, "P4")["methylation_beta__cg0001"], 0.75)

    cnv = pl.read_parquet(outputs.matrices["gene_copy_number"])
    assert math.isclose(_row(cnv, "P4")["gene_copy_number__tp53"], 2.5)

    mirna = pl.read_parquet(outputs.matrices["mirna_expression"])
    assert math.isclose(_row(mirna, "P1")["mirna_expression__hsa_let_7a"], math.log1p(10.0))

    protein = pl.read_parquet(outputs.matrices["protein_expression"])
    assert math.isclose(_row(protein, "P1")["protein_expression__tp53_p53"], -0.2)


def test_model_dataset_honors_modality_and_feature_limits(tmp_path: Path):
    cohort = _fixture_cohort(tmp_path)

    outputs = write_model_dataset(
        cohort,
        {
            "model_dataset": {
                "label_column": "lineage",
                "modalities": ["rna_expression"],
                "feature_min_samples": 2,
                "max_features_per_modality": 1,
            }
        },
    )

    assert set(outputs.matrices) == {"rna_expression"}
    rna = pl.read_parquet(outputs.matrices["rna_expression"])
    assert rna.columns == ["case_id", "submitter_id", "rna_expression__ensg1"]


def test_model_dataset_rejects_missing_label_column(tmp_path: Path):
    cohort = _fixture_cohort(tmp_path)

    with pytest.raises(ValueError, match="missing required"):
        write_model_dataset(cohort, label_column="not_a_label")


def test_model_dataset_rejects_unknown_modality(tmp_path: Path):
    cohort = _fixture_cohort(tmp_path)

    with pytest.raises(ValueError, match="unknown model dataset modalities"):
        write_model_dataset(cohort, modalities=["rna_expression", "nope"])
