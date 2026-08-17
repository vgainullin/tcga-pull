"""Offline tests for entity-aware multiomics integration."""

from __future__ import annotations

import gzip
import hashlib
import json
import math
from pathlib import Path

import polars as pl
import pytest

from tcga_pull import load_cohort
from tcga_pull.annotations import (
    METHYLATION_ANNOTATIONS,
    MethylationAnnotation,
    build_methylation_mapping,
    resolve_methylation_annotations,
)
from tcga_pull.integration import write_integrated_dataset
from tcga_pull.services import write_integration_recipe


def _cohort(tmp_path: Path) -> Path:
    cohort = tmp_path / "cohort"
    cohort.mkdir()
    pl.DataFrame(
        [
            {"case_id": "c1", "submitter_id": "P1"},
            {"case_id": "c2", "submitter_id": "P2"},
        ]
    ).write_parquet(cohort / "samples.parquet")
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
                "hugo_symbol": "KRAS",
                "primary_aliquot": True,
                "is_coding": True,
                "is_rare": True,
            },
            {
                "submitter_id": "P2",
                "hugo_symbol": "TP53",
                "primary_aliquot": False,
                "is_coding": True,
                "is_rare": True,
            },
        ]
    ).write_parquet(cohort / "variants.parquet")
    pl.DataFrame(
        [
            {"submitter_id": "P1", "gene_name": "TP53", "unstranded": 10},
            {"submitter_id": "P1", "gene_name": "PIK3CA", "unstranded": 5},
            {"submitter_id": "P2", "gene_name": "TP53", "unstranded": 20},
        ]
    ).write_parquet(cohort / "rna_expression.parquet")
    pl.DataFrame(
        [
            {
                "submitter_id": "P1",
                "probe_id": "cg1",
                "beta_value": 0.2,
                "platform": "Illumina Human Methylation 450",
            },
            {
                "submitter_id": "P1",
                "probe_id": "cg2",
                "beta_value": 0.4,
                "platform": "Illumina Human Methylation 450",
            },
            {
                "submitter_id": "P2",
                "probe_id": "cg1",
                "beta_value": 0.8,
                "platform": "Illumina Human Methylation 450",
            },
        ]
    ).write_parquet(cohort / "methylation_beta.parquet")
    pl.DataFrame([{"submitter_id": "P1", "gene_name": "TP53", "copy_number": 2.0}]).write_parquet(
        cohort / "gene_copy_number.parquet"
    )
    pl.DataFrame(
        [
            {
                "submitter_id": "P1",
                "mirna_id": "hsa-let-7a",
                "reads_per_million_mirna_mapped": 7.0,
            }
        ]
    ).write_parquet(cohort / "mirna_expression.parquet")
    pl.DataFrame(
        [
            {
                "submitter_id": "P1",
                "protein_id": "TP53|p53",
                "gene_symbol": "TP53",
                "expression_value": -0.2,
            }
        ]
    ).write_parquet(cohort / "protein_expression.parquet")
    return cohort


def _integrated_row(frame: pl.DataFrame, submitter_id: str, entity_id: str) -> dict:
    rows = frame.filter(
        pl.col("submitter_id") == submitter_id,
        pl.col("entity_id") == entity_id,
    ).to_dicts()
    assert len(rows) == 1
    return rows[0]


def _fake_annotations(*args, **kwargs):
    mapping = pl.DataFrame(
        [
            {
                "source_feature_id": "cg1",
                "source_platform": "Illumina Human Methylation 450",
                "source_platform_key": "illumina human methylation 450",
                "entity_id": "TP53",
                "gene_region": "promoter",
            },
            {
                "source_feature_id": "cg2",
                "source_platform": "Illumina Human Methylation 450",
                "source_platform_key": "illumina human methylation 450",
                "entity_id": "PIK3CA",
                "gene_region": "promoter",
            },
        ]
    )
    metadata = [
        {
            "platform": "Illumina Human Methylation 450",
            "array": "HM450",
            "release": "InfiniumAnnotation-v8.1-gencode-v41-hg38",
            "genome_build": "hg38",
            "url": "https://example.test/hm450.tsv.gz",
            "sha256": "abc123",
            "size_bytes": 123,
        }
    ]
    return mapping, metadata


def test_default_integration_merges_direct_and_annotated_gene_entities(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    cohort = _cohort(tmp_path)
    monkeypatch.setattr("tcga_pull.integration.resolve_methylation_annotations", _fake_annotations)

    outputs = write_integrated_dataset(cohort)

    integrated = pl.read_parquet(outputs.integrated)
    tp53 = _integrated_row(integrated, "P1", "TP53")
    assert tp53["entity_type"] == "gene_symbol"
    assert tp53["snv"] == 1.0
    assert tp53["rna_expression"] == 10.0
    assert tp53["gene_copy_number"] == 2.0
    assert tp53["protein_expression"] == -0.2
    assert tp53["methylation_beta"] == 0.2

    mirna = _integrated_row(integrated, "P1", "hsa-let-7a")
    assert mirna["entity_type"] == "mirna"

    values = pl.read_parquet(outputs.values)
    assert not values.filter(
        pl.col("submitter_id") == "P2",
        pl.col("entity_id") == "TP53",
        pl.col("modality") == "snv",
    ).height

    api_dataset = load_cohort(cohort).integrated
    assert api_dataset is not None
    assert api_dataset.manifest["n_values"] == values.height
    assert api_dataset.integrated.height == integrated.height
    mappings = pl.read_parquet(outputs.mappings)
    methylation = mappings.filter(pl.col("modality") == "methylation_beta")
    assert methylation["mapping_type"].unique().to_list() == ["builtin_annotation"]
    assert methylation["annotation_sha256"].unique().to_list() == ["abc123"]
    provenance = json.loads((cohort / "cohort.json").read_text())
    annotation_input = provenance["recipes"]["inputs"]["integrate"]["modalities"][
        "methylation_beta"
    ]["annotations"][0]
    assert annotation_input["sha256"] == "abc123"


def test_built_in_methylation_registry_covers_gdc_arrays():
    assert set(METHYLATION_ANNOTATIONS) == {
        "illumina human methylation 27",
        "illumina human methylation 450",
        "illumina methylation epic",
        "illumina methylation epic v2",
    }


def test_infinium_annotation_maps_promoter_probes_to_gene_symbols(tmp_path: Path):
    annotation = tmp_path / "annotation.tsv.gz"
    content = (
        "probeID\tgeneNames\tdistToTSS\ncg1\tTP53;TP53;WRONG\t100;-1400;2000\ncg2\tPIK3CA\t-200\n"
    )
    with gzip.open(annotation, "wt") as handle:
        handle.write(content)

    mapping = build_methylation_mapping(
        annotation,
        "Illumina Human Methylation 450",
    )

    assert mapping.select("source_feature_id", "entity_id", "gene_region").to_dicts() == [
        {"source_feature_id": "cg1", "entity_id": "TP53", "gene_region": "promoter"},
        {"source_feature_id": "cg2", "entity_id": "PIK3CA", "gene_region": "promoter"},
    ]
    all_regions = build_methylation_mapping(
        annotation,
        "Illumina Human Methylation 450",
        regions=("promoter", "body"),
    )
    assert (
        all_regions.filter(
            pl.col("source_feature_id") == "cg1", pl.col("entity_id") == "TP53"
        ).height
        == 1
    )


def test_annotation_resolver_downloads_verifies_and_reuses_cohort_copy(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    payload = gzip.compress(b"probeID\tgeneNames\tdistToTSS\ncg1\tTP53\t100\n")
    spec = MethylationAnnotation(
        platform="Test Array",
        array="TEST",
        url="https://example.test/annotation.tsv.gz",
        sha256=hashlib.sha256(payload).hexdigest(),
        size_bytes=len(payload),
    )
    monkeypatch.setitem(METHYLATION_ANNOTATIONS, "test array", spec)
    calls = 0

    class Response:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return None

        def raise_for_status(self):
            return None

        def iter_content(self, chunk_size: int):
            yield payload

    def fake_get(*args, **kwargs):
        nonlocal calls
        calls += 1
        return Response()

    monkeypatch.setattr("tcga_pull.annotations.requests.get", fake_get)

    mapping, metadata = resolve_methylation_annotations(tmp_path, ["Test Array"])
    assert mapping.select("source_feature_id", "entity_id").to_dicts() == [
        {"source_feature_id": "cg1", "entity_id": "TP53"}
    ]
    assert metadata[0]["source_sha256"] == spec.sha256
    assert (tmp_path / "annotations" / "test.gencode-v41.hg38.tsv.gz").exists()
    assert (tmp_path / "annotations" / "test.probe-gene.parquet").exists()

    resolve_methylation_annotations(tmp_path, ["Test Array"])
    assert calls == 1


def test_auto_methylation_rejects_unknown_platform(tmp_path: Path):
    cohort = _cohort(tmp_path)
    methylation = pl.read_parquet(cohort / "methylation_beta.parquet").with_columns(
        pl.lit("Unknown Array").alias("platform")
    )
    methylation.write_parquet(cohort / "methylation_beta.parquet")

    with pytest.raises(RuntimeError, match="unsupported GDC methylation platform"):
        write_integrated_dataset(
            cohort,
            {"integrate": {"modalities": ["methylation_beta"]}},
        )


def test_auto_methylation_rejects_multiple_platforms_per_join(tmp_path: Path):
    cohort = _cohort(tmp_path)
    methylation = pl.read_parquet(cohort / "methylation_beta.parquet")
    methylation = methylation.with_columns(
        pl.when(pl.col("probe_id") == "cg2")
        .then(pl.lit("Illumina Human Methylation 27"))
        .otherwise(pl.col("platform"))
        .alias("platform")
    )
    methylation.write_parquet(cohort / "methylation_beta.parquet")

    with pytest.raises(RuntimeError, match="multiple methylation platforms"):
        write_integrated_dataset(
            cohort,
            {"integrate": {"modalities": ["methylation_beta"]}},
        )


def test_external_mapping_merges_methylation_into_gene_namespace(tmp_path: Path):
    cohort = _cohort(tmp_path)
    mapping = tmp_path / "methylation.tsv"
    mapping.write_text(
        "IlmnID\tUCSC_RefGene_Name\tUCSC_RefGene_Group\n"
        "cg1\tTP53;WRONG\tTSS200;Body\n"
        "cg2\tTP53;PIK3CA\tTSS1500;TSS200\n"
    )
    options = {
        "integrate": {
            "modalities": {
                "methylation_beta": {
                    "entity_column": "gene_symbol",
                    "entity_type": "gene_symbol",
                    "normalize": "uppercase",
                    "aggregation": "mean",
                    "mapping": {
                        "file": str(mapping),
                        "source_column": "IlmnID",
                        "entity_column": "UCSC_RefGene_Name",
                        "separator": ";",
                        "filters": {
                            "UCSC_RefGene_Group": ["TSS200", "TSS1500"],
                        },
                    },
                }
            }
        }
    }

    outputs = write_integrated_dataset(cohort, options)

    values = pl.read_parquet(outputs.values)
    p1_tp53 = values.filter(
        pl.col("submitter_id") == "P1",
        pl.col("entity_id") == "TP53",
    ).row(0, named=True)
    assert math.isclose(p1_tp53["value"], 0.3)
    assert p1_tp53["n_source_features"] == 2
    assert values.filter(pl.col("entity_id") == "WRONG").is_empty()
    assert (
        values.filter(
            pl.col("submitter_id") == "P1",
            pl.col("entity_id") == "PIK3CA",
        ).row(0, named=True)["value"]
        == 0.4
    )

    mappings = pl.read_parquet(outputs.mappings)
    assert mappings.filter(pl.col("mapping_type") == "external").height == 3


def test_custom_modality_uses_same_entity_contract(tmp_path: Path):
    cohort = _cohort(tmp_path)
    pl.DataFrame(
        [
            {"submitter_id": "P1", "marker": "m1", "gene": "tp53", "score": 2.5},
            {"submitter_id": "P2", "marker": "m2", "gene": "kras", "score": 1.5},
        ]
    ).write_parquet(cohort / "custom.parquet")

    outputs = write_integrated_dataset(
        cohort,
        {
            "integrate": {
                "modalities": {
                    "custom_score": {
                        "source_file": "custom.parquet",
                        "feature_column": "marker",
                        "entity_column": "gene",
                        "entity_type": "gene_symbol",
                        "value_column": "score",
                        "normalize": "uppercase",
                        "aggregation": "max",
                    }
                }
            }
        },
    )

    integrated = pl.read_parquet(outputs.integrated)
    assert _integrated_row(integrated, "P1", "TP53")["custom_score"] == 2.5


def test_sample_level_join_does_not_collapse_distinct_samples(tmp_path: Path):
    cohort = _cohort(tmp_path)
    pl.DataFrame(
        [
            {
                "case_id": "c1",
                "submitter_id": "P1",
                "sample_submitter_id": "P1-TUMOR",
                "gene_name": "TP53",
                "unstranded": 10,
            },
            {
                "case_id": "c1",
                "submitter_id": "P1",
                "sample_submitter_id": "P1-NORMAL",
                "gene_name": "TP53",
                "unstranded": 20,
            },
        ]
    ).write_parquet(cohort / "rna_expression.parquet")
    pl.DataFrame(
        [
            {
                "case_id": "c1",
                "submitter_id": "P1",
                "sample_submitter_id": "P1-TUMOR",
                "gene_name": "TP53",
                "copy_number": 2.0,
            }
        ]
    ).write_parquet(cohort / "gene_copy_number.parquet")

    outputs = write_integrated_dataset(
        cohort,
        {
            "integrate": {
                "join_on": "sample_submitter_id",
                "modalities": ["rna_expression", "gene_copy_number"],
            }
        },
    )

    integrated = pl.read_parquet(outputs.integrated).sort("join_id")
    assert integrated.select("join_id", "rna_expression", "gene_copy_number").to_dicts() == [
        {
            "join_id": "P1-NORMAL",
            "rna_expression": 20.0,
            "gene_copy_number": None,
        },
        {
            "join_id": "P1-TUMOR",
            "rna_expression": 10.0,
            "gene_copy_number": 2.0,
        },
    ]
    manifest = json.loads(outputs.manifest.read_text())
    assert manifest["join_on"] == "sample_submitter_id"
    assert manifest["n_cases"] == 1
    assert manifest["n_samples"] == 2


def test_explicit_modality_rejects_missing_merge_columns(tmp_path: Path):
    cohort = _cohort(tmp_path)
    pl.DataFrame([{"submitter_id": "P1", "unstranded": 2}]).write_parquet(
        cohort / "rna_expression.parquet"
    )

    with pytest.raises(ValueError, match=r"source is missing columns.*gene_name"):
        write_integrated_dataset(
            cohort,
            {"integrate": {"modalities": ["rna_expression"]}},
        )


def test_external_mapping_requires_explicit_aggregation(tmp_path: Path):
    cohort = _cohort(tmp_path)
    mapping = tmp_path / "mapping.tsv"
    mapping.write_text("probe_id\tgene_symbol\ncg1\tTP53\n")

    with pytest.raises(ValueError, match="aggregation is required"):
        write_integrated_dataset(
            cohort,
            {
                "integrate": {
                    "modalities": {
                        "methylation_beta": {
                            "entity_column": "gene_symbol",
                            "entity_type": "gene_symbol",
                            "mapping": {"file": str(mapping)},
                        }
                    }
                }
            },
        )


def test_integration_service_records_mapping_digest(tmp_path: Path):
    cohort = _cohort(tmp_path)
    mapping = tmp_path / "mapping.tsv"
    content = "probe_id\tgene_symbol\ncg1\tTP53\n"
    mapping.write_text(content)
    options = {
        "integrate": {
            "modalities": {
                "methylation_beta": {
                    "entity_column": "gene_symbol",
                    "entity_type": "gene_symbol",
                    "aggregation": "mean",
                    "mapping": {"file": str(mapping)},
                }
            }
        }
    }

    write_integration_recipe(cohort, recipe_options=options)

    provenance = json.loads((cohort / "cohort.json").read_text())
    mapping_input = provenance["recipes"]["inputs"]["integrate"]["modalities"]["methylation_beta"][
        "mapping_file"
    ]
    assert mapping_input["source"] == str(mapping.resolve())
    assert mapping_input["sha256"] == hashlib.sha256(content.encode()).hexdigest()
    assert mapping_input["size_bytes"] == len(content)
