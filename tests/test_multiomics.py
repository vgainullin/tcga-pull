"""Offline tests for multiomics processing recipes."""

from __future__ import annotations

from pathlib import Path
from typing import cast

import pandas as pd
import polars as pl
from rich.console import Console

from tcga_pull.config import CohortSpec, ProcessingSpec
from tcga_pull.gdc import GDCClient
from tcga_pull.multiomics import (
    finalize_multiomics_parts,
    record_handled_by_multiomics,
    write_multiomics,
    write_multiomics_parts,
)
from tcga_pull.pipeline import run as pipeline_run


def _write(path: Path, text: str) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return str(path)


def test_write_multiomics_outputs_typed_parquets(tmp_path: Path):
    cohort = tmp_path / "cohort"
    data = cohort / "data" / "TCGA-XX-0001"

    rna_path = _write(
        data / "transcriptome_profiling" / "rna.tsv",
        "gene_id\tgene_name\tgene_type\tunstranded\tstranded_first\tstranded_second\t"
        "tpm_unstranded\tfpkm_unstranded\tfpkm_uq_unstranded\n"
        "N_unmapped\tN_unmapped\tN/A\t10\t0\t0\t0\t0\t0\n"
        "ENSG00000141510.18\tTP53\tprotein_coding\t42\t40\t2\t8.5\t4.1\t9.2\n",
    )
    mirna_path = _write(
        data / "transcriptome_profiling" / "mirna.tsv",
        "miRNA_ID\tread_count\treads_per_million_miRNA_mapped\tcross-mapped\n"
        "hsa-let-7a-1\t100\t25.5\tN\n",
    )
    methylation_path = _write(
        data / "dna_methylation" / "meth.tsv",
        "Composite Element REF\tBeta_value\ncg00000029\t0.75\n",
    )
    cnv_path = _write(
        data / "copy_number_variation" / "seg.tsv",
        "Sample\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\n"
        "TCGA-XX-0001-01A\t1\t100\t200\t7\t0.34\n",
    )
    gene_cnv_path = _write(
        data / "copy_number_variation" / "gene.tsv",
        "gene_id\tgene_name\tchromosome\tstart\tend\tcopy_number\n"
        "ENSG00000141510\tTP53\t17\t7668402\t7687550\t1.5\n",
    )
    protein_path = _write(
        data / "proteome_profiling" / "rppa.tsv",
        "protein_id\tgene_symbol\tantibody\texpression_value\nTP53|p53\tTP53\tp53\t-0.2\n",
    )

    pd.DataFrame(
        [
            {
                "file_id": "rna-file",
                "file_name": "rna.tsv",
                "case_id": "case-1",
                "submitter_id": "TCGA-XX-0001",
                "data_category": "Transcriptome Profiling",
                "data_type": "Gene Expression Quantification",
                "experimental_strategy": "RNA-Seq",
                "workflow_type": "STAR - Counts",
                "local_path": rna_path,
                "status": "ok",
            },
            {
                "file_id": "mirna-file",
                "file_name": "mirna.tsv",
                "case_id": "case-1",
                "submitter_id": "TCGA-XX-0001",
                "data_category": "Transcriptome Profiling",
                "data_type": "miRNA Expression Quantification",
                "experimental_strategy": "miRNA-Seq",
                "workflow_type": None,
                "local_path": mirna_path,
                "status": "ok",
            },
            {
                "file_id": "meth-file",
                "file_name": "meth.tsv",
                "case_id": "case-1",
                "submitter_id": "TCGA-XX-0001",
                "data_category": "DNA Methylation",
                "data_type": "Methylation Beta Value",
                "experimental_strategy": "Methylation Array",
                "workflow_type": None,
                "local_path": methylation_path,
                "status": "ok",
            },
            {
                "file_id": "cnv-file",
                "file_name": "seg.tsv",
                "case_id": "case-1",
                "submitter_id": "TCGA-XX-0001",
                "data_category": "Copy Number Variation",
                "data_type": "Copy Number Segment",
                "experimental_strategy": "Genotyping Array",
                "workflow_type": None,
                "local_path": cnv_path,
                "status": "ok",
            },
            {
                "file_id": "gene-cnv-file",
                "file_name": "gene.tsv",
                "case_id": "case-1",
                "submitter_id": "TCGA-XX-0001",
                "data_category": "Copy Number Variation",
                "data_type": "Gene Level Copy Number",
                "experimental_strategy": "Genotyping Array",
                "workflow_type": None,
                "local_path": gene_cnv_path,
                "status": "ok",
            },
            {
                "file_id": "protein-file",
                "file_name": "rppa.tsv",
                "case_id": "case-1",
                "submitter_id": "TCGA-XX-0001",
                "data_category": "Proteome Profiling",
                "data_type": "Protein Expression Quantification",
                "experimental_strategy": "Reverse Phase Protein Array",
                "workflow_type": None,
                "local_path": protein_path,
                "status": "ok",
            },
        ]
    ).to_parquet(cohort / "manifest.parquet", index=False)

    outputs = write_multiomics(cohort)

    assert {p.name for p in outputs} == {
        "rna_expression.parquet",
        "mirna_expression.parquet",
        "methylation_beta.parquet",
        "copy_number_segments.parquet",
        "gene_copy_number.parquet",
        "protein_expression.parquet",
    }

    rna = pl.read_parquet(cohort / "rna_expression.parquet")
    assert rna.height == 1
    assert rna.row(0, named=True)["gene_name"] == "TP53"
    assert rna.row(0, named=True)["unstranded"] == 42

    mirna = pl.read_parquet(cohort / "mirna_expression.parquet")
    assert mirna.row(0, named=True)["mirna_id"] == "hsa-let-7a-1"
    assert mirna.row(0, named=True)["reads_per_million_mirna_mapped"] == 25.5

    methylation = pl.read_parquet(cohort / "methylation_beta.parquet")
    assert methylation.row(0, named=True)["probe_id"] == "cg00000029"
    assert methylation.row(0, named=True)["beta_value"] == 0.75

    cnv = pl.read_parquet(cohort / "copy_number_segments.parquet")
    assert cnv.row(0, named=True)["num_probes"] == 7
    assert cnv.row(0, named=True)["segment_mean"] == 0.34

    gene_cnv = pl.read_parquet(cohort / "gene_copy_number.parquet")
    assert gene_cnv.row(0, named=True)["gene_name"] == "TP53"
    assert gene_cnv.row(0, named=True)["copy_number"] == 1.5

    protein = pl.read_parquet(cohort / "protein_expression.parquet")
    assert protein.row(0, named=True)["gene_symbol"] == "TP53"
    assert protein.row(0, named=True)["expression_value"] == -0.2


def test_write_multiomics_parts_finalizes_directory_parquet(tmp_path: Path):
    cohort = tmp_path / "cohort"
    data = cohort / "data" / "TCGA-XX-0001"
    rna_path = _write(
        data / "transcriptome_profiling" / "rna.tsv",
        "gene_id\tgene_name\tgene_type\tunstranded\nENSG00000141510.18\tTP53\tprotein_coding\t42\n",
    )
    records = [
        {
            "file_id": "rna-file",
            "file_name": "rna.tsv",
            "case_id": "case-1",
            "submitter_id": "TCGA-XX-0001",
            "data_category": "Transcriptome Profiling",
            "data_type": "Gene Expression Quantification",
            "experimental_strategy": "RNA-Seq",
            "workflow_type": "STAR - Counts",
            "local_path": rna_path,
            "status": "ok",
        }
    ]

    recipe_options = {"rna_expression": {"columns": ["gene_id", "gene_name", "unstranded"]}}
    part_paths = write_multiomics_parts(
        cohort,
        records,
        part_id=1,
        recipes=["rna_expression"],
        recipe_options=recipe_options,
    )
    outputs = finalize_multiomics_parts(
        cohort,
        recipes=["rna_expression"],
        recipe_options=recipe_options,
    )

    assert [p.name for p in part_paths] == ["part-000001.parquet"]
    assert [p.name for p in outputs] == ["rna_expression.parquet"]
    assert (cohort / "rna_expression.parquet").is_dir()
    df = pl.read_parquet(cohort / "rna_expression.parquet")
    assert df.height == 1
    assert df.columns == [
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
    assert df.row(0, named=True)["gene_name"] == "TP53"


def test_methylation_probe_allowlist_in_parts(tmp_path: Path):
    cohort = tmp_path / "cohort"
    data = cohort / "data" / "TCGA-XX-0001"
    methylation_path = _write(
        data / "dna_methylation" / "meth.tsv",
        "Composite Element REF\tBeta_value\ncg00000029\t0.75\ncg99999999\t0.10\n",
    )
    records = [
        {
            "file_id": "meth-file",
            "file_name": "meth.tsv",
            "case_id": "case-1",
            "submitter_id": "TCGA-XX-0001",
            "data_category": "DNA Methylation",
            "data_type": "Methylation Beta Value",
            "experimental_strategy": "Methylation Array",
            "workflow_type": None,
            "local_path": methylation_path,
            "status": "ok",
        }
    ]

    write_multiomics_parts(
        cohort,
        records,
        part_id=1,
        recipes=["methylation"],
        recipe_options={"methylation": {"probes": ["cg00000029"]}},
    )
    finalize_multiomics_parts(
        cohort,
        recipes=["methylation"],
        recipe_options={"methylation": {"probes": ["cg00000029"]}},
    )

    df = pl.read_parquet(cohort / "methylation_beta.parquet")
    assert df.height == 1
    assert df.row(0, named=True)["probe_id"] == "cg00000029"


def test_record_handled_by_multiomics_respects_selected_recipes():
    rna_record = {
        "data_category": "Transcriptome Profiling",
        "data_type": "Gene Expression Quantification",
        "experimental_strategy": "RNA-Seq",
    }
    methylation_record = {
        "data_category": "DNA Methylation",
        "data_type": "Methylation Beta Value",
    }

    assert record_handled_by_multiomics(rna_record, ["multiomics"])
    assert record_handled_by_multiomics(rna_record, ["rna_expression"])
    assert not record_handled_by_multiomics(rna_record, ["methylation"])
    assert record_handled_by_multiomics(methylation_record, ["methylation"])


def test_incremental_pipeline_processes_and_deletes_handled_raw(tmp_path: Path, monkeypatch):
    class FakeClient:
        def fetch_files(self, flt: dict) -> list[dict]:
            return [
                {
                    "file_id": "rna-file",
                    "file_name": "rna.tsv",
                    "data_category": "Transcriptome Profiling",
                    "data_type": "Gene Expression Quantification",
                    "data_format": "TSV",
                    "experimental_strategy": "RNA-Seq",
                    "analysis": {"workflow_type": "STAR - Counts"},
                    "md5sum": "md5",
                    "file_size": 128,
                    "cases": [
                        {
                            "case_id": "case-1",
                            "submitter_id": "TCGA-XX-0001",
                            "project": {"project_id": "TCGA-CHOL"},
                        }
                    ],
                }
            ]

        def fetch_clinical(self, flt: dict) -> list[dict]:
            return [
                {
                    "case_id": "case-1",
                    "submitter_id": "TCGA-XX-0001",
                    "project": {"project_id": "TCGA-CHOL"},
                }
            ]

    def fake_bulk(file_hits, download_dir, **kwargs):
        for hit in file_hits:
            path = Path(download_dir) / hit["file_id"] / hit["file_name"]
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(
                "gene_id\tgene_name\tgene_type\tunstranded\n"
                "ENSG00000141510.18\tTP53\tprotein_coding\t42\n"
            )

    monkeypatch.setattr("tcga_pull.pipeline.bulk_download_via_api", fake_bulk)

    spec = CohortSpec(
        name="incremental",
        out_dir=tmp_path,
        filters={"project": "TCGA-CHOL"},
        recipes=["multiomics"],
        processing=ProcessingSpec(
            mode="incremental",
            batch_size=1,
            delete_raw_after_processing=True,
        ),
    )
    import io

    pipeline_run(spec, client=cast(GDCClient, FakeClient()), console=Console(file=io.StringIO()))
    cohort = tmp_path / "incremental"

    raw_path = cohort / "data" / "TCGA-XX-0001" / "transcriptome_profiling" / "rna.tsv"
    assert not raw_path.exists()

    manifest = pl.read_parquet(cohort / "manifest.parquet")
    assert manifest.row(0, named=True)["status"] == "processed_deleted"
    assert manifest.row(0, named=True)["local_path"] is None

    rna = pl.read_parquet(cohort / "rna_expression.parquet")
    assert rna.height == 1
    assert rna.row(0, named=True)["gene_name"] == "TP53"
