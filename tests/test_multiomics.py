"""Offline tests for multiomics processing recipes."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import polars as pl

from tcga_pull.multiomics import write_multiomics


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
        "Composite Element REF\tBeta_value\n"
        "cg00000029\t0.75\n",
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
        "protein_id\tgene_symbol\tantibody\texpression_value\n"
        "TP53|p53\tTP53\tp53\t-0.2\n",
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
