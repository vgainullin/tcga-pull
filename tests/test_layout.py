"""Offline tests for layout / download helpers — no network required."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from tcga_pull.download import primary_case, slugify
from tcga_pull.gdc import f_and, f_in, for_cases_endpoint
from tcga_pull.layout import flatten_case, write_clinical, write_manifest


def test_for_cases_endpoint_strips_cases_prefix():
    flt = f_and(
        f_in("cases.project.project_id", ["TCGA-CHOL"]),
        f_in("files.data_type", ["Gene Expression Quantification"]),
        f_in("analysis.workflow_type", ["STAR - Counts"]),
        f_in("cases.samples.sample_type", ["Primary Tumor"]),
    )
    out = for_cases_endpoint(flt)
    fields = [c["content"]["field"] for c in out["content"]]
    assert fields == [
        "project.project_id",  # cases. stripped
        "files.data_type",  # files. preserved (join)
        "files.analysis.workflow_type",  # analysis is rooted in file → re-prefix with files.
        "samples.sample_type",  # cases. stripped
    ]


def test_for_cases_endpoint_handles_nested_or():
    flt = {
        "op": "or",
        "content": [
            f_in("cases.project.project_id", ["TCGA-CHOL"]),
            f_in("cases.project.project_id", ["TCGA-BRCA"]),
        ],
    }
    out = for_cases_endpoint(flt)
    assert out["op"] == "or"
    assert all(c["content"]["field"] == "project.project_id" for c in out["content"])


def test_for_cases_endpoint_prefixes_bare_file_fields():
    # bare names on /files → files.X on /cases (it's a join from the case root)
    flt = f_and(
        f_in("access", ["open"]),
        f_in("data_format", ["MAF"]),
        f_in("data_category", ["Simple Nucleotide Variation"]),
    )
    out = for_cases_endpoint(flt)
    fields = [c["content"]["field"] for c in out["content"]]
    assert fields == ["files.access", "files.data_format", "files.data_category"]


def test_for_cases_endpoint_passes_through_already_prefixed_files_fields():
    # If a caller already wrote `files.X` (e.g. handwritten gdc_filter), leave it.
    flt = f_in("files.access", ["open"])
    out = for_cases_endpoint(flt)
    assert out["content"]["field"] == "files.access"


def test_f_and_flattens_nested_ands():
    # GDC /files 500s on faceted queries when the filter tree has nested ANDs.
    # f_and() must produce a single top-level `and` with all leaves.
    inner = f_and(f_in("a", [1]), f_in("b", [2]))
    outer = f_and(inner, f_in("c", [3]))
    assert outer["op"] == "and"
    fields = [clause["content"]["field"] for clause in outer["content"]]
    assert fields == ["a", "b", "c"], f"expected flat tree, got nested: {outer}"


def test_f_and_single_clause_unwraps():
    only = f_and(f_in("a", [1]))
    assert only["op"] == "in"
    assert only["content"]["field"] == "a"


def test_f_and_skips_empty_clauses():
    out = f_and(None, f_in("a", [1]), {})
    assert out["op"] == "in"
    assert out["content"]["field"] == "a"


def test_slugify():
    assert slugify("Transcriptome Profiling") == "transcriptome_profiling"
    assert slugify(None) == "unknown"
    assert slugify("DNA  Methylation/Beta") == "dna_methylation_beta"
    assert slugify("") == "unknown"


def test_primary_case_single():
    h = {"cases": [{"case_id": "c1", "submitter_id": "TCGA-XX-0001"}]}
    assert primary_case(h) == ("c1", "TCGA-XX-0001")


def test_primary_case_multi_returns_none():
    h = {"cases": [{"case_id": "c1"}, {"case_id": "c2"}]}
    assert primary_case(h) is None


def test_primary_case_no_cases():
    assert primary_case({"cases": []}) is None


def test_flatten_case():
    case = {
        "case_id": "c1",
        "submitter_id": "TCGA-XX-0001",
        "project": {"project_id": "TCGA-BRCA"},
        "primary_site": "Breast",
        "demographic": {"gender": "female", "race": "white", "vital_status": "Alive"},
        "diagnoses": [
            {
                "primary_diagnosis": "Infiltrating duct carcinoma, NOS",
                "tumor_stage": "stage iia",
                "treatments": [{"treatment_type": "Pharmaceutical Therapy, NOS"}],
            }
        ],
        "exposures": [{"alcohol_history": "No"}],
    }
    row = flatten_case(case)
    assert row["case_id"] == "c1"
    assert row["project_id"] == "TCGA-BRCA"
    assert row["demographic_gender"] == "female"
    assert row["diagnosis_primary_diagnosis"] == "Infiltrating duct carcinoma, NOS"
    assert row["n_diagnoses"] == 1
    assert row["n_treatments"] == 1
    assert row["treatment_treatment_type"] == "Pharmaceutical Therapy, NOS"
    assert row["exposure_alcohol_history"] == "No"


def test_write_clinical_and_manifest_roundtrip(tmp_path: Path):
    cases = [
        {"case_id": "c1", "submitter_id": "TCGA-XX-0001", "project": {"project_id": "TCGA-BRCA"}},
        {"case_id": "c2", "submitter_id": "TCGA-XX-0002", "project": {"project_id": "TCGA-BRCA"}},
    ]
    clinical_path, raw_path = write_clinical(cases, tmp_path)
    df = pd.read_parquet(clinical_path)
    assert len(df) == 2
    assert set(df["case_id"]) == {"c1", "c2"}
    assert raw_path.exists()
    assert raw_path.read_text().count("\n") == 2

    records = [
        {
            "file_id": "f1",
            "case_id": "c1",
            "submitter_id": "TCGA-XX-0001",
            "data_category": "Transcriptome Profiling",
            "local_path": str(
                tmp_path / "data" / "TCGA-XX-0001" / "transcriptome_profiling" / "x.tsv"
            ),
        }
    ]
    mpath = write_manifest(records, tmp_path)
    mdf = pd.read_parquet(mpath)
    assert len(mdf) == 1
    assert mdf.iloc[0]["case_id"] == "c1"
