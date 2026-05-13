"""Offline tests for samples.py — pure-frame transformation, no I/O."""

from __future__ import annotations

import pandas as pd

from tcga_pull.samples import (
    _per_case_burden,
    _per_case_pair_structure,
    _project_to_program,
    build_samples_from_frames,
)


def test_project_to_program():
    assert _project_to_program("TCGA-BRCA") == "TCGA"
    assert _project_to_program("CGCI-BLGSP") == "CGCI"
    assert _project_to_program("BEATAML1.0-COHORT") == "BEATAML1.0"
    assert _project_to_program(None) is None
    assert _project_to_program("malformed") is None


# -------------------------------------------------------- pair structure


def test_per_case_pair_structure_single_aliquot():
    v = pd.DataFrame(
        {
            "submitter_id": ["P1", "P1", "P1"],
            "tumor_barcode": ["T1", "T1", "T1"],
            "normal_barcode": ["N1", "N1", "N1"],
            "normal_source": ["Blood Derived Normal"] * 3,
            "primary_aliquot": [True, True, True],
        }
    )
    out = _per_case_pair_structure(v)
    row = out.iloc[0]
    assert row["submitter_id"] == "P1"
    assert row["primary_tumor_barcode"] == "T1"
    assert row["primary_normal_barcode"] == "N1"
    assert row["normal_source"] == "Blood Derived Normal"
    assert row["n_tumor_aliquots"] == 1


def test_per_case_pair_structure_multi_aliquot():
    # P1 has two tumor aliquots; only T1 is primary
    v = pd.DataFrame(
        {
            "submitter_id": ["P1", "P1", "P1"],
            "tumor_barcode": ["T1", "T2", "T2"],
            "normal_barcode": ["N1", "N1", "N1"],
            "normal_source": ["Solid Tissue Normal"] * 3,
            "primary_aliquot": [True, False, False],
        }
    )
    out = _per_case_pair_structure(v)
    row = out.iloc[0]
    assert row["primary_tumor_barcode"] == "T1"  # the primary, not T2
    assert row["n_tumor_aliquots"] == 2


# ----------------------------------------------------------------- burden


def test_per_case_burden_counts_only_primary():
    # P1 has 5 variants on its primary aliquot, 2 on a non-primary
    v = pd.DataFrame(
        {
            "submitter_id": ["P1"] * 7,
            "tumor_barcode": ["T1"] * 5 + ["T2"] * 2,
            "primary_aliquot": [True] * 5 + [False] * 2,
            "is_coding": [True, True, False, True, False, True, True],
            "is_high_impact": [True, False, False, True, False, True, True],
        }
    )
    out = _per_case_burden(v)
    row = out.iloc[0]
    assert row["submitter_id"] == "P1"
    assert row["n_variants_total"] == 5  # primary aliquot only
    assert row["n_variants_coding"] == 3  # of the primary 5
    assert row["n_variants_high_impact"] == 2


def test_per_case_burden_empty_frame_handled():
    out = _per_case_burden(pd.DataFrame())
    assert list(out.columns) == [
        "submitter_id",
        "n_variants_total",
        "n_variants_coding",
        "n_variants_high_impact",
    ]
    assert len(out) == 0


# -------------------------------------- end-to-end build_samples_from_frames


def test_build_samples_end_to_end():
    # Two cases; one has a matched variant, one doesn't (silent tumor)
    clinical = pd.DataFrame(
        [
            {
                "case_id": "c1",
                "submitter_id": "TCGA-XX-0001",
                "project_id": "TCGA-BRCA",
                "primary_site": "Breast",
                "disease_type": "Breast Invasive Carcinoma",
                "demographic_gender": "female",
                "demographic_race": "white",
                "demographic_vital_status": "Alive",
                "diagnosis_primary_diagnosis": "Infiltrating duct carcinoma, NOS",
                "diagnosis_age_at_diagnosis": 365.25 * 50,  # 50 years in days
            },
            {
                "case_id": "c2",
                "submitter_id": "TCGA-XX-0002",
                "project_id": "TCGA-BRCA",
                "primary_site": "Breast",
                "disease_type": "Breast Invasive Carcinoma",
                "demographic_gender": "male",
                "demographic_race": "asian",
                "demographic_vital_status": "Dead",
                "demographic_days_to_death": 1200,
                "diagnosis_primary_diagnosis": "Lobular carcinoma, NOS",
                "diagnosis_age_at_diagnosis": 365.25 * 65,
            },
        ]
    )
    variants = pd.DataFrame(
        {
            # Only c1 (TCGA-XX-0001) has variants
            "submitter_id": ["TCGA-XX-0001"] * 4,
            "tumor_barcode": ["T1"] * 4,
            "normal_barcode": ["N1"] * 4,
            "normal_source": ["Blood Derived Normal"] * 4,
            "primary_aliquot": [True] * 4,
            "is_coding": [True, True, False, True],
            "is_high_impact": [True, False, False, True],
        }
    )
    s = build_samples_from_frames(clinical, variants)
    assert len(s) == 2

    c1 = s[s["case_id"] == "c1"].iloc[0]
    assert c1["program"] == "TCGA"
    assert c1["lineage"] == "breast"  # TCGA-BRCA → breast via tissue.derive_tissue
    assert c1["age_at_diagnosis_years"] == 50.00
    assert c1["primary_tumor_barcode"] == "T1"
    assert c1["primary_normal_barcode"] == "N1"
    assert c1["normal_source"] == "Blood Derived Normal"
    assert c1["n_tumor_aliquots"] == 1
    assert c1["n_variants_total"] == 4
    assert c1["n_variants_coding"] == 3
    assert c1["n_variants_high_impact"] == 2
    # OncoTree crosswalk from project_id (TCGA-BRCA → BRCA)
    assert c1["oncotree_code"] == "BRCA"
    assert c1["oncotree_name"] == "Invasive Breast Carcinoma"
    assert c1["oncotree_main_type"] == "Breast Cancer"
    assert c1["oncotree_tissue"] == "Breast"

    c2 = s[s["case_id"] == "c2"].iloc[0]
    # c2 has no variants — burden zeros, barcodes null
    assert c2["n_variants_total"] == 0
    assert c2["n_variants_coding"] == 0
    assert c2["n_variants_high_impact"] == 0
    assert pd.isna(c2["primary_tumor_barcode"])
    assert pd.isna(c2["primary_normal_barcode"])
