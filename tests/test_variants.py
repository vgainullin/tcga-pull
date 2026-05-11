"""Offline tests for the variants module — flag computations and primary-aliquot
selection. No network or downloaded MAFs required."""

from __future__ import annotations

import pandas as pd

from tcga_pull.variants import (
    HIGH_IMPACT_LEVELS,
    NON_CODING_CLASSES,
    RARE_GNOMAD_AF_THRESHOLD,
    _add_flags,
    _mark_primary_aliquot,
    _n_callers,
    _normal_source_from_barcode,
)

# ---------------------------------------------------------------- barcode parse


def test_normal_source_blood():
    # TCGA barcode position 13-14 = "10" means blood-derived normal
    assert _normal_source_from_barcode("TCGA-W5-AA39-10A-01D-A12B-09") == "Blood Derived Normal"


def test_normal_source_solid():
    assert _normal_source_from_barcode("TCGA-A7-A13G-11A-51D-A13O-09") == "Solid Tissue Normal"


def test_normal_source_tumor_returns_none():
    # "01" = primary tumor, not in the normal lookup
    assert _normal_source_from_barcode("TCGA-W5-AA39-01A-11R-A41B-07") is None


def test_normal_source_handles_garbage():
    assert _normal_source_from_barcode(None) is None
    assert _normal_source_from_barcode("") is None
    assert _normal_source_from_barcode("not-a-barcode") is None


# -------------------------------------------------------------------- callers


def test_n_callers_counts_semicolons_plus_one():
    s = pd.Series(["muse;mutect2;varscan2", "mutect2", "muse;mutect2", None, ""])
    out = _n_callers(s)
    assert out.tolist() == [3, 1, 2, 0, 0]


# ---------------------------------------------------------------------- flags


def test_add_flags_is_coding():
    df = pd.DataFrame(
        {
            "variant_class": [
                "Missense_Mutation",
                "Silent",
                "Nonsense_Mutation",
                "Intron",
                "Frame_Shift_Del",
                "5'UTR",
            ],
        }
    )
    out = _add_flags(df.copy())
    assert out["is_coding"].tolist() == [True, False, True, False, True, False]


def test_add_flags_is_high_impact():
    df = pd.DataFrame({"impact": ["HIGH", "MODERATE", "LOW", "MODIFIER", None]})
    out = _add_flags(df.copy())
    assert out["is_high_impact"].tolist() == [True, True, False, False, False]


def test_add_flags_is_rare_treats_nan_and_below_threshold_as_rare():
    df = pd.DataFrame({"gnomad_af": [None, 1e-5, 1e-3, 0.05, 0.99]})
    out = _add_flags(df.copy())
    # NaN → rare; strictly < 1e-3 → rare; >= 1e-3 → not rare
    assert out["is_rare"].tolist() == [True, True, False, False, False]


def test_add_flags_thresholds_are_consistent_with_constants():
    # The flags are computed against these constants — keep them stable so the
    # parquet schema doesn't drift silently between releases.
    assert "Silent" in NON_CODING_CLASSES
    assert "Missense_Mutation" not in NON_CODING_CLASSES
    assert {"HIGH", "MODERATE"} == HIGH_IMPACT_LEVELS
    assert RARE_GNOMAD_AF_THRESHOLD == 1e-3


# -------------------------------------------------------- primary_aliquot


def test_primary_aliquot_single_patient_single_barcode():
    df = pd.DataFrame(
        {
            "submitter_id": ["P1", "P1", "P1"],
            "tumor_barcode": ["T1", "T1", "T1"],
            "t_depth": [50, 80, 60],
        }
    )
    out = _mark_primary_aliquot(df.copy())
    assert out["primary_aliquot"].all()


def test_primary_aliquot_picks_highest_mean_depth():
    df = pd.DataFrame(
        {
            "submitter_id": ["P1", "P1", "P1", "P1"],
            "tumor_barcode": ["T1", "T1", "T2", "T2"],
            # T1 mean=20, T2 mean=100 → T2 wins
            "t_depth": [10, 30, 100, 100],
        }
    )
    out = _mark_primary_aliquot(df.copy())
    by_barcode = dict(zip(out["tumor_barcode"], out["primary_aliquot"], strict=False))
    assert by_barcode["T2"] is True or by_barcode["T2"] == True  # noqa: E712
    assert by_barcode["T1"] is False or by_barcode["T1"] == False  # noqa: E712


def test_primary_aliquot_ties_broken_by_barcode_lex():
    df = pd.DataFrame(
        {
            "submitter_id": ["P1", "P1"],
            "tumor_barcode": ["T2", "T1"],
            "t_depth": [50, 50],  # tie
        }
    )
    out = _mark_primary_aliquot(df.copy())
    by_barcode = dict(zip(out["tumor_barcode"], out["primary_aliquot"], strict=False))
    # T1 < T2 lexicographically → T1 wins the tie
    assert by_barcode["T1"]
    assert not by_barcode["T2"]


def test_primary_aliquot_independent_across_patients():
    df = pd.DataFrame(
        {
            "submitter_id": ["P1", "P1", "P2", "P2"],
            "tumor_barcode": ["T1", "T2", "T3", "T4"],
            "t_depth": [100, 10, 10, 100],
        }
    )
    out = _mark_primary_aliquot(df.copy())
    by_barcode = dict(zip(out["tumor_barcode"], out["primary_aliquot"], strict=False))
    assert by_barcode["T1"]  # highest in P1
    assert not by_barcode["T2"]
    assert not by_barcode["T3"]
    assert by_barcode["T4"]  # highest in P2
