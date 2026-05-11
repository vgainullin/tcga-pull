"""Tests for the project_id -> tissue mapping + primary_site fallback."""

from __future__ import annotations

import pytest

from tcga_pull.tissue import UNKNOWN, _normalize_primary_site, derive_tissue

# -------------------------------------------------------- project-pinned


@pytest.mark.parametrize(
    "project_id, expected",
    [
        ("TCGA-BRCA", "breast"),
        ("TCGA-LUAD", "lung"),
        ("TCGA-LUSC", "lung"),
        ("TCGA-LIHC", "liver"),
        ("TCGA-PAAD", "pancreas"),
        ("TCGA-GBM", "brain"),
        ("TCGA-LGG", "brain"),
        ("TCGA-LAML", "bone_marrow"),
        ("TARGET-AML", "bone_marrow"),
        ("MMRF-COMMPASS", "bone_marrow"),
        ("TARGET-NBL", "peripheral_nervous_system"),
        ("ALCHEMIST-ALCH", "lung"),
        ("CGCI-BLGSP", "lymph_node"),
        ("RC-PTCL", "lymph_node"),
    ],
)
def test_project_pinned_tissues(project_id: str, expected: str):
    # primary_site is None — the pinned mapping should win regardless
    assert derive_tissue(project_id, None) == expected


# ---------------------------------------------------- primary_site fallback


@pytest.mark.parametrize(
    "primary_site, expected",
    [
        ("Bronchus and lung", "lung"),
        ("Liver and intrahepatic bile ducts", "liver"),
        ("Hematopoietic and reticuloendothelial systems", "bone_marrow"),
        ("Breast", "breast"),
        ("Pancreas", "pancreas"),
        ("Kidney", "kidney"),
        ("Colon, Rectosigmoid junction", "colon"),  # comma-joined → first token
        ("BREAST", "breast"),  # case-insensitive
    ],
)
def test_primary_site_normalize(primary_site: str, expected: str):
    assert _normalize_primary_site(primary_site) == expected


def test_primary_site_unknown_strings():
    assert _normalize_primary_site(None) == UNKNOWN
    assert _normalize_primary_site("") == UNKNOWN
    assert _normalize_primary_site("Some Made Up Place") == UNKNOWN


# ----------------------------------------------------- heterogeneous projects


def test_heterogeneous_project_falls_back_to_primary_site():
    # CPTAC-3 maps to None — primary_site should drive the result
    assert derive_tissue("CPTAC-3", "Breast") == "breast"
    assert derive_tissue("CPTAC-3", "Pancreas") == "pancreas"
    assert derive_tissue("CPTAC-3", "Bronchus and lung") == "lung"


def test_unmapped_project_uses_primary_site():
    # A project we've never heard of: rely on primary_site
    assert derive_tissue("NEW-PROGRAM-XYZ", "Breast") == "breast"
    assert derive_tissue("NEW-PROGRAM-XYZ", None) == UNKNOWN


def test_pinned_project_ignores_primary_site():
    # If the project is pinned, the tissue label doesn't change even if the
    # primary_site disagrees (real data sometimes has weird values).
    assert derive_tissue("TCGA-BRCA", "Lung") == "breast"
