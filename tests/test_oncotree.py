"""Tests for the OncoTree crosswalk."""

from __future__ import annotations

import pytest

from tcga_pull.oncotree import (
    ONCOTREE_VERSION,
    PROJECT_TO_ONCOTREE,
    _load_codes,
    get_node,
    oncotree_for,
)


def test_snapshot_loads_with_expected_size():
    codes = _load_codes()
    # Snapshot has ~900 codes; tolerate small drift on future updates
    assert len(codes) > 500
    assert ONCOTREE_VERSION.startswith("oncotree_")


def test_every_project_mapping_resolves_to_a_real_node():
    """The module's import-time _validate() should have caught this, but
    pin it explicitly so a future edit to PROJECT_TO_ONCOTREE doesn't break
    silently in a way that only manifests at sample-build time."""
    codes = _load_codes()
    for project_id, code in PROJECT_TO_ONCOTREE.items():
        assert code in codes, f"{project_id} -> {code} not in OncoTree snapshot"


def test_oncotree_for_canonical_mappings():
    n = oncotree_for("TCGA-BRCA")
    assert n is not None
    assert n.code == "BRCA"
    assert n.name == "Invasive Breast Carcinoma"
    assert n.main_type == "Breast Cancer"
    assert n.tissue == "Breast"


@pytest.mark.parametrize(
    "project_id, expected_code",
    [
        ("TCGA-LUAD", "LUAD"),
        ("TCGA-GBM", "GB"),  # Glioblastoma, IDH-Wildtype
        ("TCGA-TGCT", "TESTIS"),  # heterogeneous histology → parent tissue node
        ("MMRF-COMMPASS", "PCM"),  # Plasma Cell Myeloma
        ("TARGET-AML", "AML"),
        ("CCG-CUPP", "CUP"),
    ],
)
def test_specific_known_mappings(project_id: str, expected_code: str):
    n = oncotree_for(project_id)
    assert n is not None
    assert n.code == expected_code


def test_heterogeneous_projects_map_to_other():
    for project_id in ("CPTAC-2", "CPTAC-3", "HCMI-CMDC", "EXCEPTIONAL_RESPONDERS-ER"):
        n = oncotree_for(project_id)
        assert n is not None
        assert n.code == "OTHER"
        assert n.main_type == "Other Cancer"


def test_unknown_project_returns_none():
    assert oncotree_for("MADE-UP-PROJECT") is None
    assert oncotree_for(None) is None
    assert oncotree_for("") is None


def test_get_node_lookup():
    n = get_node("LUAD")
    assert n is not None
    assert n.name == "Lung Adenocarcinoma"
    assert get_node("not-a-real-code") is None
