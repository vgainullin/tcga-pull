"""OncoTree crosswalk: GDC project_id -> OncoTree code + tree lookup.

OncoTree is MSKCC's cancer-type ontology used by cBioPortal. The live API at
oncotree.mskcc.org is dead; the canonical source is now versioned JSON in
github.com/cBioPortal/oncotree. We bundle a pinned snapshot (no network at
runtime) and validate at import time that every mapped code exists.

Heterogeneous projects (CPTAC, HCMI, EXCEPTIONAL_RESPONDERS) map to the
`OTHER` node rather than guessing per-case.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files
from typing import Any

ONCOTREE_VERSION = "oncotree_2025_10_03"
_DATA_FILE = f"{ONCOTREE_VERSION}.json"


@dataclass(frozen=True)
class OncoNode:
    """One OncoTree node, flattened to the fields we actually use."""

    code: str
    name: str
    main_type: str
    tissue: str


@lru_cache(maxsize=1)
def _load_codes() -> dict[str, OncoNode]:
    raw_text = files("tcga_pull.data").joinpath(_DATA_FILE).read_text()
    raw = json.loads(raw_text)
    out: dict[str, OncoNode] = {}

    def walk(node: dict[str, Any]) -> None:
        code = node.get("code")
        if code:
            out[code] = OncoNode(
                code=code,
                name=node.get("name", ""),
                main_type=node.get("mainType", ""),
                tissue=node.get("tissue", ""),
            )
        for child in node.get("children", {}).values():
            walk(child)

    walk(raw["TISSUE"])
    return out


# GDC project_id -> OncoTree code. Codes are validated at import (see below).
# Choices favor the most-common histology when a TCGA project is heterogeneous
# at the histology level (e.g. TCGA-ESCA is adeno-leaning -> ESCA; TCGA-TGCT
# spans seminoma + non-seminoma -> the testis tissue node TESTIS).
PROJECT_TO_ONCOTREE: dict[str, str] = {
    # TCGA -------------------------------------------------------------------
    "TCGA-ACC": "ACC",
    "TCGA-BLCA": "BLCA",
    "TCGA-BRCA": "BRCA",
    "TCGA-CESC": "CESC",
    "TCGA-CHOL": "CHOL",
    "TCGA-COAD": "COAD",
    "TCGA-DLBC": "DLBCLNOS",
    "TCGA-ESCA": "ESCA",
    "TCGA-GBM": "GB",
    "TCGA-HNSC": "HNSC",
    "TCGA-KICH": "CHRCC",
    "TCGA-KIRC": "CCRCC",
    "TCGA-KIRP": "PRCC",
    "TCGA-LAML": "AML",
    "TCGA-LGG": "DIFG",
    "TCGA-LIHC": "HCC",
    "TCGA-LUAD": "LUAD",
    "TCGA-LUSC": "LUSC",
    "TCGA-MESO": "PLMESO",
    "TCGA-OV": "HGSOC",
    "TCGA-PAAD": "PAAD",
    "TCGA-PCPG": "PHC",
    "TCGA-PRAD": "PRAD",
    "TCGA-READ": "READ",
    "TCGA-SARC": "SARCNOS",
    "TCGA-SKCM": "SKCM",
    "TCGA-STAD": "STAD",
    "TCGA-TGCT": "TESTIS",
    "TCGA-THCA": "THPA",
    "TCGA-THYM": "THYM",
    "TCGA-UCEC": "UCEC",
    "TCGA-UCS": "UCS",
    "TCGA-UVM": "UM",
    # TARGET (pediatric) -----------------------------------------------------
    "TARGET-AML": "AML",
    "TARGET-ALL-P2": "BLL",
    "TARGET-ALL-P3": "BLL",
    "TARGET-NBL": "NBL",
    "TARGET-OS": "OS",
    "TARGET-WT": "WT",
    # Other programs ---------------------------------------------------------
    "MMRF-COMMPASS": "PCM",
    "ALCHEMIST-ALCH": "NSCLC",
    "MP2PRT-ALL": "BLL",
    "BEATAML1.0-COHORT": "AML",
    "CMI-MBC": "BRCA",
    "CMI-MPC": "PRAD",
    "CMI-ASC": "SARCNOS",
    "CGCI-HTMCP-CC": "CESC",
    "CGCI-BLGSP": "DLBCLNOS",
    "RC-PTCL": "PTCL",
    "CDDP_EAGLE-1": "LUAD",
    "CCG-CUPP": "CUP",
    # Heterogeneous (multi-cancer-type) projects -> OTHER --------------------
    "HCMI-CMDC": "OTHER",
    "CPTAC-2": "OTHER",
    "CPTAC-3": "OTHER",
    "EXCEPTIONAL_RESPONDERS-ER": "OTHER",
}


def _validate() -> None:
    """Assert every code in PROJECT_TO_ONCOTREE exists in the bundled tree."""
    codes = _load_codes()
    missing = {p: c for p, c in PROJECT_TO_ONCOTREE.items() if c not in codes}
    if missing:
        raise RuntimeError(
            f"OncoTree codes in PROJECT_TO_ONCOTREE not found in {_DATA_FILE}: {missing}"
        )


_validate()


def oncotree_for(project_id: str | None) -> OncoNode | None:
    """Resolve a GDC project_id to its OncoTree node. Returns None if unmapped."""
    if not isinstance(project_id, str):
        return None
    code = PROJECT_TO_ONCOTREE.get(project_id)
    if code is None:
        return None
    return _load_codes().get(code)


def get_node(code: str) -> OncoNode | None:
    """Look up any OncoTree node by code."""
    return _load_codes().get(code)
