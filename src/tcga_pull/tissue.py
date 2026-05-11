"""Curated GDC project_id → tissue mapping plus a primary_site fallback.

The output `tissue` label is the canonical organ-of-origin used as `lineage`
on samples.parquet. Designed for analyses asking "which mutations are
enriched in breast / lung / pancreas / liver / ...".

Project-pinned values are explicit; heterogeneous projects map to `None`
here and the caller falls back to a per-case `primary_site` cleanup.
"""

from __future__ import annotations

UNKNOWN: str = "unknown"

# Project ID → tissue label (lowercase snake_case).
# None = heterogeneous project, fall back to per-case primary_site.
PROJECT_TO_TISSUE: dict[str, str | None] = {
    # TCGA
    "TCGA-ACC": "adrenal",
    "TCGA-BLCA": "bladder",
    "TCGA-BRCA": "breast",
    "TCGA-CESC": "cervix",
    "TCGA-CHOL": "bile_duct",
    "TCGA-COAD": "colon",
    "TCGA-DLBC": "lymph_node",
    "TCGA-ESCA": "esophagus",
    "TCGA-GBM": "brain",
    "TCGA-HNSC": "head_and_neck",
    "TCGA-KICH": "kidney",
    "TCGA-KIRC": "kidney",
    "TCGA-KIRP": "kidney",
    "TCGA-LAML": "bone_marrow",
    "TCGA-LGG": "brain",
    "TCGA-LIHC": "liver",
    "TCGA-LUAD": "lung",
    "TCGA-LUSC": "lung",
    "TCGA-MESO": "pleura",
    "TCGA-OV": "ovary",
    "TCGA-PAAD": "pancreas",
    "TCGA-PCPG": "adrenal",
    "TCGA-PRAD": "prostate",
    "TCGA-READ": "rectum",
    "TCGA-SARC": "soft_tissue",
    "TCGA-SKCM": "skin",
    "TCGA-STAD": "stomach",
    "TCGA-TGCT": "testis",
    "TCGA-THCA": "thyroid",
    "TCGA-THYM": "thymus",
    "TCGA-UCEC": "uterus",
    "TCGA-UCS": "uterus",
    "TCGA-UVM": "eye",
    # TARGET — pediatric
    "TARGET-AML": "bone_marrow",
    "TARGET-ALL-P2": "bone_marrow",
    "TARGET-ALL-P3": "bone_marrow",
    "TARGET-NBL": "peripheral_nervous_system",
    "TARGET-OS": "bone",
    "TARGET-WT": "kidney",
    # other programs
    "MMRF-COMMPASS": "bone_marrow",
    "ALCHEMIST-ALCH": "lung",
    "MP2PRT-ALL": "bone_marrow",
    "BEATAML1.0-COHORT": "bone_marrow",
    "CMI-MBC": "breast",
    "CMI-MPC": "prostate",
    "CMI-ASC": "soft_tissue",
    "CGCI-HTMCP-CC": "cervix",
    "CGCI-BLGSP": "lymph_node",
    "RC-PTCL": "lymph_node",
    "CDDP_EAGLE-1": "lung",
    "CCG-CUPP": "unknown_primary",
    # heterogeneous → fall back to primary_site
    "HCMI-CMDC": None,
    "CPTAC-2": None,
    "CPTAC-3": None,
    "EXCEPTIONAL_RESPONDERS-ER": None,
}

# primary_site (ICD-O topography) → tissue. Keys are matched case-insensitive
# and via 'starts-with' against the cleaned-up site string. Order matters:
# the first match wins so put more specific patterns first.
_PRIMARY_SITE_RULES: list[tuple[str, str]] = [
    ("bronchus and lung", "lung"),
    ("liver and intrahepatic bile ducts", "liver"),
    ("hematopoietic and reticuloendothelial", "bone_marrow"),
    ("bones, joints and articular cartilage", "bone"),
    ("heart, mediastinum, and pleura", "pleura"),
    ("breast", "breast"),
    ("lung", "lung"),
    ("liver", "liver"),
    ("pancreas", "pancreas"),
    ("brain", "brain"),
    ("kidney", "kidney"),
    ("colon", "colon"),
    ("rectum", "rectum"),
    ("stomach", "stomach"),
    ("esophagus", "esophagus"),
    ("ovary", "ovary"),
    ("prostate", "prostate"),
    ("uterus", "uterus"),
    ("corpus uteri", "uterus"),
    ("cervix", "cervix"),
    ("bladder", "bladder"),
    ("skin", "skin"),
    ("thyroid", "thyroid"),
    ("testis", "testis"),
    ("thymus", "thymus"),
    ("adrenal", "adrenal"),
    ("eye", "eye"),
    ("lymph nodes", "lymph_node"),
    ("connective, subcutaneous", "soft_tissue"),
    ("bile ducts", "bile_duct"),
    ("head and neck", "head_and_neck"),
    ("lip", "head_and_neck"),
    ("tonsil", "head_and_neck"),
    ("gum", "head_and_neck"),
    ("biliary tract", "bile_duct"),
    ("spinal cord", "brain"),
]


def _normalize_primary_site(primary_site: str | None) -> str:
    """Map a (possibly multi-tissue, possibly messy) ICD-O `primary_site` to a
    tissue label. Returns `UNKNOWN` if nothing matches."""
    if not isinstance(primary_site, str) or not primary_site.strip():
        return UNKNOWN
    # GDC sometimes returns comma-joined multi-tissue strings; take the first.
    first = primary_site.split(",", 1)[0].strip().lower()
    if not first:
        return UNKNOWN
    for pat, label in _PRIMARY_SITE_RULES:
        if first.startswith(pat):
            return label
    return UNKNOWN


def derive_tissue(project_id: str | None, primary_site: str | None) -> str:
    """Resolve a tissue label for one case. Project-pinned mapping wins;
    fall back to primary_site normalization for heterogeneous projects."""
    if isinstance(project_id, str) and project_id in PROJECT_TO_TISSUE:
        pinned = PROJECT_TO_TISSUE[project_id]
        if pinned is not None:
            return pinned
    return _normalize_primary_site(primary_site)
