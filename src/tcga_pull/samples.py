"""Build samples.parquet — one row per case — from clinical.parquet + variants.parquet.

Schema:

    # identifiers
    case_id, submitter_id, project_id, program
    # raw GDC labels (kept for traceability)
    primary_site, disease_type, primary_diagnosis
    # lineage (`lineage` = curated tissue label; `oncotree_*` = MSKCC OncoTree
    # crosswalk from project_id, see oncotree.py)
    lineage, oncotree_code, oncotree_name, oncotree_main_type, oncotree_tissue
    # demographics
    gender, race, ethnicity, age_at_diagnosis_days, age_at_diagnosis_years,
    vital_status, days_to_death, days_to_last_followup, ajcc_pathologic_stage
    # pair structure (derived from variants)
    n_tumor_aliquots, primary_tumor_barcode, primary_normal_barcode, normal_source
    # mutation burden (derived from variants where primary_aliquot=True)
    n_variants_total, n_variants_coding, n_variants_high_impact

Burden columns count variants on the primary tumor aliquot only, so they are
patient-level and free of multi-aliquot double-counting. Patients with no
matched variants (rare — e.g. a MAF that went to data/_multi/) get zeros for
burden and nulls for the tumor/normal barcodes.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

# Map flatten_case() output → samples.parquet column name
CLINICAL_COLUMN_MAP: dict[str, str] = {
    "case_id": "case_id",
    "submitter_id": "submitter_id",
    "project_id": "project_id",
    "primary_site": "primary_site",
    "disease_type": "disease_type",
    "demographic_gender": "gender",
    "demographic_race": "race",
    "demographic_ethnicity": "ethnicity",
    "demographic_vital_status": "vital_status",
    "demographic_days_to_death": "days_to_death",
    "diagnosis_primary_diagnosis": "primary_diagnosis",
    "diagnosis_age_at_diagnosis": "age_at_diagnosis_days",
    "diagnosis_days_to_last_follow_up": "days_to_last_followup",
    "diagnosis_ajcc_pathologic_stage": "ajcc_pathologic_stage",
}

OUTPUT_COLUMN_ORDER: list[str] = [
    # identifiers
    "case_id",
    "submitter_id",
    "project_id",
    "program",
    # raw GDC labels
    "primary_site",
    "disease_type",
    "primary_diagnosis",
    # lineage
    "lineage",
    "oncotree_code",
    "oncotree_name",
    "oncotree_main_type",
    "oncotree_tissue",
    # demographics
    "gender",
    "race",
    "ethnicity",
    "age_at_diagnosis_days",
    "age_at_diagnosis_years",
    "vital_status",
    "days_to_death",
    "days_to_last_followup",
    "ajcc_pathologic_stage",
    # pair structure
    "n_tumor_aliquots",
    "primary_tumor_barcode",
    "primary_normal_barcode",
    "normal_source",
    # burden (primary aliquot only)
    "n_variants_total",
    "n_variants_coding",
    "n_variants_high_impact",
]


def _project_to_program(project_id: str | None) -> str | None:
    """TCGA-BRCA → TCGA, CGCI-BLGSP → CGCI, BEATAML1.0-COHORT → BEATAML1.0."""
    if not isinstance(project_id, str) or "-" not in project_id:
        return None
    return project_id.split("-", 1)[0]


def _per_case_pair_structure(variants: pd.DataFrame) -> pd.DataFrame:
    """One row per submitter_id: primary tumor + normal barcode, normal source,
    total distinct tumor aliquots."""
    if variants.empty:
        return pd.DataFrame(
            columns=[
                "submitter_id",
                "primary_tumor_barcode",
                "primary_normal_barcode",
                "normal_source",
                "n_tumor_aliquots",
            ]
        )

    # Collapse to one row per (submitter, tumor_barcode) — picks any normal_barcode/source for that pair
    pairs = variants[
        ["submitter_id", "tumor_barcode", "normal_barcode", "normal_source", "primary_aliquot"]
    ].drop_duplicates(["submitter_id", "tumor_barcode"])

    n_aliquots = pairs.groupby("submitter_id").size().rename("n_tumor_aliquots")

    primary = pairs[pairs["primary_aliquot"].astype(bool)].drop_duplicates("submitter_id").copy()
    primary = primary.rename(
        columns={
            "tumor_barcode": "primary_tumor_barcode",
            "normal_barcode": "primary_normal_barcode",
        }
    )
    primary = primary[
        ["submitter_id", "primary_tumor_barcode", "primary_normal_barcode", "normal_source"]
    ]

    out: pd.DataFrame = primary.merge(n_aliquots.reset_index(), on="submitter_id", how="outer")
    return out


def _per_case_burden(variants: pd.DataFrame) -> pd.DataFrame:
    """Mutation burden on the primary tumor aliquot only. One row per submitter_id."""
    if variants.empty:
        return pd.DataFrame(
            columns=[
                "submitter_id",
                "n_variants_total",
                "n_variants_coding",
                "n_variants_high_impact",
            ]
        )
    primary = variants[variants["primary_aliquot"].astype(bool)]
    grouped = primary.groupby("submitter_id")
    out = pd.DataFrame(
        {
            "n_variants_total": grouped.size(),
            "n_variants_coding": grouped["is_coding"].sum(),
            "n_variants_high_impact": grouped["is_high_impact"].sum(),
        }
    ).reset_index()
    for col in ("n_variants_total", "n_variants_coding", "n_variants_high_impact"):
        out[col] = out[col].astype("Int64")
    burden_df: pd.DataFrame = out
    return burden_df


def build_samples_from_frames(
    clinical: pd.DataFrame,
    variants: pd.DataFrame,
) -> pd.DataFrame:
    """Pure transformation: clinical + variants → samples table.

    Left-joins variants-derived aggregates onto clinical so cases with no
    matched variants still appear (with zero burden + null barcodes).
    """
    # Project clinical to the columns we want, renaming on the way
    keep_in = [c for c in CLINICAL_COLUMN_MAP if c in clinical.columns]
    clin = clinical[keep_in].rename(columns=CLINICAL_COLUMN_MAP).copy()

    # Derived columns
    if "project_id" in clin.columns:
        from .tissue import derive_tissue

        clin["program"] = clin["project_id"].map(_project_to_program)
        # lineage = curated tissue (breast/lung/.../pancreas) with primary_site fallback
        clin["lineage"] = [
            derive_tissue(p, s)
            for p, s in zip(clin["project_id"], clin.get("primary_site", []), strict=False)
        ]

    if "age_at_diagnosis_days" in clin.columns:
        days = pd.to_numeric(clin["age_at_diagnosis_days"], errors="coerce")
        clin["age_at_diagnosis_years"] = (days / 365.25).round(2)

    # OncoTree crosswalk from project_id
    from .oncotree import oncotree_for

    onco_nodes = [oncotree_for(p) for p in clin.get("project_id", [])]
    clin["oncotree_code"] = [n.code if n else None for n in onco_nodes]
    clin["oncotree_name"] = [n.name if n else None for n in onco_nodes]
    clin["oncotree_main_type"] = [n.main_type if n else None for n in onco_nodes]
    clin["oncotree_tissue"] = [n.tissue if n else None for n in onco_nodes]

    # Join variant-derived aggregates
    pair_structure = _per_case_pair_structure(variants)
    burden = _per_case_burden(variants)
    samples = clin.merge(pair_structure, on="submitter_id", how="left").merge(
        burden, on="submitter_id", how="left"
    )

    # Fill burden zeros for cases without matched variants
    for col in ("n_variants_total", "n_variants_coding", "n_variants_high_impact"):
        if col in samples.columns:
            samples[col] = samples[col].fillna(0).astype("Int64")

    front = [c for c in OUTPUT_COLUMN_ORDER if c in samples.columns]
    rest = [c for c in samples.columns if c not in front]
    ordered: pd.DataFrame = samples[front + rest]
    return ordered


def write_samples(cohort_dir: Path) -> Path:
    """Read clinical.parquet + variants.parquet from cohort_dir, build
    samples.parquet next to them. Returns the output path."""
    cohort_dir = Path(cohort_dir)
    clin_path = cohort_dir / "clinical.parquet"
    variants_path = cohort_dir / "variants.parquet"
    if not clin_path.exists():
        raise FileNotFoundError(f"missing {clin_path} — run `tcga-pull pull` first")
    if not variants_path.exists():
        raise FileNotFoundError(
            f"missing {variants_path} — run `tcga-pull variants {cohort_dir}` first"
        )
    clinical = pd.read_parquet(clin_path)
    variants = pd.read_parquet(variants_path)
    samples = build_samples_from_frames(clinical, variants)
    out = cohort_dir / "samples.parquet"
    samples.to_parquet(out, index=False)
    return out
