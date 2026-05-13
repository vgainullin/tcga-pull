"""Polars implementation of the samples builder. Mirrors `samples.py`.

Loads clinical.parquet + variants.parquet, builds the per-case sample table
the same way the pandas version does, writes samples.parquet.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from .samples import CLINICAL_COLUMN_MAP, OUTPUT_COLUMN_ORDER, _project_to_program


def _per_case_pair_structure(variants: pl.DataFrame) -> pl.DataFrame:
    """One row per submitter_id: primary tumor + normal barcode, normal source,
    total distinct tumor aliquots."""
    if variants.is_empty():
        return pl.DataFrame(
            {
                "submitter_id": [],
                "primary_tumor_barcode": [],
                "primary_normal_barcode": [],
                "normal_source": [],
                "n_tumor_aliquots": [],
            }
        )

    pairs = variants.select(
        ["submitter_id", "tumor_barcode", "normal_barcode", "normal_source", "primary_aliquot"]
    ).unique(subset=["submitter_id", "tumor_barcode"])

    n_aliquots = pairs.group_by("submitter_id").agg(pl.len().alias("n_tumor_aliquots"))

    primary = (
        pairs.filter(pl.col("primary_aliquot"))
        .unique(subset=["submitter_id"])
        .rename(
            {
                "tumor_barcode": "primary_tumor_barcode",
                "normal_barcode": "primary_normal_barcode",
            }
        )
        .select(
            ["submitter_id", "primary_tumor_barcode", "primary_normal_barcode", "normal_source"]
        )
    )
    return primary.join(n_aliquots, on="submitter_id", how="full", coalesce=True)


def _per_case_burden(variants: pl.DataFrame) -> pl.DataFrame:
    """Mutation burden on the primary aliquot only. One row per submitter_id."""
    if variants.is_empty():
        return pl.DataFrame(
            {
                "submitter_id": [],
                "n_variants_total": [],
                "n_variants_coding": [],
                "n_variants_high_impact": [],
            }
        )
    return (
        variants.filter(pl.col("primary_aliquot"))
        .group_by("submitter_id")
        .agg(
            n_variants_total=pl.len().cast(pl.Int64),
            n_variants_coding=pl.col("is_coding").cast(pl.Int64).sum().cast(pl.Int64),
            n_variants_high_impact=pl.col("is_high_impact").cast(pl.Int64).sum().cast(pl.Int64),
        )
    )


def build_samples_from_frames(
    clinical: pl.DataFrame,
    variants: pl.DataFrame,
) -> pl.DataFrame:
    """Pure transformation: clinical + variants → samples table."""
    keep_in = [c for c in CLINICAL_COLUMN_MAP if c in clinical.columns]
    clin = clinical.select(keep_in).rename({k: CLINICAL_COLUMN_MAP[k] for k in keep_in})

    if "project_id" in clin.columns:
        from .tissue import derive_tissue

        clin = clin.with_columns(
            program=pl.col("project_id").map_elements(_project_to_program, return_dtype=pl.Utf8),
            lineage=pl.struct(["project_id", "primary_site"]).map_elements(
                lambda r: derive_tissue(r.get("project_id"), r.get("primary_site")),
                return_dtype=pl.Utf8,
            ),
        )
    if "age_at_diagnosis_days" in clin.columns:
        clin = clin.with_columns(
            age_at_diagnosis_years=(
                pl.col("age_at_diagnosis_days").cast(pl.Float64) / 365.25
            ).round(2)
        )

    # OncoTree crosswalk from project_id
    from .oncotree import oncotree_for

    project_ids = clin["project_id"].to_list() if "project_id" in clin.columns else []
    onco_nodes = [oncotree_for(p) for p in project_ids]
    clin = clin.with_columns(
        oncotree_code=pl.Series([n.code if n else None for n in onco_nodes], dtype=pl.Utf8),
        oncotree_name=pl.Series([n.name if n else None for n in onco_nodes], dtype=pl.Utf8),
        oncotree_main_type=pl.Series(
            [n.main_type if n else None for n in onco_nodes], dtype=pl.Utf8
        ),
        oncotree_tissue=pl.Series(
            [n.tissue if n else None for n in onco_nodes], dtype=pl.Utf8
        ),
    )

    pair_structure = _per_case_pair_structure(variants)
    burden = _per_case_burden(variants)
    samples = clin.join(pair_structure, on="submitter_id", how="left").join(
        burden, on="submitter_id", how="left"
    )

    # Fill burden zeros for cases without matched variants
    for col in ("n_variants_total", "n_variants_coding", "n_variants_high_impact"):
        if col in samples.columns:
            samples = samples.with_columns(pl.col(col).fill_null(0).cast(pl.Int64))

    front = [c for c in OUTPUT_COLUMN_ORDER if c in samples.columns]
    rest = [c for c in samples.columns if c not in front]
    return samples.select(front + rest)


def write_samples(cohort_dir: Path) -> Path:
    cohort_dir = Path(cohort_dir)
    clin_path = cohort_dir / "clinical.parquet"
    variants_path = cohort_dir / "variants.parquet"
    if not clin_path.exists():
        raise FileNotFoundError(f"missing {clin_path}")
    if not variants_path.exists():
        raise FileNotFoundError(f"missing {variants_path}")
    clinical = pl.read_parquet(clin_path)
    variants = pl.read_parquet(variants_path)
    samples = build_samples_from_frames(clinical, variants)
    out = cohort_dir / "samples.parquet"
    samples.write_parquet(out)
    return out
