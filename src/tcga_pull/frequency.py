"""Per-lineage variant / gene frequency aggregations with three side-by-side
controls.

Reads `variants.parquet` and `samples.parquet` from a cohort dir; writes:

  gene_frequency.parquet     — one row per (hugo_symbol, lineage)
  variant_frequency.parquet  — one row per (chrom, pos, ref, alt, lineage)

Both tables carry:
  * raw cohort frequency (mutated patients / total patients in the lineage)
  * frequency in all OTHER lineages combined (cross-cancer comparator)
  * gnomAD comparator (population non-cancer signal)

The variant table denormalises gnomad_af directly from the source MAF. The
gene table attaches the max gnomad_af we observed across that gene's
variants, as a "is this gene ever populated in gnomAD" indicator. We do
not fabricate per-gene gnomAD carrier rates — that requires the full
gnomAD VCF and is a separate scope.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

# Smoothing for log2 ratios when one side is 0
_EPSILON: float = 1e-6


# ---------------------------------------------------------------------- gene


def _attach_lineage(variants: pl.DataFrame, samples: pl.DataFrame) -> pl.DataFrame:
    """Join `lineage` from samples onto variants by submitter_id. `lineage` is a
    per-case attribute that lives on samples.parquet, not variants.parquet."""
    if "lineage" in variants.columns:
        variants = variants.drop("lineage")
    return variants.join(
        samples.select(["submitter_id", "lineage"]).unique(subset=["submitter_id"]),
        on="submitter_id",
        how="inner",
    )


def gene_frequency(
    variants: pl.DataFrame,
    samples: pl.DataFrame,
    *,
    require_rare: bool = True,
    require_coding: bool = True,
) -> pl.DataFrame:
    """One row per (hugo_symbol, lineage) with rates + enrichment columns.

    Filters applied to variants before counting:
      - primary_aliquot=True (patient-level dedup across multi-aliquot cases)
      - is_coding=True if require_coding (drop silent/intron/UTR)
      - is_rare=True if require_rare (gnomAD AF < 1e-3 or null)
    """
    flt = pl.col("primary_aliquot")
    if require_coding:
        flt = flt & pl.col("is_coding")
    if require_rare:
        flt = flt & pl.col("is_rare")
    v = _attach_lineage(variants.filter(flt), samples)

    # Per-lineage total patients (from samples)
    n_per_lineage = samples.group_by("lineage").agg(pl.len().alias("n_total_patients"))
    total_patients = samples.height

    # Mutated patients per (gene, lineage)
    mutated = (
        v.group_by(["hugo_symbol", "lineage"])
        .agg(
            pl.col("submitter_id").n_unique().alias("n_mutated_patients"),
            (
                pl.col("submitter_id")
                .filter(pl.col("is_high_impact"))
                .n_unique()
                .alias("n_high_impact_patients")
            ),
        )
        .with_columns(pl.col("n_high_impact_patients").fill_null(0))
    )

    gf = mutated.join(n_per_lineage, on="lineage", how="left").with_columns(
        freq=(
            pl.col("n_mutated_patients").cast(pl.Float64)
            / pl.col("n_total_patients").cast(pl.Float64)
        ),
        freq_high_impact=(
            pl.col("n_high_impact_patients").cast(pl.Float64)
            / pl.col("n_total_patients").cast(pl.Float64)
        ),
    )

    # Pan-cancer counts per gene (sum across all lineages)
    gene_totals = mutated.group_by("hugo_symbol").agg(
        pl.col("n_mutated_patients").sum().alias("_pancan_mutated")
    )

    gf = (
        gf.join(gene_totals, on="hugo_symbol", how="left")
        .with_columns(
            n_mutated_other=(pl.col("_pancan_mutated") - pl.col("n_mutated_patients")),
            n_total_other=(pl.lit(total_patients) - pl.col("n_total_patients")),
        )
        .with_columns(
            freq_other_lineages=(
                pl.col("n_mutated_other").cast(pl.Float64)
                / pl.col("n_total_other").cast(pl.Float64)
            ),
        )
        .with_columns(
            log2_enrichment_vs_other=(
                ((pl.col("freq") + _EPSILON).log(base=2))
                - ((pl.col("freq_other_lineages") + _EPSILON).log(base=2))
            ),
        )
        .drop("_pancan_mutated")
    )

    # gnomAD signal: max population AF across observed variants in the gene.
    # "How common in non-cancer controls is the most polymorphic site we saw?"
    gnomad_max = (
        variants.filter(pl.col("gnomad_af").is_not_null())
        .group_by("hugo_symbol")
        .agg(pl.col("gnomad_af").max().alias("gnomad_max_af_in_gene"))
    )
    gf = gf.join(gnomad_max, on="hugo_symbol", how="left").with_columns(
        log2_enrichment_vs_gnomad=(
            ((pl.col("freq") + _EPSILON).log(base=2))
            - ((pl.col("gnomad_max_af_in_gene").fill_null(0.0) + _EPSILON).log(base=2))
        )
    )

    cols = [
        "hugo_symbol",
        "lineage",
        "n_mutated_patients",
        "n_high_impact_patients",
        "n_total_patients",
        "freq",
        "freq_high_impact",
        "n_mutated_other",
        "n_total_other",
        "freq_other_lineages",
        "log2_enrichment_vs_other",
        "gnomad_max_af_in_gene",
        "log2_enrichment_vs_gnomad",
    ]
    return gf.select([c for c in cols if c in gf.columns]).sort(
        ["lineage", "log2_enrichment_vs_other"], descending=[False, True]
    )


# ---------------------------------------------------------------------- variant


def variant_frequency(
    variants: pl.DataFrame,
    samples: pl.DataFrame,
    *,
    require_coding: bool = True,
) -> pl.DataFrame:
    """One row per (chrom, pos, ref, alt, lineage) with rates + enrichments.

    No is_rare filter here — the gnomad_af column is the comparison axis,
    so we keep population-common variants too. They show up with high
    gnomad_af and low log2 enrichment, which is the correct signal.
    """
    flt = pl.col("primary_aliquot")
    if require_coding:
        flt = flt & pl.col("is_coding")
    v = _attach_lineage(variants.filter(flt), samples)

    n_per_lineage = samples.group_by("lineage").agg(pl.len().alias("n_total_patients"))
    total_patients = samples.height

    grouped = (
        v.group_by(
            [
                "chrom",
                "pos",
                "ref",
                "alt",
                "hugo_symbol",
                "hgvsp_short",
                "variant_class",
                "consequence",
                "impact",
                "lineage",
            ]
        )
        .agg(
            pl.col("submitter_id").n_unique().alias("n_patients_with_variant"),
            pl.col("gnomad_af").max().alias("gnomad_af"),  # max across rows; same variant
        )
        .join(n_per_lineage, on="lineage", how="left")
        .with_columns(
            cohort_freq=(
                pl.col("n_patients_with_variant").cast(pl.Float64)
                / pl.col("n_total_patients").cast(pl.Float64)
            )
        )
    )

    # Vs other lineages: how often does THIS variant show up across all other lineages?
    variant_totals = grouped.group_by(["chrom", "pos", "ref", "alt"]).agg(
        pl.col("n_patients_with_variant").sum().alias("_pancan_with_variant")
    )
    grouped = (
        grouped.join(variant_totals, on=["chrom", "pos", "ref", "alt"], how="left")
        .with_columns(
            n_with_variant_other=(
                pl.col("_pancan_with_variant") - pl.col("n_patients_with_variant")
            ),
            n_total_other=(pl.lit(total_patients) - pl.col("n_total_patients")),
        )
        .with_columns(
            freq_other_lineages=(
                pl.col("n_with_variant_other").cast(pl.Float64)
                / pl.col("n_total_other").cast(pl.Float64)
            ),
        )
        .with_columns(
            log2_enrichment_vs_other=(
                ((pl.col("cohort_freq") + _EPSILON).log(base=2))
                - ((pl.col("freq_other_lineages") + _EPSILON).log(base=2))
            ),
            log2_enrichment_vs_gnomad=(
                ((pl.col("cohort_freq") + _EPSILON).log(base=2))
                - ((pl.col("gnomad_af").fill_null(0.0) + _EPSILON).log(base=2))
            ),
        )
        .drop("_pancan_with_variant")
    )

    cols = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "hugo_symbol",
        "hgvsp_short",
        "variant_class",
        "consequence",
        "impact",
        "lineage",
        "n_patients_with_variant",
        "n_total_patients",
        "cohort_freq",
        "n_with_variant_other",
        "n_total_other",
        "freq_other_lineages",
        "log2_enrichment_vs_other",
        "gnomad_af",
        "log2_enrichment_vs_gnomad",
    ]
    return grouped.select([c for c in cols if c in grouped.columns]).sort(
        ["lineage", "log2_enrichment_vs_other"], descending=[False, True]
    )


# ---------------------------------------------------------------- write helpers


def write_gene_frequency(cohort_dir: Path) -> Path:
    cohort_dir = Path(cohort_dir)
    variants = pl.read_parquet(cohort_dir / "variants.parquet")
    samples = pl.read_parquet(cohort_dir / "samples.parquet")
    out = cohort_dir / "gene_frequency.parquet"
    gene_frequency(variants, samples).write_parquet(out)
    return out


def write_variant_frequency(cohort_dir: Path) -> Path:
    cohort_dir = Path(cohort_dir)
    variants = pl.read_parquet(cohort_dir / "variants.parquet")
    samples = pl.read_parquet(cohort_dir / "samples.parquet")
    out = cohort_dir / "variant_frequency.parquet"
    variant_frequency(variants, samples).write_parquet(out)
    return out
