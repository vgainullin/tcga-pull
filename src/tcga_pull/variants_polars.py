"""Polars implementation of the variants aggregation pipeline.

Mirrors the public API of `variants.py` (read_maf, aggregate_mafs,
aggregate_cohort, write_variants) so the two can be benchmarked head-to-head
on the same cohort. The parquet schema is identical to the pandas version
within nullable-type equivalence (Int64/Float64/Utf8/Boolean).

Key differences vs pandas:
  * Lazy by default — column pushdown means we only read 25 cols from disk
    instead of reading all 140 and dropping later.
  * Multi-threaded under the hood.
  * Strict typing — no silent coercion between Int64 and float NaN.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from .variants import (
    _TCGA_NORMAL_TYPES,
    HIGH_IMPACT_LEVELS,
    MAF_KEEP,
    MAF_RENAME,
    NON_CODING_CLASSES,
    OUTPUT_COLUMN_ORDER,
    RARE_GNOMAD_AF_THRESHOLD,
)

# Per-column schema hints; everything else is left to inference.
# polars accepts type classes here at runtime; mypy thinks the dict value
# should be DataType instances. Suppress the noise.
_SCAN_SCHEMA: dict[str, type[pl.DataType]] = {
    "Chromosome": pl.Utf8,
    "Start_Position": pl.Int64,
    "End_Position": pl.Int64,
    "t_depth": pl.Int64,
    "t_alt_count": pl.Int64,
    "gnomAD_AF": pl.Float64,
}

# MAFs use "" for missing; some columns also use ".".
_NULL_VALUES: list[str] = ["", "."]


def _scan_maf(path: Path) -> pl.LazyFrame:
    """Lazy scan of one MAF, projected to MAF_KEEP and renamed to snake_case.

    `scan_csv` pushes the column selection down into the reader, so only the
    columns we care about are decompressed from the gzip.
    """
    lf = pl.scan_csv(
        path,
        separator="\t",
        comment_prefix="#",
        schema_overrides=_SCAN_SCHEMA,
        null_values=_NULL_VALUES,
        infer_schema_length=2000,
        ignore_errors=False,
        truncate_ragged_lines=False,
    )
    # Take only columns we actually keep — polars pushes this into scan
    lf = lf.select([c for c in MAF_KEEP])
    # Rename to snake_case
    lf = lf.rename({k: v for k, v in MAF_RENAME.items() if k in MAF_KEEP})
    # vaf = t_alt_count / t_depth (Float64; null where t_depth is null or 0)
    lf = lf.with_columns(
        vaf=(pl.col("t_alt_count").cast(pl.Float64) / pl.col("t_depth").cast(pl.Float64))
    )
    return lf.with_columns(source_file=pl.lit(path.name))


def read_maf(path: Path) -> pl.DataFrame:
    """Eager read of one MAF; convenience wrapper around `_scan_maf`."""
    return _scan_maf(path).collect()


def _add_flags(lf: pl.LazyFrame) -> pl.LazyFrame:
    """Add is_coding / is_high_impact / is_rare / n_callers / normal_source.

    All vectorised; nothing per-row.
    """
    # is_coding: variant_class not in NON_CODING_CLASSES
    lf = lf.with_columns(
        is_coding=~pl.col("variant_class").is_in(list(NON_CODING_CLASSES)),
        is_high_impact=pl.col("impact").is_in(list(HIGH_IMPACT_LEVELS)),
        is_rare=(pl.col("gnomad_af").is_null() | (pl.col("gnomad_af") < RARE_GNOMAD_AF_THRESHOLD)),
    )
    # n_callers: count semicolons + 1 if non-empty, else 0
    lf = lf.with_columns(
        n_callers=pl.when(pl.col("callers").is_null() | (pl.col("callers") == ""))
        .then(0)
        .otherwise(pl.col("callers").str.count_matches(";") + 1)
        .cast(pl.Int64)
    )
    # normal_source: parse TCGA barcode position 13-14 -> known kind
    nn = pl.col("normal_barcode").str.split("-").list.get(3, null_on_oob=True).str.slice(0, 2)
    lf = lf.with_columns(
        normal_source=nn.replace_strict(_TCGA_NORMAL_TYPES, default=None, return_dtype=pl.Utf8)
    )
    return lf


def _mark_primary_aliquot(df: pl.DataFrame) -> pl.DataFrame:
    """For each submitter_id, mark one tumor_barcode as primary
    (highest mean t_depth; ties broken by barcode lex)."""
    if "submitter_id" not in df.columns or "tumor_barcode" not in df.columns:
        return df.with_columns(primary_aliquot=pl.lit(False))
    if "t_depth" not in df.columns:
        df = df.with_columns(_mean_depth=pl.lit(0.0))
    else:
        agg = (
            df.group_by(["submitter_id", "tumor_barcode"], maintain_order=False)
            .agg(pl.col("t_depth").mean().alias("_mean_depth"))
            .sort(
                ["submitter_id", "_mean_depth", "tumor_barcode"],
                descending=[False, True, False],
            )
        )
        primary = agg.unique(subset=["submitter_id"], keep="first")[["tumor_barcode"]]
    # Pass a plain list to is_in() so polars treats it as a value set, not
    # a column expression (avoids the "ambiguous collection" deprecation).
    primary_set = primary["tumor_barcode"].to_list()
    return df.with_columns(primary_aliquot=pl.col("tumor_barcode").is_in(primary_set))


def aggregate_mafs(paths: list[Path]) -> pl.DataFrame:
    """Read + project + flag + concat. Stays lazy until the final collect."""
    if not paths:
        return pl.DataFrame()
    scans = [_scan_maf(p) for p in paths]
    lf = pl.concat(scans, how="vertical_relaxed")
    lf = _add_flags(lf)
    return lf.collect()


def aggregate_cohort(cohort_dir: Path) -> pl.DataFrame:
    """End-to-end: walk cohort, read MAFs, add flags, join clinical, mark
    primary aliquot, order columns."""
    cohort_dir = Path(cohort_dir)
    mafs = sorted(cohort_dir.glob("data/*/simple_nucleotide_variation/*.maf.gz"))
    if not mafs:
        raise FileNotFoundError(
            f"no .maf.gz files under {cohort_dir}/data/*/simple_nucleotide_variation/"
        )

    variants = aggregate_mafs(mafs)

    # Join clinical (project_id, submitter_id, primary_diagnosis)
    clin_path = cohort_dir / "clinical.parquet"
    if clin_path.exists():
        clin = pl.read_parquet(clin_path)
        keep = [
            c
            for c in ("case_id", "submitter_id", "project_id", "diagnosis_primary_diagnosis")
            if c in clin.columns
        ]
        clin = clin.select(keep)
        if "diagnosis_primary_diagnosis" in clin.columns:
            clin = clin.rename({"diagnosis_primary_diagnosis": "primary_diagnosis"})
        variants = variants.join(clin, on="case_id", how="left")

    variants = _mark_primary_aliquot(variants)

    # Final column order — keep only those that exist
    front = [c for c in OUTPUT_COLUMN_ORDER if c in variants.columns]
    rest = [c for c in variants.columns if c not in front]
    return variants.select(front + rest)


def write_variants(cohort_dir: Path) -> Path:
    df = aggregate_cohort(cohort_dir)
    out = Path(cohort_dir) / "variants.parquet"
    df.write_parquet(out)
    return out
