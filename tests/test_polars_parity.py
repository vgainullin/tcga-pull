"""Parity tests: the polars implementation must produce the same parquet
contents as the pandas reference implementation, cell for cell.

We build a tiny synthetic cohort on disk (a few MAFs + clinical.parquet),
run both pipelines, and compare the resulting DataFrames after normalising
nullable types (pandas NaN vs polars Null) and column order.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd
import polars as pl

from tcga_pull import samples as samples_pandas
from tcga_pull import samples_polars, variants_polars
from tcga_pull import variants as variants_pandas


def _maf_row(
    *,
    chrom: str = "chr1",
    pos: int = 1000,
    hugo: str = "GENE1",
    case: str = "case-uuid-1",
    tumor: str = "TCGA-AA-AAAA-01A-11D-0000-09",
    normal: str = "TCGA-AA-AAAA-10A-01D-0000-09",
    variant_class: str = "Missense_Mutation",
    impact: str = "MODERATE",
    callers: str = "muse;mutect2;varscan2",
    t_depth: int = 100,
    t_alt: int = 25,
    gnomad_af: float | None = None,
) -> dict:
    return {
        "Hugo_Symbol": hugo,
        "Chromosome": chrom,
        "Start_Position": pos,
        "End_Position": pos,
        "Variant_Type": "SNP",
        "Reference_Allele": "C",
        "Tumor_Seq_Allele2": "T",
        "Variant_Classification": variant_class,
        "Consequence": "missense_variant",
        "IMPACT": impact,
        "t_depth": t_depth,
        "t_alt_count": t_alt,
        "Tumor_Sample_Barcode": tumor,
        "Matched_Norm_Sample_Barcode": normal,
        "case_id": case,
        "Transcript_ID": "ENST00000000001",
        "EXON": "5/12",
        "HGVSc": "c.100C>T",
        "HGVSp_Short": "p.P34L",
        "callers": callers,
        "gnomAD_AF": gnomad_af if gnomad_af is not None else "",
        "COSMIC": "",
        "SIFT": "",
        "PolyPhen": "",
        "CONTEXT": "ACGTACGT",
        "hotspot": "N",
    }


def _write_maf(path: Path, rows: list[dict]) -> None:
    cols = list(rows[0].keys())
    df = pd.DataFrame(rows, columns=cols)
    body = "#version gdc-1.0.0\n" + df.to_csv(sep="\t", index=False)
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wb") as fh:
        fh.write(body.encode("utf-8"))


def _build_fixture_cohort(tmp_path: Path) -> Path:
    cohort = tmp_path / "fixture"
    data = cohort / "data"

    # Patient 1 — single aliquot, three variants (one silent, one missense, one nonsense)
    _write_maf(
        data / "TCGA-AA-AAAA" / "simple_nucleotide_variation" / "P1.maf.gz",
        [
            _maf_row(
                pos=100,
                case="c1",
                tumor="TCGA-AA-AAAA-01A-11D-0000-09",
                normal="TCGA-AA-AAAA-10A-01D-0000-09",
                variant_class="Missense_Mutation",
                impact="MODERATE",
            ),
            _maf_row(
                pos=200,
                case="c1",
                tumor="TCGA-AA-AAAA-01A-11D-0000-09",
                normal="TCGA-AA-AAAA-10A-01D-0000-09",
                variant_class="Silent",
                impact="LOW",
                callers="mutect2",
            ),
            _maf_row(
                pos=300,
                case="c1",
                tumor="TCGA-AA-AAAA-01A-11D-0000-09",
                normal="TCGA-AA-AAAA-10A-01D-0000-09",
                variant_class="Nonsense_Mutation",
                impact="HIGH",
            ),
        ],
    )
    # Patient 2 — two aliquots so primary_aliquot logic runs
    _write_maf(
        data / "TCGA-BB-BBBB" / "simple_nucleotide_variation" / "P2a.maf.gz",
        [
            _maf_row(
                pos=400,
                case="c2",
                tumor="TCGA-BB-BBBB-01A-11D-0000-09",
                normal="TCGA-BB-BBBB-10A-01D-0000-09",
                impact="MODERATE",
                t_depth=50,
            ),
        ],
    )
    _write_maf(
        data / "TCGA-BB-BBBB" / "simple_nucleotide_variation" / "P2b.maf.gz",
        [
            _maf_row(
                pos=500,
                case="c2",
                tumor="TCGA-BB-BBBB-01B-04D-0000-09",
                normal="TCGA-BB-BBBB-10A-01D-0000-09",
                impact="MODERATE",
                t_depth=200,  # deeper coverage → wins primary_aliquot
            ),
        ],
    )

    # Minimal clinical.parquet — what flatten_case produces (subset of cols)
    clinical = pd.DataFrame(
        [
            {
                "case_id": "c1",
                "submitter_id": "TCGA-AA-AAAA",
                "project_id": "TCGA-XX",
                "primary_site": "Site",
                "disease_type": "Disease",
                "demographic_gender": "female",
                "demographic_race": "white",
                "demographic_vital_status": "Alive",
                "diagnosis_primary_diagnosis": "Diag One",
                "diagnosis_age_at_diagnosis": 365.25 * 50,
            },
            {
                "case_id": "c2",
                "submitter_id": "TCGA-BB-BBBB",
                "project_id": "TCGA-XX",
                "primary_site": "Site",
                "disease_type": "Disease",
                "demographic_gender": "male",
                "demographic_race": "asian",
                "demographic_vital_status": "Dead",
                "demographic_days_to_death": 1200,
                "diagnosis_primary_diagnosis": "Diag Two",
                "diagnosis_age_at_diagnosis": 365.25 * 65,
            },
        ]
    )
    clinical.to_parquet(cohort / "clinical.parquet", index=False)
    return cohort


# -------------------------------------------------------- helpers


def _assert_frames_match(pd_path: Path, pl_path: Path, sort_cols: list[str]) -> None:
    """Read both parquets back via polars (uniform null handling) and compare.
    Cells where both sides are null are considered equal; numerics use rtol.
    """
    a = pl.read_parquet(pd_path).sort(sort_cols)
    b = pl.read_parquet(pl_path).sort(sort_cols)

    assert a.shape == b.shape, f"shape mismatch: pandas={a.shape} polars={b.shape}"
    assert set(a.columns) == set(b.columns), (
        f"columns differ:\n  only in pandas: {set(a.columns) - set(b.columns)}\n"
        f"  only in polars: {set(b.columns) - set(a.columns)}"
    )

    # Align column order
    a = a.select(sorted(a.columns))
    b = b.select(sorted(b.columns))

    # Compare column-by-column so we get useful error messages
    diffs: list[str] = []
    for col in a.columns:
        ca, cb = a[col], b[col]
        # Both null → equal. Otherwise compare values.
        both_null = ca.is_null() & cb.is_null()
        if ca.dtype.is_numeric() and cb.dtype.is_numeric():
            # Cast both to Float64 for tolerant numeric compare
            ca_f = ca.cast(pl.Float64, strict=False)
            cb_f = cb.cast(pl.Float64, strict=False)
            close = ((ca_f - cb_f).abs() <= 1e-9 * cb_f.abs().fill_null(1.0)).fill_null(False)
            eq = both_null | close
        else:
            # Coerce to Utf8 for tolerant string-or-object compare
            ca_s = ca.cast(pl.Utf8, strict=False)
            cb_s = cb.cast(pl.Utf8, strict=False)
            eq = both_null | (ca_s == cb_s).fill_null(False)
        if not eq.all():
            n_diff = (~eq).sum()
            diffs.append(f"  {col}: {n_diff}/{len(eq)} rows differ")
    if diffs:
        raise AssertionError("Column diffs:\n" + "\n".join(diffs))


# -------------------------------------------------------- the actual tests


def test_variants_parity(tmp_path: Path):
    cohort = _build_fixture_cohort(tmp_path)

    # Pandas pipeline
    pd_path = variants_pandas.write_variants(cohort)
    pd_path = pd_path.rename(pd_path.with_suffix(".pandas.parquet"))

    # Polars pipeline (overwrites variants.parquet)
    pl_out = variants_polars.write_variants(cohort)
    pl_path = pl_out.rename(pl_out.with_suffix(".polars.parquet"))

    _assert_frames_match(
        pd_path,
        pl_path,
        sort_cols=["submitter_id", "tumor_barcode", "chrom", "pos"],
    )


def test_samples_parity(tmp_path: Path):
    cohort = _build_fixture_cohort(tmp_path)

    # Build variants.parquet via pandas (so samples has its input)
    variants_pandas.write_variants(cohort)

    pd_out = samples_pandas.write_samples(cohort)
    pd_path = pd_out.rename(pd_out.with_suffix(".pandas.parquet"))

    pl_out = samples_polars.write_samples(cohort)
    pl_path = pl_out.rename(pl_out.with_suffix(".polars.parquet"))

    _assert_frames_match(pd_path, pl_path, sort_cols=["submitter_id"])


# Quick standalone sanity that polars at least produces non-empty output
def test_polars_variants_writes_something(tmp_path: Path):
    cohort = _build_fixture_cohort(tmp_path)
    out = variants_polars.write_variants(cohort)
    df = pl.read_parquet(out)
    assert len(df) == 5  # 3 + 1 + 1 variants
    assert "primary_aliquot" in df.columns
    # Patient 2's deeper-coverage aliquot should be the primary
    p2 = df.filter(pl.col("submitter_id") == "TCGA-BB-BBBB")
    primaries = p2.filter(pl.col("primary_aliquot"))["tumor_barcode"].unique().to_list()
    assert primaries == ["TCGA-BB-BBBB-01B-04D-0000-09"]


def test_pandas_marks_one_primary_per_patient(tmp_path: Path):
    cohort = _build_fixture_cohort(tmp_path)
    out = variants_pandas.write_variants(cohort)
    df_pd = pd.read_parquet(out)
    per_patient = df_pd[df_pd["primary_aliquot"]].groupby("submitter_id")["tumor_barcode"].nunique()
    assert (per_patient == 1).all()


def test_polars_marks_one_primary_per_patient(tmp_path: Path):
    cohort = _build_fixture_cohort(tmp_path)
    out = variants_polars.write_variants(cohort)
    df_pl = pl.read_parquet(out)
    per_patient = (
        df_pl.filter(pl.col("primary_aliquot"))
        .group_by("submitter_id")
        .agg(pl.col("tumor_barcode").n_unique().alias("n"))
    )
    assert (per_patient["n"] == 1).all()
