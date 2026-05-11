"""Walk a cohort produced by `tcga-pull pull` (Simple Nucleotide Variation MAFs)
and produce a tidy per-variant parquet.

Schema (one row per (variant x tumor aliquot)):

    # identifiers
    project_id, case_id, submitter_id, primary_diagnosis
    # locus
    chrom, pos, end_pos, ref, alt, variant_type
    # gene / transcript / consequence
    hugo_symbol, transcript_id, exon, hgvsc, hgvsp_short,
    variant_class, consequence, impact
    # support
    t_depth, t_alt_count, vaf
    # caller agreement
    callers, n_callers
    # population / functional annotations
    gnomad_af, cosmic_id, sift, polyphen, context, hotspot
    # convenience flags (cheap to recompute, costly to forget)
    is_coding, is_high_impact, is_rare
    # pair information
    tumor_barcode, normal_barcode, normal_source, primary_aliquot
    # provenance
    source_file

`vaf` = t_alt_count / t_depth (mutant allele frequency in tumor sample).
`primary_aliquot` is True for one tumor barcode per patient (picked by
highest mean t_depth) — filter on this for patient-level analysis to avoid
double-counting multi-aliquot patients.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

# Source columns we keep from the MAF (GDC's "Aliquot Ensemble Masked" schema).
# Order is just for readability; final column order is set at the end.
MAF_KEEP: list[str] = [
    # locus
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    # gene / transcript / consequence
    "Hugo_Symbol",
    "Transcript_ID",
    "EXON",
    "HGVSc",
    "HGVSp_Short",
    "Variant_Classification",
    "Consequence",
    "IMPACT",
    # support
    "t_depth",
    "t_alt_count",
    # caller agreement
    "callers",
    # annotations
    "gnomAD_AF",
    "COSMIC",
    "SIFT",
    "PolyPhen",
    "CONTEXT",
    "hotspot",
    # case / pair
    "case_id",
    "Tumor_Sample_Barcode",
    "Matched_Norm_Sample_Barcode",
]

MAF_RENAME: dict[str, str] = {
    "Chromosome": "chrom",
    "Start_Position": "pos",
    "End_Position": "end_pos",
    "Variant_Type": "variant_type",
    "Reference_Allele": "ref",
    "Tumor_Seq_Allele2": "alt",
    "Hugo_Symbol": "hugo_symbol",
    "Transcript_ID": "transcript_id",
    "EXON": "exon",
    "HGVSc": "hgvsc",
    "HGVSp_Short": "hgvsp_short",
    "Variant_Classification": "variant_class",
    "Consequence": "consequence",
    "IMPACT": "impact",
    "gnomAD_AF": "gnomad_af",
    "COSMIC": "cosmic_id",
    "SIFT": "sift",
    "PolyPhen": "polyphen",
    "CONTEXT": "context",
    "Tumor_Sample_Barcode": "tumor_barcode",
    "Matched_Norm_Sample_Barcode": "normal_barcode",
}

# variant_class values that we consider NON-coding (silent / regulatory / RNA).
# Everything else (Missense, Nonsense, Frame_Shift_*, Splice_Site, In_Frame_*,
# Translation_Start_Site, Nonstop) is_coding=True.
NON_CODING_CLASSES: frozenset[str] = frozenset(
    {
        "Silent",
        "Intron",
        "IGR",
        "3'UTR",
        "5'UTR",
        "3'Flank",
        "5'Flank",
        "RNA",
        "Targeted_Region",
    }
)

HIGH_IMPACT_LEVELS: frozenset[str] = frozenset({"HIGH", "MODERATE"})

RARE_GNOMAD_AF_THRESHOLD: float = 1e-3


# TCGA barcode position 13-14 encodes the sample type (NN in TCGA-XX-YYYY-NN[ABC]-...).
# Normals are 10-19; tumors are 01-09. We only enumerate normals here.
_TCGA_NORMAL_TYPES: dict[str, str] = {
    "10": "Blood Derived Normal",
    "11": "Solid Tissue Normal",
    "12": "Buccal Cell Normal",
    "13": "EBV Immortalized Normal",
    "14": "Bone Marrow Normal",
    "15": "Lymphoid Normal",
}


def _normal_source_from_barcode(barcode: str | None) -> str | None:
    """Parse a TCGA aliquot barcode and return the matched-normal source kind."""
    if not isinstance(barcode, str):
        return None
    parts = barcode.split("-")
    if len(parts) < 4:
        return None
    nn = parts[3][:2]
    return _TCGA_NORMAL_TYPES.get(nn)


def _n_callers(callers_field: pd.Series) -> pd.Series:
    """Count distinct callers in the semicolon-joined `callers` field."""
    s = callers_field.fillna("").astype(str)
    counts = s.str.count(";") + 1
    counts = counts.where(s != "", 0)
    out: pd.Series = counts.astype("Int64")
    return out


def _add_flags(df: pd.DataFrame) -> pd.DataFrame:
    """Add is_coding / is_high_impact / is_rare / n_callers / normal_source."""
    if "variant_class" in df.columns:
        df["is_coding"] = ~df["variant_class"].isin(NON_CODING_CLASSES)
    if "impact" in df.columns:
        df["is_high_impact"] = df["impact"].isin(HIGH_IMPACT_LEVELS)
    if "gnomad_af" in df.columns:
        af = pd.to_numeric(df["gnomad_af"], errors="coerce")
        df["is_rare"] = af.isna() | (af < RARE_GNOMAD_AF_THRESHOLD)
    if "callers" in df.columns:
        df["n_callers"] = _n_callers(df["callers"])
    if "normal_barcode" in df.columns:
        df["normal_source"] = df["normal_barcode"].map(_normal_source_from_barcode)
    return df


def _mark_primary_aliquot(df: pd.DataFrame) -> pd.DataFrame:
    """For each submitter_id with multiple tumor aliquots, mark exactly one as
    primary (highest mean t_depth → best coverage). Ties broken by barcode string."""
    if "submitter_id" not in df.columns or "tumor_barcode" not in df.columns:
        df["primary_aliquot"] = False
        return df
    if "t_depth" in df.columns:
        agg = df.groupby(["submitter_id", "tumor_barcode"], dropna=False)["t_depth"].mean()
        agg = agg.reset_index().rename(columns={"t_depth": "_mean_depth"})
    else:
        agg = df[["submitter_id", "tumor_barcode"]].drop_duplicates()
        agg["_mean_depth"] = 0
    # one tumor_barcode per submitter: highest mean depth, then lex barcode for stability
    agg = agg.sort_values(
        ["submitter_id", "_mean_depth", "tumor_barcode"],
        ascending=[True, False, True],
    )
    primary = set(agg.drop_duplicates("submitter_id")["tumor_barcode"])
    df["primary_aliquot"] = df["tumor_barcode"].isin(primary)
    return df


def read_maf(path: Path) -> pd.DataFrame:
    """Read one GDC MAF (gzipped, tab-separated, ## comment headers).

    Projects to MAF_KEEP, renames to snake_case, computes VAF. Flags + primary
    aliquot are computed in `aggregate_cohort` because they require either
    cohort-wide context or are cheap enough to defer.
    """
    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        low_memory=False,
        dtype={
            "Chromosome": str,
            "Start_Position": "Int64",
            "End_Position": "Int64",
            "t_depth": "Int64",
            "t_alt_count": "Int64",
        },
    )
    cols = [c for c in MAF_KEEP if c in df.columns]
    out: pd.DataFrame = df[cols].rename(columns=MAF_RENAME).copy()
    if "t_depth" in out and "t_alt_count" in out:
        # avoid divide-by-zero; pandas turns 0/0 into NaN naturally with float division
        out["vaf"] = out["t_alt_count"].astype("Float64") / out["t_depth"].astype("Float64")
    return out


# Column order in the output (only those that actually exist are emitted).
OUTPUT_COLUMN_ORDER: list[str] = [
    # identifiers
    "project_id",
    "case_id",
    "submitter_id",
    "primary_diagnosis",
    # locus
    "chrom",
    "pos",
    "end_pos",
    "ref",
    "alt",
    "variant_type",
    # gene / transcript / consequence
    "hugo_symbol",
    "transcript_id",
    "exon",
    "hgvsc",
    "hgvsp_short",
    "variant_class",
    "consequence",
    "impact",
    # support
    "t_depth",
    "t_alt_count",
    "vaf",
    # caller agreement
    "callers",
    "n_callers",
    # population / functional
    "gnomad_af",
    "cosmic_id",
    "sift",
    "polyphen",
    "context",
    "hotspot",
    # flags
    "is_coding",
    "is_high_impact",
    "is_rare",
    # pair
    "tumor_barcode",
    "normal_barcode",
    "normal_source",
    "primary_aliquot",
    # provenance
    "source_file",
]


def aggregate_cohort(cohort_dir: Path) -> pd.DataFrame:
    """Walk cohort_dir/data/<submitter>/simple_nucleotide_variation/*.maf.gz,
    concatenate, join clinical lineage cols, and add convenience flags.
    """
    cohort_dir = Path(cohort_dir)
    mafs = sorted(cohort_dir.glob("data/*/simple_nucleotide_variation/*.maf.gz"))
    if not mafs:
        raise FileNotFoundError(
            f"no .maf.gz files under {cohort_dir}/data/*/simple_nucleotide_variation/"
        )

    frames: list[pd.DataFrame] = []
    for p in mafs:
        df = read_maf(p)
        df["source_file"] = p.name
        frames.append(df)
    variants = pd.concat(frames, ignore_index=True)

    # Join in project_id + primary_diagnosis from clinical.parquet
    clin_path = cohort_dir / "clinical.parquet"
    if clin_path.exists():
        clin = pd.read_parquet(clin_path)
        keep = ["case_id", "submitter_id", "project_id", "diagnosis_primary_diagnosis"]
        clin = clin[[c for c in keep if c in clin.columns]].rename(
            columns={"diagnosis_primary_diagnosis": "primary_diagnosis"}
        )
        variants = variants.merge(clin, on="case_id", how="left")

    variants = _add_flags(variants)
    variants = _mark_primary_aliquot(variants)

    front = [c for c in OUTPUT_COLUMN_ORDER if c in variants.columns]
    rest = [c for c in variants.columns if c not in front]
    ordered: pd.DataFrame = variants[front + rest]
    return ordered


def write_variants(cohort_dir: Path) -> Path:
    df = aggregate_cohort(cohort_dir)
    out = Path(cohort_dir) / "variants.parquet"
    df.to_parquet(out, index=False)
    return out
