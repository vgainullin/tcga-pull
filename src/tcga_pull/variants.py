"""Walk a cohort produced by `tcga-pull pull` (Simple Nucleotide Variation MAFs)
and produce a tidy per-variant parquet:

    project_id, case_id, submitter_id, primary_diagnosis,
    hugo_symbol, chrom, pos, ref, alt,
    variant_class, consequence, impact, hgvsp_short,
    t_depth, t_alt_count, vaf, callers, gnomad_af

vaf = t_alt_count / t_depth (mutant allele frequency in tumor).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

# Source columns we want from the MAF, in source order
MAF_KEEP: list[str] = [
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "Variant_Classification",
    "Consequence",
    "IMPACT",
    "HGVSp_Short",
    "t_depth",
    "t_alt_count",
    "callers",
    "gnomAD_AF",
    "case_id",
    "Tumor_Sample_Barcode",
    "Matched_Norm_Sample_Barcode",
]

# Renamed in the output for tidy snake_case
MAF_RENAME: dict[str, str] = {
    "Hugo_Symbol": "hugo_symbol",
    "Chromosome": "chrom",
    "Start_Position": "pos",
    "Reference_Allele": "ref",
    "Tumor_Seq_Allele2": "alt",
    "Variant_Classification": "variant_class",
    "Consequence": "consequence",
    "IMPACT": "impact",
    "HGVSp_Short": "hgvsp_short",
    "Tumor_Sample_Barcode": "tumor_barcode",
    "Matched_Norm_Sample_Barcode": "normal_barcode",
    "gnomAD_AF": "gnomad_af",
}


def read_maf(path: Path) -> pd.DataFrame:
    """Read one GDC MAF (gzipped, tab-separated, ## comment headers). Return projected df."""
    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        low_memory=False,
        dtype={
            "Chromosome": str,
            "Start_Position": "Int64",
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


def aggregate_cohort(cohort_dir: Path) -> pd.DataFrame:
    """Walk cohort_dir/data/<submitter>/simple_nucleotide_variation/*.maf.gz,
    concatenate, and join project_id / primary_diagnosis from clinical.parquet."""
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

    # column order
    front = [
        "project_id",
        "case_id",
        "submitter_id",
        "primary_diagnosis",
        "hugo_symbol",
        "chrom",
        "pos",
        "ref",
        "alt",
        "variant_class",
        "consequence",
        "impact",
        "hgvsp_short",
        "t_depth",
        "t_alt_count",
        "vaf",
        "callers",
        "gnomad_af",
    ]
    front = [c for c in front if c in variants.columns]
    rest = [c for c in variants.columns if c not in front]
    ordered: pd.DataFrame = variants[front + rest]
    return ordered


def write_variants(cohort_dir: Path) -> Path:
    df = aggregate_cohort(cohort_dir)
    out = Path(cohort_dir) / "variants.parquet"
    df.to_parquet(out, index=False)
    return out
