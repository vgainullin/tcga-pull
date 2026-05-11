"""Offline tests for the MAF validator. Builds tiny synthetic MAFs on disk."""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

from tcga_pull.validate import (
    REQUIRED_OUTPUT_COLUMNS,
    _program_from_barcode,
    find_mafs,
    validate_mafs,
)


def _write_maf(path: Path, rows: list[dict], extra_cols: list[str] | None = None) -> None:
    """Synthesize a minimal GDC-style MAF.gz with a leading comment header.

    Always includes the columns required by `read_maf`'s projection so the
    downstream parsing path is exercised end-to-end.
    """
    base_cols = [
        "Hugo_Symbol",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Variant_Type",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
        "Variant_Classification",
        "Consequence",
        "IMPACT",
        "t_depth",
        "t_alt_count",
        "Tumor_Sample_Barcode",
        "Matched_Norm_Sample_Barcode",
        "case_id",
        "Transcript_ID",
        "EXON",
        "HGVSc",
        "HGVSp_Short",
        "callers",
        "gnomAD_AF",
        "COSMIC",
        "SIFT",
        "PolyPhen",
        "CONTEXT",
        "hotspot",
    ]
    cols = base_cols + (extra_cols or [])
    df = pd.DataFrame(rows, columns=cols)
    body = "#version gdc-1.0.0\n" + df.to_csv(sep="\t", index=False)
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wb") as fh:
        fh.write(body.encode("utf-8"))


def _row(
    *,
    chrom: str = "chr1",
    pos: int = 1000,
    hugo: str = "GENE1",
    case: str = "case-uuid-1",
    tumor: str = "TCGA-AA-AAAA-01A-11D-0000-09",
    normal: str = "TCGA-AA-AAAA-10A-01D-0000-09",
    variant_class: str = "Missense_Mutation",
    callers: str = "muse;mutect2;varscan2",
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
        "IMPACT": "MODERATE",
        "t_depth": 100,
        "t_alt_count": 25,
        "Tumor_Sample_Barcode": tumor,
        "Matched_Norm_Sample_Barcode": normal,
        "case_id": case,
        "Transcript_ID": "ENST00000000001",
        "EXON": "5/12",
        "HGVSc": "c.100C>T",
        "HGVSp_Short": "p.P34L",
        "callers": callers,
        "gnomAD_AF": None,
        "COSMIC": None,
        "SIFT": None,
        "PolyPhen": None,
        "CONTEXT": "ACGTACGT",
        "hotspot": "N",
    }


def test_program_from_barcode():
    assert _program_from_barcode("TCGA-AA-AAAA-01A-11D-A123-09") == "TCGA"
    assert _program_from_barcode("TARGET-20-PADDXA-09A-01D") == "TARGET"
    assert _program_from_barcode("CPTAC-C3-00001") == "CPTAC"
    assert _program_from_barcode("") == "unknown"
    assert _program_from_barcode(None) == "unknown"
    assert _program_from_barcode("no-hyphen-at-start") == "no"


def test_find_mafs_recurses_any_layout(tmp_path: Path):
    # mimic gdc-client _downloads layout
    _write_maf(tmp_path / "_downloads" / "abc-123" / "a.maf.gz", [_row()])
    # mimic restructured cohort layout
    _write_maf(
        tmp_path / "data" / "TCGA-XX-0001" / "simple_nucleotide_variation" / "b.maf.gz",
        [_row()],
    )
    # noise: parquet file shouldn't be picked up
    (tmp_path / "x.parquet").write_bytes(b"PAR1")
    found = find_mafs(tmp_path)
    assert len(found) == 2
    assert all(p.suffix == ".gz" for p in found)


def test_validate_mafs_clean_run(tmp_path: Path):
    paths = []
    for i in range(3):
        p = tmp_path / f"f{i}.maf.gz"
        _write_maf(p, [_row(pos=1000 + j, case=f"case{i}") for j in range(4)])
        paths.append(p)
    report = validate_mafs(paths)
    assert report.files_total == 3
    assert report.files_ok == 3
    assert report.files_failed == 0
    assert report.rows_total == 12
    assert report.distinct_case_ids == 3
    assert report.distinct_tumor_barcodes == 1
    assert report.required_columns_complete
    assert report.missing_required_columns == []
    assert report.rows_by_program == {"TCGA": 12}
    assert report.rows_with_unknown_barcode == 0


def test_validate_mafs_collects_parse_errors(tmp_path: Path):
    good = tmp_path / "good.maf.gz"
    _write_maf(good, [_row()])
    bad = tmp_path / "bad.maf.gz"
    # invalid gzip
    bad.write_bytes(b"NOT GZIP DATA")
    report = validate_mafs([good, bad])
    assert report.files_ok == 1
    assert report.files_failed == 1
    assert len(report.parse_errors) == 1
    name, etype, _msg = report.parse_errors[0]
    assert name == "bad.maf.gz"
    assert "Error" in etype or "BadGzip" in etype


def test_validate_mafs_flags_non_tcga_barcodes(tmp_path: Path):
    p = tmp_path / "mixed.maf.gz"
    rows = [
        _row(tumor="TCGA-A1-A0SB-01A-11D-A099-09"),
        _row(tumor="C3L-00001-01"),  # CPTAC-style
        _row(tumor=""),  # missing barcode
    ]
    _write_maf(p, rows)
    report = validate_mafs([p])
    assert report.rows_by_program == {"TCGA": 1, "C3L": 1, "unknown": 1}
    assert report.rows_with_unknown_barcode == 1


def test_required_columns_constant_aligns_with_read_maf_output():
    # The validator's REQUIRED set should be a subset of what read_maf surfaces.
    # If we add a required column that read_maf doesn't produce, this fails fast.
    from tcga_pull.variants import MAF_KEEP, MAF_RENAME

    # Source columns are either renamed (MAF_RENAME) or pass through unchanged.
    output_cols = set(MAF_RENAME.values()) | (set(MAF_KEEP) - set(MAF_RENAME))
    missing = [c for c in REQUIRED_OUTPUT_COLUMNS if c not in output_cols]
    assert not missing, f"REQUIRED columns not produced by read_maf: {missing}"
