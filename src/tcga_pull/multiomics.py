"""Processing recipes for non-SNV GDC omics tables.

The writers in this module are intentionally conservative: they process one
downloaded file at a time from manifest.parquet, project common GDC tabular
columns into typed long-form parquet outputs, and tolerate absent file types by
writing an empty parquet with the expected schema.
"""

from __future__ import annotations

import shutil
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Any, cast

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

COMMON_FIELDS: tuple[pa.Field, ...] = (
    pa.field("case_id", pa.string()),
    pa.field("submitter_id", pa.string()),
    pa.field("file_id", pa.string()),
    pa.field("file_name", pa.string()),
    pa.field("data_type", pa.string()),
    pa.field("experimental_strategy", pa.string()),
    pa.field("workflow_type", pa.string()),
)

RNA_SCHEMA = pa.schema(
    [
        *COMMON_FIELDS,
        pa.field("gene_id", pa.string()),
        pa.field("gene_name", pa.string()),
        pa.field("gene_type", pa.string()),
        pa.field("unstranded", pa.int64()),
        pa.field("stranded_first", pa.int64()),
        pa.field("stranded_second", pa.int64()),
        pa.field("tpm_unstranded", pa.float64()),
        pa.field("fpkm_unstranded", pa.float64()),
        pa.field("fpkm_uq_unstranded", pa.float64()),
    ]
)

MIRNA_SCHEMA = pa.schema(
    [
        *COMMON_FIELDS,
        pa.field("mirna_id", pa.string()),
        pa.field("read_count", pa.int64()),
        pa.field("reads_per_million_mirna_mapped", pa.float64()),
        pa.field("cross_mapped", pa.string()),
    ]
)

METHYLATION_SCHEMA = pa.schema(
    [
        *COMMON_FIELDS,
        pa.field("probe_id", pa.string()),
        pa.field("beta_value", pa.float64()),
    ]
)

CNV_SEGMENT_SCHEMA = pa.schema(
    [
        *COMMON_FIELDS,
        pa.field("cnv_type", pa.string()),
        pa.field("sample", pa.string()),
        pa.field("chrom", pa.string()),
        pa.field("start", pa.int64()),
        pa.field("end", pa.int64()),
        pa.field("num_probes", pa.int64()),
        pa.field("segment_mean", pa.float64()),
    ]
)

GENE_CNV_SCHEMA = pa.schema(
    [
        *COMMON_FIELDS,
        pa.field("gene_id", pa.string()),
        pa.field("gene_name", pa.string()),
        pa.field("chrom", pa.string()),
        pa.field("start", pa.int64()),
        pa.field("end", pa.int64()),
        pa.field("copy_number", pa.float64()),
    ]
)

PROTEIN_SCHEMA = pa.schema(
    [
        *COMMON_FIELDS,
        pa.field("protein_id", pa.string()),
        pa.field("gene_symbol", pa.string()),
        pa.field("antibody", pa.string()),
        pa.field("expression_value", pa.float64()),
    ]
)

MULTIOMICS_RECIPE_NAMES = {
    "rna_expression",
    "mirna_expression",
    "methylation",
    "copy_number",
    "protein_expression",
    "multiomics",
}


def write_rna_expression(cohort_dir: Path) -> Path:
    cohort_dir = Path(cohort_dir)
    rows = _manifest_rows(
        cohort_dir,
        data_category="Transcriptome Profiling",
        data_type="Gene Expression Quantification",
        experimental_strategy="RNA-Seq",
    )
    out = cohort_dir / "rna_expression.parquet"
    _write_stream(out, (_rna_frame(row) for row in rows), RNA_SCHEMA)
    return out


def write_mirna_expression(cohort_dir: Path) -> Path:
    cohort_dir = Path(cohort_dir)
    rows = _manifest_rows(
        cohort_dir,
        data_category="Transcriptome Profiling",
        data_type="miRNA Expression Quantification",
        experimental_strategy="miRNA-Seq",
    )
    out = cohort_dir / "mirna_expression.parquet"
    _write_stream(out, (_mirna_frame(row) for row in rows), MIRNA_SCHEMA)
    return out


def write_methylation_beta(cohort_dir: Path) -> Path:
    cohort_dir = Path(cohort_dir)
    rows = _manifest_rows(
        cohort_dir,
        data_category="DNA Methylation",
        data_type="Methylation Beta Value",
    )
    out = cohort_dir / "methylation_beta.parquet"
    _write_stream(out, (_methylation_frame(row) for row in rows), METHYLATION_SCHEMA)
    return out


def write_copy_number(cohort_dir: Path) -> tuple[Path, Path]:
    cohort_dir = Path(cohort_dir)
    segment_rows = list(
        _manifest_rows(
            cohort_dir,
            data_category="Copy Number Variation",
            data_type=["Copy Number Segment", "Masked Copy Number Segment"],
        )
    )
    gene_rows = _manifest_rows(
        cohort_dir,
        data_category="Copy Number Variation",
        data_type="Gene Level Copy Number",
    )

    segments_out = cohort_dir / "copy_number_segments.parquet"
    gene_out = cohort_dir / "gene_copy_number.parquet"
    _write_stream(
        segments_out, (_cnv_segment_frame(row) for row in segment_rows), CNV_SEGMENT_SCHEMA
    )
    _write_stream(gene_out, (_gene_cnv_frame(row) for row in gene_rows), GENE_CNV_SCHEMA)
    return segments_out, gene_out


def write_protein_expression(cohort_dir: Path) -> Path:
    cohort_dir = Path(cohort_dir)
    rows = _manifest_rows(
        cohort_dir,
        data_category="Proteome Profiling",
        data_type="Protein Expression Quantification",
    )
    out = cohort_dir / "protein_expression.parquet"
    _write_stream(out, (_protein_frame(row) for row in rows), PROTEIN_SCHEMA)
    return out


def write_multiomics(cohort_dir: Path) -> list[Path]:
    cnv_segments, gene_cnv = write_copy_number(cohort_dir)
    return [
        write_rna_expression(cohort_dir),
        write_mirna_expression(cohort_dir),
        write_methylation_beta(cohort_dir),
        cnv_segments,
        gene_cnv,
        write_protein_expression(cohort_dir),
    ]


def write_multiomics_parts(
    cohort_dir: Path,
    records: list[dict[str, Any]],
    *,
    part_id: int,
    recipes: Iterable[str],
) -> list[Path]:
    """Write selected multiomics recipe outputs for one downloaded batch.

    Parts land under <cohort>/_parts/<output-name>/part-XXXXXX.parquet and are
    finalized into top-level parquet files after all batches finish.
    """
    selected = _expanded_recipes(recipes)
    out: list[Path] = []
    for output_name, schema, frames in _recipe_frame_iterators(records, selected):
        part_path = _parts_dir(cohort_dir, output_name) / f"part-{part_id:06d}.parquet"
        if _write_stream(part_path, frames, schema, write_empty=False):
            out.append(part_path)
    return out


def finalize_multiomics_parts(cohort_dir: Path, *, recipes: Iterable[str]) -> list[Path]:
    selected = _expanded_recipes(recipes)
    out: list[Path] = []
    for output_name, schema, _ in _recipe_frame_iterators([], selected):
        out_path = Path(cohort_dir) / f"{output_name}.parquet"
        _combine_parts(_parts_dir(cohort_dir, output_name), out_path, schema)
        out.append(out_path)
    return out


def record_handled_by_multiomics(record: dict[str, Any], recipes: Iterable[str]) -> bool:
    selected = _expanded_recipes(recipes)
    return any(_matches_output(record, output_name) for output_name in selected)


def _manifest_rows(
    cohort_dir: Path,
    *,
    data_category: str,
    data_type: str | Iterable[str],
    experimental_strategy: str | None = None,
) -> Iterator[dict[str, Any]]:
    manifest_path = cohort_dir / "manifest.parquet"
    if not manifest_path.exists():
        raise FileNotFoundError(f"missing {manifest_path} - run `tcga-pull pull` first")

    df = pd.read_parquet(manifest_path)
    if df.empty:
        return

    data_types = {data_type} if isinstance(data_type, str) else set(data_type)
    mask = pd.Series(True, index=df.index)
    if "status" in df.columns:
        mask &= df["status"] == "ok"
    mask &= df["local_path"].notna()
    mask &= df["data_category"] == data_category
    mask &= df["data_type"].isin(data_types)
    if experimental_strategy is not None and "experimental_strategy" in df.columns:
        mask &= df["experimental_strategy"] == experimental_strategy

    for row in df[mask].to_dict("records"):
        path = row.get("local_path")
        if path and Path(path).exists():
            yield row


def _expanded_recipes(recipes: Iterable[str]) -> set[str]:
    selected = set(recipes) & MULTIOMICS_RECIPE_NAMES
    if "multiomics" in selected:
        return {
            "rna_expression",
            "mirna_expression",
            "methylation_beta",
            "copy_number_segments",
            "gene_copy_number",
            "protein_expression",
        }
    out: set[str] = set()
    if "rna_expression" in selected:
        out.add("rna_expression")
    if "mirna_expression" in selected:
        out.add("mirna_expression")
    if "methylation" in selected:
        out.add("methylation_beta")
    if "copy_number" in selected:
        out.update({"copy_number_segments", "gene_copy_number"})
    if "protein_expression" in selected:
        out.add("protein_expression")
    return out


def _recipe_frame_iterators(
    records: list[dict[str, Any]],
    selected: set[str],
) -> Iterator[tuple[str, pa.Schema, Iterable[pd.DataFrame]]]:
    if "rna_expression" in selected:
        yield (
            "rna_expression",
            RNA_SCHEMA,
            (_rna_frame(row) for row in records if _matches_output(row, "rna_expression")),
        )
    if "mirna_expression" in selected:
        yield (
            "mirna_expression",
            MIRNA_SCHEMA,
            (_mirna_frame(row) for row in records if _matches_output(row, "mirna_expression")),
        )
    if "methylation_beta" in selected:
        yield (
            "methylation_beta",
            METHYLATION_SCHEMA,
            (
                _methylation_frame(row)
                for row in records
                if _matches_output(row, "methylation_beta")
            ),
        )
    if "copy_number_segments" in selected:
        yield (
            "copy_number_segments",
            CNV_SEGMENT_SCHEMA,
            (
                _cnv_segment_frame(row)
                for row in records
                if _matches_output(row, "copy_number_segments")
            ),
        )
    if "gene_copy_number" in selected:
        yield (
            "gene_copy_number",
            GENE_CNV_SCHEMA,
            (_gene_cnv_frame(row) for row in records if _matches_output(row, "gene_copy_number")),
        )
    if "protein_expression" in selected:
        yield (
            "protein_expression",
            PROTEIN_SCHEMA,
            (_protein_frame(row) for row in records if _matches_output(row, "protein_expression")),
        )


def _matches_output(record: dict[str, Any], output_name: str) -> bool:
    category = record.get("data_category")
    data_type = record.get("data_type")
    strategy = record.get("experimental_strategy")
    if output_name == "rna_expression":
        return (
            category == "Transcriptome Profiling"
            and data_type == "Gene Expression Quantification"
            and strategy == "RNA-Seq"
        )
    if output_name == "mirna_expression":
        return (
            category == "Transcriptome Profiling"
            and data_type == "miRNA Expression Quantification"
            and strategy == "miRNA-Seq"
        )
    if output_name == "methylation_beta":
        return category == "DNA Methylation" and data_type == "Methylation Beta Value"
    if output_name == "copy_number_segments":
        return category == "Copy Number Variation" and data_type in {
            "Copy Number Segment",
            "Masked Copy Number Segment",
        }
    if output_name == "gene_copy_number":
        return category == "Copy Number Variation" and data_type == "Gene Level Copy Number"
    if output_name == "protein_expression":
        return category == "Proteome Profiling" and data_type == "Protein Expression Quantification"
    return False


def _parts_dir(cohort_dir: Path, output_name: str) -> Path:
    return Path(cohort_dir) / "_parts" / output_name


def _read_table(path: Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep="\t", comment="#", dtype=str, low_memory=False)
    except pd.errors.ParserError:
        df = pd.read_csv(path, sep="\t", comment="#", dtype=str, engine="python")
    df = df.dropna(axis=1, how="all")
    df.columns = [_normalize_col(c) for c in df.columns]
    return cast(pd.DataFrame, df)


def _normalize_col(value: object) -> str:
    out = str(value).strip().lower()
    for old, new in ((" ", "_"), ("-", "_"), (".", "_"), ("/", "_"), ("(", ""), (")", "")):
        out = out.replace(old, new)
    while "__" in out:
        out = out.replace("__", "_")
    return out.strip("_")


def _rna_frame(row: dict[str, Any]) -> pd.DataFrame:
    df = _read_table(Path(row["local_path"]))
    if "gene_id" not in df.columns and len(df.columns) > 0:
        df = df.rename(columns={df.columns[0]: "gene_id"})
    if "gene_id" not in df.columns:
        return _empty_df(RNA_SCHEMA)
    df = df[df["gene_id"].astype(str).str.startswith("ENSG", na=False)].copy()
    out = _with_common(df, row)
    out["gene_id"] = _text(df, "gene_id")
    out["gene_name"] = _text(df, "gene_name", "gene")
    out["gene_type"] = _text(df, "gene_type")
    for col in ("unstranded", "stranded_first", "stranded_second"):
        out[col] = _int(df, col)
    for col in ("tpm_unstranded", "fpkm_unstranded", "fpkm_uq_unstranded"):
        out[col] = _float(df, col)
    return _select_schema(out, RNA_SCHEMA)


def _mirna_frame(row: dict[str, Any]) -> pd.DataFrame:
    df = _read_table(Path(row["local_path"]))
    out = _with_common(df, row)
    out["mirna_id"] = _text(df, "mirna_id", "mirna", "mi_rna_id")
    out["read_count"] = _int(df, "read_count", "reads")
    out["reads_per_million_mirna_mapped"] = _float(
        df,
        "reads_per_million_mirna_mapped",
        "reads_per_million_mi_rna_mapped",
        "rpm",
    )
    out["cross_mapped"] = _text(df, "cross_mapped", "cross_mapped_mirna")
    return _select_schema(out, MIRNA_SCHEMA)


def _methylation_frame(row: dict[str, Any]) -> pd.DataFrame:
    df = _read_table(Path(row["local_path"]))
    out = _with_common(df, row)
    out["probe_id"] = _text(df, "probe_id", "composite_element_ref", "id_ref", "name")
    out["beta_value"] = _float(df, "beta_value", "beta", "value")
    return _select_schema(out, METHYLATION_SCHEMA)


def _cnv_segment_frame(row: dict[str, Any]) -> pd.DataFrame:
    df = _read_table(Path(row["local_path"]))
    out = _with_common(df, row)
    out["cnv_type"] = row.get("data_type")
    out["sample"] = _text(df, "sample", "sample_id", "sample_name")
    out["chrom"] = _text(df, "chrom", "chromosome", "chr")
    out["start"] = _int(df, "start", "start_position", "loc_start")
    out["end"] = _int(df, "end", "end_position", "loc_end")
    out["num_probes"] = _int(df, "num_probes", "num_mark", "num_markers", "num_marked")
    out["segment_mean"] = _float(df, "segment_mean", "seg_mean", "segment_mean_value")
    return _select_schema(out, CNV_SEGMENT_SCHEMA)


def _gene_cnv_frame(row: dict[str, Any]) -> pd.DataFrame:
    df = _read_table(Path(row["local_path"]))
    out = _with_common(df, row)
    out["gene_id"] = _text(df, "gene_id", "ensembl_gene_id")
    out["gene_name"] = _text(df, "gene_name", "gene_symbol", "symbol", "gene")
    out["chrom"] = _text(df, "chrom", "chromosome", "chr")
    out["start"] = _int(df, "start", "start_position", "loc_start")
    out["end"] = _int(df, "end", "end_position", "loc_end")
    out["copy_number"] = _float(
        df,
        "copy_number",
        "copy_number_value",
        "copy_number_estimate",
        "value",
        "log2_copy_ratio",
    )
    return _select_schema(out, GENE_CNV_SCHEMA)


def _protein_frame(row: dict[str, Any]) -> pd.DataFrame:
    df = _read_table(Path(row["local_path"]))
    out = _with_common(df, row)
    out["protein_id"] = _text(df, "protein_id", "id_ref", "protein")
    out["gene_symbol"] = _text(df, "gene_symbol", "gene_name", "gene", "symbol")
    out["antibody"] = _text(df, "antibody", "antibody_name", "peptide_target")
    out["expression_value"] = _float(
        df,
        "expression_value",
        "protein_expression",
        "normalized_expression",
        "value",
        "rppa_value",
    )
    return _select_schema(out, PROTEIN_SCHEMA)


def _with_common(df: pd.DataFrame, row: dict[str, Any]) -> pd.DataFrame:
    out = pd.DataFrame(index=df.index)
    for col in (
        "case_id",
        "submitter_id",
        "file_id",
        "file_name",
        "data_type",
        "experimental_strategy",
        "workflow_type",
    ):
        out[col] = row.get(col)
    return out


def _text(df: pd.DataFrame, *aliases: str) -> pd.Series:
    col = _first_col(df, aliases)
    if col is None:
        return pd.Series([None] * len(df), index=df.index, dtype="object")
    return cast(pd.Series, df[col].where(df[col].notna(), None))


def _int(df: pd.DataFrame, *aliases: str) -> pd.Series:
    col = _first_col(df, aliases)
    if col is None:
        return pd.Series([pd.NA] * len(df), index=df.index, dtype="Int64")
    return pd.to_numeric(df[col], errors="coerce").round().astype("Int64")


def _float(df: pd.DataFrame, *aliases: str) -> pd.Series:
    col = _first_col(df, aliases)
    if col is None:
        return pd.Series([pd.NA] * len(df), index=df.index, dtype="Float64")
    return pd.to_numeric(df[col], errors="coerce").astype("Float64")


def _first_col(df: pd.DataFrame, aliases: Iterable[str]) -> str | None:
    for alias in aliases:
        normalized = _normalize_col(alias)
        if normalized in df.columns:
            return normalized
    return None


def _select_schema(df: pd.DataFrame, schema: pa.Schema) -> pd.DataFrame:
    for field in schema:
        if field.name not in df.columns:
            df[field.name] = pd.NA
    return cast(pd.DataFrame, df[[field.name for field in schema]])


def _empty_df(schema: pa.Schema) -> pd.DataFrame:
    return pd.DataFrame({field.name: pd.Series(dtype="object") for field in schema})


def _combine_parts(parts_dir: Path, out_path: Path, schema: pa.Schema) -> None:
    part_paths = sorted(parts_dir.glob("part-*.parquet")) if parts_dir.exists() else []
    if not part_paths:
        _write_stream(out_path, [], schema)
        return

    if out_path.exists():
        if out_path.is_dir():
            shutil.rmtree(out_path)
        else:
            out_path.unlink()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    parts_dir.rename(out_path)


def _write_stream(
    path: Path,
    frames: Iterable[pd.DataFrame],
    schema: pa.Schema,
    *,
    write_empty: bool = True,
) -> bool:
    path.parent.mkdir(parents=True, exist_ok=True)
    writer: pq.ParquetWriter | None = None
    wrote_any = False
    try:
        for frame in frames:
            if frame.empty:
                continue
            table = pa.Table.from_pandas(frame, schema=schema, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(path, schema=schema)
            writer.write_table(table)
            wrote_any = True
        if not wrote_any and write_empty:
            empty = pa.Table.from_pandas(_empty_df(schema), schema=schema, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(path, schema=schema)
            writer.write_table(empty)
    finally:
        if writer is not None:
            writer.close()
    return wrote_any
