"""Validate a collection of MAFs against our expected schema.

Useful for:
  * QC of a freshly downloaded cohort (catch schema drift before parquet write).
  * Sanity-checking the `_downloads/` staging area after a partial pull.
  * Stress-testing `read_maf` against new data programs (TARGET, CPTAC, MMRF, ...).

The validator never raises on a single bad MAF — it collects failures so you
see the full picture of what's wrong in one pass.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from .variants import MAF_RENAME, OUTPUT_COLUMN_ORDER, read_maf

if TYPE_CHECKING:
    from rich.console import Console

# Columns we expect read_maf to surface for every MAF (after rename to snake_case).
# These are the ones load-bearing for the downstream pipeline.
REQUIRED_OUTPUT_COLUMNS: tuple[str, ...] = (
    "chrom",
    "pos",
    "ref",
    "alt",
    "hugo_symbol",
    "variant_class",
    "consequence",
    "impact",
    "t_depth",
    "t_alt_count",
    "tumor_barcode",
    "normal_barcode",
    "case_id",
)


@dataclass
class ValidationReport:
    files_total: int = 0
    files_ok: int = 0
    files_failed: int = 0
    rows_total: int = 0

    # one entry per failed file: (path_name, exception_type, message)
    parse_errors: list[tuple[str, str, str]] = field(default_factory=list)

    # column -> number of files in which it appeared (post-rename)
    column_presence: dict[str, int] = field(default_factory=dict)

    # rows grouped by inferred program (from tumor_barcode prefix, "TCGA-XX-..." → "TCGA")
    rows_by_program: dict[str, int] = field(default_factory=dict)

    # rows whose tumor_barcode does not match a recognisable program prefix
    rows_with_unknown_barcode: int = 0

    distinct_case_ids: int = 0
    distinct_tumor_barcodes: int = 0
    distinct_normal_barcodes: int = 0

    @property
    def required_columns_complete(self) -> bool:
        return all(self.column_presence.get(c, 0) == self.files_ok for c in REQUIRED_OUTPUT_COLUMNS)

    @property
    def missing_required_columns(self) -> list[str]:
        """Columns we expected but that weren't present in every OK file."""
        return [
            c for c in REQUIRED_OUTPUT_COLUMNS if self.column_presence.get(c, 0) != self.files_ok
        ]


def _program_from_barcode(barcode: object) -> str:
    """Extract the program prefix (text before the first '-') from a barcode.
    Returns 'unknown' for empty/non-string/no-hyphen values."""
    if not isinstance(barcode, str) or "-" not in barcode:
        return "unknown"
    prefix = barcode.split("-", 1)[0]
    return prefix or "unknown"


def validate_mafs(paths: list[Path]) -> ValidationReport:
    """Read each MAF; collect parse errors, column coverage, row totals."""
    report = ValidationReport(files_total=len(paths))
    frames: list[pd.DataFrame] = []
    for p in paths:
        try:
            df = read_maf(p)
        except Exception as e:
            report.files_failed += 1
            report.parse_errors.append((p.name, type(e).__name__, str(e)[:200]))
            continue
        report.files_ok += 1
        for col in df.columns:
            report.column_presence[col] = report.column_presence.get(col, 0) + 1
        frames.append(df)

    if frames:
        all_rows = pd.concat(frames, ignore_index=True)
        report.rows_total = len(all_rows)
        if "case_id" in all_rows.columns:
            report.distinct_case_ids = int(all_rows["case_id"].nunique(dropna=True))
        if "tumor_barcode" in all_rows.columns:
            report.distinct_tumor_barcodes = int(all_rows["tumor_barcode"].nunique(dropna=True))
            programs = all_rows["tumor_barcode"].map(_program_from_barcode)
            counts = programs.value_counts(dropna=False).to_dict()
            report.rows_by_program = {str(k): int(v) for k, v in counts.items()}
            report.rows_with_unknown_barcode = int(counts.get("unknown", 0))
        if "normal_barcode" in all_rows.columns:
            report.distinct_normal_barcodes = int(all_rows["normal_barcode"].nunique(dropna=True))

    return report


def render_report(report: ValidationReport, console: Console | None = None) -> None:
    """Pretty-print a ValidationReport to the console."""
    from rich.console import Console as _Console
    from rich.table import Table

    console = console or _Console()
    summary = Table(title="MAF validation", show_header=False)
    summary.add_column("k", style="dim")
    summary.add_column("v")
    summary.add_row("Files scanned", f"{report.files_total:,}")
    summary.add_row("Files OK", f"[green]{report.files_ok:,}[/green]")
    failed_style = "red" if report.files_failed else "dim"
    summary.add_row("Files failed", f"[{failed_style}]{report.files_failed:,}[/{failed_style}]")
    summary.add_row("Variant rows", f"{report.rows_total:,}")
    summary.add_row("Distinct case_ids", f"{report.distinct_case_ids:,}")
    summary.add_row("Distinct tumor barcodes", f"{report.distinct_tumor_barcodes:,}")
    summary.add_row("Distinct normal barcodes", f"{report.distinct_normal_barcodes:,}")
    if report.rows_with_unknown_barcode:
        summary.add_row(
            "Rows with non-standard barcode",
            f"[yellow]{report.rows_with_unknown_barcode:,}[/yellow]",
        )
    console.print(summary)

    if report.parse_errors:
        err_tbl = Table(title=f"Parse errors ({len(report.parse_errors)})")
        err_tbl.add_column("file")
        err_tbl.add_column("error type")
        err_tbl.add_column("message", overflow="fold")
        for name, etype, msg in report.parse_errors[:15]:
            err_tbl.add_row(name, etype, msg)
        if len(report.parse_errors) > 15:
            err_tbl.add_row("...", f"+{len(report.parse_errors) - 15} more", "")
        console.print(err_tbl)

    if report.missing_required_columns:
        miss_tbl = Table(title="Missing required columns")
        miss_tbl.add_column("column")
        miss_tbl.add_column("present_in", justify="right")
        miss_tbl.add_column("of", justify="right")
        for c in report.missing_required_columns:
            miss_tbl.add_row(c, str(report.column_presence.get(c, 0)), str(report.files_ok))
        console.print(miss_tbl)
    elif report.files_ok:
        console.print(
            f"[green]all {len(REQUIRED_OUTPUT_COLUMNS)} required columns present "
            f"in every OK file[/green]"
        )

    # Column-completeness (any non-required column that's not in every file)
    optional_missing = sorted(
        [
            (c, n)
            for c, n in report.column_presence.items()
            if c in MAF_RENAME.values() or c in OUTPUT_COLUMN_ORDER
            if n != report.files_ok
            if c not in REQUIRED_OUTPUT_COLUMNS
        ],
        key=lambda x: x[1],
    )
    if optional_missing:
        tbl = Table(title="Optional columns with partial coverage")
        tbl.add_column("column")
        tbl.add_column("present_in", justify="right")
        tbl.add_column("of", justify="right")
        for c, n in optional_missing:
            tbl.add_row(c, str(n), str(report.files_ok))
        console.print(tbl)

    if report.rows_by_program:
        prog_tbl = Table(title="Variant rows by inferred program")
        prog_tbl.add_column("program")
        prog_tbl.add_column("rows", justify="right")
        for k, v in sorted(report.rows_by_program.items(), key=lambda x: -x[1])[:20]:
            prog_tbl.add_row(k, f"{v:,}")
        console.print(prog_tbl)


def find_mafs(path: Path) -> list[Path]:
    """Recursively find *.maf.gz under `path`. Works on either layout:
    flat `_downloads/<file_id>/<file>.maf.gz` or restructured
    `data/<submitter>/simple_nucleotide_variation/<file>.maf.gz`."""
    return sorted(Path(path).rglob("*.maf.gz"))
