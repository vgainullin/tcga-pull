"""Side-by-side benchmark: pandas vs polars on the variants + samples pipelines.

Runs each engine end-to-end on a cohort dir, captures wall time + peak
resident memory, then diffs the parquet outputs cell-by-cell so we can see
whether they actually agree on real data.
"""

from __future__ import annotations

import json
import resource
import time
from collections.abc import Callable
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import polars as pl

from . import samples as samples_pandas
from . import samples_polars, variants_polars
from . import variants as variants_pandas


def _peak_rss_mb() -> float:
    """Return peak RSS in MB for this process. On macOS getrusage returns bytes;
    on Linux it returns KB. Handle both."""
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # macOS reports in bytes; Linux in kilobytes
    if rss > 1 << 30:  # > 1 GB sample, assume bytes
        return rss / (1 << 20)
    return rss / 1024.0


@dataclass
class StageResult:
    engine: str
    stage: str
    wall_seconds: float
    peak_rss_mb_after: float
    rows: int
    cols: int
    parquet_bytes: int


@dataclass
class BenchResult:
    cohort_dir: str
    stages: list[StageResult]
    diffs: dict[str, dict[str, int]]  # {"variants": {"col_a": n_diff, ...}, "samples": {...}}


def _time_call(fn: Callable[[], Path]) -> tuple[Path, float]:
    t0 = time.perf_counter()
    path = fn()
    return path, time.perf_counter() - t0


def _summary(path: Path) -> tuple[int, int, int]:
    df = pl.read_parquet(path)
    return len(df), df.width, path.stat().st_size


def _diff_parquets(pd_path: Path, pl_path: Path, sort_cols: list[str]) -> dict[str, int]:
    """Return {col: n_diffs} for cells that aren't equal (treating nulls as equal)."""
    a = pl.read_parquet(pd_path).sort(sort_cols)
    b = pl.read_parquet(pl_path).sort(sort_cols)

    if a.shape != b.shape:
        return {"_shape": 1, "pandas_rows": a.height, "polars_rows": b.height}
    if set(a.columns) != set(b.columns):
        return {"_columns_differ": 1}

    a = a.select(sorted(a.columns))
    b = b.select(sorted(b.columns))

    out: dict[str, int] = {}
    for col in a.columns:
        ca, cb = a[col], b[col]
        both_null = ca.is_null() & cb.is_null()
        if ca.dtype.is_numeric() and cb.dtype.is_numeric():
            ca_f = ca.cast(pl.Float64, strict=False)
            cb_f = cb.cast(pl.Float64, strict=False)
            close = ((ca_f - cb_f).abs() <= 1e-9 * cb_f.abs().fill_null(1.0)).fill_null(False)
            eq = both_null | close
        else:
            ca_s = ca.cast(pl.Utf8, strict=False)
            cb_s = cb.cast(pl.Utf8, strict=False)
            eq = both_null | (ca_s == cb_s).fill_null(False)
        n_diff = int((~eq).sum())
        if n_diff:
            out[col] = n_diff
    return out


def run_bench(cohort_dir: Path) -> BenchResult:
    """Run variants + samples through pandas then polars; diff outputs."""
    cohort_dir = Path(cohort_dir)
    stages: list[StageResult] = []

    # ---- pandas: variants
    out_pd_variants, dt = _time_call(lambda: variants_pandas.write_variants(cohort_dir))
    pd_v_kept = out_pd_variants.with_suffix(".pandas.parquet")
    out_pd_variants.rename(pd_v_kept)
    rows, cols, bytes_ = _summary(pd_v_kept)
    stages.append(StageResult("pandas", "variants", dt, _peak_rss_mb(), rows, cols, bytes_))

    # ---- polars: variants
    out_pl_variants, dt = _time_call(lambda: variants_polars.write_variants(cohort_dir))
    pl_v_kept = out_pl_variants.with_suffix(".polars.parquet")
    out_pl_variants.rename(pl_v_kept)
    rows, cols, bytes_ = _summary(pl_v_kept)
    stages.append(StageResult("polars", "variants", dt, _peak_rss_mb(), rows, cols, bytes_))

    variants_diff = _diff_parquets(
        pd_v_kept, pl_v_kept, sort_cols=["submitter_id", "tumor_barcode", "chrom", "pos"]
    )

    # samples needs variants.parquet; restore from the pandas copy and rerun each
    pd_v_kept.replace(cohort_dir / "variants.parquet")

    out_pd_samples, dt = _time_call(lambda: samples_pandas.write_samples(cohort_dir))
    pd_s_kept = out_pd_samples.with_suffix(".pandas.parquet")
    out_pd_samples.rename(pd_s_kept)
    rows, cols, bytes_ = _summary(pd_s_kept)
    stages.append(StageResult("pandas", "samples", dt, _peak_rss_mb(), rows, cols, bytes_))

    # swap in polars variants for the polars samples pass
    (cohort_dir / "variants.parquet").replace(pd_v_kept)  # rename back to .pandas.parquet
    pl_v_kept.replace(cohort_dir / "variants.parquet")

    out_pl_samples, dt = _time_call(lambda: samples_polars.write_samples(cohort_dir))
    pl_s_kept = out_pl_samples.with_suffix(".polars.parquet")
    out_pl_samples.rename(pl_s_kept)
    rows, cols, bytes_ = _summary(pl_s_kept)
    stages.append(StageResult("polars", "samples", dt, _peak_rss_mb(), rows, cols, bytes_))

    samples_diff = _diff_parquets(pd_s_kept, pl_s_kept, sort_cols=["submitter_id"])

    # leave variants.parquet pointing at the polars version; user can pick
    return BenchResult(
        cohort_dir=str(cohort_dir),
        stages=stages,
        diffs={"variants": variants_diff, "samples": samples_diff},
    )


def to_dict(result: BenchResult) -> dict[str, Any]:
    return {
        "cohort_dir": result.cohort_dir,
        "stages": [asdict(s) for s in result.stages],
        "diffs": result.diffs,
    }


def write_json(result: BenchResult, path: Path) -> Path:
    path.write_text(json.dumps(to_dict(result), indent=2))
    return path
