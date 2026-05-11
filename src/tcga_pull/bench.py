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
    """Run variants + samples through pandas then polars; diff outputs.

    Bench artefacts go to `bench/` inside the cohort dir, so the canonical
    `variants.parquet` / `samples.parquet` are restored from their pre-bench
    state on success. Downstream commands (frequency, etc.) keep working.
    """
    cohort_dir = Path(cohort_dir)
    stages: list[StageResult] = []
    bench_dir = cohort_dir / "bench"
    bench_dir.mkdir(exist_ok=True)

    # Snapshot anything that already exists so we can put it back
    variants_canonical = cohort_dir / "variants.parquet"
    samples_canonical = cohort_dir / "samples.parquet"
    variants_was = variants_canonical.read_bytes() if variants_canonical.exists() else None
    samples_was = samples_canonical.read_bytes() if samples_canonical.exists() else None

    try:
        # ---- variants: pandas
        _, dt_pd_v = _time_call(lambda: variants_pandas.write_variants(cohort_dir))
        pd_v_kept = bench_dir / "variants.pandas.parquet"
        variants_canonical.replace(pd_v_kept)
        rows, cols, bytes_ = _summary(pd_v_kept)
        stages.append(
            StageResult("pandas", "variants", dt_pd_v, _peak_rss_mb(), rows, cols, bytes_)
        )

        # ---- variants: polars
        _, dt_pl_v = _time_call(lambda: variants_polars.write_variants(cohort_dir))
        pl_v_kept = bench_dir / "variants.polars.parquet"
        variants_canonical.replace(pl_v_kept)
        rows, cols, bytes_ = _summary(pl_v_kept)
        stages.append(
            StageResult("polars", "variants", dt_pl_v, _peak_rss_mb(), rows, cols, bytes_)
        )

        variants_diff = _diff_parquets(
            pd_v_kept, pl_v_kept, sort_cols=["submitter_id", "tumor_barcode", "chrom", "pos"]
        )

        # ---- samples: needs variants.parquet; use the pandas variants for both passes
        # for a fair comparison (samples math is the only thing varying here).
        pd_v_kept.replace(variants_canonical)

        _, dt_pd_s = _time_call(lambda: samples_pandas.write_samples(cohort_dir))
        pd_s_kept = bench_dir / "samples.pandas.parquet"
        samples_canonical.replace(pd_s_kept)
        rows, cols, bytes_ = _summary(pd_s_kept)
        stages.append(StageResult("pandas", "samples", dt_pd_s, _peak_rss_mb(), rows, cols, bytes_))

        _, dt_pl_s = _time_call(lambda: samples_polars.write_samples(cohort_dir))
        pl_s_kept = bench_dir / "samples.polars.parquet"
        samples_canonical.replace(pl_s_kept)
        rows, cols, bytes_ = _summary(pl_s_kept)
        stages.append(StageResult("polars", "samples", dt_pl_s, _peak_rss_mb(), rows, cols, bytes_))

        samples_diff = _diff_parquets(pd_s_kept, pl_s_kept, sort_cols=["submitter_id"])

        # Move the canonical artefacts back to bench-side names so we can restore
        # the snapshot.
        variants_canonical.replace(bench_dir / "_variants.last.parquet")
        # samples_canonical was just moved to pl_s_kept above; nothing to clear.
    finally:
        # Always restore the canonical filenames
        if variants_was is not None:
            variants_canonical.write_bytes(variants_was)
        if samples_was is not None:
            samples_canonical.write_bytes(samples_was)
        # Clean up the throwaway last-variants copy
        last = bench_dir / "_variants.last.parquet"
        if last.exists():
            last.unlink()

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
