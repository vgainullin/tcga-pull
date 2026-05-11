"""Tests for the public load_cohort API."""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl
import pytest

from tcga_pull import Cohort, load_cohort


def _make_cohort_dir(tmp_path: Path, *, with_optional: bool = False) -> Path:
    cohort = tmp_path / "fixture"
    cohort.mkdir()

    pl.DataFrame({"case_id": ["c1", "c2"], "submitter_id": ["s1", "s2"]}).write_parquet(
        cohort / "clinical.parquet"
    )
    pl.DataFrame({"file_id": ["f1"], "case_id": ["c1"]}).write_parquet(cohort / "manifest.parquet")
    pl.DataFrame({"chrom": ["chr1"], "pos": [1000]}).write_parquet(cohort / "variants.parquet")
    pl.DataFrame({"case_id": ["c1", "c2"], "lineage": ["breast", "lung"]}).write_parquet(
        cohort / "samples.parquet"
    )
    (cohort / "cohort.json").write_text(
        json.dumps({"name": "fixture", "n_files": 1, "filter": {"op": "and"}})
    )
    if with_optional:
        pl.DataFrame({"hugo_symbol": ["TP53"], "lineage": ["breast"]}).write_parquet(
            cohort / "gene_frequency.parquet"
        )
        pl.DataFrame({"chrom": ["chr1"], "pos": [1000], "lineage": ["breast"]}).write_parquet(
            cohort / "variant_frequency.parquet"
        )
    return cohort


def test_load_cohort_returns_cohort(tmp_path: Path):
    path = _make_cohort_dir(tmp_path)
    cohort = load_cohort(path)
    assert isinstance(cohort, Cohort)
    assert cohort.name == "fixture"


def test_load_cohort_raises_for_missing_dir(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        load_cohort(tmp_path / "does-not-exist")


def test_cohort_required_frames_load_lazily(tmp_path: Path):
    cohort = load_cohort(_make_cohort_dir(tmp_path))
    assert len(cohort.clinical) == 2
    assert len(cohort.manifest) == 1
    assert len(cohort.variants) == 1
    assert len(cohort.samples) == 2


def test_cohort_required_frame_raises_when_missing(tmp_path: Path):
    cohort_dir = _make_cohort_dir(tmp_path)
    (cohort_dir / "variants.parquet").unlink()
    cohort = load_cohort(cohort_dir)
    with pytest.raises(FileNotFoundError, match=r"variants\.parquet"):
        _ = cohort.variants


def test_cohort_optional_frames_return_none_when_missing(tmp_path: Path):
    cohort = load_cohort(_make_cohort_dir(tmp_path, with_optional=False))
    assert cohort.gene_frequency is None
    assert cohort.variant_frequency is None


def test_cohort_optional_frames_load_when_present(tmp_path: Path):
    cohort = load_cohort(_make_cohort_dir(tmp_path, with_optional=True))
    assert cohort.gene_frequency is not None
    assert len(cohort.gene_frequency) == 1
    assert cohort.variant_frequency is not None


def test_cohort_provenance_parses_json(tmp_path: Path):
    cohort = load_cohort(_make_cohort_dir(tmp_path))
    p = cohort.provenance
    assert p["name"] == "fixture"
    assert p["n_files"] == 1


def test_cohort_provenance_empty_when_missing(tmp_path: Path):
    cohort_dir = _make_cohort_dir(tmp_path)
    (cohort_dir / "cohort.json").unlink()
    cohort = load_cohort(cohort_dir)
    assert cohort.provenance == {}


def test_cohort_summary_counts_rows(tmp_path: Path):
    cohort = load_cohort(_make_cohort_dir(tmp_path, with_optional=True))
    s = cohort.summary()
    assert s["n_clinical"] == 2
    assert s["n_variants"] == 1
    assert s["n_samples"] == 2
    assert s["n_gene_frequency"] == 1
    assert s["n_variant_frequency"] == 1


def test_repeated_access_is_cached(tmp_path: Path):
    cohort = load_cohort(_make_cohort_dir(tmp_path))
    df1 = cohort.variants
    df2 = cohort.variants
    assert df1 is df2  # cached_property
