"""Stable Python API for consuming a tcga-pull cohort from another project.

Usage:

    from tcga_pull import load_cohort
    cohort = load_cohort("/path/to/cohort")
    variants = cohort.variants          # polars DataFrame
    samples = cohort.samples
    gene_freq = cohort.gene_frequency   # None if `tcga-pull frequency` not run
    print(cohort.summary())

Each frame is lazy-loaded on first access and cached on the instance. Create a
new `Cohort` if you want to force a re-read from disk. The parquet column
schemas are documented in `SCHEMAS.md`.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import Any

import polars as pl


@dataclass
class Cohort:
    """Read-only view over a tcga-pull cohort directory.

    Required files (raised as `FileNotFoundError` if missing):
      - clinical.parquet
      - manifest.parquet
      - variants.parquet
      - samples.parquet

    Optional files (returned as `None` if missing):
      - gene_frequency.parquet
      - variant_frequency.parquet
      - cohort.json    (provenance sidecar)
    """

    path: Path

    @property
    def name(self) -> str:
        return self.path.name

    def _read_required(self, fname: str) -> pl.DataFrame:
        p = self.path / fname
        if not p.exists():
            raise FileNotFoundError(f"missing {p} — did you run the matching tcga-pull step?")
        return pl.read_parquet(p)

    def _read_optional(self, fname: str) -> pl.DataFrame | None:
        p = self.path / fname
        return pl.read_parquet(p) if p.exists() else None

    # --- required parquets -----------------------------------------------------

    @cached_property
    def clinical(self) -> pl.DataFrame:
        """One row per case. Demographics + diagnosis fields flattened from GDC."""
        return self._read_required("clinical.parquet")

    @cached_property
    def manifest(self) -> pl.DataFrame:
        """One row per file. Includes local_path after a successful download."""
        return self._read_required("manifest.parquet")

    @cached_property
    def variants(self) -> pl.DataFrame:
        """One row per (variant x tumor aliquot). ~37 columns. See SCHEMAS.md."""
        return self._read_required("variants.parquet")

    @cached_property
    def samples(self) -> pl.DataFrame:
        """One row per case. Includes the curated `lineage` (tissue) column."""
        return self._read_required("samples.parquet")

    # --- optional parquets -----------------------------------------------------

    @cached_property
    def gene_frequency(self) -> pl.DataFrame | None:
        """One row per (gene, lineage). Produced by `tcga-pull frequency`."""
        return self._read_optional("gene_frequency.parquet")

    @cached_property
    def variant_frequency(self) -> pl.DataFrame | None:
        """One row per (variant, lineage). Produced by `tcga-pull frequency`."""
        return self._read_optional("variant_frequency.parquet")

    # --- provenance ------------------------------------------------------------

    @cached_property
    def provenance(self) -> dict[str, Any]:
        """Resolved filter + counts + timestamp from cohort.json. {} if missing."""
        p = self.path / "cohort.json"
        return json.loads(p.read_text()) if p.exists() else {}

    # --- summary ---------------------------------------------------------------

    def summary(self) -> dict[str, int | str]:
        """Cheap shape summary — number of rows in each available parquet."""
        out: dict[str, int | str] = {"path": str(self.path), "name": self.name}
        for attr in ("clinical", "manifest", "samples", "variants"):
            try:
                df = getattr(self, attr)
                out[f"n_{attr}"] = len(df)
            except FileNotFoundError:
                out[f"n_{attr}"] = 0
        for opt in ("gene_frequency", "variant_frequency"):
            df = getattr(self, opt)
            out[f"n_{opt}"] = len(df) if df is not None else 0
        return out


def load_cohort(path: str | Path) -> Cohort:
    """Open a tcga-pull cohort directory produced by `tcga-pull pull`.

    Raises `FileNotFoundError` if the path isn't a directory. Individual parquet
    files are checked on first access, not here — so an empty cohort dir is
    still a valid handle; reading `cohort.variants` later will raise.
    """
    p = Path(path).expanduser()
    if not p.is_dir():
        raise FileNotFoundError(f"not a directory: {p}")
    return Cohort(path=p)
