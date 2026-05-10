"""tcga-pull: query GDC and produce a per-case structured data product."""

from .config import CohortSpec
from .gdc import GDCClient
from .pipeline import fetch_preview, run

__all__ = ["CohortSpec", "GDCClient", "fetch_preview", "run"]
