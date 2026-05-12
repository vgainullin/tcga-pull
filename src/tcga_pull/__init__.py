"""tcga-pull: query GDC and produce a per-case structured data product.

Public API:
    load_cohort(path) -> Cohort  — open a cohort directory; lazy parquet views
    CohortSpec                   — declarative cohort spec (used by the CLI)
    GDCClient                    — low-level GDC API client
"""

__version__ = "0.1.2"

from .api import Cohort, load_cohort
from .config import CohortSpec
from .gdc import GDCClient
from .pipeline import fetch_preview, run

__all__ = [
    "Cohort",
    "CohortSpec",
    "GDCClient",
    "__version__",
    "fetch_preview",
    "load_cohort",
    "run",
]
