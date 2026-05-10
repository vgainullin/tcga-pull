"""GDC schema introspection — cached `_mapping` for files/cases/projects.

Used both to validate filter fields before issuing a query and to ground the
LLM in actual GDC field names rather than hallucinated ones.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

from .gdc import GDCClient

CACHE_DIR = Path.home() / ".cache" / "tcga-pull"
ENDPOINTS = ("files", "cases", "projects")


@dataclass
class Schema:
    endpoint: str
    raw: dict

    @property
    def fields(self) -> list[str]:
        """All queryable fields for this endpoint."""
        return list(self.raw.get("fields", []))

    @property
    def expand(self) -> list[str]:
        return list(self.raw.get("expand", []))

    @property
    def descriptions(self) -> dict[str, str]:
        """field -> description (shape comes from /_mapping `_mapping` key)."""
        out: dict[str, str] = {}
        mapping = self.raw.get("_mapping", {})
        for fname, meta in mapping.items():
            if isinstance(meta, dict):
                desc = meta.get("description") or ""
                if desc:
                    out[fname] = desc
        return out


class SchemaCache:
    def __init__(self, client: GDCClient | None = None, cache_dir: Path = CACHE_DIR):
        self.client = client or GDCClient()
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _path(self, endpoint: str) -> Path:
        return self.cache_dir / f"{endpoint}_mapping.json"

    def get(self, endpoint: str, refresh: bool = False) -> Schema:
        if endpoint not in ENDPOINTS:
            raise ValueError(f"unknown endpoint: {endpoint}")
        path = self._path(endpoint)
        if not refresh and path.exists():
            raw = json.loads(path.read_text())
        else:
            raw = self.client.mapping(endpoint)
            path.write_text(json.dumps(raw))
        return Schema(endpoint=endpoint, raw=raw)

    @cached_property
    def files(self) -> Schema:
        return self.get("files")

    @cached_property
    def cases(self) -> Schema:
        return self.get("cases")

    def validate_field(self, field: str) -> bool:
        """field may live on any endpoint (files.X, cases.X, etc.) — accept if any matches."""
        return any(field in self.get(ep).fields for ep in ENDPOINTS)
