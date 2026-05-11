"""GDC API client. Open-access TCGA queries only by default.

Filter grammar is GDC's native JSON tree:
    {"op": "and", "content": [<clause>, <clause>, ...]}
    {"op": "in",  "content": {"field": "...", "value": [...]}}
    {"op": "=",   "content": {"field": "...", "value": "..."}}

Every query is wrapped with files.access = "open" before being sent.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field
from typing import Any
from urllib.parse import urljoin

import requests

GDC_API = "https://api.gdc.cancer.gov/"
DEFAULT_TIMEOUT = 60
PAGE_SIZE = 1000

FILE_FIELDS: list[str] = [
    "file_id",
    "file_name",
    "data_category",
    "data_type",
    "data_format",
    "experimental_strategy",
    "analysis.workflow_type",
    "access",
    "md5sum",
    "file_size",
    "cases.case_id",
    "cases.submitter_id",
    "cases.project.project_id",
    "cases.samples.sample_type",
    "cases.samples.tissue_type",
]

CASE_FIELDS: list[str] = [
    "case_id",
    "submitter_id",
    "project.project_id",
    "disease_type",
    "primary_site",
]

CLINICAL_EXPAND: list[str] = [
    "demographic",
    "diagnoses",
    "diagnoses.treatments",
    "exposures",
    "family_histories",
]


def f_in(field_name: str, values: Any) -> dict:
    if isinstance(values, (str, int, float)):
        values = [values]
    return {"op": "in", "content": {"field": field_name, "value": list(values)}}


def f_eq(field_name: str, value: Any) -> dict:
    return {"op": "=", "content": {"field": field_name, "value": value}}


def f_and(*clauses: dict | None) -> dict:
    """Combine clauses under AND. Flattens nested AND children so the result
    is always a single-level `and` (or just the lone clause if one input).
    Falsy clauses (None, {}) are skipped — caller can pass conditional clauses.

    GDC's /files endpoint silently 500s on faceted queries when the filter
    tree contains a nested AND. Keeping the tree flat dodges that and reads
    more naturally too.
    """
    flat: list[dict] = []
    for c in clauses:
        if not c:
            continue
        if c.get("op") == "and" and isinstance(c.get("content"), list):
            flat.extend(c["content"])
        else:
            flat.append(c)
    if len(flat) == 1:
        return flat[0]
    return {"op": "and", "content": flat}


def f_or(*clauses: dict | None) -> dict:
    """Combine clauses under OR, flattening nested OR children for symmetry
    with f_and (and the same GDC quirk likely applies)."""
    flat: list[dict] = []
    for c in clauses:
        if not c:
            continue
        if c.get("op") == "or" and isinstance(c.get("content"), list):
            flat.extend(c["content"])
        else:
            flat.append(c)
    if len(flat) == 1:
        return flat[0]
    return {"op": "or", "content": flat}


def open_access(user_filter: dict | None) -> dict:
    """Wrap a user filter with the mandatory access = open clause.

    Note: bare `access` (no `files.` prefix). GDC's /files endpoint silently
    500s on faceted queries when file-rooted fields carry the `files.` prefix;
    using bare names avoids that. `for_cases_endpoint` re-prefixes to
    `files.access` when the filter is translated to /cases.
    """
    access = f_in("access", ["open"])
    return f_and(user_filter, access) if user_filter else access


# Field-name prefixes that survive translation unchanged when going /files → /cases
# (they're either case-rooted or already-prefixed file joins).
_FILES_ROOTED_PREFIXES_ON_CASES: tuple[str, ...] = ("files.",)


def for_cases_endpoint(flt: dict) -> dict:
    """Translate a /files-rooted filter for use on /cases.

    Our convention on /files:
      * bare names (e.g. "access", "data_format") = file's own fields
      * "cases.X" = case-side join
      * "analysis.X" = analysis-side join (analysis belongs to file)

    On /cases the root flips, so:
      * bare names → "files.X" (file's own fields become a join from case)
      * "cases.X" → "X" (case-rooted, strip prefix)
      * "analysis.X" → "files.analysis.X" (analysis joins via file)
      * Already-prefixed "files.X" passes through unchanged.
    """

    def walk(node: Any) -> Any:
        if not isinstance(node, dict):
            return node
        c = node.get("content")
        if isinstance(c, list):
            return {**node, "content": [walk(x) for x in c]}
        if isinstance(c, dict) and "field" in c:
            f = c["field"]
            if f.startswith("cases."):
                f = f[len("cases.") :]
            elif f.startswith("analysis."):
                f = "files." + f
            elif not f.startswith(_FILES_ROOTED_PREFIXES_ON_CASES):
                # bare file field on /files → `files.X` join on /cases
                f = "files." + f
            return {**node, "content": {**c, "field": f}}
        return node

    result: dict = walk(flt)
    return result


@dataclass
class GDCClient:
    base_url: str = GDC_API
    timeout: int = DEFAULT_TIMEOUT
    session: requests.Session = field(default_factory=requests.Session)

    def _post(self, endpoint: str, body: dict) -> dict[str, Any]:
        url = urljoin(self.base_url, endpoint)
        resp = self.session.post(url, json=body, timeout=self.timeout)
        resp.raise_for_status()
        data: dict[str, Any] = resp.json()
        return data

    def _get(self, endpoint: str, params: dict | None = None) -> dict[str, Any]:
        url = urljoin(self.base_url, endpoint)
        resp = self.session.get(url, params=params, timeout=self.timeout)
        resp.raise_for_status()
        data: dict[str, Any] = resp.json()
        return data

    # --- counts (cheap dry-run) -------------------------------------------------

    def count(self, endpoint: str, filters: dict) -> int:
        """Hit pagination metadata only. Use with /files or /cases."""
        body = {"filters": filters, "size": 0}
        data = self._post(endpoint, body)
        return int(data["data"]["pagination"]["total"])

    def count_files(self, filters: dict) -> int:
        return self.count("files", filters)

    def count_cases(self, filters: dict) -> int:
        return self.count("cases", for_cases_endpoint(filters))

    # --- projects ---------------------------------------------------------------

    def list_projects(self, program: str = "TCGA") -> list[dict]:
        body = {
            "filters": f_in("program.name", [program]),
            "fields": "project_id,name,disease_type,primary_site",
            "size": 1000,
        }
        data = self._post("projects", body)
        hits: list[dict] = data["data"]["hits"]
        return hits

    # --- files ------------------------------------------------------------------

    def iter_files(
        self,
        filters: dict,
        fields: list[str] | None = None,
        page_size: int = PAGE_SIZE,
    ) -> Iterator[dict]:
        """Yield file hits page by page."""
        fields = fields or FILE_FIELDS
        body_base = {
            "filters": filters,
            "fields": ",".join(fields),
            "size": page_size,
        }
        offset = 0
        while True:
            body = {**body_base, "from": offset}
            data = self._post("files", body)
            hits = data["data"]["hits"]
            yield from hits
            pagination = data["data"]["pagination"]
            if offset + len(hits) >= pagination["total"] or not hits:
                return
            offset += len(hits)

    def fetch_files(self, filters: dict, fields: list[str] | None = None) -> list[dict]:
        return list(self.iter_files(filters, fields=fields))

    # --- cases / clinical -------------------------------------------------------

    def iter_cases(
        self,
        filters: dict,
        fields: list[str] | None = None,
        expand: list[str] | None = None,
        page_size: int = PAGE_SIZE,
    ) -> Iterator[dict]:
        fields = fields or CASE_FIELDS
        expand = expand if expand is not None else CLINICAL_EXPAND
        body_base = {
            "filters": for_cases_endpoint(filters),
            "fields": ",".join(fields),
            "expand": ",".join(expand),
            "size": page_size,
        }
        offset = 0
        while True:
            body = {**body_base, "from": offset}
            data = self._post("cases", body)
            hits = data["data"]["hits"]
            yield from hits
            pagination = data["data"]["pagination"]
            if offset + len(hits) >= pagination["total"] or not hits:
                return
            offset += len(hits)

    def fetch_clinical(self, filters: dict) -> list[dict]:
        return list(self.iter_cases(filters))

    # --- schema -----------------------------------------------------------------

    def mapping(self, endpoint: str) -> dict:
        """GET /<endpoint>/_mapping. Returns schema for filters/fields/expand."""
        return self._get(f"{endpoint}/_mapping")
