"""Live regression test for the GDC facet-with-filter bug.

A filter built from open_access() + sugar mapping should NOT 500 when combined
with a nested facet. Bug history: before fix, file-rooted fields carried the
`files.` prefix, which trips a server-side 500 on /files when paired with any
facet (e.g. cases.project.project_id).

Marked `network` — only runs when you actually want to hit GDC.
"""

from __future__ import annotations

import json

import pytest
import requests

from tcga_pull.config import CohortSpec


@pytest.mark.network
def test_open_access_filter_does_not_500_with_facets():
    spec = CohortSpec(
        name="regression",
        out_dir=__import__("pathlib").Path("/tmp"),
        filters={
            "project": "TCGA-CHOL",
            "data_category": "Simple Nucleotide Variation",
            "data_format": "MAF",
        },
    )
    flt = spec.resolve_filter()

    # Same shape as the call that used to 500
    r = requests.post(
        "https://api.gdc.cancer.gov/files",
        json={
            "filters": flt,
            "facets": "cases.project.project_id",
            "size": 0,
        },
        timeout=30,
    )
    # Surface the response body when it fails — we want to know if GDC changes
    # behavior on us in the future.
    assert r.status_code == 200, f"GDC returned {r.status_code}: {r.text[:300]}"

    data = r.json()
    total = data["data"]["pagination"]["total"]
    buckets = data["data"]["aggregations"]["cases.project.project_id"]["buckets"]
    assert total > 0, f"expected files for TCGA-CHOL SNV MAF, got {total}"
    chol_bucket = next((b for b in buckets if b["key"] == "TCGA-CHOL"), None)
    assert chol_bucket is not None, (
        f"TCGA-CHOL missing from facet buckets: {json.dumps([b['key'] for b in buckets[:5]])}"
    )
    assert chol_bucket["doc_count"] == total
