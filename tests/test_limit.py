"""Tests for the per-project limit applied at the file-list stage."""

from __future__ import annotations

from tcga_pull.config import LimitSpec
from tcga_pull.pipeline import apply_limit


def _hit(file_id: str, submitter: str, project: str, *, case_id: str | None = None) -> dict:
    return {
        "file_id": file_id,
        "file_name": f"{file_id}.maf.gz",
        "file_size": 1000,
        "cases": [
            {
                "case_id": case_id or f"case-{submitter}",
                "submitter_id": submitter,
                "project": {"project_id": project},
            }
        ],
    }


def test_apply_limit_disabled_returns_input_unchanged():
    hits = [_hit("f1", "TCGA-AA-0001", "TCGA-BRCA")]
    assert apply_limit(hits, LimitSpec()) is hits  # short-circuit


def test_apply_limit_per_project_caps_each_project_at_n():
    # BRCA has 5 patients, LUAD has 3. limit=2 should keep first 2 of each.
    brca = [_hit(f"f{i}", f"TCGA-BR-{i:04d}", "TCGA-BRCA") for i in range(5)]
    luad = [_hit(f"l{i}", f"TCGA-LU-{i:04d}", "TCGA-LUAD") for i in range(3)]
    kept = apply_limit(brca + luad, LimitSpec(per_project=2))

    # 2 patients x 2 projects = 4 files (one per patient here)
    assert len(kept) == 4
    kept_submitters = sorted(h["cases"][0]["submitter_id"] for h in kept)
    assert kept_submitters == [
        "TCGA-BR-0000",
        "TCGA-BR-0001",
        "TCGA-LU-0000",
        "TCGA-LU-0001",
    ]


def test_apply_limit_is_deterministic_by_submitter_sort():
    # Input order shouldn't matter — sorted submitter_id picks the same cases.
    hits_a = [_hit(f"f{i}", f"S{i:02d}", "TCGA-X") for i in [3, 1, 4, 0, 2]]
    hits_b = [_hit(f"g{i}", f"S{i:02d}", "TCGA-X") for i in [0, 1, 2, 3, 4]]
    kept_a = apply_limit(hits_a, LimitSpec(per_project=3))
    kept_b = apply_limit(hits_b, LimitSpec(per_project=3))
    subs_a = sorted(h["cases"][0]["submitter_id"] for h in kept_a)
    subs_b = sorted(h["cases"][0]["submitter_id"] for h in kept_b)
    assert subs_a == subs_b == ["S00", "S01", "S02"]


def test_apply_limit_keeps_all_when_project_smaller_than_n():
    hits = [_hit(f"f{i}", f"S{i}", "TCGA-X") for i in range(3)]
    kept = apply_limit(hits, LimitSpec(per_project=10))
    assert len(kept) == 3


def test_apply_limit_keeps_multi_aliquot_files_for_kept_cases():
    # One patient has 3 aliquots; limit=1 should keep that 1 patient with all 3 files
    hits = [_hit(f"a{i}", "TCGA-XX-0001", "TCGA-X", case_id="case-A") for i in range(3)] + [
        _hit(f"b{i}", "TCGA-XX-0002", "TCGA-X", case_id="case-B") for i in range(2)
    ]
    kept = apply_limit(hits, LimitSpec(per_project=1))
    # Should keep the 3 aliquots of patient A, drop patient B's 2 files
    assert len(kept) == 3
    assert {h["cases"][0]["submitter_id"] for h in kept} == {"TCGA-XX-0001"}


def test_apply_limit_passes_multi_case_files_through():
    # A file mapped to >1 case (project-wide MAF). limit shouldn't drop it.
    multi_case_hit = {
        "file_id": "m1",
        "file_name": "m1.maf.gz",
        "file_size": 1000,
        "cases": [
            {"case_id": "c1", "submitter_id": "S1", "project": {"project_id": "TCGA-X"}},
            {"case_id": "c2", "submitter_id": "S2", "project": {"project_id": "TCGA-X"}},
        ],
    }
    hits = [_hit(f"f{i}", f"S{i + 10}", "TCGA-X") for i in range(5)] + [multi_case_hit]
    kept = apply_limit(hits, LimitSpec(per_project=2))
    # 2 single-case + 1 multi-case = 3 files
    assert multi_case_hit in kept
    assert sum(1 for h in kept if len(h["cases"]) == 1) == 2
