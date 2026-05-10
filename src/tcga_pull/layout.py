"""Flatten clinical / manifest into parquet outputs.

clinical.parquet  : one row per case (first-of for nested lists, plus counts)
clinical_raw.jsonl: full nested case records, one JSON object per line
manifest.parquet  : one row per file with case mapping + local_path
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


def _flatten_scalar(prefix: str, obj: dict | None) -> dict[str, Any]:
    """Promote scalar leaves of `obj` to <prefix>_<key>."""
    if not obj:
        return {}
    out: dict[str, Any] = {}
    for k, v in obj.items():
        if isinstance(v, (dict, list)):
            continue
        out[f"{prefix}_{k}"] = v
    return out


def flatten_case(case: dict) -> dict[str, Any]:
    """One row per case. First-of for nested lists; counts for awareness."""
    row: dict[str, Any] = {
        "case_id": case.get("case_id"),
        "submitter_id": case.get("submitter_id"),
        "project_id": (case.get("project") or {}).get("project_id"),
        "disease_type": case.get("disease_type"),
        "primary_site": case.get("primary_site"),
    }
    row.update(_flatten_scalar("demographic", case.get("demographic")))

    diagnoses = case.get("diagnoses") or []
    row["n_diagnoses"] = len(diagnoses)
    if diagnoses:
        row.update(_flatten_scalar("diagnosis", diagnoses[0]))
        treatments = (diagnoses[0].get("treatments")) or []
        row["n_treatments"] = sum(len(d.get("treatments") or []) for d in diagnoses)
        if treatments:
            row.update(_flatten_scalar("treatment", treatments[0]))
    else:
        row["n_treatments"] = 0

    exposures = case.get("exposures") or []
    row["n_exposures"] = len(exposures)
    if exposures:
        row.update(_flatten_scalar("exposure", exposures[0]))

    family = case.get("family_histories") or []
    row["n_family_histories"] = len(family)
    if family:
        row.update(_flatten_scalar("family_history", family[0]))

    return row


def write_clinical(cases: list[dict], cohort_dir: Path) -> tuple[Path, Path]:
    cohort_dir.mkdir(parents=True, exist_ok=True)
    raw_path = cohort_dir / "clinical_raw.jsonl"
    with raw_path.open("w") as fh:
        for c in cases:
            fh.write(json.dumps(c) + "\n")

    df = pd.DataFrame(flatten_case(c) for c in cases)
    parquet_path = cohort_dir / "clinical.parquet"
    df.to_parquet(parquet_path, index=False)
    return parquet_path, raw_path


def write_manifest(records: list[dict], cohort_dir: Path) -> Path:
    cohort_dir.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(records)
    path = cohort_dir / "manifest.parquet"
    df.to_parquet(path, index=False)
    return path


def write_provenance(cohort_dir: Path, payload: dict) -> Path:
    """Resolved filter + counts + timestamps. Lets users re-run a cohort."""
    import time

    cohort_dir.mkdir(parents=True, exist_ok=True)
    payload = {**payload, "created_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")}
    path = cohort_dir / "cohort.json"
    path.write_text(json.dumps(payload, indent=2))
    return path
