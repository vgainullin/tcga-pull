"""Write GDC manifest, fetch files (bulk API or gdc-client), restructure into per-case folders."""

from __future__ import annotations

import contextlib
import csv
import re
import shutil
import subprocess
import tarfile
from collections.abc import Iterable
from pathlib import Path

import requests
from rich.console import Console

from .gdc import GDC_API

GDC_CLIENT_BIN = "gdc-client"
BULK_BATCH_SIZE = 200
BULK_FILE_SIZE_LIMIT = 100 * 1024 * 1024  # use bulk only when every file is <100MB


def slugify(s: str | None) -> str:
    out = (s or "unknown").strip().lower()
    out = re.sub(r"[^\w]+", "_", out)
    return out.strip("_") or "unknown"


def primary_case(file_hit: dict) -> tuple[str, str] | None:
    """Return (case_id, submitter_id) for files mapped to exactly one case."""
    cases = file_hit.get("cases") or []
    if len(cases) != 1:
        return None
    c = cases[0]
    case_id = c.get("case_id")
    submitter = c.get("submitter_id") or case_id
    if not case_id:
        return None
    return case_id, submitter


def write_manifest_tsv(file_hits: list[dict], path: Path) -> int:
    """Write the gdc-client manifest TSV. Returns row count."""
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["id", "filename", "md5", "size", "state"]
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for h in file_hits:
            w.writerow(
                [
                    h["file_id"],
                    h["file_name"],
                    h.get("md5sum", ""),
                    h.get("file_size", ""),
                    "validated",
                ]
            )
    return len(file_hits)


def check_gdc_client() -> str:
    """Return the resolved gdc-client path or raise if missing."""
    found = shutil.which(GDC_CLIENT_BIN)
    if not found:
        raise FileNotFoundError(
            "gdc-client not found on PATH. "
            "Install from https://gdc.cancer.gov/access-data/gdc-data-transfer-tool"
        )
    return found


def should_use_bulk(file_hits: list[dict]) -> bool:
    """Bulk API beats gdc-client by orders of magnitude for many small files
    (no per-file handshake overhead). Skip it for any file >100MB — gdc-client's
    segmented/resumable download earns its keep there."""
    return all(int(h.get("file_size") or 0) < BULK_FILE_SIZE_LIMIT for h in file_hits)


def bulk_download_via_api(
    file_hits: list[dict],
    download_dir: Path,
    *,
    base_url: str = GDC_API,
    batch_size: int = BULK_BATCH_SIZE,
    console: Console | None = None,
) -> None:
    """POST file IDs to GDC `/data` and stream-extract the tar.gz response into
    `download_dir/<file_id>/<file_name>` — the same layout gdc-client produces,
    so the existing `restructure()` step works unchanged.

    Open-access only (no token). For controlled files use gdc-client.
    """
    download_dir.mkdir(parents=True, exist_ok=True)
    console = console or Console()
    url = base_url.rstrip("/") + "/data"
    ids = [h["file_id"] for h in file_hits]
    name_by_id = {h["file_id"]: h["file_name"] for h in file_hits}
    total = len(ids)
    n_batches = (total + batch_size - 1) // batch_size
    done = 0

    for bi in range(n_batches):
        batch = ids[bi * batch_size : (bi + 1) * batch_size]
        with (
            console.status(
                f"[cyan]POST /data  batch {bi + 1}/{n_batches}  ({len(batch)} files)[/cyan]"
            ),
            requests.post(url, json={"ids": batch}, stream=True, timeout=600) as r,
        ):
            r.raise_for_status()
            # GDC behavior: 1 id → raw file; ≥2 ids → tar.gz of <id>/<name>
            if len(batch) > 1:
                with tarfile.open(fileobj=r.raw, mode="r|gz") as tar:
                    for m in tar:
                        if m.isfile():
                            tar.extract(m, path=download_dir)
            else:
                fid = batch[0]
                fname = name_by_id[fid]
                dest = download_dir / fid / fname
                dest.parent.mkdir(parents=True, exist_ok=True)
                with dest.open("wb") as fh:
                    for chunk in r.iter_content(chunk_size=1 << 16):
                        if chunk:
                            fh.write(chunk)
        done += len(batch)
        console.log(f"  bulk: {done}/{total}")


def run_gdc_client(
    manifest: Path,
    download_dir: Path,
    *,
    n_processes: int = 4,
    extra_args: Iterable[str] = (),
    console: Console | None = None,
) -> None:
    """Invoke gdc-client, streaming stdout to console."""
    bin_path = check_gdc_client()
    download_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        bin_path,
        "download",
        "-m",
        str(manifest),
        "-d",
        str(download_dir),
        "-n",
        str(n_processes),
        *extra_args,
    ]
    console = console or Console()
    console.log(f"[dim]$ {' '.join(cmd)}[/dim]")
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    assert proc.stdout is not None
    for line in proc.stdout:
        console.out(line.rstrip())
    rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"gdc-client exited with code {rc}")


def restructure(
    file_hits: list[dict],
    download_dir: Path,
    cohort_data_dir: Path,
    *,
    move: bool = True,
) -> list[dict]:
    """Move/copy gdc-client output (download_dir/<file_id>/<file_name>) into
    cohort_data_dir/<submitter_id>/<data_category_slug>/<file_name>.

    Returns a list of post-move records: file_id, case_id, submitter_id, data_category, local_path.
    Files mapped to !=1 case go into cohort_data_dir/_multi/<data_category_slug>/.
    """
    cohort_data_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict] = []

    for h in file_hits:
        file_id = h["file_id"]
        file_name = h["file_name"]
        src = download_dir / file_id / file_name
        if not src.exists():
            # gdc-client sometimes nests differently if state is partial; skip with note
            records.append({**_record_base(h), "local_path": None, "status": "missing"})
            continue

        cat = slugify(h.get("data_category"))
        case = primary_case(h)
        if case is None:
            dest_dir = cohort_data_dir / "_multi" / cat
            submitter = None
            case_id = None
        else:
            case_id, submitter = case
            dest_dir = cohort_data_dir / submitter / cat
        dest_dir.mkdir(parents=True, exist_ok=True)
        dest = dest_dir / file_name

        if move:
            shutil.move(str(src), str(dest))
            # remove the now-empty <file_id> dir
            with contextlib.suppress(OSError):
                src.parent.rmdir()
        else:
            shutil.copy2(str(src), str(dest))

        records.append(
            {
                **_record_base(h),
                "case_id": case_id,
                "submitter_id": submitter,
                "data_category": h.get("data_category"),
                "local_path": str(dest),
                "status": "ok",
            }
        )
    return records


def _record_base(h: dict) -> dict:
    workflow = (h.get("analysis") or {}).get("workflow_type")
    return {
        "file_id": h["file_id"],
        "file_name": h["file_name"],
        "data_type": h.get("data_type"),
        "data_format": h.get("data_format"),
        "experimental_strategy": h.get("experimental_strategy"),
        "workflow_type": workflow,
        "md5sum": h.get("md5sum"),
        "file_size": h.get("file_size"),
    }
