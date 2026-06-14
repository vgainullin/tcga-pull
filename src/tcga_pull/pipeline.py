"""End-to-end orchestration: filter -> preview -> download -> restructure -> parquets."""

from __future__ import annotations

from collections.abc import Callable, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from rich.console import Console
from rich.table import Table

from .config import CohortSpec, LimitSpec
from .download import (
    bulk_download_via_api,
    primary_case,
    restructure,
    run_gdc_client,
    should_use_bulk,
    write_manifest_tsv,
)
from .gdc import GDCClient
from .layout import write_clinical, write_manifest, write_provenance


def _chunks(items: list[dict], size: int) -> Iterator[list[dict]]:
    for i in range(0, len(items), size):
        yield items[i : i + size]


@dataclass
class Preview:
    spec: CohortSpec
    resolved_filter: dict
    n_files: int
    n_cases: int
    file_hits: list[dict]
    total_size: int

    @property
    def by_data_category(self) -> dict[str, int]:
        out: dict[str, int] = {}
        for h in self.file_hits:
            cat = h.get("data_category") or "unknown"
            out[cat] = out.get(cat, 0) + 1
        return out


def apply_limit(file_hits: list[dict], limit: LimitSpec) -> list[dict]:
    """Trim file_hits according to a LimitSpec. Deterministic: cases are
    sorted by submitter_id within each project, and the first N are kept.

    Files mapped to a single case participate in the limit. Multi-case files
    (rare; project-level outputs) pass through unchanged — they belong to
    every case in their group and we don't want to silently drop them.
    """
    n = limit.per_project
    if n is None:
        return file_hits

    # Build {project_id: {submitter_id: case_id}} from single-case files
    by_project: dict[str, dict[str, str]] = {}
    for h in file_hits:
        cases = h.get("cases") or []
        if len(cases) != 1:
            continue
        c = cases[0]
        project_id = (c.get("project") or {}).get("project_id")
        case_id = c.get("case_id")
        submitter = c.get("submitter_id") or case_id
        if project_id and case_id and submitter:
            by_project.setdefault(project_id, {})[submitter] = case_id

    kept_case_ids: set[str] = set()
    for submitters in by_project.values():
        for submitter in sorted(submitters)[:n]:
            kept_case_ids.add(submitters[submitter])

    def keep(h: dict) -> bool:
        cases = h.get("cases") or []
        if len(cases) != 1:
            return True  # multi-case → pass through
        return cases[0].get("case_id") in kept_case_ids

    return [h for h in file_hits if keep(h)]


def fetch_preview(spec: CohortSpec, client: GDCClient | None = None) -> Preview:
    client = client or GDCClient()
    flt = spec.resolve_filter()
    file_hits = client.fetch_files(flt)
    file_hits = apply_limit(file_hits, spec.limit)
    n_files = len(file_hits)
    case_ids = {
        c.get("case_id") for h in file_hits for c in (h.get("cases") or []) if c.get("case_id")
    }
    total_size = sum(int(h.get("file_size") or 0) for h in file_hits)
    return Preview(
        spec=spec,
        resolved_filter=flt,
        n_files=n_files,
        n_cases=len(case_ids),
        file_hits=file_hits,
        total_size=total_size,
    )


def render_preview(p: Preview, console: Console | None = None) -> None:
    console = console or Console()
    spec = p.spec
    t = Table(title=f"Cohort preview: [bold]{spec.name}[/bold]", show_header=False)
    t.add_column("k", style="dim")
    t.add_column("v")
    t.add_row("Output", str(spec.cohort_dir))
    t.add_row("Files", f"{p.n_files:,}")
    t.add_row("Cases", f"{p.n_cases:,}")
    t.add_row("Total size", _human_size(p.total_size))
    multi = sum(1 for h in p.file_hits if primary_case(h) is None)
    if multi:
        t.add_row("Multi-case files", f"{multi:,} (will go to data/_multi/)")
    console.print(t)

    cat = p.by_data_category
    if cat:
        c = Table(title="By data_category", show_header=True)
        c.add_column("category")
        c.add_column("files", justify="right")
        for k, v in sorted(cat.items(), key=lambda x: -x[1]):
            c.add_row(k, f"{v:,}")
        console.print(c)


def run(
    spec: CohortSpec,
    *,
    console: Console | None = None,
    client: GDCClient | None = None,
) -> Path:
    """Full pipeline. Returns the cohort directory.

    The preview table is rendered before any bytes move; user can Ctrl-C if
    the resolved filter is wrong. Use `tcga-pull preview <yaml>` if you want
    a true dry-run.
    """
    console = console or Console()
    client = client or GDCClient()

    with console.status("[cyan]Querying GDC…[/cyan]"):
        preview = fetch_preview(spec, client=client)
    render_preview(preview, console=console)

    if preview.n_files == 0:
        console.print("[yellow]No files match. Refine the filter.[/yellow]")
        return spec.cohort_dir

    if spec.processing.mode == "incremental":
        return _run_incremental(spec, preview=preview, console=console, client=client)

    cohort_dir = spec.cohort_dir
    cohort_dir.mkdir(parents=True, exist_ok=True)
    download_tmp = cohort_dir / "_downloads"
    data_dir = cohort_dir / "data"

    manifest_tsv = cohort_dir / "manifest.tsv"
    write_manifest_tsv(preview.file_hits, manifest_tsv)
    console.log(f"[green]wrote[/green] {manifest_tsv}")

    if should_use_bulk(preview.file_hits):
        console.log("[dim]download: bulk /data endpoint (small files)[/dim]")
        bulk_download_via_api(
            preview.file_hits,
            download_tmp,
            console=console,
        )
    else:
        console.log("[dim]download: gdc-client (some files >100MB)[/dim]")
        run_gdc_client(
            manifest=manifest_tsv,
            download_dir=download_tmp,
            n_processes=spec.n_processes,
            console=console,
        )

    with console.status("[cyan]Restructuring into per-case folders…[/cyan]"):
        records = restructure(preview.file_hits, download_tmp, data_dir, move=True)

    with console.status("[cyan]Fetching clinical metadata…[/cyan]"):
        cases = client.fetch_clinical(preview.resolved_filter)
    clinical_path, raw_path = write_clinical(cases, cohort_dir)
    manifest_path = write_manifest(records, cohort_dir)
    prov_path = write_provenance(
        cohort_dir,
        {
            "name": spec.name,
            "filter": preview.resolved_filter,
            "n_files": preview.n_files,
            "n_cases": preview.n_cases,
            "total_size": preview.total_size,
        },
    )

    # cleanup gdc-client's working dir (per-file logs etc.)
    import shutil

    shutil.rmtree(download_tmp, ignore_errors=True)

    console.print(f"[green]done[/green] -> {cohort_dir}")
    console.print(f"  clinical : {clinical_path}")
    console.print(f"  manifest : {manifest_path}")
    console.print(f"  provenance: {prov_path}")
    console.print(f"  raw      : {raw_path}")

    # Run any post-processing recipes declared in the YAML
    for recipe_name in spec.recipes:
        recipe = RECIPE_REGISTRY.get(recipe_name)
        if recipe is None:  # already validated in CohortSpec.__post_init__
            continue
        console.print(f"\n[cyan]==> recipe: {recipe_name}[/cyan]")
        recipe(cohort_dir, spec.recipe_options)

    return cohort_dir


def _run_incremental(
    spec: CohortSpec,
    *,
    preview: Preview,
    console: Console,
    client: GDCClient,
) -> Path:
    """Batch-oriented pull for large cohorts.

    Batch-aware multiomics recipes process each batch immediately. When
    configured, raw files handled by those recipes are deleted before the next
    batch downloads. Non-batch-aware recipes still run at the end against any
    raw files intentionally retained, e.g. SNV MAFs for variants/samples.
    """
    from .multiomics import (
        MULTIOMICS_RECIPE_NAMES,
        finalize_multiomics_parts,
        write_multiomics_parts,
    )

    if not should_use_bulk(preview.file_hits):
        raise ValueError("incremental processing currently supports only bulk-API sized files")

    cohort_dir = spec.cohort_dir
    cohort_dir.mkdir(parents=True, exist_ok=True)
    download_tmp = cohort_dir / "_downloads"
    data_dir = cohort_dir / "data"

    manifest_tsv = cohort_dir / "manifest.tsv"
    write_manifest_tsv(preview.file_hits, manifest_tsv)
    console.log(f"[green]wrote[/green] {manifest_tsv}")

    batch_recipes = [r for r in spec.recipes if r in MULTIOMICS_RECIPE_NAMES]
    records: list[dict] = []
    total_batches = (len(preview.file_hits) + spec.processing.batch_size - 1) // (
        spec.processing.batch_size
    )
    for bi, file_batch in enumerate(
        _chunks(preview.file_hits, spec.processing.batch_size), start=1
    ):
        console.print(
            f"\n[cyan]==> incremental batch {bi}/{total_batches}[/cyan] ({len(file_batch):,} files)"
        )
        bulk_download_via_api(
            file_batch,
            download_tmp,
            batch_size=spec.processing.batch_size,
            console=console,
        )
        batch_records = restructure(file_batch, download_tmp, data_dir, move=True)
        if batch_recipes:
            write_multiomics_parts(
                cohort_dir,
                batch_records,
                part_id=bi,
                recipes=batch_recipes,
                recipe_options=spec.recipe_options,
            )
            if spec.processing.delete_raw_after_processing:
                _delete_batch_raw(batch_records, recipes=batch_recipes)
        records.extend(batch_records)

    with console.status("[cyan]Fetching clinical metadata…[/cyan]"):
        cases = client.fetch_clinical(preview.resolved_filter)
    clinical_path, raw_path = write_clinical(cases, cohort_dir)
    manifest_path = write_manifest(records, cohort_dir)
    prov_path = write_provenance(
        cohort_dir,
        {
            "name": spec.name,
            "filter": preview.resolved_filter,
            "n_files": preview.n_files,
            "n_cases": preview.n_cases,
            "total_size": preview.total_size,
            "processing": {
                "mode": spec.processing.mode,
                "batch_size": spec.processing.batch_size,
                "delete_raw_after_processing": spec.processing.delete_raw_after_processing,
            },
        },
    )

    import shutil

    shutil.rmtree(download_tmp, ignore_errors=True)

    console.print(f"[green]done[/green] -> {cohort_dir}")
    console.print(f"  clinical : {clinical_path}")
    console.print(f"  manifest : {manifest_path}")
    console.print(f"  provenance: {prov_path}")
    console.print(f"  raw      : {raw_path}")

    if batch_recipes:
        console.print("\n[cyan]==> finalize incremental multiomics[/cyan]")
        finalize_multiomics_parts(
            cohort_dir,
            recipes=batch_recipes,
            recipe_options=spec.recipe_options,
        )

    for recipe_name in spec.recipes:
        if recipe_name in MULTIOMICS_RECIPE_NAMES:
            continue
        recipe = RECIPE_REGISTRY.get(recipe_name)
        if recipe is None:
            continue
        console.print(f"\n[cyan]==> recipe: {recipe_name}[/cyan]")
        recipe(cohort_dir, spec.recipe_options)

    return cohort_dir


def _delete_batch_raw(records: list[dict], *, recipes: list[str]) -> None:
    from .multiomics import record_handled_by_multiomics

    for record in records:
        local_path = record.get("local_path")
        if not local_path or not record_handled_by_multiomics(record, recipes):
            continue
        path = Path(local_path)
        if path.exists():
            path.unlink()
        record["local_path"] = None
        record["status"] = "processed_deleted"


def _recipe_variants(cohort_dir: Path, recipe_options: dict[str, Any] | None = None) -> None:
    from .variants_polars import write_variants

    write_variants(cohort_dir)


def _recipe_samples(cohort_dir: Path, recipe_options: dict[str, Any] | None = None) -> None:
    from .samples_polars import write_samples

    write_samples(cohort_dir)


def _recipe_frequency(cohort_dir: Path, recipe_options: dict[str, Any] | None = None) -> None:
    from .frequency import write_gene_frequency, write_variant_frequency

    write_gene_frequency(cohort_dir)
    write_variant_frequency(cohort_dir)


def _recipe_rna_expression(cohort_dir: Path, recipe_options: dict[str, Any] | None = None) -> None:
    from .multiomics import write_rna_expression

    write_rna_expression(cohort_dir, recipe_options)


def _recipe_mirna_expression(
    cohort_dir: Path, recipe_options: dict[str, Any] | None = None
) -> None:
    from .multiomics import write_mirna_expression

    write_mirna_expression(cohort_dir, recipe_options)


def _recipe_methylation(cohort_dir: Path, recipe_options: dict[str, Any] | None = None) -> None:
    from .multiomics import write_methylation_beta

    write_methylation_beta(cohort_dir, recipe_options)


def _recipe_copy_number(cohort_dir: Path, recipe_options: dict[str, Any] | None = None) -> None:
    from .multiomics import write_copy_number

    write_copy_number(cohort_dir, recipe_options)


def _recipe_protein_expression(
    cohort_dir: Path, recipe_options: dict[str, Any] | None = None
) -> None:
    from .multiomics import write_protein_expression

    write_protein_expression(cohort_dir, recipe_options)


def _recipe_multiomics(cohort_dir: Path, recipe_options: dict[str, Any] | None = None) -> None:
    from .multiomics import write_multiomics

    write_multiomics(cohort_dir, recipe_options)


# Name → function. Names must match KNOWN_RECIPES in config.py.
RECIPE_REGISTRY: dict[str, Callable[[Path, dict[str, Any] | None], None]] = {
    "variants": _recipe_variants,
    "samples": _recipe_samples,
    "frequency": _recipe_frequency,
    "rna_expression": _recipe_rna_expression,
    "mirna_expression": _recipe_mirna_expression,
    "methylation": _recipe_methylation,
    "copy_number": _recipe_copy_number,
    "protein_expression": _recipe_protein_expression,
    "multiomics": _recipe_multiomics,
}


def _human_size(n: int) -> str:
    units = ["B", "KB", "MB", "GB", "TB"]
    s = float(n)
    for u in units:
        if s < 1024:
            return f"{s:.1f} {u}"
        s /= 1024
    return f"{s:.1f} PB"
