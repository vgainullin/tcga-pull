"""End-to-end orchestration: filter -> preview -> download -> restructure -> parquets."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from rich.console import Console
from rich.table import Table

from .config import CohortSpec
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


def fetch_preview(spec: CohortSpec, client: GDCClient | None = None) -> Preview:
    client = client or GDCClient()
    flt = spec.resolve_filter()
    file_hits = client.fetch_files(flt)
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
    yes: bool = False,
    console: Console | None = None,
    client: GDCClient | None = None,
) -> Path:
    """Full pipeline. Returns the cohort directory."""
    console = console or Console()
    client = client or GDCClient()

    with console.status("[cyan]Querying GDC…[/cyan]"):
        preview = fetch_preview(spec, client=client)
    render_preview(preview, console=console)

    if preview.n_files == 0:
        console.print("[yellow]No files match. Refine the filter.[/yellow]")
        return spec.cohort_dir

    if not yes:
        import questionary

        if not questionary.confirm("Proceed with download?", default=True).ask():
            console.print("[yellow]Aborted.[/yellow]")
            return spec.cohort_dir

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
        recipe(cohort_dir)

    return cohort_dir


def _recipe_variants(cohort_dir: Path) -> None:
    from .variants_polars import write_variants

    write_variants(cohort_dir)


def _recipe_samples(cohort_dir: Path) -> None:
    from .samples_polars import write_samples

    write_samples(cohort_dir)


def _recipe_frequency(cohort_dir: Path) -> None:
    from .frequency import write_gene_frequency, write_variant_frequency

    write_gene_frequency(cohort_dir)
    write_variant_frequency(cohort_dir)


# Name → function. Names must match KNOWN_RECIPES in config.py.
RECIPE_REGISTRY: dict[str, Callable[[Path], None]] = {
    "variants": _recipe_variants,
    "samples": _recipe_samples,
    "frequency": _recipe_frequency,
}


def _human_size(n: int) -> str:
    units = ["B", "KB", "MB", "GB", "TB"]
    s = float(n)
    for u in units:
        if s < 1024:
            return f"{s:.1f} {u}"
        s /= 1024
    return f"{s:.1f} PB"
