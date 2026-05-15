"""tcga-pull CLI."""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from . import pipeline, services

app = typer.Typer(
    name="tcga-pull",
    help="Query the NCI GDC and produce a per-case structured TCGA data product.",
    no_args_is_help=True,
    add_completion=False,
)
console = Console()


def _spec_from_args(
    config: Path | None,
    *,
    name: str | None,
    out: Path | None,
    project: list[str] | None,
    projects_file: Path | None,
    data_type: list[str] | None,
    data_category: list[str] | None,
    data_format: list[str] | None,
    experimental_strategy: list[str] | None,
    workflow: list[str] | None,
    sample_type: list[str] | None,
    n_processes: int,
    limit_per_project: int | None = None,
):
    try:
        return services.build_cohort_spec(
            services.CohortBuildOptions(
                config=config,
                name=name,
                out=out,
                project=project,
                projects_file=projects_file,
                data_type=data_type,
                data_category=data_category,
                data_format=data_format,
                experimental_strategy=experimental_strategy,
                workflow=workflow,
                sample_type=sample_type,
                n_processes=n_processes,
                limit_per_project=limit_per_project,
            )
        )
    except services.SpecBuildError as e:
        console.print(f"[yellow]{e}[/yellow]")
        raise typer.Exit(2) from None


@app.command()
def pull(
    config: Path | None = typer.Argument(None, exists=True, dir_okay=False, help="Cohort YAML."),
    name: str | None = typer.Option(None, "--name"),
    out: Path | None = typer.Option(
        None,
        "--out",
        "-o",
        help="Output base dir. Overrides YAML out_dir. Defaults to ./cohorts in flag mode.",
    ),
    project: list[str] | None = typer.Option(None, "--project", help="e.g. TCGA-BRCA. Repeatable."),
    projects_file: Path | None = typer.Option(
        None,
        "--projects-file",
        exists=True,
        dir_okay=False,
        help="File with one project_id per line. # comments + blank lines allowed.",
    ),
    data_type: list[str] | None = typer.Option(None, "--data-type"),
    data_category: list[str] | None = typer.Option(None, "--data-category"),
    data_format: list[str] | None = typer.Option(None, "--data-format"),
    experimental_strategy: list[str] | None = typer.Option(None, "--strategy"),
    workflow: list[str] | None = typer.Option(None, "--workflow"),
    sample_type: list[str] | None = typer.Option(None, "--sample-type"),
    n_processes: int = typer.Option(4, "--n-processes", "-n"),
    limit_per_project: int | None = typer.Option(
        None,
        "--limit-per-project",
        help="Cap at N cases per GDC project (deterministic by submitter_id sort).",
    ),
) -> None:
    """Build a cohort from a YAML config or flags. Run `preview` first to dry-run."""
    spec = _spec_from_args(
        config,
        name=name,
        out=out,
        project=project,
        projects_file=projects_file,
        data_type=data_type,
        data_category=data_category,
        data_format=data_format,
        experimental_strategy=experimental_strategy,
        workflow=workflow,
        sample_type=sample_type,
        n_processes=n_processes,
        limit_per_project=limit_per_project,
    )
    services.run_cohort(spec, console=console)


@app.command()
def preview(
    config: Path | None = typer.Argument(None, exists=True, dir_okay=False),
    project: list[str] | None = typer.Option(None, "--project"),
    projects_file: Path | None = typer.Option(None, "--projects-file", exists=True, dir_okay=False),
    data_type: list[str] | None = typer.Option(None, "--data-type"),
    data_category: list[str] | None = typer.Option(None, "--data-category"),
    data_format: list[str] | None = typer.Option(None, "--data-format"),
    experimental_strategy: list[str] | None = typer.Option(None, "--strategy"),
    workflow: list[str] | None = typer.Option(None, "--workflow"),
    sample_type: list[str] | None = typer.Option(None, "--sample-type"),
) -> None:
    """Show what a filter would pull, without downloading."""
    spec = _spec_from_args(
        config,
        name="preview",
        out=Path("/tmp"),
        project=project,
        projects_file=projects_file,
        data_type=data_type,
        data_category=data_category,
        data_format=data_format,
        experimental_strategy=experimental_strategy,
        workflow=workflow,
        sample_type=sample_type,
        n_processes=4,
    )
    p = services.preview_cohort(spec)
    pipeline.render_preview(p, console=console)


@app.command()
def projects(program: str = typer.Option("TCGA", "--program")) -> None:
    """List GDC projects (default: TCGA)."""
    rows = services.list_projects(program=program)
    t = Table(title=f"Projects in {program}")
    t.add_column("project_id")
    t.add_column("name")
    t.add_column("primary_site")
    for r in sorted(rows, key=lambda x: x.get("project_id") or ""):
        site = r.get("primary_site") or []
        if isinstance(site, list):
            site = ", ".join(site)
        t.add_row(r.get("project_id", ""), r.get("name", ""), site)
    console.print(t)


@app.command()
def agent(
    query: str | None = typer.Option(None, "--query", "-q", help="One-shot NL query."),
    out: Path = typer.Option(Path("./cohorts"), "--out", "-o"),
) -> None:
    """Conversational cohort builder over OpenRouter (needs OPENROUTER_API_KEY)."""
    from .agent import run_agent

    run_agent(query=query, out=out, console=console)


def _resolve_variants_engine(engine: str):
    try:
        return services.resolve_variants_writer(engine)
    except ValueError:
        console.print(f"[red]unknown engine: {engine!r} (use 'polars' or 'pandas')[/red]")
        raise typer.Exit(2) from None


def _resolve_samples_engine(engine: str):
    try:
        return services.resolve_samples_writer(engine)
    except ValueError:
        console.print(f"[red]unknown engine: {engine!r} (use 'polars' or 'pandas')[/red]")
        raise typer.Exit(2) from None


@app.command()
def variants(
    cohort_dir: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=False,
        help="A cohort dir produced by `tcga-pull pull` containing SNV MAFs.",
    ),
    engine: str = typer.Option(
        "polars", "--engine", help="polars (default) or pandas. Both produce identical parquet."
    ),
) -> None:
    """Aggregate per-case MAFs into <cohort>/variants.parquet (one row per variant)."""
    import polars as pl

    _resolve_variants_engine(engine)
    out = services.write_variants_recipe(cohort_dir, engine=engine).path
    df = pl.read_parquet(out)
    console.print(
        f"[green]wrote[/green] {out}  ({len(df):,} variants x {df.width} cols, engine={engine})"
    )
    if "primary_diagnosis" in df.columns:
        console.print("\n[dim]top primary_diagnosis:[/dim]")
        console.print(df["primary_diagnosis"].value_counts(sort=True).head(5))
    if "variant_class" in df.columns:
        console.print("\n[dim]variant_class:[/dim]")
        console.print(df["variant_class"].value_counts(sort=True).head(8))


@app.command()
def samples(
    cohort_dir: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=False,
        help="A cohort dir with clinical.parquet + variants.parquet.",
    ),
    engine: str = typer.Option("polars", "--engine", help="polars (default) or pandas."),
) -> None:
    """Build per-case samples.parquet (demographics + lineage + burden)."""
    import polars as pl

    _resolve_samples_engine(engine)
    out = services.write_samples_recipe(cohort_dir, engine=engine).path
    df = pl.read_parquet(out)
    console.print(
        f"[green]wrote[/green] {out}  ({len(df):,} cases x {df.width} cols, engine={engine})"
    )
    if "lineage" in df.columns:
        console.print("\n[dim]lineage (top 8):[/dim]")
        console.print(df["lineage"].value_counts(sort=True).head(8))
    if "n_variants_total" in df.columns:
        stats = df.select(
            median=pl.col("n_variants_total").median(),
            q25=pl.col("n_variants_total").quantile(0.25),
            q75=pl.col("n_variants_total").quantile(0.75),
            smax=pl.col("n_variants_total").max(),
        ).row(0, named=True)
        console.print(
            f"\n[dim]burden (n_variants_total, primary aliquot only):[/dim] "
            f"median={int(stats['median'] or 0)}  "
            f"IQR={int(stats['q25'] or 0)}-{int(stats['q75'] or 0)}  "
            f"max={int(stats['smax'] or 0)}"
        )


@app.command("frequency")
def frequency_cmd(
    cohort_dir: Path = typer.Argument(
        ..., exists=True, file_okay=False, help="A cohort dir with variants + samples parquets."
    ),
    top_per_lineage: int = typer.Option(
        10, "--top", help="Print top N genes per lineage to stdout."
    ),
) -> None:
    """Build gene_frequency.parquet + variant_frequency.parquet with raw rate,
    vs-other-lineages, and vs-gnomAD comparators side by side."""
    import polars as pl

    outputs = services.write_frequency_recipe(cohort_dir)
    gp = outputs.gene_frequency
    vp = outputs.variant_frequency
    gf = pl.read_parquet(gp)
    vf = pl.read_parquet(vp)
    console.print(
        f"[green]wrote[/green] {gp}  ({len(gf):,} (gene, lineage) rows)\n"
        f"[green]wrote[/green] {vp}  ({len(vf):,} (variant, lineage) rows)"
    )

    # Top genes per lineage by raw frequency, restricted to lineages with ≥30 patients.
    if "lineage" not in gf.columns:
        return
    big_lineages = (
        gf.group_by("lineage")
        .agg(pl.col("n_total_patients").max().alias("n"))
        .filter(pl.col("n") >= 30)
        .sort("n", descending=True)["lineage"]
        .to_list()
    )
    for lin in big_lineages[:8]:
        top = (
            gf.filter(pl.col("lineage") == lin)
            .sort("freq", descending=True)
            .head(top_per_lineage)
            .select(["hugo_symbol", "n_mutated_patients", "n_total_patients", "freq"])
        )
        if top.is_empty():
            continue
        from rich.table import Table

        t = Table(title=f"{lin} — top {top_per_lineage} mutated genes (rare+coding)")
        t.add_column("gene")
        t.add_column("mutated", justify="right")
        t.add_column("of", justify="right")
        t.add_column("freq", justify="right")
        for row in top.to_dicts():
            t.add_row(
                row["hugo_symbol"],
                str(row["n_mutated_patients"]),
                str(row["n_total_patients"]),
                f"{row['freq']:.1%}",
            )
        console.print(t)


@app.command("bench")
def bench_cmd(
    cohort_dir: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=False,
        help="A cohort dir with raw MAFs already pulled.",
    ),
    save: Path | None = typer.Option(None, "--save", help="Write bench.json to this path."),
) -> None:
    """Benchmark pandas vs polars on the variants + samples pipelines.
    Reports wall time, peak RSS, and a cell-by-cell diff of the two outputs."""
    from rich.table import Table

    from .bench import run_bench, write_json

    result = run_bench(cohort_dir)

    timing = Table(title="bench results")
    timing.add_column("engine")
    timing.add_column("stage")
    timing.add_column("wall (s)", justify="right")
    timing.add_column("peak RSS (MB)", justify="right")
    timing.add_column("rows", justify="right")
    timing.add_column("cols", justify="right")
    timing.add_column("parquet (MB)", justify="right")
    for s in result.stages:
        timing.add_row(
            s.engine,
            s.stage,
            f"{s.wall_seconds:.2f}",
            f"{s.peak_rss_mb_after:,.0f}",
            f"{s.rows:,}",
            str(s.cols),
            f"{s.parquet_bytes / (1 << 20):.1f}",
        )
    console.print(timing)

    # Speedups
    by_key: dict[tuple[str, str], float] = {
        (s.engine, s.stage): s.wall_seconds for s in result.stages
    }
    for stage in ("variants", "samples"):
        if ("pandas", stage) in by_key and ("polars", stage) in by_key:
            t_pd = by_key[("pandas", stage)]
            t_pl = by_key[("polars", stage)]
            console.print(
                f"  [bold]{stage}[/bold]: polars is [green]{t_pd / t_pl:.1f}x[/green] faster"
            )

    # Diff summary
    for stage in ("variants", "samples"):
        d = result.diffs.get(stage, {})
        if d:
            console.print(f"\n[red]{stage}: {len(d)} columns differ[/red]")
            for col, n in sorted(d.items(), key=lambda x: -x[1])[:10]:
                console.print(f"  {col}: {n:,}")
        else:
            console.print(f"  [bold]{stage}[/bold] diff: [green]identical[/green]")

    if save:
        write_json(result, save)
        console.print(f"  wrote {save}")


@app.command("install-gdc-client")
def install_gdc_client_cmd(
    dest: Path = typer.Option(
        Path.home() / ".local" / "bin",
        "--dest",
        help="Install directory. Make sure it's on your PATH.",
    ),
    force: bool = typer.Option(False, "--force", help="Overwrite an existing gdc-client."),
) -> None:
    """Download, verify, and install the official NCI gdc-client binary.

    Only needed for cohorts pulling files larger than ~100 MB (BAMs, segmented
    archives). MAFs and clinical data go through the bulk API and don't need
    this.
    """
    from .install_gdc import GDC_CLIENT_VERSION, InstallError, install_gdc_client

    console.print(f"[dim]installing gdc-client v{GDC_CLIENT_VERSION} into {dest}...[/dim]")
    try:
        path = install_gdc_client(dest, force=force)
    except InstallError as e:
        console.print(f"[red]{e}[/red]")
        raise typer.Exit(1) from e

    console.print(f"[green]installed:[/green] {path}")
    import os

    on_path = any(Path(p) == path.parent for p in os.environ.get("PATH", "").split(os.pathsep))
    if not on_path:
        console.print(
            f"[yellow]warning:[/yellow] {path.parent} is not on your PATH. "
            f"Add it, or move/symlink {path.name} somewhere that is."
        )


@app.command("validate-mafs")
def validate_mafs_cmd(
    path: Path = typer.Argument(
        ...,
        exists=True,
        help="Directory to recurse through, or a single .maf.gz file.",
    ),
) -> None:
    """Read every .maf.gz under PATH; report parse errors, schema coverage,
    row totals, and program breakdown. Does not write anything."""
    from .validate import find_mafs, render_report, validate_mafs

    mafs = [path] if path.is_file() else find_mafs(path)
    if not mafs:
        console.print(f"[yellow]no .maf.gz files under {path}[/yellow]")
        raise typer.Exit(2)
    console.print(f"[dim]scanning {len(mafs):,} MAFs under {path}...[/dim]")
    report = validate_mafs(mafs)
    render_report(report, console=console)
    if report.files_failed or report.missing_required_columns:
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
