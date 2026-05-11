"""tcga-pull CLI."""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from . import pipeline
from .config import from_flags, load_yaml, read_projects_file
from .gdc import GDCClient

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
    out: Path,
    project: list[str] | None,
    projects_file: Path | None,
    data_type: list[str] | None,
    data_category: list[str] | None,
    data_format: list[str] | None,
    experimental_strategy: list[str] | None,
    workflow: list[str] | None,
    sample_type: list[str] | None,
    n_processes: int,
):
    if config:
        return load_yaml(config)

    # Merge --project (repeatable) with --projects-file (newline-separated)
    projects: list[str] = list(project or [])
    if projects_file is not None:
        projects.extend(read_projects_file(projects_file))
    # de-dupe while preserving order
    projects = list(dict.fromkeys(projects))

    if not projects:
        console.print(
            "[yellow]Need either a YAML config, --project, or --projects-file. "
            "Try `tcga-pull projects` to list available projects.[/yellow]"
        )
        raise typer.Exit(2)
    cohort_name = name or _default_name(projects, data_type)
    return from_flags(
        name=cohort_name,
        out_dir=out.expanduser().resolve(),
        project=projects,
        data_type=data_type,
        data_category=data_category,
        data_format=data_format,
        experimental_strategy=experimental_strategy,
        workflow=workflow,
        sample_type=sample_type,
        n_processes=n_processes,
    )


@app.command()
def pull(
    config: Path | None = typer.Argument(None, exists=True, dir_okay=False, help="Cohort YAML."),
    name: str | None = typer.Option(None, "--name"),
    out: Path = typer.Option(Path("./cohorts"), "--out", "-o"),
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
    yes: bool = typer.Option(False, "--yes", "-y", help="Skip download confirmation."),
) -> None:
    """Build a cohort from a YAML config or flags. Prompts before downloading."""
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
    )
    pipeline.run(spec, yes=yes, console=console)


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
    p = pipeline.fetch_preview(spec)
    pipeline.render_preview(p, console=console)


@app.command()
def projects(program: str = typer.Option("TCGA", "--program")) -> None:
    """List GDC projects (default: TCGA)."""
    client = GDCClient()
    rows = client.list_projects(program=program)
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


@app.command()
def variants(
    cohort_dir: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=False,
        help="A cohort dir produced by `tcga-pull pull` containing SNV MAFs.",
    ),
) -> None:
    """Aggregate per-case MAFs into <cohort>/variants.parquet (one row per variant)."""
    from .variants import write_variants

    out = write_variants(cohort_dir)
    import pandas as pd

    df = pd.read_parquet(out)
    console.print(f"[green]wrote[/green] {out}  ({len(df):,} variants x {df.shape[1]} cols)")
    if "primary_diagnosis" in df.columns:
        console.print("\n[dim]top primary_diagnosis:[/dim]")
        console.print(df["primary_diagnosis"].value_counts().head(5).to_string())
    if "variant_class" in df.columns:
        console.print("\n[dim]variant_class:[/dim]")
        console.print(df["variant_class"].value_counts().head(8).to_string())


@app.command()
def samples(
    cohort_dir: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=False,
        help="A cohort dir with clinical.parquet + variants.parquet.",
    ),
) -> None:
    """Build per-case samples.parquet (demographics + lineage + burden)."""
    from .samples import write_samples

    out = write_samples(cohort_dir)
    import pandas as pd

    df = pd.read_parquet(out)
    console.print(f"[green]wrote[/green] {out}  ({len(df):,} cases x {df.shape[1]} cols)")
    if "lineage" in df.columns:
        console.print("\n[dim]lineage (top 8):[/dim]")
        console.print(df["lineage"].value_counts().head(8).to_string())
    if "n_variants_total" in df.columns:
        bs = df["n_variants_total"].describe()
        console.print(
            f"\n[dim]burden (n_variants_total, primary aliquot only):[/dim] "
            f"median={int(bs['50%'])}  IQR={int(bs['25%'])}-{int(bs['75%'])}  "
            f"max={int(bs['max'])}"
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


def _default_name(project: list[str], data_type: list[str] | None) -> str:
    if len(project) == 1:
        p = project[0].lower().replace("tcga-", "")
    else:
        p = f"multi_{len(project)}_projects"
    if data_type:
        dt = data_type[0].lower().replace(" ", "_")
        return f"{p}_{dt}"
    return p


if __name__ == "__main__":
    app()
