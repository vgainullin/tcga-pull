# tcga-pull

Query the NCI Genomic Data Commons (GDC), download cohorts, and produce
per-case structured TCGA data products.

> Status: **alpha** — interface is unstable, breaking changes likely.

## What it does

Given a filter over the GDC catalogue, `tcga-pull`:

1. Resolves the filter against the live GDC API and shows you counts + size
   before any bytes move.
2. Downloads the matching files via the bulk `/data` endpoint (fast for many
   small files) or `gdc-client` (for files > 100 MB; segmented + resumable).
3. Restructures the flat output into `data/<case_submitter_id>/<data_category>/`
   so files are grouped per individual rather than by opaque file UUID.
4. Flattens nested clinical records (`demographic`, `diagnoses`,
   `exposures`, …) into a tidy `clinical.parquet` (one row per case).
5. Writes `manifest.parquet` (one row per file with `local_path` + provenance).

A `variants` post-processing recipe walks the cohort's MAFs and produces a
unified `variants.parquet` with `project_id, case_id, primary_diagnosis,
chrom, pos, ref, alt, variant_class, consequence, vaf, …` — ready for
analysis with pandas / DuckDB / Polars without any further wrangling.

## Install

```sh
git clone <repo>
cd tcga-pull
uv sync                # base
uv sync --extra agent  # base + conversational agent over OpenRouter
```

You also need the GDC Data Transfer Tool for files > 100 MB; small-file
cohorts use the bulk API and don't need it.

```sh
# Easiest on macOS with micromamba:
brew install micromamba
micromamba create -y -n gdc -c bioconda gdc-client
ln -sf "$(micromamba env list | awk '/gdc / {print $NF}')/bin/gdc-client" \
       ~/.local/bin/gdc-client
gdc-client --version   # 2.3+
```

## Usage

```sh
# Browse projects
tcga-pull projects

# Dry-run: counts and breakdown, no download
tcga-pull preview --project TCGA-CHOL --data-type "Gene Expression Quantification" \
                  --workflow "STAR - Counts"

# Pull a cohort (with confirmation)
tcga-pull pull --project TCGA-CHOL --data-type "Gene Expression Quantification" \
               --workflow "STAR - Counts" --name chol_rnaseq

# Or from a YAML config
tcga-pull pull cohort.yaml

# Post-process: aggregate MAFs into a single variants table
tcga-pull variants /path/to/cohort_dir

# Conversational mode (needs OPENROUTER_API_KEY)
tcga-pull agent
tcga-pull agent -q "BRCA primary tumor RNA-seq STAR counts"
```

### YAML cohort spec

```yaml
name: brca_snv
out_dir: ./cohorts
filters:
  project: TCGA-BRCA
  data_category: Simple Nucleotide Variation
# Or pass-through raw GDC filter:
# gdc_filter: { op: "and", content: [ ... ] }
download:
  n_processes: 4
```

## Output layout

```
cohorts/<name>/
  clinical.parquet      # one row per case, demographic + diagnoses flattened
  clinical_raw.jsonl    # full nested clinical records (lossless)
  manifest.parquet      # one row per file: local_path, md5, size, data_*
  manifest.tsv          # the input sent to gdc-client (provenance)
  cohort.json           # resolved filter + counts + timestamps
  variants.parquet      # (after `tcga-pull variants`) one row per somatic mutation
  data/
    <TCGA-XX-XXXX>/
      <data_category>/
        <file>
```

## Scope

- **Open-access only.** Controlled-access data (raw BAMs, per-caller VCFs)
  needs dbGaP authorization and a token — not currently supported.
- **GDC only.** ICGC / cBioPortal / depmap support is on the roadmap; the
  filter and download paths are designed to be source-agnostic.

## Development

```sh
uv sync --extra agent
uv run ruff check
uv run ruff format --check
uv run mypy
uv run pytest

# Full integration ladder (some rungs hit the GDC API, one does a real download)
./scripts/verify.sh        # rungs 1-3 (no real download)
./scripts/verify.sh all    # rungs 1-5 (needs OPENROUTER_API_KEY for rung 5)
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for the dev workflow.

## License

Apache License 2.0 — see [LICENSE](LICENSE).
