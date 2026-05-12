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

Post-processing **recipes** can be declared in the YAML and run as part of
`pull`, so one command produces the whole cohort data product. Currently
three recipes ship:

- **variants** — walks the cohort's MAFs and writes `variants.parquet`
  (~37 columns: locus, gene/transcript, consequence, support, caller
  agreement, gnomAD/COSMIC/SIFT/PolyPhen, and flags `is_coding` /
  `is_high_impact` / `is_rare` / `primary_aliquot`).
- **samples** — joins clinical + variants into `samples.parquet` (one row
  per case, with curated `lineage` tissue label and mutation burden).
- **frequency** — emits `gene_frequency.parquet` and
  `variant_frequency.parquet` with raw rate, vs-other-lineages, and
  vs-gnomAD comparators side by side.

Each recipe is also runnable standalone (`tcga-pull variants <dir>`,
`tcga-pull samples <dir>`, `tcga-pull frequency <dir>`) for cohorts that
weren't built with recipes declared.

## Reproducing the pancancer SNV cohort

A cohort is fully specified by a YAML file. `examples/pancancer_snv.yaml`
captures the 55-project pancancer somatic-SNV cohort (~19,552 cases, ~21,300
MAFs, ~4.1 M variants) plus the post-processing recipes that build
`variants.parquet`, `samples.parquet`, and the frequency tables. One command
rebuilds the whole thing:

```sh
tcga-pull pull examples/pancancer_snv.yaml
```

Resumable: a killed `pull` picks up from `_downloads/` on restart. Expected
wall time on the live GDC: ~60–75 min for the download, a few minutes for
each recipe.

For any other cohort, write your own YAML — see
[examples/pancancer_snv.yaml](examples/pancancer_snv.yaml) for the format.

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

# Multi-project pulls: repeat --project, point at a file, or use YAML
tcga-pull pull --project TCGA-BRCA --project TCGA-LUAD \
               --data-category "Simple Nucleotide Variation" --data-format MAF
tcga-pull pull --projects-file projects.txt \
               --data-category "Simple Nucleotide Variation" --data-format MAF

# Cap per project for quick prototypes (also available in YAML as `limit.per_project`)
tcga-pull pull examples/pancancer_snv.yaml --limit-per-project 20

# Or from a YAML config
tcga-pull pull cohort.yaml

# Post-processing recipes (also auto-runnable via the YAML `recipes:` section)
tcga-pull variants  /path/to/cohort_dir            # one row per (variant x tumor aliquot)
tcga-pull samples   /path/to/cohort_dir            # one row per case (lineage + burden)
tcga-pull frequency /path/to/cohort_dir            # per-(gene|variant, tissue) tables

# All recipes default to polars; --engine pandas keeps the reference path alive
tcga-pull variants /path/to/cohort_dir --engine pandas

# QC + bench
tcga-pull validate-mafs /path/to/cohort_dir        # parse every .maf.gz, report drift
tcga-pull bench         /path/to/cohort_dir        # pandas vs polars head-to-head

# Conversational mode (needs OPENROUTER_API_KEY)
tcga-pull agent
tcga-pull agent -q "BRCA primary tumor RNA-seq STAR counts"
```

### Python API

For downstream projects that consume the parquets:

```python
from tcga_pull import load_cohort

cohort = load_cohort("/path/to/cohort")
cohort.variants            # polars DataFrame — one row per (variant x tumor aliquot)
cohort.samples             # one row per case, with curated `lineage` (tissue)
cohort.gene_frequency      # one row per (gene, lineage); None if frequency wasn't run
cohort.provenance          # dict from cohort.json (resolved filter, counts, timestamps)
cohort.summary()           # cheap shape readout
```

Column schemas, semantics, and stability promises live in [SCHEMAS.md](SCHEMAS.md).

### YAML cohort spec

```yaml
name: pancancer_snv
# out_dir is optional — defaults to ./cohorts; CLI --out overrides.
filters:
  project:                 # any sugar field can be a list
    - TCGA-BRCA
    - TCGA-LUAD
    - TARGET-AML
  data_category: Simple Nucleotide Variation
  data_format: MAF
# Or pass-through raw GDC filter:
# gdc_filter: { op: "and", content: [ ... ] }

# Optional: cap each project at N cases (deterministic submitter_id sort).
# Useful for quick prototypes — the pancancer YAML at full size pulls 1.2 GB.
limit:
  per_project: 20

download:
  n_processes: 4

# Post-processing recipes that run after the pull completes, in order.
# Omit to do a pull-only run; the canonical pancancer spec runs all three.
recipes:
  - variants
  - samples
  - frequency
```

See [examples/pancancer_snv.yaml](examples/pancancer_snv.yaml) for the full
55-project pancancer spec.

## Output layout

```
cohorts/<name>/
  clinical.parquet            # one row per case, demographic + diagnoses flattened
  clinical_raw.jsonl          # full nested clinical records (lossless)
  manifest.parquet            # one row per file: local_path, md5, size, data_*
  manifest.tsv                # the input sent to gdc-client (provenance)
  cohort.json                 # resolved filter + counts + timestamps
  variants.parquet            # one row per somatic mutation       (variants recipe)
  samples.parquet             # one row per case, with burden       (samples recipe)
  gene_frequency.parquet      # one row per (gene, tissue)          (frequency recipe)
  variant_frequency.parquet   # one row per (variant, tissue)       (frequency recipe)
  data/
    <TCGA-XX-XXXX>/
      <data_category>/
        <file>
```

Column-by-column schemas in [SCHEMAS.md](SCHEMAS.md).

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
uv run pytest                                # offline tests only (default)
uv run pytest -m network                     # + live GDC API queries
uv run pytest -m "network or download"       # + live download (~30 s, ~3 MB)
```

CI runs the offline suite on every push across a 3.10–3.13 matrix, and runs
the network + download suite once on Ubuntu 3.13 after the offline suite
passes.

See [CONTRIBUTING.md](CONTRIBUTING.md) for the dev workflow.

## License

Apache License 2.0 — see [LICENSE](LICENSE).
