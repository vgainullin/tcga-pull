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
unified `variants.parquet` (~37 columns) with locus, gene/transcript,
consequence, support, caller agreement, gnomAD/COSMIC/SIFT/PolyPhen
annotations, and convenience flags (`is_coding`, `is_high_impact`,
`is_rare`, `primary_aliquot`). Ready for analysis with pandas / DuckDB /
Polars without any further wrangling.

## Reproducing the pancancer SNV cohort

A cohort is fully specified by a YAML file. `examples/pancancer_snv.yaml`
captures the 55-project pancancer somatic-SNV cohort (~19,552 cases, ~21,300
MAFs, ~4.1 M variants) plus the post-processing recipes that build
`variants.parquet`, `samples.parquet`, and the frequency tables. One command
rebuilds the whole thing:

```sh
tcga-pull pull examples/pancancer_snv.yaml --yes
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

# Or from a YAML config
tcga-pull pull cohort.yaml

# Post-process: aggregate MAFs into a single variants table
tcga-pull variants /path/to/cohort_dir

# And then build a per-case table (demographics + lineage + burden)
tcga-pull samples /path/to/cohort_dir

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
out_dir: ./cohorts
filters:
  project:                 # any sugar field can be a list
    - TCGA-BRCA
    - TCGA-LUAD
    - TARGET-AML
  data_category: Simple Nucleotide Variation
  data_format: MAF
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
  samples.parquet       # (after `tcga-pull samples`) one row per case, with burden
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
