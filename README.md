# tcga-pull

CLI that pulls cohorts from the NCI Genomic Data Commons into parquet tables.

## Install

```sh
git clone https://github.com/vgainullin/tcga-pull && cd tcga-pull
uv sync
```

For downloads where individual files exceed 100 MB:

```sh
brew install micromamba
micromamba create -y -n gdc -c bioconda gdc-client
ln -sf "$(micromamba env list | awk '/gdc / {print $NF}')/bin/gdc-client" \
       ~/.local/bin/gdc-client
```

## Usage

Define a cohort:

```yaml
# brca_snv.yaml
name: brca_snv
filters:
  project: TCGA-BRCA
  data_category: Simple Nucleotide Variation
  data_format: MAF
recipes: [variants, samples, frequency]
```

Run it:

```sh
tcga-pull pull brca_snv.yaml
```

Output goes to `./cohorts/brca_snv/`:

```
clinical.parquet            manifest.parquet            cohort.json
clinical_raw.jsonl          manifest.tsv                data/<patient>/<data_category>/<file>
variants.parquet            samples.parquet
gene_frequency.parquet      variant_frequency.parquet
```

The full GDC filter grammar, the YAML schema, and a 55-project pancancer
worked example are in [`examples/pancancer_snv.yaml`](examples/pancancer_snv.yaml).

## Python API

```python
from tcga_pull import load_cohort

cohort = load_cohort("./cohorts/brca_snv")
cohort.clinical
cohort.manifest
cohort.variants
cohort.samples
cohort.gene_frequency
cohort.provenance
```

Recipe outputs (`variants`, `samples`, `gene_frequency`, `variant_frequency`)
return `None` if the recipe didn't run. Column-by-column schemas in
[SCHEMAS.md](SCHEMAS.md).

## Shipped recipes

| recipe | output | rows |
|---|---|---|
| `variants` | `variants.parquet` | one per (variant × tumor aliquot) |
| `samples` | `samples.parquet` | one per case (clinical + tissue + burden) |
| `frequency` | `gene_frequency.parquet`, `variant_frequency.parquet` | per (gene or variant, tissue) |

All three are somatic-SNV-specific. The pull / restructure / clinical /
manifest layer is data-type-agnostic; new recipes plug in through
`pipeline.RECIPE_REGISTRY`.

## Scope

- Open-access GDC.
- Controlled-access data (raw BAMs, per-caller VCFs) and cross-source merges
  (ICGC, cBioPortal, depmap) are not implemented.

## Development

See [CONTRIBUTING.md](CONTRIBUTING.md).

## License

MIT — [LICENSE](LICENSE).
