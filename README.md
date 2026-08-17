# tcga-pull

CLI that pulls cohorts from the NCI Genomic Data Commons into parquet tables.

## Install

```sh
uv tool install git+https://github.com/vgainullin/tcga-pull
```

Puts `tcga-pull` on your PATH. To hack on the code instead, see
[CONTRIBUTING.md](CONTRIBUTING.md).

Most open-access pulls (MAFs, clinical, small derived files) work through the
GDC bulk API and need nothing else. If your filter resolves to files over
100 MB (BAMs, segmented archives), install the official NCI binary:

```sh
tcga-pull install-gdc-client
```

Detects your platform, downloads the prebuilt binary from gdc.cancer.gov,
verifies it against a pinned MD5, and drops it in `~/.local/bin/`. The pull
pipeline picks it up automatically.

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
variant_summary.parquet
gene_frequency.parquet      variant_frequency.parquet
```

The full GDC filter grammar, the YAML schema, and a 55-project pancancer
worked example are in [`examples/pancancer_snv.yaml`](examples/pancancer_snv.yaml).
That pancancer example also declares optional omics add-ons that can be pulled
separately and merged downstream by `case_id` or `submitter_id`:

```sh
tcga-pull omics examples/pancancer_snv.yaml
tcga-pull preview examples/pancancer_snv.yaml --omics rna_gene_expression_star_counts
tcga-pull overlap examples/pancancer_snv.yaml \
  --omics dna_methylation_beta --omics protein_expression_rppa \
  --json overlap.json --parquet overlap.parquet
tcga-pull pull examples/pancancer_snv.yaml --omics rna_gene_expression_star_counts
```

`overlap` queries metadata only. It reports files, unique cases, bytes, project
and sample-type breakdowns for each selection, plus pairwise and all-selected
case intersections. Machine-readable outputs include the resolved GDC filters
and UTC query timestamp.

To cap every modality to the same cases, create a shared case-set artifact from
the uncapped N-way intersection, then reuse it for each preview or pull:

```sh
tcga-pull case-set examples/pancancer_snv.yaml \
  --omics dna_methylation_beta --omics protein_expression_rppa \
  --limit-per-project 100 --out shared-cases.json
tcga-pull pull examples/pancancer_snv.yaml --case-set shared-cases.json
tcga-pull pull examples/pancancer_snv.yaml \
  --omics dna_methylation_beta --case-set shared-cases.json
tcga-pull pull examples/pancancer_snv.yaml \
  --omics protein_expression_rppa --case-set shared-cases.json
```

The artifact records selected modalities, resolved filters, candidate and
selected counts, per-project shortfalls, deterministic ordering, and exact case
IDs. Each resulting `cohort.json` records the artifact path and SHA-256. Without
`--case-set`, `limit.per_project` keeps its existing independent-per-cohort
behavior.

For one processed download containing SNVs plus all declared omics file types, use
[`examples/pancancer_multiomics.yaml`](examples/pancancer_multiomics.yaml):

```sh
tcga-pull preview examples/pancancer_multiomics.yaml
tcga-pull pull examples/pancancer_multiomics.yaml
```

That config uses incremental processing: each download batch is converted into
recipe-specific parquet parts, handled non-SNV raw files are deleted, and the
final top-level omics parquets are assembled from those parts. For custom
cohorts, the same behavior can be enabled with:

```yaml
processing:
  mode: incremental
  batch_size: 200
  delete_raw_after_processing: true
```

Recipe-specific reduction options can shrink processed outputs further:

```yaml
recipe_options:
  variants:
    genes_file: ./panels/cancer_genes.txt
  rna_expression:
    columns: [gene_id, gene_name, gene_type, unstranded]
    genes_file: ./panels/cancer_genes.txt
  methylation:
    probes_file: ./probe_panels/immune_probes.txt
  copy_number:
    outputs: [segments]  # choose from: segments, gene
```

RNA gene allowlists match gene symbols, versioned Ensembl IDs, or unversioned
Ensembl IDs. Variant allowlists match `Hugo_Symbol`. These options reduce rows
during processing; GDC still transfers complete per-sample source files.
`variant_summary.parquet` retains unfiltered per-case aliquot structure and
burden, so `samples.parquet` keeps its canonical whole-cohort semantics.
File-backed panels are resolved once before processing. Their effective sorted
values, resolved source path, and SHA-256 digest are recorded in `cohort.json`.

The same mode can be enabled from flags for ad hoc pulls:

```sh
tcga-pull pull cohort.yaml --incremental --processing-batch-size 200 --delete-raw-after-processing
```

Normalized omics rows retain their GDC file, case, sample, sample type, tissue
type, project, workflow, and platform identifiers. `manifest.parquet` also keeps
the complete case/sample associations as JSON when a file maps ambiguously.
This allows downstream tools to construct sample-level tumor/normal cohorts
without inferring sample identity from filenames or patient-level joins.

If the raw files are already downloaded, build just the non-SNV omics parquets
with:

```sh
tcga-pull multiomics ./cohorts/pancancer_multiomics
```

To merge processed modalities by both sample and biological entity, run:

```sh
tcga-pull integrate ./cohorts/pancancer_multiomics --config cohort.yaml
```

Direct gene-level sources (SNV, RNA, gene CNV, and RPPA) share the
`gene_symbol` namespace. Methylation is mapped into that namespace
automatically using the `platform` retained from GDC. tcga-pull supports HM27,
HM450, EPIC, and EPIC v2 and downloads the matching hg38/GENCODE v41
InfiniumAnnotation v8.1 table into the cohort's `annotations/` directory. The
source and normalized mapping are checksum-verified and recorded in both the
integration manifest and `cohort.json`. The annotation source follows the
[GDC improved probe annotation guidance](https://gdc.cancer.gov/content/improved-dna-methylation-array-probe-annotation)
and is pinned by content hash because its upstream download URL is mutable.

By default, methylation uses probes within 1,500 bases of a transcription start
site and averages their beta values per gene. Include annotated gene-body probes
when required:

```yaml
recipe_options:
  integrate:
    modalities:
      methylation_beta:
        annotation_regions: [promoter, body]
```

Each methylation join identifier must resolve to one array platform; filter the
cohort by `platform` or use an exact sample join if a case contains multiple
arrays. Unknown or missing platforms fail with an actionable error rather than
dropping probes.

For custom assays or annotations, `mapping.file` remains available as an
explicit override. Mapping tables may be TSV, CSV, or parquet. One-to-many
mappings are expanded; many source features mapping to the same sample/entity
are combined using the declared aggregation. Custom modalities use the same
general contract by declaring
`source_file`, `feature_column`, `entity_column`, `entity_type`, `value_column`,
and `aggregation`. The default join is case-level (`submitter_id`). Set
`integrate.join_on` to `sample_id` or `sample_submitter_id` when exact assay
sample alignment is required. No cross-namespace biological relationship is
inferred.

Apply a variant panel to an existing download with the same cohort config:

```sh
tcga-pull variants ./cohorts/pancancer_multiomics --config cohort.yaml
tcga-pull samples ./cohorts/pancancer_multiomics
```

To export case-aligned matrices for model training:

```sh
tcga-pull dataset ./cohorts/pancancer_multiomics \
  --label-column oncotree_tissue \
  --modality snv \
  --modality rna_expression \
  --modality methylation_beta \
  --min-class-count 20 \
  --max-features-per-modality 5000
```

This writes `model_dataset/` under the cohort with `samples.parquet`,
`feature_index.parquet`, one matrix parquet per modality, and a JSON manifest.
For strict tumor-only model inputs, include `sample_type: [Primary Tumor]` in
the cohort filter before download. Use `platform` when selections must use the
same assay, for example `platform: Illumina Human Methylation 450`.

## Python API

```python
from tcga_pull import load_cohort

cohort = load_cohort("./cohorts/brca_snv")
cohort.clinical
cohort.manifest
cohort.variants
cohort.variant_summary
cohort.samples
cohort.gene_frequency
cohort.rna_expression
cohort.methylation_beta
cohort.integrated
cohort.model_dataset
cohort.provenance
```

Optional recipe outputs return `None` if the recipe didn't run. Column-by-column
schemas in [SCHEMAS.md](SCHEMAS.md).

## Shipped recipes

| recipe | output | rows |
|---|---|---|
| `variants` | `variants.parquet`, `variant_summary.parquet` | filtered variant rows plus unfiltered per-case summary |
| `samples` | `samples.parquet` | one per case (clinical + tissue + burden) |
| `frequency` | `gene_frequency.parquet`, `variant_frequency.parquet` | per (gene or variant, tissue) |
| `rna_expression` | `rna_expression.parquet` | one per (sample file × gene) |
| `mirna_expression` | `mirna_expression.parquet` | one per (sample file × miRNA) |
| `methylation` | `methylation_beta.parquet` | one per (sample file × methylation probe) |
| `copy_number` | `copy_number_segments.parquet`, `gene_copy_number.parquet` | sample-file segment-level and gene-level CNV |
| `protein_expression` | `protein_expression.parquet` | one per (sample file × RPPA target) |
| `multiomics` | all non-SNV omics parquets above | batch processor |
| `integrate` | `integrated/*.parquet`, `integrated/manifest.json` | sample/entity-aligned values |
| `model_dataset` | `model_dataset/*.parquet`, `model_dataset/manifest.json` | case-aligned training matrices |

The pull / restructure / clinical / manifest layer is data-type-agnostic; new
recipes plug in through `pipeline.RECIPE_REGISTRY`.

## Coverage matrix

Audit current open-access GDC coverage before adding or validating recipes:

```sh
tcga-pull coverage --program TCGA --out-dir ./coverage
```

This writes `tcga_open_access_coverage_matrix.parquet` and
`tcga_open_access_coverage_matrix.md`. Each row is a project-level file group
by `data_category`, `data_type`, `data_format`, `experimental_strategy`, and
`workflow_type`, classified as:

- `supported` — tcga-pull has a normalizing parquet recipe.
- `raw_only` — files can be inventoried/downloaded, but no normalizing recipe
  exists yet.

## Scope

- Open-access GDC.
- Controlled-access data (raw BAMs, per-caller VCFs) and cross-source merges
  (ICGC, cBioPortal, depmap) are not implemented.

## Development

See [CONTRIBUTING.md](CONTRIBUTING.md).

## License

MIT — [LICENSE](LICENSE).
