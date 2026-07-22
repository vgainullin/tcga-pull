# Cohort parquet schemas

This is the stable contract between `tcga-pull` and downstream consumers. All
parquets live at the top of the cohort directory. Column types are the values
polars writes (Int64, Float64, Utf8, Boolean). Pandas reads them back with the
expected nullable dtypes.

```
cohort/
  clinical.parquet           # one row per case
  clinical_raw.jsonl         # full nested GDC clinical (lossless)
  manifest.parquet           # one row per file
  manifest.tsv               # the manifest sent to gdc-client (provenance)
  samples.parquet            # one row per case, with derived tissue + burden
  variants.parquet           # one row per (variant × tumor aliquot)
  gene_frequency.parquet     # optional: one row per (gene × tissue)
  variant_frequency.parquet  # optional: one row per (variant × tissue)
  rna_expression.parquet     # optional: one row per (case × gene)
  mirna_expression.parquet   # optional: one row per (case × miRNA)
  methylation_beta.parquet   # optional: one row per (case × methylation probe)
  copy_number_segments.parquet # optional: one row per CNV segment
  gene_copy_number.parquet   # optional: one row per (case × gene CNV)
  protein_expression.parquet # optional: one row per (case × RPPA target)
  model_dataset/             # optional: case-aligned training matrices
  cohort.json                # resolved filter + counts + timestamp
  data/<submitter_id>/<data_category>/<file>   # raw downloaded files
```

For the Python API: `from tcga_pull import load_cohort`.

Multiomics recipe outputs all include these provenance columns:

```
case_id, submitter_id, file_id, file_name, data_type,
experimental_strategy, workflow_type
```

Additional columns:

| parquet | additional columns |
|---|---|
| `rna_expression.parquet` | `gene_id`, `gene_name`, `gene_type`, `unstranded`, `stranded_first`, `stranded_second`, `tpm_unstranded`, `fpkm_unstranded`, `fpkm_uq_unstranded` |
| `mirna_expression.parquet` | `mirna_id`, `read_count`, `reads_per_million_mirna_mapped`, `cross_mapped` |
| `methylation_beta.parquet` | `probe_id`, `beta_value` |
| `copy_number_segments.parquet` | `cnv_type`, `sample`, `chrom`, `start`, `end`, `num_probes`, `segment_mean` |
| `gene_copy_number.parquet` | `gene_id`, `gene_name`, `chrom`, `start`, `end`, `copy_number` |
| `protein_expression.parquet` | `protein_id`, `gene_symbol`, `antibody`, `expression_value` |

When `processing.mode: incremental` is used, these parquet outputs may be
directory-style parquet datasets with `part-*.parquet` files under the named
`*.parquet` path. The Python API reads both file and directory parquet outputs.
Recipe-specific options may reduce columns or rows, for example RNA column
selection, methylation probe allowlists, or choosing only segment-level CNV.

---

## `clinical.parquet`

One row per case. Demographics + first-of-list flatten from the GDC clinical
nested structure (full nested form preserved in `clinical_raw.jsonl`).

Columns produced by `flatten_case()`:

| column | type | notes |
|---|---|---|
| `case_id` | Utf8 | GDC UUID, primary key |
| `submitter_id` | Utf8 | TCGA-style barcode or program-specific ID |
| `project_id` | Utf8 | e.g. `TCGA-BRCA` |
| `disease_type` | Utf8 | GDC broad histology |
| `primary_site` | Utf8 | ICD-O-3 topographical site |
| `demographic_*` | mixed | one column per scalar field in GDC's `demographic` block (gender, race, ethnicity, vital_status, days_to_death, year_of_birth, …) |
| `diagnosis_*` | mixed | one column per scalar field in the **first** diagnosis (primary_diagnosis, age_at_diagnosis, ajcc_pathologic_stage, days_to_last_follow_up, …) |
| `n_diagnoses` | Int64 | total count of diagnoses (multi-primary cases) |
| `treatment_*` | mixed | first treatment from the first diagnosis |
| `n_treatments` | Int64 | total treatments across diagnoses |
| `exposure_*` | mixed | first exposure record |
| `n_exposures` | Int64 |  |
| `family_history_*` | mixed | first family history record |
| `n_family_histories` | Int64 |  |

Non-scalar fields beyond the first list element are not flattened here. Use
`clinical_raw.jsonl` if you need them.

---

## `manifest.parquet`

One row per downloaded file.

| column | type | notes |
|---|---|---|
| `file_id` | Utf8 | GDC UUID |
| `file_name` | Utf8 | original filename |
| `case_id` | Utf8 | the case this file belongs to (null for multi-case files) |
| `submitter_id` | Utf8 | matches the per-case folder name on disk |
| `data_type` | Utf8 | e.g. `Masked Somatic Mutation` |
| `data_format` | Utf8 | e.g. `MAF`, `VCF`, `TSV` |
| `experimental_strategy` | Utf8 | e.g. `WXS`, `RNA-Seq` |
| `workflow_type` | Utf8 | analysis workflow name |
| `md5sum` | Utf8 | from GDC |
| `file_size` | Int64 | bytes |
| `local_path` | Utf8 | absolute path on disk after `restructure()`; null if file was routed to `_multi/` |
| `status` | Utf8 | `ok`, `missing` (rare; failed restructure) |

---

## `samples.parquet`

One row per case (matches `clinical.parquet` cardinality). The canonical
per-patient view, with curated tissue + per-patient mutation burden.

```
identifiers + raw GDC labels:
  case_id, submitter_id, project_id, program
  primary_site, disease_type, primary_diagnosis

curated lineage:
  lineage             — tissue label (breast, lung, pancreas, …); see src/tcga_pull/tissue.py
  oncotree_code       — MSKCC OncoTree code, e.g. BRCA, LUAD, AML; see src/tcga_pull/oncotree.py
  oncotree_name       — OncoTree diagnosis name, e.g. "Invasive Breast Carcinoma"
  oncotree_main_type  — OncoTree mainType, e.g. "Breast Cancer", "Non-Small Cell Lung Cancer"
  oncotree_tissue     — OncoTree tissue group, e.g. "Breast", "Lung"
                        Heterogeneous projects (CPTAC, HCMI, EXCEPTIONAL_RESPONDERS) → "OTHER".
                        Projects outside the curated map → all four columns null.

demographics:
  gender, race, ethnicity
  age_at_diagnosis_days, age_at_diagnosis_years
  vital_status, days_to_death, days_to_last_followup, ajcc_pathologic_stage

pair structure (derived from variants.parquet):
  n_tumor_aliquots, primary_tumor_barcode, primary_normal_barcode, normal_source

mutation burden (counted on primary aliquot only):
  n_variants_total, n_variants_coding, n_variants_high_impact
```

Burden columns count the **primary tumor aliquot only**, so a case with three
sequenced aliquots is still one row with one burden total. Use `n_variants_coding`
as the canonical "TMB-like" count.

---

## `variants.parquet`

One row per (variant × tumor aliquot). 37 columns. The big one — ~4 M rows
on a pancancer cohort.

```
identifiers (left-joined from clinical):
  project_id, case_id, submitter_id, primary_diagnosis

locus:
  chrom (Utf8, e.g. "chr1"), pos, end_pos, ref, alt, variant_type (SNP/INS/DEL)

gene / transcript / consequence:
  hugo_symbol, transcript_id, exon, hgvsc, hgvsp_short
  variant_class    — TCGA category (Missense_Mutation, Silent, Nonsense_Mutation, …)
  consequence      — VEP term  (missense_variant, stop_gained, splice_donor_variant, …)
  impact           — HIGH / MODERATE / LOW / MODIFIER

support:
  t_depth, t_alt_count, vaf  (= t_alt_count / t_depth)

caller agreement:
  callers          — semicolon-joined list (muse;mutect2;varscan2)
  n_callers        — count derived from `callers`

population / functional annotations:
  gnomad_af, cosmic_id, sift, polyphen, context (96-channel trinucleotide), hotspot

convenience flags:
  is_coding        — variant_class ∉ {Silent, Intron, IGR, 3'/5'UTR, Flank, RNA, Targeted_Region}
  is_high_impact   — impact ∈ {HIGH, MODERATE}
  is_rare          — gnomad_af.is_null() | (gnomad_af < 1e-3)

pair info:
  tumor_barcode, normal_barcode, normal_source     — Blood Derived Normal / Solid Tissue Normal / …
  primary_aliquot  — True for one tumor barcode per patient (highest mean t_depth; ties: lex)

provenance:
  source_file      — the MAF.gz this row came from
```

Filter to `primary_aliquot=True` for patient-level analysis to avoid
double-counting multi-aliquot cases. Filter to `is_coding & is_rare` for the
canonical "likely-somatic, gene-affecting" subset.

---

## `gene_frequency.parquet` (optional, produced by `tcga-pull frequency`)

One row per (gene × tissue). Filtering: `primary_aliquot ∧ is_coding ∧ is_rare`.

| column | type | notes |
|---|---|---|
| `hugo_symbol` | Utf8 | |
| `lineage` | Utf8 | tissue |
| `n_mutated_patients` | Int64 | unique patients in this tissue with ≥1 hit |
| `n_high_impact_patients` | Int64 | subset with HIGH/MODERATE impact |
| `n_total_patients` | Int64 | patients in the tissue |
| `freq` | Float64 | `n_mutated / n_total` |
| `freq_high_impact` | Float64 | same numerator restricted to high impact |
| `n_mutated_other`, `n_total_other`, `freq_other_lineages` | | counts collapsed across all other tissues |
| `log2_enrichment_vs_other` | Float64 | log2 ratio with smoothing ε = 1e-6 |
| `gnomad_max_af_in_gene` | Float64 | max gnomAD AF across observed variants in the gene (sparse) |
| `log2_enrichment_vs_gnomad` | Float64 | log2 (cohort freq / gnomad_max_af). Treat as advisory, not rigorous. |

---

## `variant_frequency.parquet` (optional, produced by `tcga-pull frequency`)

One row per (chrom × pos × ref × alt × tissue). No `is_rare` filter — keeps
population-common variants too so `gnomad_af` can drive the comparison.

```
locus + annotation:
  chrom, pos, ref, alt, hugo_symbol, hgvsp_short, variant_class, consequence, impact

lineage:
  lineage           — tissue label

cohort signal:
  n_patients_with_variant, n_total_patients, cohort_freq

vs other tissues:
  n_with_variant_other, n_total_other, freq_other_lineages
  log2_enrichment_vs_other

vs gnomAD population:
  gnomad_af
  log2_enrichment_vs_gnomad
```

---

## `model_dataset/` (optional, produced by `tcga-pull dataset` or recipe `model_dataset`)

Directory layout:

```
model_dataset/
  manifest.json
  samples.parquet
  feature_index.parquet
  snv.parquet                  # optional
  rna_expression.parquet       # optional
  methylation_beta.parquet     # optional
  gene_copy_number.parquet     # optional
  mirna_expression.parquet     # optional
  protein_expression.parquet   # optional
```

`samples.parquet` contains one row per selected case:

```
case_id, submitter_id, <label_column>, split, has_<modality>...
```

The default label is `oncotree_tissue`. Splits are deterministic, stratified by
label, and controlled by `seed`, `train_fraction`, `val_fraction`, and
`test_fraction`.

Each modality matrix contains `case_id`, `submitter_id`, then feature columns
named `<modality>__<feature>`. SNV features are binary rare-coding gene hits
from the primary aliquot. RNA and miRNA values are `log1p` transformed from
counts/RPM columns. Methylation, gene CNV, and protein values keep the processed
parquet values.

`feature_index.parquet` maps generated matrix columns back to source features:

| column | type | notes |
|---|---|---|
| `modality` | Utf8 | e.g. `snv`, `rna_expression` |
| `feature_id` | Utf8 | source gene/probe/miRNA/protein id |
| `feature_name` | Utf8 | currently same as `feature_id` |
| `column_name` | Utf8 | matrix column name |
| `value_column` | Utf8 | source value column or `value` for derived SNV |
| `transform` | Utf8 | `binary`, `log1p`, or `none` |
| `n_samples` | Int64 | selected cases with non-null source values |

`manifest.json` records options, label counts, split counts, emitted modality
paths, feature counts, and skipped modalities.

---

## Shared case-set JSON

`tcga-pull case-set ... --out shared-cases.json` writes a reusable selection
artifact. `cases` is the exact deterministic subset applied by
`preview` / `pull --case-set`.

```json
{
  "schema_version": 1,
  "cohort": "<source cohort name>",
  "created_at": "<ISO-8601 timestamp>",
  "selections": ["primary", "<optional_omics name>", "..."],
  "resolved_filters": {"<selection>": {<GDC filter JSON>}},
  "requested_per_project": <int>,
  "candidate_count": <int>,
  "selected_count": <int>,
  "candidate_counts_by_project": {"TCGA-LUAD": <int>},
  "selected_counts_by_project": {"TCGA-LUAD": <int>},
  "shortfalls_by_project": {
    "TCGA-LUAD": {"requested": <int>, "available": <int>}
  },
  "ordering": "project_id, submitter_id, case_id ascending",
  "cases": [
    {"case_id": "<GDC UUID>", "submitter_id": "<barcode>", "project_id": "TCGA-LUAD"}
  ]
}
```

---

## `cohort.json`

```json
{
  "name": "<cohort_name>",
  "filter": {<resolved GDC filter JSON, open-access wrapping applied>},
  "n_files": <int>,
  "n_cases": <int>,
  "total_size": <bytes>,
  "case_set": {                         // present only with --case-set
    "source": "<resolved artifact path>",
    "sha256": "<artifact digest>",
    "selection": "primary|<optional_omics name>",
    "selections": ["primary", "<optional_omics name>", "..."],
    "candidate_count": <int>,
    "selected_count": <int>,
    "requested_per_project": <int>,
    "ordering": "project_id, submitter_id, case_id ascending"
  },
  "created_at": "<ISO-8601 timestamp>"
}
```

The `filter` value is the exact JSON sent to the GDC `/files` endpoint, so this
file is enough to re-pull the same cohort (idempotent up to GDC data releases).
For shared selections, `case_set` identifies the reusable JSON artifact and
records its digest so downstream cohorts can prove they used the same case set.

---

## Stability promises

- **`load_cohort()` API** is part of the v0.x semver contract.
- **Column names + types** in the tables above are part of the contract too.
  Additive changes (new columns) are minor-version; removals/renames are
  major-version.
- **Lineage values** in `samples.parquet` (the strings `breast`, `lung`, …)
  may grow when OncoTree integration lands but won't be renamed.
- **`oncotree_*` columns** are populated from a pinned OncoTree snapshot
  (`src/tcga_pull/data/oncotree_2025_10_03.json`). Bumping the snapshot is a
  minor-version change; the column set is stable. Codes may change between
  OncoTree versions — treat the column type as Optional[Utf8].
