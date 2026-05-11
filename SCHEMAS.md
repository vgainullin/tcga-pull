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
  cohort.json                # resolved filter + counts + timestamp
  data/<submitter_id>/<data_category>/<file>   # raw downloaded files
```

For the Python API: `from tcga_pull import load_cohort`.

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
  oncotree_code       — reserved null; populated when OncoTree crosswalk lands
  oncotree_main_type  — reserved null
  oncotree_tissue     — reserved null

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

## `cohort.json`

```json
{
  "name": "<cohort_name>",
  "filter": {<resolved GDC filter JSON, open-access wrapping applied>},
  "n_files": <int>,
  "n_cases": <int>,
  "total_size": <bytes>,
  "created_at": "<ISO-8601 timestamp>"
}
```

The `filter` value is the exact JSON sent to the GDC `/files` endpoint, so this
file is enough to re-pull the same cohort (idempotent up to GDC data releases).

---

## Stability promises

- **`load_cohort()` API** is part of the v0.x semver contract.
- **Column names + types** in the tables above are part of the contract too.
  Additive changes (new columns) are minor-version; removals/renames are
  major-version.
- **Lineage values** in `samples.parquet` (the strings `breast`, `lung`, …)
  may grow when OncoTree integration lands but won't be renamed.
- **`oncotree_*` columns** are currently all null, reserved for the eventual
  OncoTree crosswalk. Downstream code should treat them as Optional[Utf8].
