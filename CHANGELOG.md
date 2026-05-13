# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `tcga-pull install-gdc-client` — downloads the official NCI binary, verifies
  against pinned per-platform MD5, drops it in `~/.local/bin/`. Replaces the
  micromamba-based install instructions in the README.
- OncoTree crosswalk (`tcga_pull.oncotree`): GDC project_id → MSKCC OncoTree
  code/name/mainType/tissue. Bundles the `oncotree_2025_10_03` snapshot
  (898 codes) and a curated 55-project mapping. Validates at import that every
  mapped code exists in the snapshot. The `oncotree_*` columns on
  `samples.parquet` are no longer null placeholders — they're populated from
  the crosswalk. Heterogeneous projects (CPTAC, HCMI) map to `OTHER`.

## [0.1.2] — 2026-05-12

### Changed
- Relicensed from Apache-2.0 to MIT.
- README rewritten in terse, factual style.

### Fixed
- `pyproject.toml` version now matches `__init__.py` (the v0.1.1 release bumped
  only the latter).

## [0.1.1] — 2026-05-12

### Added
- YAML `recipes:` section runs variants / samples / frequency after `pull`.
- YAML `limit.per_project` and CLI `--limit-per-project` for cohort sampling.
- `load_cohort()` public API + `SCHEMAS.md`.
- CI `e2e` job: live `@network` + `@download` pytest tests on every push.

### Changed
- `polars` promoted from optional `[engines]` extra to core dep; `duckdb`
  stays in (renamed) `[queries]` extra.
- `tcga-pull pull` no longer prompts; the preview is the consent gate.
  `--yes` flag removed (use `tcga-pull preview` for dry-run).
- CLI `--out` now overrides the YAML's `out_dir` when explicitly passed.

### Removed
- `scripts/verify.sh` — replaced by pytest `@network` / `@download` markers.

## [0.1.0] — 2026-05-11

First tagged release. Cohort artefacts and the `load_cohort` Python API are
covered by semver from here.

### Added
- Initial release.
- GDC API client with filter helpers (`f_in`, `f_and`, `f_or`, `open_access`)
  and a `for_cases_endpoint` filter translator (`/files` paths → `/cases`).
- Two download backends: bulk `/data` endpoint (fast for many small files) and
  `gdc-client` subprocess (segmented + resumable for files > 100 MB). The
  pipeline picks automatically based on per-file size.
- Per-case restructure: flat `_downloads/<file_id>/<file>` →
  `data/<case_submitter_id>/<data_category_slug>/<file>`.
- Tidy clinical flattening into `clinical.parquet` (one row per case) plus
  lossless `clinical_raw.jsonl`.
- `manifest.parquet` with `local_path`, `md5`, `file_size`, data-type metadata.
- `tcga-pull variants` recipe: aggregates per-case MAFs into a unified
  `variants.parquet`. Schema covers locus (chrom, pos, end_pos, ref, alt,
  variant_type), gene/transcript (hugo_symbol, transcript_id, exon, hgvsc,
  hgvsp_short), consequence (variant_class, consequence, impact), support
  (t_depth, t_alt_count, vaf, callers, n_callers), and annotations (gnomad_af,
  cosmic_id, sift, polyphen, context, hotspot).
- Convenience flag columns on variants: `is_coding`, `is_high_impact`,
  `is_rare` (gnomAD AF < 1e-3 or null), `normal_source` (parsed from TCGA
  barcode position 13-14), and `primary_aliquot` — True for one tumor
  barcode per patient, picked by highest mean t_depth; lets you filter to
  patient-level analysis without double-counting multi-aliquot cases.
- `tcga-pull samples` recipe: builds `samples.parquet` (one row per case)
  from `clinical.parquet` + `variants.parquet`. Carries identifiers,
  raw GDC labels (primary_site, disease_type, primary_diagnosis),
  `lineage` (defaults to project_id; OncoTree integration deferred),
  demographics, pair structure (n_tumor_aliquots, primary tumor/normal
  barcode, normal_source), and patient-level mutation burden
  (n_variants_total/coding/high_impact, computed on the primary aliquot
  only to avoid double-counting).
- Multi-project pulls: `--project` is repeatable, `--projects-file` reads
  newline-separated project IDs (# comments + blank lines ignored), and
  YAML accepts a list for any sugar filter (e.g. `project: [TCGA-BRCA, TCGA-LUAD]`).
- `--data-format` flag added to `pull` and `preview` (was missing from the
  CLI surface even though the sugar mapping supported it).

### Fixed
- Faceted queries on the GDC `/files` endpoint no longer 500. Root cause was
  nested `AND` clauses in our filter tree — `f_and` now auto-flattens, so a
  filter combined with any nested facet works reliably. Added an offline
  unit test for the flattening and a network-tagged regression test that
  asserts the GDC API returns 200 for a faceted open-access SNV query.
- Internal convention: file-rooted filter fields now use bare names
  (`access`, `data_type`, `data_format`, …) instead of the `files.` prefix.
  `for_cases_endpoint` re-prefixes when translating to /cases. Handwritten
  `files.X` filters still pass through unchanged.
- `bulk_download_via_api` is now resilient: transient 5xx / connection /
  timeout / tar errors retry with exponential backoff (5 tries, capped at
  60s). 4xx errors fail fast. Restarting a partial pull skips files already
  on disk so an interrupted run resumes from where it stopped.

### Added (continued)
- `tcga-pull validate-mafs <path>` subcommand: walks a directory recursively
  (or single file), reads every `*.maf.gz` through `read_maf`, and reports
  parse errors, schema completeness against `REQUIRED_OUTPUT_COLUMNS`,
  optional-column coverage, row totals, and a program breakdown inferred
  from `tumor_barcode` prefixes. Designed to catch schema drift between
  TCGA / TARGET / CPTAC / MMRF / CMI / etc. before the parquet write.
  Backed by a new `validate.py` module + `aggregate_mafs(paths)` building
  block factored out of `aggregate_cohort`.
- Conversational agent over OpenRouter (`tcga-pull agent`) with five tools
  (`list_projects`, `search_fields`, `count_files`, `preview_clinical`,
  `download`) and a hard `questionary.confirm` gate before any download.
- Five-rung verification ladder (`scripts/verify.sh`).
