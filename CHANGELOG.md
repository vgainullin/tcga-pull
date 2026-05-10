# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
  `variants.parquet` with project, diagnosis, chrom/pos/ref/alt, consequence,
  impact, VAF, caller agreement.
- Conversational agent over OpenRouter (`tcga-pull agent`) with five tools
  (`list_projects`, `search_fields`, `count_files`, `preview_clinical`,
  `download`) and a hard `questionary.confirm` gate before any download.
- Five-rung verification ladder (`scripts/verify.sh`).
