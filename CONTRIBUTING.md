# Contributing

## Dev setup

```sh
uv sync --extra agent
uv run pre-commit install  # installs the lint+type hook
```

Python 3.10+ is supported; the prototype targets 3.13. CI runs on 3.10–3.13.

## Local checks

The same checks CI runs:

```sh
uv run ruff check
uv run ruff format --check
uv run mypy
uv run pytest                # offline tests
./scripts/verify.sh          # rungs 1–3 (hits the GDC API, no real download)
./scripts/verify.sh all      # rungs 1–5 (real download + OpenRouter agent)
```

`pre-commit` runs `ruff check --fix` and `ruff format` on staged files; CI
will reject anything that disagrees with what `ruff format` produces.

## Tests

- **Offline tests** (`tests/`) — run in CI, no network. Cover pure helpers
  (filters, slugify, MAF projection, clinical flattening, parquet round-trip).
- **Integration ladder** (`scripts/verify.sh`) — hits the live GDC API.
  Rungs are independent; gated by CLI args so a single failing rung doesn't
  mask later coverage.

Mark tests that need network or downloads:

```python
import pytest

@pytest.mark.network
def test_live_gdc(): ...

@pytest.mark.download
def test_full_pipeline(): ...
```

## Commits

- Small, focused commits. Description-first commit messages, no rambling
  bodies.
- Don't push commits that fail CI locally — run `uv run pytest` and
  `uv run ruff check` first.
- Don't commit large local artifacts (`cohorts/`, `*.parquet`, downloaded
  files). The `.gitignore` covers the common cases; double-check `git status`.

## PRs

- Reference an issue if one exists.
- If you add a feature, update `CHANGELOG.md` under `## [Unreleased]`.
- If you change a CLI surface, update `README.md`.
- If you add a new module, add tests in the same PR.
