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
uv run pytest                              # offline tests only (default)
uv run pytest -m network                   # + live GDC API queries
uv run pytest -m "network or download"     # + live download (~30 s, ~3 MB)
```

`pre-commit` runs `ruff check --fix` and `ruff format` on staged files; CI
will reject anything that disagrees with what `ruff format` produces.

## Tests

- **Offline tests** (`tests/test_*.py` without `@pytest.mark.network`) — run
  in CI on every push across the 3.10–3.13 matrix. Cover pure helpers,
  schema parity (pandas vs polars), and synthetic-cohort recipe math.
- **Network tests** (`@pytest.mark.network`) — hit the live GDC API for
  metadata. Fast (<5 s), no downloads. Run by CI's `e2e` job.
- **Download tests** (`@pytest.mark.download`) — pull a tiny live cohort
  through the full pull → recipes pipeline and assert artefacts.
  ~30 s wall, ~3 MB on the wire. Also run by CI's `e2e` job.

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
