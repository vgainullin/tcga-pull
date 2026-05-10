#!/usr/bin/env bash
# tcga-pull verification ladder. Each rung must pass before the next.
#
# Usage:
#   scripts/verify.sh           # rungs 1..3 (no real downloads, no API key needed)
#   scripts/verify.sh all       # rungs 1..5 (downloads ~10MB, needs OPENROUTER_API_KEY)
#   scripts/verify.sh 1 2       # specific rungs
#
# Rung 4 actually downloads files via gdc-client. Rung 5 calls OpenRouter.
set -euo pipefail

cd "$(dirname "$0")/.."

GREEN=$'\033[32m'; RED=$'\033[31m'; DIM=$'\033[2m'; BOLD=$'\033[1m'; RESET=$'\033[0m'
ok()     { echo "  ${GREEN}✓${RESET} $*"; }
fail()   { echo "  ${RED}✗${RESET} $*"; exit 1; }
header() { echo; echo "${BOLD}═══ $* ═══${RESET}"; }
note()   { echo "  ${DIM}$*${RESET}"; }

WANT="${*:-1 2 3}"
if [[ "$WANT" == "all" ]]; then WANT="1 2 3 4 5"; fi
want() { [[ " $WANT " == *" $1 "* ]]; }

SMOKE_OUT=/tmp/tcga_pull_smoke
COHORT_NAME=smoke

# ------------------------------------------------------------------------------
# Rung 1 — GDC API client: filters, pagination, expected fields
# ------------------------------------------------------------------------------
if want 1; then
header "1/5  GDC API client"
uv run python - <<'PY' || fail "GDC client smoke test failed (see traceback above)"
import sys
from tcga_pull.gdc import GDCClient, f_in, f_and, open_access, for_cases_endpoint

c = GDCClient()
flt = open_access(f_and(
    f_in("cases.project.project_id", ["TCGA-CHOL"]),
    f_in("files.data_type", ["Gene Expression Quantification"]),
    f_in("analysis.workflow_type", ["STAR - Counts"]),
))
n_files = c.count_files(flt)
n_cases = c.count_cases(flt)
print(f"    files={n_files} cases={n_cases}")
assert n_files > 0, "expected files > 0 for TCGA-CHOL STAR Counts"
if n_cases == 0:
    # Diagnostic: dump the translated filter and the canonical /cases field schema
    # so we can see exactly why the join didn't match.
    import json
    print("\n--- /cases filter we sent (after translation): ---")
    print(json.dumps(for_cases_endpoint(flt), indent=2))
    mapping = c.mapping("cases")
    fields = mapping.get("fields", [])
    relevant = [f for f in fields if any(k in f for k in ("project", "workflow", "data_type", "sample_type", "access"))]
    print(f"\n--- relevant /cases fields ({len(relevant)} of {len(fields)}): ---")
    for f in sorted(relevant):
        print(f"  {f}")
    sys.exit(2)
assert n_cases > 0, "expected cases > 0"

hits = c.fetch_files(flt)
print(f"    fetched={len(hits)} (pagination terminated)")
assert len(hits) == n_files, f"pagination mismatch: fetched {len(hits)} != count {n_files}"

h = hits[0]
required = ["file_id", "file_name", "data_type", "md5sum", "file_size", "cases"]
missing = [k for k in required if not h.get(k)]
assert not missing, f"missing fields on first hit: {missing}"
case = (h.get("cases") or [{}])[0]
assert case.get("submitter_id", "").startswith("TCGA-"), f"unexpected submitter_id: {case.get('submitter_id')}"
print(f"    sample: {h['file_id'][:8]}… {h['file_name']} (case {case['submitter_id']})")
PY
ok "GDC client: count, fetch, pagination, expected fields"
fi

# ------------------------------------------------------------------------------
# Rung 2 — CLI wiring + sugar→GDC mapping
# ------------------------------------------------------------------------------
if want 2; then
header "2/5  CLI commands"
PROJECTS_OUT=$(uv run tcga-pull projects 2>&1) || fail "tcga-pull projects errored"
echo "$PROJECTS_OUT" | grep -q "TCGA-CHOL" || fail "tcga-pull projects didn't list TCGA-CHOL"
PROJ_COUNT=$(echo "$PROJECTS_OUT" | grep -c "TCGA-" || true)
ok "tcga-pull projects (found $PROJ_COUNT TCGA projects)"

PREVIEW_OUT=$(uv run tcga-pull preview --project TCGA-CHOL --data-type "Gene Expression Quantification" --workflow "STAR - Counts" 2>&1) \
    || fail "tcga-pull preview errored"
echo "$PREVIEW_OUT" | grep -qE "Files.*[1-9]" || fail "preview didn't show non-zero file count"
echo "$PREVIEW_OUT" | grep -qE "Cases.*[1-9]" || fail "preview didn't show non-zero case count"
ok "tcga-pull preview shows non-zero counts"
fi

# ------------------------------------------------------------------------------
# Rung 3 — gdc-client present and runnable
# ------------------------------------------------------------------------------
if want 3; then
header "3/5  gdc-client binary"
if ! command -v gdc-client >/dev/null 2>&1; then
    cat <<EOF
  ${RED}✗${RESET} gdc-client not on PATH.
    Install from: https://gdc.cancer.gov/access-data/gdc-data-transfer-tool
    Or try:       uv tool install gdc-client
EOF
    exit 1
fi
GDC_VER=$(gdc-client --version 2>&1 | head -1)
ok "gdc-client present: $GDC_VER"
fi

# ------------------------------------------------------------------------------
# Rung 4 — real end-to-end download (~10MB into $SMOKE_OUT)
# ------------------------------------------------------------------------------
if want 4; then
header "4/5  end-to-end download (~10MB)"
note "downloading TCGA-CHOL STAR-Counts into $SMOKE_OUT/$COHORT_NAME/"
rm -rf "$SMOKE_OUT"
uv run tcga-pull pull \
    --project TCGA-CHOL \
    --data-type "Gene Expression Quantification" \
    --workflow "STAR - Counts" \
    --name "$COHORT_NAME" \
    -o "$SMOKE_OUT" \
    --yes \
    || fail "download pipeline errored"

COHORT="$SMOKE_OUT/$COHORT_NAME"
test -f "$COHORT/clinical.parquet"  || fail "missing clinical.parquet"
test -f "$COHORT/manifest.parquet"  || fail "missing manifest.parquet"
test -f "$COHORT/cohort.json"       || fail "missing cohort.json"
test -d "$COHORT/data"              || fail "missing data/ dir"

CASE_DIRS=$(find "$COHORT/data" -mindepth 1 -maxdepth 1 -type d | wc -l | tr -d ' ')
[[ "$CASE_DIRS" -gt 0 ]] || fail "no per-case folders under data/"
ok "produced clinical.parquet, manifest.parquet, $CASE_DIRS per-case folders"

uv run python - <<PY || fail "parquet contents check failed"
import pandas as pd, sys
clin = pd.read_parquet("$COHORT/clinical.parquet")
manif = pd.read_parquet("$COHORT/manifest.parquet")
print(f"    clinical: {len(clin)} cases, {len(clin.columns)} columns")
print(f"    manifest: {len(manif)} files")
assert len(clin) > 0 and len(manif) > 0, "empty parquet outputs"
assert "case_id" in clin.columns and "submitter_id" in clin.columns, "clinical missing key cols"
assert "local_path" in manif.columns, "manifest missing local_path"
ok_paths = manif["local_path"].dropna().astype(str)
import os
assert all(os.path.exists(p) for p in ok_paths), "some local_path entries don't exist on disk"
PY
ok "parquet outputs round-trip and reference real files"
fi

# ------------------------------------------------------------------------------
# Rung 5 — agent over OpenRouter (one-shot, no download)
# ------------------------------------------------------------------------------
if want 5; then
header "5/5  agent (OpenRouter)"
if [[ -z "${OPENROUTER_API_KEY:-}" ]]; then
    fail "OPENROUTER_API_KEY not set — export it or skip rung 5"
fi
note "model=${TCGA_PULL_MODEL:-anthropic/claude-haiku-4-5}"
AGENT_OUT=$(echo "" | uv run tcga-pull agent \
    -q "Use count_files to count TCGA-CHOL Gene Expression Quantification STAR Counts files. Report the number. Do NOT call download." \
    2>&1) \
    || fail "agent run errored"
echo "$AGENT_OUT" | tail -20
echo "$AGENT_OUT" | grep -qE "→ count_files" \
    || fail "agent didn't invoke count_files (no '→ count_files' line)"
ok "agent invoked count_files via OpenRouter tool-use"
fi

echo
echo "${GREEN}${BOLD}all requested rungs passed${RESET}"
