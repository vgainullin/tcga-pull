"""Conversational cohort-building agent over OpenRouter (OpenAI-compatible).

Default model: anthropic/claude-haiku-4-5. Override with TCGA_PULL_MODEL.
Requires OPENROUTER_API_KEY in env.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

from openai import OpenAI
from rich.console import Console
from rich.panel import Panel

from .config import CohortSpec
from .gdc import GDCClient
from .pipeline import fetch_preview, render_preview
from .pipeline import run as run_pipeline
from .schema import SchemaCache

DEFAULT_MODEL = "anthropic/claude-haiku-4-5"
FALLBACK_MODEL = "anthropic/claude-sonnet-4-6"
OPENROUTER_BASE = "https://openrouter.ai/api/v1"

SYSTEM_PROMPT = """\
You are a TCGA cohort-building assistant. You help the user construct a cohort \
of files from the NCI Genomic Data Commons (GDC) and produce a structured per-case \
data product on disk.

Workflow:
1. Understand what the user wants (cancer type, data modality, sample type, etc.).
2. Use list_projects / search_fields to discover correct GDC field names if unsure.
3. Build a GDC filter and call count_files to show file/case counts and download size.
4. Optionally call preview_clinical for a sanity check.
5. ALWAYS get explicit user confirmation before calling the download tool.
6. Open access only — every filter is automatically wrapped with files.access = "open".

GDC filter grammar (JSON):
- {"op": "and", "content": [<clauses>]}
- {"op": "or",  "content": [<clauses>]}
- {"op": "in", "content": {"field": "<path>", "value": [<list>]}}
- {"op": "=",  "content": {"field": "<path>", "value": <single>}}

Common fields:
- cases.project.project_id     e.g. "TCGA-BRCA", "TCGA-LUAD"
- cases.primary_site            e.g. "Breast", "Lung"
- cases.disease_type
- cases.samples.sample_type     "Primary Tumor", "Solid Tissue Normal", ...
- files.data_category           "Transcriptome Profiling", "DNA Methylation", ...
- files.data_type               "Gene Expression Quantification", ...
- files.experimental_strategy   "RNA-Seq", "WXS", "Methylation Array", ...
- analysis.workflow_type        "STAR - Counts", "BWA-aligned reads", ...

Example filter (BRCA primary tumor RNA-seq STAR counts):
{
  "op": "and",
  "content": [
    {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-BRCA"]}},
    {"op": "in", "content": {"field": "cases.samples.sample_type", "value": ["Primary Tumor"]}},
    {"op": "in", "content": {"field": "files.data_type", "value": ["Gene Expression Quantification"]}},
    {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}}
  ]
}

Be concise. After count_files, summarise in plain language and ask whether to refine, \
preview, or download. Do not invent field names — use search_fields if uncertain.
"""

TOOLS: list[dict] = [
    {
        "type": "function",
        "function": {
            "name": "list_projects",
            "description": "List GDC projects under a program (default TCGA). Returns project_id, name, primary_site for each.",
            "parameters": {
                "type": "object",
                "properties": {
                    "program": {"type": "string", "default": "TCGA"},
                },
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "search_fields",
            "description": "Search GDC schema for queryable fields whose path or description matches a keyword. Use to discover field names like 'workflow', 'sample_type', 'stage'.",
            "parameters": {
                "type": "object",
                "properties": {
                    "keyword": {"type": "string"},
                    "limit": {"type": "integer", "default": 20},
                },
                "required": ["keyword"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "count_files",
            "description": "Cheap dry-run: counts files and unique cases for a GDC filter, plus total estimated download size. Always call before downloading.",
            "parameters": {
                "type": "object",
                "properties": {
                    "filter": {
                        "type": "object",
                        "description": "GDC filter JSON tree (see system prompt).",
                        "additionalProperties": True,
                    },
                },
                "required": ["filter"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "preview_clinical",
            "description": "Fetch a few example clinical records matching the filter — sanity check before downloading.",
            "parameters": {
                "type": "object",
                "properties": {
                    "filter": {"type": "object", "additionalProperties": True},
                    "n": {"type": "integer", "default": 3},
                },
                "required": ["filter"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "download",
            "description": "Download files matching the filter and produce the structured cohort. Only call after the user has explicitly confirmed.",
            "parameters": {
                "type": "object",
                "properties": {
                    "name": {"type": "string", "description": "Cohort folder name."},
                    "filter": {"type": "object", "additionalProperties": True},
                },
                "required": ["name", "filter"],
            },
        },
    },
]


class AgentTools:
    def __init__(self, gdc: GDCClient, schema: SchemaCache, out_dir: Path, console: Console):
        self.gdc = gdc
        self.schema = schema
        self.out_dir = out_dir
        self.console = console

    def list_projects(self, program: str = "TCGA") -> list[dict]:
        rows = self.gdc.list_projects(program=program)
        return [
            {
                "project_id": r.get("project_id"),
                "name": r.get("name"),
                "primary_site": r.get("primary_site"),
            }
            for r in rows
        ]

    def search_fields(self, keyword: str, limit: int = 20) -> list[dict]:
        kw = keyword.lower()
        out: list[dict] = []
        for ep in ("files", "cases"):
            sch = self.schema.get(ep)
            descs = sch.descriptions
            for f in sch.fields:
                if kw in f.lower() or kw in descs.get(f, "").lower():
                    out.append({"endpoint": ep, "field": f, "description": descs.get(f, "")})
                    if len(out) >= limit:
                        return out
        return out

    def count_files(self, filter: dict) -> dict:
        from .gdc import open_access

        flt = open_access(filter)
        n_files = self.gdc.count_files(flt)
        n_cases = self.gdc.count_cases(flt)
        return {"files": n_files, "cases": n_cases, "filter_with_open_access": flt}

    def preview_clinical(self, filter: dict, n: int = 3) -> list[dict]:
        from .gdc import open_access

        flt = open_access(filter)
        cases: list[dict] = []
        for case in self.gdc.iter_cases(flt, page_size=n):
            cases.append(
                {
                    "case_id": case.get("case_id"),
                    "submitter_id": case.get("submitter_id"),
                    "project_id": (case.get("project") or {}).get("project_id"),
                    "primary_site": case.get("primary_site"),
                    "demographic": case.get("demographic"),
                    "first_diagnosis": (case.get("diagnoses") or [None])[0],
                }
            )
            if len(cases) >= n:
                break
        return cases

    def download(self, name: str, filter: dict) -> dict:
        # hard confirmation gate, regardless of agent assurances
        import questionary

        spec = CohortSpec(
            name=name,
            out_dir=self.out_dir,
            gdc_filter=filter,
        )
        preview = fetch_preview(spec, client=self.gdc)
        render_preview(preview, console=self.console)
        if preview.n_files == 0:
            return {"status": "no_match"}
        if not questionary.confirm(
            f"Download {preview.n_files} files into {spec.cohort_dir}?", default=False
        ).ask():
            return {"status": "user_declined"}
        cohort_dir = run_pipeline(spec, console=self.console, client=self.gdc)
        return {"status": "ok", "cohort_dir": str(cohort_dir)}


def _execute(tools: AgentTools, name: str, args: dict) -> Any:
    fn = getattr(tools, name, None)
    if fn is None:
        return {"error": f"unknown tool: {name}"}
    try:
        return fn(**args)
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}


def _client() -> OpenAI:
    api_key = os.environ.get("OPENROUTER_API_KEY")
    if not api_key:
        raise RuntimeError(
            "OPENROUTER_API_KEY not set. Get a key at https://openrouter.ai/ "
            "and export it: `export OPENROUTER_API_KEY=...`"
        )
    return OpenAI(
        base_url=OPENROUTER_BASE,
        api_key=api_key,
        default_headers={
            "HTTP-Referer": "https://github.com/tcga-pull",
            "X-Title": "tcga-pull",
        },
    )


def _stream_turn(client: OpenAI, model: str, messages: list[dict], console: Console) -> dict:
    """One streamed turn. Returns the final assistant message dict (with tool_calls if any)."""
    text_chunks: list[str] = []
    tool_calls: dict[int, dict] = {}

    # OpenAI SDK has TypedDict-based message/tool params; we build them as plain
    # dicts (which work at runtime). Narrow ignores at the SDK boundary.
    stream = client.chat.completions.create(
        model=model,
        messages=messages,  # type: ignore[arg-type]
        tools=TOOLS,  # type: ignore[arg-type]
        stream=True,
    )
    printed_any_text = False
    for chunk in stream:
        if not chunk.choices:  # type: ignore[union-attr]
            continue
        delta = chunk.choices[0].delta  # type: ignore[union-attr]
        if delta.content:
            console.out(delta.content, end="", style="cyan", highlight=False)
            text_chunks.append(delta.content)
            printed_any_text = True
        if delta.tool_calls:
            for tc in delta.tool_calls:
                slot = tool_calls.setdefault(
                    tc.index,
                    {"id": None, "type": "function", "function": {"name": "", "arguments": ""}},
                )
                if tc.id:
                    slot["id"] = tc.id
                if tc.function and tc.function.name:
                    slot["function"]["name"] = tc.function.name
                if tc.function and tc.function.arguments:
                    slot["function"]["arguments"] += tc.function.arguments

    if printed_any_text:
        console.print()

    msg: dict = {"role": "assistant"}
    if text_chunks:
        msg["content"] = "".join(text_chunks)
    if tool_calls:
        msg["tool_calls"] = [tool_calls[i] for i in sorted(tool_calls)]
    return msg


def run_agent(query: str | None, out: Path, console: Console) -> None:
    client = _client()
    model = os.environ.get("TCGA_PULL_MODEL", DEFAULT_MODEL)
    gdc = GDCClient()
    schema = SchemaCache(gdc)
    tools = AgentTools(gdc=gdc, schema=schema, out_dir=out.expanduser().resolve(), console=console)

    console.print(
        Panel(
            f"[bold]tcga-pull[/bold] · model=[cyan]{model}[/cyan]\n"
            "Describe a cohort. I'll show counts before anything downloads.\n"
            "[dim]Empty line / Ctrl-D to exit.[/dim]",
            border_style="dim",
        )
    )

    messages: list[dict] = [{"role": "system", "content": SYSTEM_PROMPT}]

    def get_user_input() -> str | None:
        try:
            return console.input("[bold green]>[/bold green] ").strip()
        except EOFError:
            return None

    first = query
    while True:
        user_text: str
        if first is not None:
            user_text = first
            first = None
            console.print(f"[bold green]>[/bold green] {user_text}")
        else:
            line = get_user_input()
            if not line:
                console.print("[dim]bye.[/dim]")
                return
            user_text = line

        messages.append({"role": "user", "content": user_text})

        # inner tool-use loop
        while True:
            assistant_msg = _stream_turn(client, model, messages, console)
            messages.append(assistant_msg)
            tool_calls = assistant_msg.get("tool_calls") or []
            if not tool_calls:
                break  # back to user
            for tc in tool_calls:
                name = tc["function"]["name"]
                try:
                    args = json.loads(tc["function"]["arguments"] or "{}")
                except json.JSONDecodeError as e:
                    args = {}
                    console.log(f"[yellow]bad tool args from model[/yellow]: {e}")
                console.log(f"[dim]→ {name}({_short(args)})[/dim]")
                result = _execute(tools, name, args)
                messages.append(
                    {
                        "role": "tool",
                        "tool_call_id": tc["id"],
                        "name": name,
                        "content": json.dumps(result, default=str)[:20000],
                    }
                )


def _short(args: dict) -> str:
    s = json.dumps(args, default=str)
    return s if len(s) < 120 else s[:117] + "..."
