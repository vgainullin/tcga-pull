"""Entity-aware integration of processed cohort modalities.

The integration layer keeps sample identity and biological feature identity
separate. Every source is aligned to the cohort sample index, then normalized
to an explicit entity namespace. Modalities only share a row when both the
sample and normalized entity match.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import polars as pl


@dataclass(frozen=True)
class IntegrationOutputs:
    path: Path
    values: Path
    integrated: Path
    mappings: Path
    manifest: Path


@dataclass(frozen=True)
class ModalitySpec:
    name: str
    source_file: str
    feature_column: str
    value_column: str | None
    entity_column: str
    entity_type: str
    aggregation: str = "mean"
    normalize: str = "none"
    value_kind: str = "numeric"
    mapping: dict[str, Any] | None = None


BUILTIN_MODALITIES: dict[str, dict[str, Any]] = {
    "snv": {
        "source_file": "variants.parquet",
        "feature_column": "hugo_symbol",
        "entity_column": "hugo_symbol",
        "entity_type": "gene_symbol",
        "value_kind": "presence",
        "aggregation": "binary_any",
        "normalize": "uppercase",
    },
    "rna_expression": {
        "source_file": "rna_expression.parquet",
        "feature_column": "gene_name",
        "entity_column": "gene_name",
        "entity_type": "gene_symbol",
        "value_column": "unstranded",
        "aggregation": "mean",
        "normalize": "uppercase",
    },
    "methylation_beta": {
        "source_file": "methylation_beta.parquet",
        "feature_column": "probe_id",
        "entity_column": "probe_id",
        "entity_type": "cpg_probe",
        "value_column": "beta_value",
        "aggregation": "mean",
    },
    "gene_copy_number": {
        "source_file": "gene_copy_number.parquet",
        "feature_column": "gene_name",
        "entity_column": "gene_name",
        "entity_type": "gene_symbol",
        "value_column": "copy_number",
        "aggregation": "mean",
        "normalize": "uppercase",
    },
    "mirna_expression": {
        "source_file": "mirna_expression.parquet",
        "feature_column": "mirna_id",
        "entity_column": "mirna_id",
        "entity_type": "mirna",
        "value_column": "reads_per_million_mirna_mapped",
        "aggregation": "mean",
    },
    "protein_expression": {
        "source_file": "protein_expression.parquet",
        "feature_column": "protein_id",
        "entity_column": "gene_symbol",
        "entity_type": "gene_symbol",
        "value_column": "expression_value",
        "aggregation": "mean",
        "normalize": "uppercase",
    },
}

AGGREGATIONS = {"mean", "median", "min", "max", "sum", "binary_any"}
NORMALIZATIONS = {"none", "uppercase", "lowercase"}
JOIN_COLUMNS = {"submitter_id", "sample_id", "sample_submitter_id"}
RESERVED_MODALITY_NAMES = {
    "case_id",
    "submitter_id",
    "join_id",
    "entity_type",
    "entity_id",
}


def write_integrated_dataset(
    cohort_dir: Path,
    recipe_options: dict[str, Any] | None = None,
    *,
    out_dir: Path | None = None,
) -> IntegrationOutputs:
    """Build entity-aware long and wide integration outputs.

    ``recipe_options.integrate.modalities`` may be a list of built-in modality
    names or a mapping of modality names to overrides/custom specifications.
    External feature mappings are configured per modality with a ``mapping``
    block and are never inferred.
    """
    cohort_dir = Path(cohort_dir)
    options = _integration_options(recipe_options)
    join_on = str(options.get("join_on", "submitter_id"))
    if join_on not in JOIN_COLUMNS:
        raise ValueError(f"integrate.join_on must be one of {sorted(JOIN_COLUMNS)}")
    selected, explicitly_selected = _selected_modalities(options)
    sample_index = _sample_index(cohort_dir)

    frames: list[pl.DataFrame] = []
    mapping_frames: list[pl.DataFrame] = []
    modality_manifest: dict[str, Any] = {}
    skipped: dict[str, str] = {}

    for name, raw_spec in selected.items():
        try:
            spec = _modality_spec(name, raw_spec)
            frame, mappings, details = _integrate_modality(
                cohort_dir, sample_index, spec, join_on=join_on
            )
        except (FileNotFoundError, ValueError) as exc:
            if explicitly_selected:
                raise
            skipped[name] = str(exc)
            continue
        frames.append(frame)
        mapping_frames.append(mappings)
        modality_manifest[name] = details

    if not frames:
        reason = "; ".join(f"{name}: {why}" for name, why in sorted(skipped.items()))
        suffix = f" ({reason})" if reason else ""
        raise ValueError(f"no mergeable modality rows were found{suffix}")

    values = pl.concat(frames, how="vertical").sort(
        ["submitter_id", "join_id", "entity_type", "entity_id", "modality"]
    )
    identity_conflicts = (
        values.select("case_id", "submitter_id", "join_id")
        .unique()
        .group_by("join_id")
        .len()
        .filter(pl.col("len") > 1)
    )
    if not identity_conflicts.is_empty():
        conflicts = identity_conflicts["join_id"].sort().to_list()
        raise ValueError(f"join identifiers map to multiple cases: {conflicts[:5]}")
    mappings = (
        pl.concat(mapping_frames, how="vertical")
        .unique()
        .sort(["modality", "source_feature_id", "entity_type", "entity_id"])
    )
    integrated = values.pivot(
        index=["case_id", "submitter_id", "join_id", "entity_type", "entity_id"],
        on="modality",
        values="value",
        aggregate_function="first",
    ).sort(["submitter_id", "join_id", "entity_type", "entity_id"])

    target = Path(out_dir) if out_dir is not None else cohort_dir / "integrated"
    target.mkdir(parents=True, exist_ok=True)
    values_path = target / "values.parquet"
    integrated_path = target / "integrated.parquet"
    mappings_path = target / "mappings.parquet"
    manifest_path = target / "manifest.json"
    values.write_parquet(values_path)
    integrated.write_parquet(integrated_path)
    mappings.write_parquet(mappings_path)

    manifest = {
        "cohort_dir": str(cohort_dir.resolve()),
        "options": options,
        "join_on": join_on,
        "n_cases": values["submitter_id"].n_unique(),
        "n_samples": values["join_id"].n_unique(),
        "n_entities": values.select(["entity_type", "entity_id"]).unique().height,
        "n_values": values.height,
        "modalities": modality_manifest,
        "skipped_modalities": skipped,
        "files": {
            "values": str(values_path),
            "integrated": str(integrated_path),
            "mappings": str(mappings_path),
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    return IntegrationOutputs(
        path=target,
        values=values_path,
        integrated=integrated_path,
        mappings=mappings_path,
        manifest=manifest_path,
    )


def _integration_options(recipe_options: dict[str, Any] | None) -> dict[str, Any]:
    raw = recipe_options or {}
    if "integrate" in raw:
        value = raw["integrate"]
        if not isinstance(value, dict):
            raise ValueError("recipe_options.integrate must be a mapping")
        return dict(value)
    return dict(raw)


def _selected_modalities(
    options: dict[str, Any],
) -> tuple[dict[str, dict[str, Any]], bool]:
    configured = options.get("modalities")
    if configured is None:
        return {name: {} for name in BUILTIN_MODALITIES}, False
    if isinstance(configured, list):
        names = [str(name) for name in configured]
        if len(names) != len(set(names)):
            raise ValueError("integrate.modalities contains duplicate names")
        if any(not name for name in names):
            raise ValueError("integration modality names must not be empty")
        reserved = [name for name in names if name in RESERVED_MODALITY_NAMES]
        if reserved:
            raise ValueError(f"reserved integration modality name(s): {reserved}")
        return {name: {} for name in names}, True
    if not isinstance(configured, dict):
        raise ValueError("integrate.modalities must be a list or mapping")
    out: dict[str, dict[str, Any]] = {}
    for name, value in configured.items():
        if not str(name):
            raise ValueError("integration modality names must not be empty")
        if str(name) in RESERVED_MODALITY_NAMES:
            raise ValueError(f"reserved integration modality name: {name!r}")
        if value is None:
            out[str(name)] = {}
        elif isinstance(value, dict):
            out[str(name)] = dict(value)
        else:
            raise ValueError(f"integrate.modalities.{name} must be a mapping")
    return out, True


def _modality_spec(name: str, overrides: dict[str, Any]) -> ModalitySpec:
    raw = {**BUILTIN_MODALITIES.get(name, {}), **overrides}
    required = ("source_file", "feature_column", "entity_column", "entity_type")
    missing = [key for key in required if not raw.get(key)]
    value_kind = str(raw.get("value_kind", "numeric"))
    if value_kind not in {"numeric", "presence"}:
        raise ValueError(f"{name}: value_kind must be 'numeric' or 'presence'")
    if value_kind == "numeric" and not raw.get("value_column"):
        missing.append("value_column")
    if missing:
        raise ValueError(f"{name}: missing integration option(s): {sorted(set(missing))}")
    aggregation = str(raw.get("aggregation", "mean"))
    if aggregation not in AGGREGATIONS:
        raise ValueError(f"{name}: unknown aggregation {aggregation!r}")
    normalize = str(raw.get("normalize", "none"))
    if normalize not in NORMALIZATIONS:
        raise ValueError(f"{name}: unknown normalization {normalize!r}")
    mapping = raw.get("mapping")
    if mapping is not None and not isinstance(mapping, dict):
        raise ValueError(f"{name}: mapping must be a mapping")
    if mapping is not None and "aggregation" not in overrides:
        raise ValueError(f"{name}: aggregation is required when mapping is configured")
    return ModalitySpec(
        name=name,
        source_file=str(raw["source_file"]),
        feature_column=str(raw["feature_column"]),
        value_column=str(raw["value_column"]) if raw.get("value_column") else None,
        entity_column=str(raw["entity_column"]),
        entity_type=str(raw["entity_type"]),
        aggregation=aggregation,
        normalize=normalize,
        value_kind=value_kind,
        mapping=dict(mapping) if mapping is not None else None,
    )


def _sample_index(cohort_dir: Path) -> pl.DataFrame:
    for file_name in ("samples.parquet", "clinical.parquet"):
        path = cohort_dir / file_name
        if not path.exists():
            continue
        df = pl.read_parquet(path)
        required = {"case_id", "submitter_id"}
        if required <= set(df.columns):
            index = (
                df.select(
                    pl.col("case_id").cast(pl.Utf8),
                    pl.col("submitter_id").cast(pl.Utf8),
                )
                .drop_nulls()
                .unique()
            )
            conflicts = index.group_by("submitter_id").len().filter(pl.col("len") > 1)
            if not conflicts.is_empty():
                identifiers = conflicts["submitter_id"].sort().to_list()
                raise ValueError(f"submitter identifiers map to multiple cases: {identifiers[:5]}")
            return index
    raise FileNotFoundError(
        "samples.parquet or clinical.parquet must contain case_id and submitter_id"
    )


def _integrate_modality(
    cohort_dir: Path,
    sample_index: pl.DataFrame,
    spec: ModalitySpec,
    *,
    join_on: str,
) -> tuple[pl.DataFrame, pl.DataFrame, dict[str, Any]]:
    source_path = Path(spec.source_file)
    if not source_path.is_absolute():
        source_path = cohort_dir / source_path
    if not source_path.exists():
        raise FileNotFoundError(f"missing {source_path}")
    source = _read_frame(source_path)
    source = _apply_builtin_filters(source, spec.name)
    required = {"submitter_id", spec.feature_column}
    if join_on != "submitter_id":
        required.update({"case_id", join_on})
    if spec.value_column is not None:
        required.add(spec.value_column)
    if spec.mapping is None:
        required.add(spec.entity_column)
    missing = sorted(required - set(source.columns))
    if missing:
        raise ValueError(f"{spec.name}: source is missing columns {missing}")

    if spec.value_kind == "presence":
        value_expression = pl.lit(1.0)
    else:
        assert spec.value_column is not None
        value_expression = pl.col(spec.value_column).cast(pl.Float64, strict=False)

    identity_columns = [
        pl.col("submitter_id").cast(pl.Utf8),
        pl.col(join_on).cast(pl.Utf8).alias("join_id"),
    ]
    if join_on != "submitter_id":
        identity_columns.insert(0, pl.col("case_id").cast(pl.Utf8))

    selected = source.select(
        *identity_columns,
        pl.col(spec.feature_column).cast(pl.Utf8).str.strip_chars().alias("source_feature_id"),
        value_expression.alias("value"),
        *(
            [pl.col(spec.entity_column).cast(pl.Utf8).alias("entity_id")]
            if spec.mapping is None
            else []
        ),
    ).filter(
        pl.col("submitter_id").is_not_null(),
        pl.col("join_id").is_not_null(),
        pl.col("source_feature_id").is_not_null(),
    )

    if spec.mapping is not None:
        mapping = _mapping_frame(spec)
        selected = selected.join(mapping, on="source_feature_id", how="inner")
        mapping_type = "external"
    else:
        mapping_type = "direct"

    selected = selected.with_columns(
        _normalize_entity(pl.col("entity_id"), spec.normalize).alias("entity_id")
    ).filter(
        pl.col("entity_id").is_not_null(),
        pl.col("entity_id").str.len_chars() > 0,
        ~pl.col("entity_id").str.to_lowercase().is_in(["nan", "none", "null"]),
        pl.col("value").is_not_null(),
    )
    if join_on == "submitter_id":
        selected = selected.join(sample_index, on="submitter_id", how="inner")
    else:
        selected = selected.join(sample_index, on=["case_id", "submitter_id"], how="inner")
    if selected.is_empty():
        raise ValueError(f"{spec.name}: no rows matched the cohort samples and entity mapping")
    observed_mappings = selected.select("source_feature_id", "entity_id").unique()

    grouped = (
        _aggregate_values(selected, spec)
        .with_columns(
            pl.lit(spec.name).alias("modality"),
            pl.lit(spec.entity_type).alias("entity_type"),
            pl.lit(spec.aggregation).alias("aggregation"),
        )
        .select(
            "case_id",
            "submitter_id",
            "join_id",
            "entity_type",
            "entity_id",
            "modality",
            "value",
            "aggregation",
            "n_source_features",
            "n_source_rows",
        )
    )
    mappings = (
        observed_mappings.with_columns(
            _normalize_entity(pl.col("entity_id"), spec.normalize).alias("entity_id"),
            pl.lit(spec.name).alias("modality"),
            pl.lit(spec.entity_type).alias("entity_type"),
            pl.lit(mapping_type).alias("mapping_type"),
        )
        .filter(pl.col("entity_id").is_not_null(), pl.col("entity_id").str.len_chars() > 0)
        .select(
            "modality",
            "source_feature_id",
            "entity_type",
            "entity_id",
            "mapping_type",
        )
        .unique()
    )
    details = {
        **asdict(spec),
        "source_file": str(source_path.resolve()),
        "mapping": _manifest_mapping(spec.mapping),
        "join_on": join_on,
        "n_cases": grouped["submitter_id"].n_unique(),
        "n_samples": grouped["join_id"].n_unique(),
        "n_entities": grouped["entity_id"].n_unique(),
        "n_values": grouped.height,
    }
    return grouped, mappings, details


def _read_frame(path: Path) -> pl.DataFrame:
    suffixes = path.suffixes
    if suffixes and suffixes[-1] == ".parquet":
        return pl.read_parquet(path)
    separator = "," if ".csv" in suffixes else "\t"
    return pl.read_csv(path, separator=separator, infer_schema_length=10000, null_values=["", "NA"])


def _apply_builtin_filters(source: pl.DataFrame, modality: str) -> pl.DataFrame:
    if modality != "snv":
        return source
    filters = []
    for column in ("primary_aliquot", "is_coding", "is_rare"):
        if column in source.columns:
            filters.append(pl.col(column).fill_null(False))
    return source.filter(*filters) if filters else source


def _mapping_frame(spec: ModalitySpec) -> pl.DataFrame:
    mapping = spec.mapping or {}
    raw_path = mapping.get("file")
    if not raw_path:
        raise ValueError(f"{spec.name}: mapping.file is required")
    path = Path(str(raw_path)).expanduser()
    if not path.exists():
        raise FileNotFoundError(f"missing mapping file {path}")
    table = _read_frame(path)
    source_column = _mapping_column(
        table,
        mapping.get("source_column"),
        (spec.feature_column, "probe_id", "IlmnID", "Name", "ID"),
        "source",
    )
    entity_column = _mapping_column(
        table,
        mapping.get("entity_column"),
        (spec.entity_column, "gene_symbol", "UCSC_RefGene_Name", "Gene_Symbol"),
        "entity",
    )
    separator = mapping.get("separator")
    filters = mapping.get("filters") or {}
    if not isinstance(filters, dict):
        raise ValueError(f"{spec.name}: mapping.filters must be a mapping")
    missing_filters = [str(column) for column in filters if str(column) not in table.columns]
    if missing_filters:
        raise ValueError(f"{spec.name}: mapping is missing filter columns {missing_filters}")

    mapping_columns = list(
        dict.fromkeys([source_column, entity_column, *(str(column) for column in filters)])
    )
    rows: list[dict[str, str]] = []
    for row in table.select(mapping_columns).iter_rows(named=True):
        source_id = _clean_text(row[source_column])
        entities = _split_values(row[entity_column], separator)
        if source_id is None or not entities:
            continue
        allowed = [True] * len(entities)
        for column, accepted in filters.items():
            accepted_values = {str(value) for value in _as_list(accepted)}
            tokens = _split_values(row[str(column)], separator)
            if not tokens:
                matches = [False] * len(entities)
            elif len(tokens) == 1:
                matches = [tokens[0] in accepted_values] * len(entities)
            elif len(tokens) == len(entities):
                matches = [token in accepted_values for token in tokens]
            else:
                raise ValueError(
                    f"{spec.name}: mapping filter {column!r} has {len(tokens)} values "
                    f"for {len(entities)} entities on source feature {source_id!r}"
                )
            allowed = [keep and match for keep, match in zip(allowed, matches, strict=True)]
        rows.extend(
            {"source_feature_id": source_id, "entity_id": entity}
            for entity, keep in zip(entities, allowed, strict=True)
            if keep
        )
    if not rows:
        raise ValueError(f"{spec.name}: mapping produced no feature-to-entity relationships")
    return pl.DataFrame(rows).unique()


def _mapping_column(
    table: pl.DataFrame,
    requested: Any,
    candidates: tuple[str, ...],
    kind: str,
) -> str:
    if requested is not None:
        column = str(requested)
        if column not in table.columns:
            raise ValueError(f"mapping {kind}_column {column!r} is missing")
        return column
    normalized = {_normalized_column(column): column for column in table.columns}
    for candidate in candidates:
        match = normalized.get(_normalized_column(candidate))
        if match is not None:
            return match
    raise ValueError(f"mapping {kind}_column is required; available columns: {table.columns}")


def _normalized_column(value: str) -> str:
    return "".join(character.lower() for character in value if character.isalnum())


def _split_values(value: Any, separator: Any) -> list[str]:
    text = _clean_text(value)
    if text is None:
        return []
    if separator is None:
        return [text]
    return [item for raw in text.split(str(separator)) if (item := raw.strip())]


def _clean_text(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text if text and text.lower() not in {"nan", "none", "null"} else None


def _as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else [value]


def _normalize_entity(expression: pl.Expr, normalization: str) -> pl.Expr:
    expression = expression.cast(pl.Utf8).str.strip_chars()
    if normalization == "uppercase":
        return expression.str.to_uppercase()
    if normalization == "lowercase":
        return expression.str.to_lowercase()
    return expression


def _aggregate_values(source: pl.DataFrame, spec: ModalitySpec) -> pl.DataFrame:
    groups = ["case_id", "submitter_id", "join_id", "entity_id"]
    if spec.aggregation == "mean":
        value = pl.col("value").mean()
    elif spec.aggregation == "median":
        value = pl.col("value").median()
    elif spec.aggregation == "min":
        value = pl.col("value").min()
    elif spec.aggregation == "max":
        value = pl.col("value").max()
    elif spec.aggregation == "sum":
        value = pl.col("value").sum()
    else:
        value = (pl.col("value").max() > 0).cast(pl.Float64)
    return source.group_by(groups).agg(
        value.cast(pl.Float64).alias("value"),
        pl.col("source_feature_id").n_unique().alias("n_source_features"),
        pl.len().alias("n_source_rows"),
    )


def _manifest_mapping(mapping: dict[str, Any] | None) -> dict[str, Any] | None:
    if mapping is None:
        return None
    out = dict(mapping)
    if out.get("file"):
        out["file"] = str(Path(str(out["file"])).expanduser().resolve())
    return out
