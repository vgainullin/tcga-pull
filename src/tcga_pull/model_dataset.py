"""Build case-aligned matrices for downstream multi-omic modeling."""

from __future__ import annotations

import hashlib
import json
import random
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import polars as pl

KNOWN_MODALITIES = {
    "snv",
    "rna_expression",
    "methylation_beta",
    "gene_copy_number",
    "mirna_expression",
    "protein_expression",
}

DEFAULT_MODALITIES = (
    "snv",
    "rna_expression",
    "methylation_beta",
    "gene_copy_number",
    "mirna_expression",
    "protein_expression",
)


@dataclass(frozen=True)
class ModelDatasetOptions:
    label_column: str = "oncotree_tissue"
    modalities: tuple[str, ...] = DEFAULT_MODALITIES
    train_fraction: float = 0.70
    val_fraction: float = 0.15
    test_fraction: float = 0.15
    min_class_count: int = 2
    feature_min_samples: int = 1
    max_features_per_modality: int | None = None
    seed: int = 13


@dataclass(frozen=True)
class ModelDatasetOutputs:
    path: Path
    manifest: Path
    samples: Path
    feature_index: Path
    matrices: dict[str, Path]


@dataclass(frozen=True)
class _MatrixResult:
    modality: str
    path: Path
    feature_rows: list[dict[str, Any]]
    available_submitters: set[str]
    n_features: int


def write_model_dataset(
    cohort_dir: Path,
    recipe_options: dict[str, Any] | None = None,
    *,
    out_dir: Path | None = None,
    label_column: str | None = None,
    modalities: list[str] | tuple[str, ...] | None = None,
    train_fraction: float | None = None,
    val_fraction: float | None = None,
    test_fraction: float | None = None,
    min_class_count: int | None = None,
    feature_min_samples: int | None = None,
    max_features_per_modality: int | None = None,
    seed: int | None = None,
) -> ModelDatasetOutputs:
    cohort_dir = Path(cohort_dir)
    options = _options_from(
        recipe_options,
        label_column=label_column,
        modalities=modalities,
        train_fraction=train_fraction,
        val_fraction=val_fraction,
        test_fraction=test_fraction,
        min_class_count=min_class_count,
        feature_min_samples=feature_min_samples,
        max_features_per_modality=max_features_per_modality,
        seed=seed,
    )
    _validate_options(options)

    dataset_dir = _dataset_dir(cohort_dir, recipe_options, out_dir)
    dataset_dir.mkdir(parents=True, exist_ok=True)

    samples = _selected_samples(cohort_dir, options)
    sample_ids = samples.select(["case_id", "submitter_id"])

    matrix_results: list[_MatrixResult] = []
    skipped: dict[str, str] = {}
    for modality in options.modalities:
        result = _write_modality_matrix(cohort_dir, dataset_dir, sample_ids, modality, options)
        if result is None:
            skipped[modality] = "source parquet missing or has no usable feature rows"
            continue
        matrix_results.append(result)

    samples = _attach_modality_flags(samples, matrix_results)
    samples_path = dataset_dir / "samples.parquet"
    samples.write_parquet(samples_path)

    feature_rows = [row for result in matrix_results for row in result.feature_rows]
    feature_index = pl.DataFrame(feature_rows) if feature_rows else _empty_feature_index()
    feature_index_path = dataset_dir / "feature_index.parquet"
    feature_index.write_parquet(feature_index_path)

    manifest_path = dataset_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "cohort_dir": str(cohort_dir),
                "options": asdict(options),
                "n_samples": samples.height,
                "label_column": options.label_column,
                "label_counts": _value_counts(samples, options.label_column),
                "split_counts": _value_counts(samples, "split"),
                "modalities": {
                    result.modality: {
                        "path": str(result.path),
                        "n_features": result.n_features,
                        "n_available_samples": len(result.available_submitters),
                    }
                    for result in matrix_results
                },
                "skipped_modalities": skipped,
            },
            indent=2,
            sort_keys=True,
        )
    )

    return ModelDatasetOutputs(
        path=dataset_dir,
        manifest=manifest_path,
        samples=samples_path,
        feature_index=feature_index_path,
        matrices={result.modality: result.path for result in matrix_results},
    )


def _options_from(
    recipe_options: dict[str, Any] | None,
    *,
    label_column: str | None,
    modalities: list[str] | tuple[str, ...] | None,
    train_fraction: float | None,
    val_fraction: float | None,
    test_fraction: float | None,
    min_class_count: int | None,
    feature_min_samples: int | None,
    max_features_per_modality: int | None,
    seed: int | None,
) -> ModelDatasetOptions:
    raw = recipe_options or {}
    if "model_dataset" in raw and isinstance(raw["model_dataset"], dict):
        raw = raw["model_dataset"]
    option_modalities = modalities
    if option_modalities is None and raw.get("modalities") is not None:
        option_modalities = tuple(str(item) for item in raw["modalities"])
    return ModelDatasetOptions(
        label_column=str(label_column or raw.get("label_column") or "oncotree_tissue"),
        modalities=tuple(option_modalities or DEFAULT_MODALITIES),
        train_fraction=float(
            train_fraction if train_fraction is not None else raw.get("train_fraction", 0.70)
        ),
        val_fraction=float(
            val_fraction if val_fraction is not None else raw.get("val_fraction", 0.15)
        ),
        test_fraction=float(
            test_fraction if test_fraction is not None else raw.get("test_fraction", 0.15)
        ),
        min_class_count=int(
            min_class_count if min_class_count is not None else raw.get("min_class_count", 2)
        ),
        feature_min_samples=int(
            feature_min_samples
            if feature_min_samples is not None
            else raw.get("feature_min_samples", 1)
        ),
        max_features_per_modality=(
            max_features_per_modality
            if max_features_per_modality is not None
            else raw.get("max_features_per_modality")
        ),
        seed=int(seed if seed is not None else raw.get("seed", 13)),
    )


def _dataset_dir(
    cohort_dir: Path,
    recipe_options: dict[str, Any] | None,
    out_dir: Path | None,
) -> Path:
    if out_dir is not None:
        return out_dir.expanduser().resolve()
    raw = recipe_options or {}
    if "model_dataset" in raw and isinstance(raw["model_dataset"], dict):
        raw = raw["model_dataset"]
    configured = raw.get("out_dir")
    if configured:
        path = Path(str(configured)).expanduser()
        return path if path.is_absolute() else cohort_dir / path
    return cohort_dir / "model_dataset"


def _validate_options(options: ModelDatasetOptions) -> None:
    unknown = sorted(set(options.modalities) - KNOWN_MODALITIES)
    if unknown:
        raise ValueError(f"unknown model dataset modalities: {unknown}")
    if options.min_class_count <= 0:
        raise ValueError("min_class_count must be > 0")
    if options.feature_min_samples <= 0:
        raise ValueError("feature_min_samples must be > 0")
    if options.max_features_per_modality is not None and options.max_features_per_modality <= 0:
        raise ValueError("max_features_per_modality must be > 0")
    total = options.train_fraction + options.val_fraction + options.test_fraction
    if abs(total - 1.0) > 1e-6:
        raise ValueError("train_fraction + val_fraction + test_fraction must equal 1.0")
    for name, value in (
        ("train_fraction", options.train_fraction),
        ("val_fraction", options.val_fraction),
        ("test_fraction", options.test_fraction),
    ):
        if value < 0:
            raise ValueError(f"{name} must be >= 0")


def _selected_samples(cohort_dir: Path, options: ModelDatasetOptions) -> pl.DataFrame:
    path = cohort_dir / "samples.parquet"
    if not path.exists():
        raise FileNotFoundError(f"missing {path} - run `tcga-pull samples` first")
    samples = pl.read_parquet(path)
    required = {"case_id", "submitter_id", options.label_column}
    missing = sorted(required - set(samples.columns))
    if missing:
        raise ValueError(f"samples.parquet missing required model dataset columns: {missing}")

    samples = samples.filter(pl.col(options.label_column).is_not_null())
    labels = (
        samples.group_by(options.label_column)
        .agg(pl.len().alias("n"))
        .filter(pl.col("n") >= options.min_class_count)
        .select(options.label_column)
    )
    samples = samples.join(labels, on=options.label_column, how="inner")
    if samples.is_empty():
        raise ValueError("no samples remain after label and min_class_count filtering")

    splits = _split_rows(samples, options)
    return samples.join(splits, on="submitter_id", how="inner").sort("submitter_id")


def _split_rows(samples: pl.DataFrame, options: ModelDatasetOptions) -> pl.DataFrame:
    rows: list[dict[str, str]] = []
    for label_value in sorted(samples[options.label_column].unique().to_list()):
        group = (
            samples.filter(pl.col(options.label_column) == label_value)
            .select("submitter_id")
            .sort("submitter_id")["submitter_id"]
            .to_list()
        )
        rng = random.Random(options.seed + _stable_int(str(label_value)))
        rng.shuffle(group)
        split_by_submitter = _assign_group_splits(
            group,
            train_fraction=options.train_fraction,
            val_fraction=options.val_fraction,
            test_fraction=options.test_fraction,
        )
        rows.extend(
            {"submitter_id": submitter, "split": split}
            for submitter, split in split_by_submitter.items()
        )
    return pl.DataFrame(rows)


def _assign_group_splits(
    submitters: list[str],
    *,
    train_fraction: float,
    val_fraction: float,
    test_fraction: float,
) -> dict[str, str]:
    n = len(submitters)
    if n == 1:
        counts = {"train": 1, "val": 0, "test": 0}
    elif n == 2:
        counts = {"train": 1, "val": 0, "test": 1 if test_fraction > 0 else 0}
        if counts["test"] == 0:
            counts["train"] = 2
    else:
        n_val = max(1, round(n * val_fraction)) if val_fraction > 0 else 0
        n_test = max(1, round(n * test_fraction)) if test_fraction > 0 else 0
        while n_val + n_test >= n:
            if n_val >= n_test and n_val > 0:
                n_val -= 1
            elif n_test > 0:
                n_test -= 1
        counts = {"train": n - n_val - n_test, "val": n_val, "test": n_test}

    out: dict[str, str] = {}
    idx = 0
    for split in ("train", "val", "test"):
        for submitter in submitters[idx : idx + counts[split]]:
            out[submitter] = split
        idx += counts[split]
    return out


def _write_modality_matrix(
    cohort_dir: Path,
    dataset_dir: Path,
    sample_ids: pl.DataFrame,
    modality: str,
    options: ModelDatasetOptions,
) -> _MatrixResult | None:
    if modality == "snv":
        return _write_snv_matrix(cohort_dir, dataset_dir, sample_ids, options)

    specs = {
        "rna_expression": ("rna_expression.parquet", "gene_id", "unstranded", "log1p"),
        "methylation_beta": ("methylation_beta.parquet", "probe_id", "beta_value", "none"),
        "gene_copy_number": ("gene_copy_number.parquet", "gene_name", "copy_number", "none"),
        "mirna_expression": (
            "mirna_expression.parquet",
            "mirna_id",
            "reads_per_million_mirna_mapped",
            "log1p",
        ),
        "protein_expression": (
            "protein_expression.parquet",
            "protein_id",
            "expression_value",
            "none",
        ),
    }
    file_name, feature_col, preferred_value_col, transform = specs[modality]
    path = cohort_dir / file_name
    if not path.exists():
        return None
    df = pl.read_parquet(path)
    if preferred_value_col not in df.columns and modality == "mirna_expression":
        preferred_value_col = "read_count"
    return _write_long_matrix(
        dataset_dir=dataset_dir,
        sample_ids=sample_ids,
        source=df,
        modality=modality,
        feature_col=feature_col,
        value_col=preferred_value_col,
        transform=transform,
        options=options,
    )


def _write_snv_matrix(
    cohort_dir: Path,
    dataset_dir: Path,
    sample_ids: pl.DataFrame,
    options: ModelDatasetOptions,
) -> _MatrixResult | None:
    path = cohort_dir / "variants.parquet"
    if not path.exists():
        return None
    variants = pl.read_parquet(path)
    required = {"submitter_id", "hugo_symbol"}
    if variants.is_empty() or not required <= set(variants.columns):
        return None
    filters = [pl.col("hugo_symbol").is_not_null()]
    if "primary_aliquot" in variants.columns:
        filters.append(pl.col("primary_aliquot").fill_null(False))
    if "is_coding" in variants.columns:
        filters.append(pl.col("is_coding").fill_null(False))
    if "is_rare" in variants.columns:
        filters.append(pl.col("is_rare").fill_null(False))
    filtered = variants.filter(*filters).select(["submitter_id", "hugo_symbol"])
    if filtered.is_empty():
        return None

    available = _snv_available_submitters(cohort_dir, variants)
    grouped = (
        filtered.group_by(["submitter_id", "hugo_symbol"])
        .agg(pl.len().alias("value"))
        .with_columns((pl.col("value") > 0).cast(pl.Int8))
    )
    return _pivot_and_write(
        dataset_dir=dataset_dir,
        sample_ids=sample_ids,
        source=grouped,
        modality="snv",
        feature_col="hugo_symbol",
        value_col="value",
        transform="binary",
        fill_null=0,
        available_submitters=available,
        options=options,
    )


def _write_long_matrix(
    *,
    dataset_dir: Path,
    sample_ids: pl.DataFrame,
    source: pl.DataFrame,
    modality: str,
    feature_col: str,
    value_col: str,
    transform: str,
    options: ModelDatasetOptions,
) -> _MatrixResult | None:
    required = {"submitter_id", feature_col, value_col}
    if source.is_empty() or not required <= set(source.columns):
        return None
    value_expr = pl.col(value_col).cast(pl.Float64, strict=False)
    if transform == "log1p":
        value_expr = pl.when(value_expr >= 0).then(value_expr.log1p()).otherwise(None)
    filtered = source.select(
        "submitter_id",
        pl.col(feature_col).cast(pl.Utf8).alias(feature_col),
        value_expr.alias("value"),
    ).filter(pl.col(feature_col).is_not_null(), pl.col("value").is_not_null())
    if filtered.is_empty():
        return None
    available = set(filtered["submitter_id"].drop_nulls().cast(pl.Utf8).unique().to_list())
    grouped = filtered.group_by(["submitter_id", feature_col]).agg(pl.col("value").mean())
    return _pivot_and_write(
        dataset_dir=dataset_dir,
        sample_ids=sample_ids,
        source=grouped,
        modality=modality,
        feature_col=feature_col,
        value_col="value",
        transform=transform,
        fill_null=None,
        available_submitters=available,
        options=options,
    )


def _pivot_and_write(
    *,
    dataset_dir: Path,
    sample_ids: pl.DataFrame,
    source: pl.DataFrame,
    modality: str,
    feature_col: str,
    value_col: str,
    transform: str,
    fill_null: int | float | None,
    available_submitters: set[str],
    options: ModelDatasetOptions,
) -> _MatrixResult | None:
    feature_stats = _selected_feature_stats(source, feature_col, options)
    if feature_stats.is_empty():
        return None
    features = feature_stats[feature_col].to_list()
    source = source.filter(pl.col(feature_col).is_in(features))
    column_map = _feature_column_map(modality, features)
    wide = source.pivot(
        index="submitter_id",
        on=feature_col,
        values=value_col,
        aggregate_function="mean",
    )
    wide = wide.rename(
        {feature: column_map[str(feature)] for feature in features if feature in wide.columns}
    )
    matrix = sample_ids.join(wide, on="submitter_id", how="left")
    feature_columns = [column_map[str(feature)] for feature in features]
    if fill_null is not None and feature_columns:
        matrix = matrix.with_columns(pl.col(feature_columns).fill_null(fill_null))

    path = dataset_dir / f"{modality}.parquet"
    matrix.write_parquet(path)
    feature_rows = [
        {
            "modality": modality,
            "feature_id": str(feature),
            "feature_name": str(feature),
            "column_name": column_map[str(feature)],
            "value_column": value_col,
            "transform": transform,
            "n_samples": int(n_samples),
        }
        for feature, n_samples in feature_stats.select([feature_col, "n_samples"]).iter_rows()
    ]
    return _MatrixResult(
        modality=modality,
        path=path,
        feature_rows=feature_rows,
        available_submitters=available_submitters,
        n_features=len(feature_columns),
    )


def _selected_feature_stats(
    source: pl.DataFrame,
    feature_col: str,
    options: ModelDatasetOptions,
) -> pl.DataFrame:
    stats = (
        source.filter(pl.col("value").is_not_null())
        .group_by(feature_col)
        .agg(pl.col("submitter_id").n_unique().alias("n_samples"))
        .filter(pl.col("n_samples") >= options.feature_min_samples)
        .sort(["n_samples", feature_col], descending=[True, False])
    )
    if options.max_features_per_modality is not None:
        stats = stats.head(options.max_features_per_modality)
    return stats.sort(feature_col)


def _feature_column_map(modality: str, features: list[Any]) -> dict[str, str]:
    out: dict[str, str] = {}
    used: dict[str, int] = {}
    for feature in features:
        key = str(feature)
        base = f"{modality}__{_slug(key)}"
        count = used.get(base, 0) + 1
        used[base] = count
        out[key] = base if count == 1 else f"{base}__{count}"
    return out


def _slug(value: str) -> str:
    chars = [ch.lower() if ch.isalnum() else "_" for ch in value.strip()]
    slug = "".join(chars)
    while "__" in slug:
        slug = slug.replace("__", "_")
    slug = slug.strip("_")
    return slug or "feature"


def _attach_modality_flags(
    samples: pl.DataFrame,
    matrix_results: list[_MatrixResult],
) -> pl.DataFrame:
    out = samples
    for result in matrix_results:
        out = out.with_columns(
            pl.col("submitter_id")
            .is_in(sorted(result.available_submitters))
            .alias(f"has_{result.modality}")
        )
    return out


def _snv_available_submitters(cohort_dir: Path, variants: pl.DataFrame) -> set[str]:
    manifest_path = cohort_dir / "manifest.parquet"
    if manifest_path.exists():
        manifest = pl.read_parquet(manifest_path)
        if {"submitter_id", "data_category", "data_format"} <= set(manifest.columns):
            rows = manifest.filter(
                (pl.col("data_category") == "Simple Nucleotide Variation")
                | (pl.col("data_format") == "MAF")
            )
            return set(rows["submitter_id"].drop_nulls().cast(pl.Utf8).unique().to_list())
    return set(variants["submitter_id"].drop_nulls().cast(pl.Utf8).unique().to_list())


def _value_counts(df: pl.DataFrame, column: str) -> dict[str, int]:
    if column not in df.columns:
        return {}
    rows = df.group_by(column).agg(pl.len().alias("n")).sort(column).to_dicts()
    return {str(row[column]): int(row["n"]) for row in rows}


def _empty_feature_index() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "modality": pl.Series([], dtype=pl.Utf8),
            "feature_id": pl.Series([], dtype=pl.Utf8),
            "feature_name": pl.Series([], dtype=pl.Utf8),
            "column_name": pl.Series([], dtype=pl.Utf8),
            "value_column": pl.Series([], dtype=pl.Utf8),
            "transform": pl.Series([], dtype=pl.Utf8),
            "n_samples": pl.Series([], dtype=pl.Int64),
        }
    )


def _stable_int(value: str) -> int:
    digest = hashlib.sha256(value.encode("utf-8")).hexdigest()
    return int(digest[:12], 16)
