"""Pinned, cohort-local annotation resources used by analysis recipes."""

from __future__ import annotations

import hashlib
import json
import os
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import polars as pl
import requests

ANNOTATION_RELEASE = "InfiniumAnnotation-v8.1-gencode-v41-hg38"
ANNOTATION_NORMALIZER_VERSION = 1
ANNOTATION_BASE_URL = "https://raw.githubusercontent.com/zhou-lab/InfiniumAnnotationData/main/Anno"


class AnnotationResolutionError(RuntimeError):
    """Raised when a built-in annotation cannot be resolved safely."""


@dataclass(frozen=True)
class MethylationAnnotation:
    platform: str
    array: str
    url: str
    sha256: str
    size_bytes: int
    release: str = ANNOTATION_RELEASE
    genome_build: str = "hg38"

    @property
    def slug(self) -> str:
        return self.array.lower()


METHYLATION_ANNOTATIONS: dict[str, MethylationAnnotation] = {
    "illumina human methylation 27": MethylationAnnotation(
        platform="Illumina Human Methylation 27",
        array="HM27",
        url=f"{ANNOTATION_BASE_URL}/HM27/HM27.hg38.manifest.gencode.v41.tsv.gz",
        sha256="3bb525a91f12ace233e35d12222aab3af5cc215c94a8ad3b16335c6848c9490e",
        size_bytes=2_473_291,
    ),
    "illumina human methylation 450": MethylationAnnotation(
        platform="Illumina Human Methylation 450",
        array="HM450",
        url=f"{ANNOTATION_BASE_URL}/HM450/HM450.hg38.manifest.gencode.v41.tsv.gz",
        sha256="4ed7385d8f487c24665a87b8de5d1c4cd147187e60747dd1f2701101c3935066",
        size_bytes=35_907_305,
    ),
    "illumina methylation epic": MethylationAnnotation(
        platform="Illumina Methylation Epic",
        array="EPIC",
        url=f"{ANNOTATION_BASE_URL}/EPIC/EPIC.hg38.manifest.gencode.v41.tsv.gz",
        sha256="cb12bf4d1cd1a1ea8994a0e7234a27a9e5cab9e0986a4b01e4b4c8c566203f68",
        size_bytes=62_477_035,
    ),
    "illumina methylation epic v2": MethylationAnnotation(
        platform="Illumina Methylation Epic v2",
        array="EPICv2",
        url=f"{ANNOTATION_BASE_URL}/EPICv2/EPICv2.hg38.manifest.gencode.v41.tsv.gz",
        sha256="92b6bd7b99e374b82601d5bf97b5b0e7b4a3de574e6c4c7dfc5caf2bb1f13878",
        size_bytes=63_652_988,
    ),
}


def normalize_platform(value: str) -> str:
    return " ".join(value.strip().lower().split())


def resolve_methylation_annotations(
    cohort_dir: Path,
    platforms: list[str],
    *,
    regions: tuple[str, ...] = ("promoter",),
) -> tuple[pl.DataFrame, list[dict[str, Any]]]:
    """Resolve platform-specific CpG-to-gene mappings inside the cohort."""
    requested_regions = tuple(dict.fromkeys(region.lower() for region in regions))
    unknown_regions = sorted(set(requested_regions) - {"promoter", "body"})
    if not requested_regions or unknown_regions:
        raise AnnotationResolutionError(
            "methylation annotation regions must contain promoter and/or body; "
            f"unknown values: {unknown_regions}"
        )

    normalized_platforms = sorted({normalize_platform(value) for value in platforms if value})
    unknown = [value for value in normalized_platforms if value not in METHYLATION_ANNOTATIONS]
    if unknown:
        supported = [spec.platform for spec in METHYLATION_ANNOTATIONS.values()]
        raise AnnotationResolutionError(
            f"unsupported GDC methylation platform(s): {unknown}; supported: {supported}"
        )
    if not normalized_platforms:
        raise AnnotationResolutionError(
            "methylation rows have no GDC platform; rebuild methylation_beta.parquet "
            "from a manifest that retains platform metadata"
        )

    frames: list[pl.DataFrame] = []
    metadata: list[dict[str, Any]] = []
    for platform_key in normalized_platforms:
        spec = METHYLATION_ANNOTATIONS[platform_key]
        frame, details = _resolve_one(cohort_dir, spec, requested_regions)
        frames.append(frame)
        metadata.append(details)
    return pl.concat(frames, how="vertical"), metadata


def _resolve_one(
    cohort_dir: Path,
    spec: MethylationAnnotation,
    regions: tuple[str, ...],
) -> tuple[pl.DataFrame, dict[str, Any]]:
    annotation_dir = Path(cohort_dir) / "annotations"
    annotation_dir.mkdir(parents=True, exist_ok=True)
    source_path = annotation_dir / f"{spec.slug}.gencode-v41.hg38.tsv.gz"
    parquet_path = annotation_dir / f"{spec.slug}.probe-gene.parquet"
    metadata_path = annotation_dir / f"{spec.slug}.probe-gene.json"

    if source_path.exists():
        actual = _sha256(source_path)
        if actual != spec.sha256:
            raise AnnotationResolutionError(
                f"annotation checksum mismatch for {source_path}: expected {spec.sha256}, "
                f"got {actual}; refusing to use or overwrite it"
            )
    else:
        _download_annotation(spec, source_path)

    cache_key = {
        "source_sha256": spec.sha256,
        "regions": list(regions),
        "normalizer_version": ANNOTATION_NORMALIZER_VERSION,
    }
    if parquet_path.exists() and metadata_path.exists():
        cached = json.loads(metadata_path.read_text())
        if all(cached.get(key) == value for key, value in cache_key.items()) and cached.get(
            "mapping_sha256"
        ) == _sha256(parquet_path):
            return pl.read_parquet(parquet_path), cached

    mapping = build_methylation_mapping(source_path, spec.platform, regions=regions)
    mapping.write_parquet(parquet_path)
    details = {
        **asdict(spec),
        **cache_key,
        "source_file": str(source_path.resolve()),
        "mapping_file": str(parquet_path.resolve()),
        "mapping_sha256": _sha256(parquet_path),
        "regions": list(regions),
        "n_mappings": mapping.height,
        "n_probes": mapping["source_feature_id"].n_unique(),
        "n_genes": mapping["entity_id"].n_unique(),
    }
    metadata_path.write_text(json.dumps(details, indent=2, sort_keys=True))
    return mapping, details


def build_methylation_mapping(
    path: Path,
    platform: str,
    *,
    regions: tuple[str, ...] = ("promoter",),
) -> pl.DataFrame:
    """Normalize one InfiniumAnnotation release into aligned probe/gene rows."""
    table = pl.read_csv(
        path,
        separator="\t",
        columns=["probeID", "geneNames", "distToTSS"],
        infer_schema_length=10_000,
        schema_overrides={
            "probeID": pl.Utf8,
            "geneNames": pl.Utf8,
            "distToTSS": pl.Utf8,
        },
        null_values=["", "NA"],
    ).drop_nulls(["probeID", "geneNames", "distToTSS"])
    split = table.with_columns(
        pl.col("geneNames").str.split(";").alias("genes"),
        pl.col("distToTSS").str.split(";").alias("distances"),
    )
    invalid = split.filter(pl.col("genes").list.len() != pl.col("distances").list.len())
    if not invalid.is_empty():
        probe = invalid["probeID"][0]
        raise AnnotationResolutionError(
            f"annotation has unaligned gene/distance values for probe {probe!r}"
        )
    return (
        split.explode("genes", "distances")
        .with_columns(
            pl.col("probeID").cast(pl.Utf8).str.strip_chars().alias("source_feature_id"),
            pl.col("genes").cast(pl.Utf8).str.strip_chars().str.to_uppercase().alias("entity_id"),
            pl.col("distances").cast(pl.Int64, strict=False).alias("distance_to_tss"),
        )
        .drop_nulls(["source_feature_id", "entity_id", "distance_to_tss"])
        .with_columns(
            pl.when(pl.col("distance_to_tss").abs() <= 1500)
            .then(pl.lit("promoter"))
            .otherwise(pl.lit("body"))
            .alias("gene_region")
        )
        .filter(
            pl.col("gene_region").is_in(regions),
            pl.col("entity_id").str.len_chars() > 0,
        )
        .select(
            "source_feature_id",
            pl.lit(platform).alias("source_platform"),
            pl.lit(normalize_platform(platform)).alias("source_platform_key"),
            "entity_id",
            "gene_region",
        )
        .group_by(
            "source_feature_id",
            "source_platform",
            "source_platform_key",
            "entity_id",
        )
        .agg(
            pl.when((pl.col("gene_region") == "promoter").any())
            .then(pl.lit("promoter"))
            .otherwise(pl.lit("body"))
            .alias("gene_region")
        )
        .sort(["source_feature_id", "entity_id", "gene_region"])
    )


def _download_annotation(spec: MethylationAnnotation, destination: Path) -> None:
    partial = destination.with_suffix(destination.suffix + ".part")
    last_error: requests.RequestException | None = None
    try:
        for attempt in range(3):
            try:
                with requests.get(spec.url, stream=True, timeout=(10, 120)) as response:
                    response.raise_for_status()
                    digest = hashlib.sha256()
                    size = 0
                    with partial.open("wb") as handle:
                        for chunk in response.iter_content(chunk_size=1 << 20):
                            if not chunk:
                                continue
                            handle.write(chunk)
                            digest.update(chunk)
                            size += len(chunk)
                break
            except requests.RequestException as exc:
                last_error = exc
                if partial.exists():
                    partial.unlink()
                if attempt < 2:
                    time.sleep(2**attempt)
        else:
            assert last_error is not None
            raise last_error
        actual = digest.hexdigest()
        if actual != spec.sha256 or size != spec.size_bytes:
            raise AnnotationResolutionError(
                f"downloaded annotation failed verification for {spec.platform}: "
                f"expected sha256={spec.sha256} size={spec.size_bytes}, "
                f"got sha256={actual} size={size}"
            )
        os.replace(partial, destination)
    except (OSError, requests.RequestException) as exc:
        raise AnnotationResolutionError(
            f"failed to download annotation for {spec.platform} from {spec.url}: {exc}"
        ) from exc
    finally:
        if partial.exists():
            partial.unlink()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()
