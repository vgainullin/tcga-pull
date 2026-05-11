"""Tests for the per-lineage frequency aggregations.

Built around a tiny hand-crafted cohort where the expected counts are
trivial to verify by eye:

  Lineage "BRCA": 4 patients
    P1: PIK3CA H1047R, TP53 R175H
    P2: PIK3CA H1047R, TTN p.X (silent — should be filtered out)
    P3: PIK3CA E545K
    P4: (no mutations)
  Lineage "LUAD": 3 patients
    P5: KRAS G12C
    P6: TP53 R175H
    P7: (no mutations)

Total: 7 patients, 2 lineages.

Gene rates we expect:
  PIK3CA in BRCA: 3/4 = 0.75, in LUAD: 0/3 = 0.0
  TP53 in BRCA:   1/4 = 0.25, in LUAD: 1/3 ≈ 0.333
  KRAS in LUAD:   1/3 ≈ 0.333, in BRCA: 0/4 = 0.0
  TTN: filtered (silent only) → not in coding output
"""

from __future__ import annotations

import math

import polars as pl

from tcga_pull.frequency import gene_frequency, variant_frequency


def _fixture() -> tuple[pl.DataFrame, pl.DataFrame]:
    rows = []

    def variant(*, submitter, lineage, hugo, pos, variant_class, impact, gnomad_af=None):
        # `lineage` is bound here only so the test can reason per-lineage; it is
        # NOT placed on the variants row because real variants.parquet doesn't
        # have lineage — the frequency code joins it from samples.parquet.
        del lineage  # silence unused-arg lint
        return {
            "submitter_id": submitter,
            "hugo_symbol": hugo,
            "chrom": "chr1",
            "pos": pos,
            "ref": "C",
            "alt": "T",
            "hgvsp_short": f"p.X{pos}Y",
            "variant_class": variant_class,
            "consequence": "missense_variant",
            "impact": impact,
            "primary_aliquot": True,
            "is_coding": variant_class != "Silent",
            "is_high_impact": impact in ("HIGH", "MODERATE"),
            "is_rare": gnomad_af is None or gnomad_af < 1e-3,
            "gnomad_af": gnomad_af,
        }

    # BRCA — 4 patients
    rows.append(
        variant(
            submitter="P1",
            lineage="BRCA",
            hugo="PIK3CA",
            pos=100,
            variant_class="Missense_Mutation",
            impact="MODERATE",
        )
    )
    rows.append(
        variant(
            submitter="P1",
            lineage="BRCA",
            hugo="TP53",
            pos=200,
            variant_class="Missense_Mutation",
            impact="HIGH",
        )
    )
    rows.append(
        variant(
            submitter="P2",
            lineage="BRCA",
            hugo="PIK3CA",
            pos=100,
            variant_class="Missense_Mutation",
            impact="MODERATE",
        )
    )
    rows.append(
        variant(
            submitter="P2",
            lineage="BRCA",
            hugo="TTN",
            pos=300,
            variant_class="Silent",
            impact="LOW",
        )
    )
    rows.append(
        variant(
            submitter="P3",
            lineage="BRCA",
            hugo="PIK3CA",
            pos=400,
            variant_class="Missense_Mutation",
            impact="MODERATE",
        )
    )
    # P4 has no variants

    # LUAD — 3 patients
    rows.append(
        variant(
            submitter="P5",
            lineage="LUAD",
            hugo="KRAS",
            pos=500,
            variant_class="Missense_Mutation",
            impact="MODERATE",
        )
    )
    rows.append(
        variant(
            submitter="P6",
            lineage="LUAD",
            hugo="TP53",
            pos=200,
            variant_class="Missense_Mutation",
            impact="HIGH",
        )
    )
    # P7 has no variants

    variants = pl.DataFrame(rows)

    samples = pl.DataFrame(
        [
            {"submitter_id": "P1", "lineage": "BRCA"},
            {"submitter_id": "P2", "lineage": "BRCA"},
            {"submitter_id": "P3", "lineage": "BRCA"},
            {"submitter_id": "P4", "lineage": "BRCA"},
            {"submitter_id": "P5", "lineage": "LUAD"},
            {"submitter_id": "P6", "lineage": "LUAD"},
            {"submitter_id": "P7", "lineage": "LUAD"},
        ]
    )
    return variants, samples


def _row(df: pl.DataFrame, gene: str, lineage: str) -> dict:
    rows = df.filter((pl.col("hugo_symbol") == gene) & (pl.col("lineage") == lineage)).to_dicts()
    assert len(rows) == 1, f"expected 1 row for ({gene}, {lineage}); got {len(rows)}"
    return rows[0]


# ------------------------------------------------------------------- gene


def test_gene_frequency_pik3ca_brca():
    variants, samples = _fixture()
    gf = gene_frequency(variants, samples)
    r = _row(gf, "PIK3CA", "BRCA")
    assert r["n_mutated_patients"] == 3
    assert r["n_total_patients"] == 4
    assert math.isclose(r["freq"], 3 / 4)
    # PIK3CA is BRCA-only; rate in "other lineages" = 0
    assert r["n_mutated_other"] == 0
    assert math.isclose(r["freq_other_lineages"], 0.0)
    # log2((0.75 + ε) / (0 + ε)) → very large positive number
    assert r["log2_enrichment_vs_other"] > 10


def test_gene_frequency_tp53_pancancer():
    variants, samples = _fixture()
    gf = gene_frequency(variants, samples)
    # BRCA side
    rb = _row(gf, "TP53", "BRCA")
    assert rb["n_mutated_patients"] == 1 and rb["n_total_patients"] == 4
    assert math.isclose(rb["freq"], 1 / 4)
    assert rb["n_mutated_other"] == 1  # the LUAD TP53 carrier
    assert rb["n_total_other"] == 3
    assert math.isclose(rb["freq_other_lineages"], 1 / 3)
    # LUAD side
    rl = _row(gf, "TP53", "LUAD")
    assert math.isclose(rl["freq"], 1 / 3)
    assert math.isclose(rl["freq_other_lineages"], 1 / 4)


def test_gene_frequency_silent_filtered_out():
    variants, samples = _fixture()
    gf = gene_frequency(variants, samples)
    # TTN only has a Silent mutation in the fixture; must not appear.
    assert gf.filter(pl.col("hugo_symbol") == "TTN").is_empty()


def test_gene_frequency_high_impact_counts():
    variants, samples = _fixture()
    gf = gene_frequency(variants, samples)
    # PIK3CA in BRCA — all 3 patients are MODERATE = high-impact
    r = _row(gf, "PIK3CA", "BRCA")
    assert r["n_high_impact_patients"] == 3
    # TP53 in BRCA — 1 patient with HIGH = high-impact
    r = _row(gf, "TP53", "BRCA")
    assert r["n_high_impact_patients"] == 1


# ----------------------------------------------------------------- variant


def test_variant_frequency_basic_counts():
    variants, samples = _fixture()
    vf = variant_frequency(variants, samples)
    # PIK3CA pos 100 in BRCA — 2 patients out of 4
    pik = vf.filter(
        (pl.col("hugo_symbol") == "PIK3CA") & (pl.col("pos") == 100) & (pl.col("lineage") == "BRCA")
    ).to_dicts()
    assert len(pik) == 1
    assert pik[0]["n_patients_with_variant"] == 2
    assert pik[0]["n_total_patients"] == 4
    assert math.isclose(pik[0]["cohort_freq"], 2 / 4)


def test_variant_frequency_no_is_rare_filter():
    """variant_frequency keeps common variants too (gnomad_af is the comparison axis)."""
    variants, samples = _fixture()
    # add a definitely-common variant
    rows = [
        *variants.to_dicts(),
        {
            "submitter_id": "P1",
            "hugo_symbol": "COMMON",
            "chrom": "chr2",
            "pos": 999,
            "ref": "C",
            "alt": "T",
            "hgvsp_short": "p.X1Y",
            "variant_class": "Missense_Mutation",
            "consequence": "missense_variant",
            "impact": "MODERATE",
            "primary_aliquot": True,
            "is_coding": True,
            "is_high_impact": True,
            "is_rare": False,  # common variant
            "gnomad_af": 0.05,
        },
    ]
    vf = variant_frequency(pl.DataFrame(rows), samples)
    common = vf.filter(pl.col("hugo_symbol") == "COMMON")
    assert len(common) == 1, "common (non-rare) variant should still appear"
    assert math.isclose(common["gnomad_af"][0], 0.05)
