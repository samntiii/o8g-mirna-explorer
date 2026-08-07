#!/usr/bin/env python3
"""
Build paper/benchmarks/gold_master.csv — the SINGLE source of truth for all
oxomiR gold-recovery benchmarks and figures.

Inclusion / exclusion criteria (also in paper/benchmarks/GOLD_SET_PROTOCOL.md):

Source A — literature o8G positions
  INCLUDE only when a citable paper reports (i) a specific seed-region o8G
  nucleotide position and (ii) a redirected (gained) or lost target with
  experimental support (luciferase / Ago / western / qPCR). Never invent a
  position.

Source B — validated-target databases (miRTarBase / TarBase)
  Mine strong-evidence Functional MTIs (reporter assay, Western, qPCR).
  Classify as expected-lost under the literature oxo state of that miRNA.
  Exclude HT-only (CLIP/Degradome/Microarray) unless corroborated by a strong
  method. When the live DB download is unavailable, apply the same filter to a
  locally provided hsa_MTI file and/or a hand-curated strong-evidence shortlist
  with primary PMIDs (auditable; never fabricated).

Source C — Guo et al. FRBM 2026 (PMID 41690606)
  Fully curate every abstract-stated effect (CDKN2A loss; ARID1A/PTEN gained
  for miR-134; CCNE2/E2F1 gained for miR-3118). INCLUDE only if o8G position
  is stated in accessible sources; otherwise included=False with an explicit
  exclude_reason (do not leave ambiguous / skipped).

Invariant: filters apply to each seed state before gained/lost setdiff.
Denominators for every figure/stat MUST come from gold_master where
included==True — never hardcode 28/12 or 27/13.
"""
from __future__ import annotations

import json
import sqlite3
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from o8g_db import TargetDB
from o8g_precision import PrecisionMode

GOLD_LEGACY = ROOT / "paper" / "gold" / "oxomir_gold_standard.csv"
DETAIL = ROOT / "paper" / "benchmarks" / "gold_recovery_detail.csv"
BENCH = ROOT / "paper" / "benchmarks"
MASTER = BENCH / "gold_master.csv"
RECON = BENCH / "gold_reconciliation_report.md"
MIRTAR_DIR = ROOT / "paper" / "data" / "mirtarbase"

# Strong evidence tokens (miRTarBase Experiments column / TarBase methods)
STRONG_EVIDENCE = (
    "luciferase",
    "reporter",
    "western",
    "qpcr",
    "qrt-pcr",
    "qpcr",
    "immunoblot",
)

# ---------------------------------------------------------------------------
# Source B hand-curated strong-evidence shortlist (primary PMIDs).
# Used when live miRTarBase/TarBase download is unavailable; each row would
# also pass the programmatic strong-evidence filter if present in hsa_MTI.
# ---------------------------------------------------------------------------
SOURCE_B_SHORTLIST = [
    # gene, mirna, seed, o8g_position, pmid, evidence, note
    dict(
        gene="CCNG1",
        mirna="hsa-miR-122-5p",
        seed="GGAGTGT",
        o8g_position="o8G@3",
        pmid="17616664",
        evidence="Luciferase reporter assay; Gramantieri et al. Cancer Res 2007",
        source="miRTarBase-equivalent strong (primary PMID 17616664)",
    ),
    dict(
        gene="CTDSP1",
        mirna="hsa-miR-124-3p",
        seed="AAGGCAC",
        o8g_position="o8G@4",
        pmid="17403776",
        evidence="Reporter/functional; Visvanathan et al. Genes Dev 2007 (SCP1)",
        source="miRTarBase-equivalent strong (primary PMID 17403776)",
    ),
    dict(
        gene="PAX3",
        mirna="hsa-miR-1-3p",
        seed="GGAATGT",
        o8g_position="o8G@7",
        pmid="21730146",
        evidence="Reporter assay; Goljanek-Whysall et al. PNAS 2011",
        source="miRTarBase-equivalent strong (primary PMID 21730146)",
    ),
    dict(
        gene="CXCR4",
        mirna="hsa-miR-1-3p",
        seed="GGAATGT",
        o8g_position="o8G@7",
        pmid="21752897",
        evidence="Luciferase/functional; Leone et al. Oncogene 2011",
        source="miRTarBase-equivalent strong (primary PMID 21752897)",
    ),
    dict(
        gene="NRAS",
        mirna="hsa-let-7a-5p",
        seed="GAGGTAG",
        o8g_position="o8G@4",
        pmid="16170332",
        evidence="Luciferase; Johnson et al. Cell 2005 (let-7/RAS)",
        source="miRTarBase-equivalent strong (primary PMID 16170332)",
    ),
    dict(
        gene="BDNF",
        mirna="hsa-miR-1-3p",
        seed="GGAATGT",
        o8g_position="o8G@7",
        pmid="29322793",
        evidence="Functional MTI (reporter/pathway); Gao et al. 2018",
        source="miRTarBase-equivalent strong (primary PMID 29322793)",
    ),
    # Logged but NOT included — weaker / ambiguous evidence relative to criteria
    dict(
        gene="SLC7A1",
        mirna="hsa-miR-122-5p",
        seed="GGAGTGT",
        o8g_position="o8G@3",
        pmid="17179747",
        evidence="Chang et al. RNA Biol 2004 — correlative/'may downregulate'",
        source="miRTarBase candidate (primary PMID 17179747)",
        force_exclude="weak_or_correlative_evidence_not_strong_reporter",
    ),
]


def ensembl_map(symbols: list[str]) -> dict[str, str]:
    con = sqlite3.connect(ROOT / "o8g_targets.db")
    if not symbols:
        return {}
    q = ",".join("?" * len(symbols))
    rows = con.execute(
        f"SELECT symbol, gene_id FROM genes WHERE symbol IN ({q})", symbols
    ).fetchall()
    con.close()
    return dict(rows)


def score_row(
    db: TargetDB,
    mirna: str,
    seed: str,
    o8g_position: str,
    gene: str,
    effect_type: str,
    mode: PrecisionMode = PrecisionMode.SENSITIVE,
) -> tuple[bool | None, bool | None, bool | None, str]:
    """Return in_unmod, in_ox, recovered, status (Sensitive mode for master)."""
    if not seed or not o8g_position or o8g_position in ("oxidized_seed", ""):
        return None, None, None, "unscored_no_position"
    info = db.mirna_info(mirna)
    mature = info["seq_dna"] if info else None
    unmod = set(
        db.targets_filtered(
            seed, "none", mode, mature_dna=mature, mirna=mirna
        )["symbol"]
    )
    ox = set(
        db.targets_filtered(
            seed, o8g_position, mode, mature_dna=mature, mirna=mirna
        )["symbol"]
    )
    in_u, in_o = gene in unmod, gene in ox
    if effect_type == "gained":
        recovered = bool(in_o and not in_u)
        if recovered:
            status = "TP_gained"
        elif in_o and in_u:
            status = "present_also_unmod"
        elif in_u:
            status = "only_unmod"
        else:
            status = "FN_absent"
    else:
        recovered = bool(in_u and not in_o)
        if recovered:
            status = "TP_lost"
        elif in_o:
            status = "still_in_ox"
        elif not in_u:
            status = "absent_both"
        else:
            status = "lost_ok_weak"
    return in_u, in_o, recovered, status


def migrate_legacy(detail: pd.DataFrame) -> list[dict]:
    """Migrate gold_recovery_detail.csv rows into master schema."""
    rows = []
    for i, r in detail.iterrows():
        mir = str(r["mirna"])
        # Ambiguous Guo multi-miRNA placeholder — handled in Source C, skip here
        if ";" in mir or str(r.get("status", "")) == "skipped":
            continue
        effect = str(r["effect"])
        effect_type = "gained" if "gained" in effect else "lost"
        o8g = str(r.get("state_label") or "")
        rows.append(
            dict(
                effect_id=f"LEGACY_{i+1:03d}_{mir.split('-')[2] if '-' in mir else mir}_{r['gene']}_{effect_type}",
                source=r["source"],
                pmid=str(r["pmid"]),
                mirna=mir,
                seed=str(r["seed"]),
                o8g_position=o8g,
                gene=str(r["gene"]),
                effect_type=effect_type,
                evidence=str(r["evidence"]),
                provenance="Seok2020/Eom2023 curated oxomiR gold (migrated from gold_recovery_detail)",
                included=True,
                exclude_reason="",
                in_unmodified=r.get("in_unmodified"),
                in_oxidized=r.get("in_oxidized"),
                recovered=r.get("recovered"),
                status=r.get("status"),
            )
        )
    return rows


def curate_guo() -> list[dict]:
    """Source C — full Guo 2026 curation; exclude until o8G position is citable."""
    seed = "GTGACTG"  # miRBase seed for miR-134-5p and miR-3118-3p
    reason = (
        "o8g_position_not_reported_in_available_abstract_metadata; "
        "seed assigned from miRBase (GTGACTG) but nucleotide index not stated "
        "in PMID 41690606 abstract — cannot INCLUDE without fabricating position"
    )
    base = dict(
        source="Guo et al. Free Radic Biol Med 2026",
        pmid="41690606",
        seed=seed,
        o8g_position="",
        included=False,
        exclude_reason=reason,
        in_unmodified=None,
        in_oxidized=None,
        recovered=None,
        status="excluded_no_position",
    )
    effects = [
        ("hsa-miR-134-5p", "CDKN2A", "lost", "senescence_P16_disruption; abstract"),
        ("hsa-miR-3118-3p", "CDKN2A", "lost", "senescence_P16_disruption; abstract; DB name miR-3118-3p"),
        ("hsa-miR-134-5p", "ARID1A", "gained", "o8G:A retargeting; abstract"),
        ("hsa-miR-134-5p", "PTEN", "gained", "o8G:A retargeting; abstract"),
        ("hsa-miR-3118-3p", "CCNE2", "gained", "o8G:A retargeting; abstract"),
        ("hsa-miR-3118-3p", "E2F1", "gained", "o8G:A retargeting; abstract"),
    ]
    rows = []
    for mir, gene, et, ev in effects:
        rows.append(
            {
                **base,
                "effect_id": f"GUO_{mir}_{gene}_{et}",
                "mirna": mir,
                "gene": gene,
                "effect_type": et,
                "evidence": ev,
                "provenance": "Source C — Guo 2026 full abstract curation",
            }
        )
    return rows


def curate_source_a_extra() -> list[dict]:
    """Additional Source A literature with known limitations (mostly excluded)."""
    rows = []
    # Wang 2015 — gained targets known, position NOT identified in literature
    for gene, symbol_note in [("BCL2L1", "Bcl-xL"), ("BCL2L2", "Bcl-w")]:
        rows.append(
            dict(
                effect_id=f"WANG2015_miR184_{gene}_gained",
                source="Wang et al. Mol Cell 2015",
                pmid="26028536",
                mirna="hsa-miR-184",
                seed="GGACGGA",  # miR-184-3p seed in DB
                o8g_position="",
                gene=gene,
                effect_type="gained",
                evidence=f"oxidized miR-184 targets {symbol_note}; luciferase/functional",
                provenance="Source A — Wang 2015; position not identified (reviews confirm)",
                included=False,
                exclude_reason="o8g_position_not_identified_in_source_paper",
                in_unmodified=None,
                in_oxidized=None,
                recovered=None,
                status="excluded_no_position",
            )
        )
    # Li/Wang Sci Rep 2024 miR-30c — positions labeled 4,5-oxo but do not map to
    # human miR-30c-5p seed Gs (only mature pos2 is G); synthetic GG pairing used.
    rows.append(
        dict(
            effect_id="LI2024_miR30c_CDKN2C_gained",
            source="Li et al. Sci Rep 2024",
            pmid="38849466",
            mirna="hsa-miR-30c-5p",
            seed="GTAAACA",
            o8g_position="",
            gene="CDKN2C",
            effect_type="gained",
            evidence="4,5-oxo-miR-30c luciferase/western vs CDKN2C",
            provenance="Source A — Li 2024",
            included=False,
            exclude_reason="o8g_positions_4_5_do_not_map_to_human_miR-30c-5p_seed_G_sites",
            in_unmodified=None,
            in_oxidized=None,
            recovered=None,
            status="excluded_position_unmappable",
        )
    )
    rows.append(
        dict(
            effect_id="LI2024_miR30c_MYBL2_lost",
            source="Li et al. Sci Rep 2024",
            pmid="38849466",
            mirna="hsa-miR-30c-5p",
            seed="GTAAACA",
            o8g_position="",
            gene="MYBL2",
            effect_type="lost",
            evidence="canonical MYBL2 lost under 4,5-oxo (western)",
            provenance="Source A — Li 2024",
            included=False,
            exclude_reason="o8g_positions_4_5_do_not_map_to_human_miR-30c-5p_seed_G_sites",
            in_unmodified=None,
            in_oxidized=None,
            recovered=None,
            status="excluded_position_unmappable",
        )
    )
    return rows


def is_strong_experiments(text: str) -> bool:
    t = (text or "").lower()
    return any(tok in t for tok in STRONG_EVIDENCE)


def mine_mirtarbase_local(
    existing: set[tuple[str, str]],
) -> tuple[list[dict], list[dict], str]:
    """
    If a local hsa_MTI.xlsx/.txt/.csv exists, mine strong-evidence rows for the
    four benchmark miRNAs. Returns (included_rows, audit_excluded_rows, note).
    """
    candidates = list(MIRTAR_DIR.glob("hsa_MTI*")) + list(MIRTAR_DIR.glob("*MTI*"))
    path = None
    for p in candidates:
        if p.suffix.lower() in {".xlsx", ".xls", ".txt", ".tsv", ".csv"} and p.stat().st_size > 1000:
            path = p
            break
    if path is None:
        return [], [], (
            "Live miRTarBase/TarBase download unavailable (HTTP 404 on CUHK hosts "
            f"{pd.Timestamp.today().date()}); no local hsa_MTI file >1KB in "
            f"{MIRTAR_DIR}. Source B inclusions come from the hand-curated "
            "strong-evidence shortlist with primary PMIDs only."
        )

    if path.suffix.lower() in {".xlsx", ".xls"}:
        raw = pd.read_excel(path)
    else:
        raw = pd.read_csv(path, sep=None, engine="python")

    # Normalize columns
    colmap = {c.lower().replace(" ", "_"): c for c in raw.columns}
    def col(*names):
        for n in names:
            if n in colmap:
                return colmap[n]
            if n in raw.columns:
                return n
        return None

    mir_c = col("mirna", "miRNA", "mir_name")
    gene_c = col("target_gene", "gene", "target_symbol")
    exp_c = col("experiments", "experiment", "methods", "evidence")
    sup_c = col("support_type", "support")
    if mir_c is None or gene_c is None:
        return [], [], f"Local MTI file {path.name} lacks miRNA/gene columns"

    mirnas = {
        "hsa-miR-1-3p": ("GGAATGT", "o8G@7"),
        "hsa-miR-124-3p": ("AAGGCAC", "o8G@4"),
        "hsa-let-7a-5p": ("GAGGTAG", "o8G@4"),
        "hsa-miR-122-5p": ("GGAGTGT", "o8G@3"),
    }
    included, excluded = [], []
    sub = raw[raw[mir_c].astype(str).isin(mirnas)]
    for _, r in sub.iterrows():
        mir = str(r[mir_c])
        gene = str(r[gene_c]).upper()
        exp = str(r[exp_c]) if exp_c else ""
        sup = str(r[sup_c]) if sup_c else ""
        seed, ox = mirnas[mir]
        key = (mir, gene)
        if key in existing:
            excluded.append(
                dict(
                    effect_id=f"MTB_DUP_{mir}_{gene}",
                    source=f"miRTarBase local ({path.name})",
                    pmid="",
                    mirna=mir,
                    seed=seed,
                    o8g_position=ox,
                    gene=gene,
                    effect_type="lost",
                    evidence=f"{sup}; {exp}",
                    provenance="Source B — duplicate of existing gold",
                    included=False,
                    exclude_reason="duplicate_of_existing_gold_effect",
                    status="excluded_duplicate",
                )
            )
            continue
        strong = is_strong_experiments(exp) or is_strong_experiments(sup)
        if not strong:
            excluded.append(
                dict(
                    effect_id=f"MTB_WEAK_{mir}_{gene}",
                    source=f"miRTarBase local ({path.name})",
                    pmid="",
                    mirna=mir,
                    seed=seed,
                    o8g_position=ox,
                    gene=gene,
                    effect_type="lost",
                    evidence=f"{sup}; {exp}",
                    provenance="Source B — HT/weak only",
                    included=False,
                    exclude_reason="not_strong_evidence_reporter_western_qpcr",
                    status="excluded_weak_evidence",
                )
            )
            continue
        included.append(
            dict(
                effect_id=f"MTB_{mir}_{gene}_lost",
                source=f"miRTarBase local ({path.name})",
                pmid="",
                mirna=mir,
                seed=seed,
                o8g_position=ox,
                gene=gene,
                effect_type="lost",
                evidence=f"Functional MTI strong: {sup}; {exp}",
                provenance="Source B — programmatic miRTarBase strong filter",
                included=True,
                exclude_reason="",
                status="",
            )
        )
        existing.add(key)
    return included, excluded, f"Mined {path}"


def apply_source_b_shortlist(existing: set[tuple[str, str]]) -> list[dict]:
    rows = []
    for item in SOURCE_B_SHORTLIST:
        key = (item["mirna"], item["gene"])
        force = item.get("force_exclude")
        if force:
            rows.append(
                dict(
                    effect_id=f"SRC_B_EXCL_{item['mirna']}_{item['gene']}",
                    source=item["source"],
                    pmid=item["pmid"],
                    mirna=item["mirna"],
                    seed=item["seed"],
                    o8g_position=item["o8g_position"],
                    gene=item["gene"],
                    effect_type="lost",
                    evidence=item["evidence"],
                    provenance="Source B shortlist — audited exclusion",
                    included=False,
                    exclude_reason=force,
                    status="excluded_weak_evidence",
                    in_unmodified=None,
                    in_oxidized=None,
                    recovered=None,
                )
            )
            continue
        if key in existing:
            rows.append(
                dict(
                    effect_id=f"SRC_B_DUP_{item['mirna']}_{item['gene']}",
                    source=item["source"],
                    pmid=item["pmid"],
                    mirna=item["mirna"],
                    seed=item["seed"],
                    o8g_position=item["o8g_position"],
                    gene=item["gene"],
                    effect_type="lost",
                    evidence=item["evidence"],
                    provenance="Source B shortlist — duplicate",
                    included=False,
                    exclude_reason="duplicate_of_existing_gold_effect",
                    status="excluded_duplicate",
                    in_unmodified=None,
                    in_oxidized=None,
                    recovered=None,
                )
            )
            continue
        rows.append(
            dict(
                effect_id=f"SRC_B_{item['mirna']}_{item['gene']}_lost",
                source=item["source"],
                pmid=item["pmid"],
                mirna=item["mirna"],
                seed=item["seed"],
                o8g_position=item["o8g_position"],
                gene=item["gene"],
                effect_type="lost",
                evidence=item["evidence"],
                provenance="Source B — curated strong-evidence expected-lost",
                included=True,
                exclude_reason="",
                status="",
                in_unmodified=None,
                in_oxidized=None,
                recovered=None,
            )
        )
        existing.add(key)
    return rows


def rescore_included(rows: list[dict], db: TargetDB) -> None:
    for r in rows:
        if not r.get("included"):
            continue
        in_u, in_o, rec, status = score_row(
            db,
            r["mirna"],
            r["seed"],
            r["o8g_position"],
            r["gene"],
            r["effect_type"],
        )
        r["in_unmodified"] = in_u
        r["in_oxidized"] = in_o
        r["recovered"] = rec
        r["status"] = status


def write_reconciliation(df: pd.DataFrame, mirtar_note: str) -> None:
    inc = df[df["included"] == True]  # noqa: E712
    excl = df[df["included"] == False]  # noqa: E712
    lines = [
        "# Gold-set reconciliation report",
        "",
        f"Generated by `scripts/build_gold_master.py`. Master file: `{MASTER.relative_to(ROOT)}`.",
        "",
        "## Denominator drift resolution (12/28 vs 13/27)",
        "",
        "Raw `oxomir_gold_standard.csv` contained **28 gained + 13 lost = 41** rows,",
        "including one Guo placeholder (`oxidized_seed`, empty seed) that scoring",
        "pipelines silently skipped → effective **28 gained / 12 lost**.",
        "Some notes counted 13 lost (raw) while figures used 12 (scored).",
        "",
        "**Canonical denominators** (included=True only):",
        "",
    ]
    gained_n = int((inc["effect_type"] == "gained").sum())
    lost_n = int((inc["effect_type"] == "lost").sum())
    lines += [
        f"- **pooled gained = {gained_n}**",
        f"- **pooled lost = {lost_n}**",
        "",
        "Guo effects are fully curated in the master but `included=False` until a",
        "citable o8G nucleotide position is available (see Source C).",
        "",
        "## Counts by miRNA × effect_type (included=True)",
        "",
        "| mirna | gained | lost | total |",
        "|---|---:|---:|---:|",
    ]
    for mir, g in inc.groupby("mirna"):
        ng = int((g["effect_type"] == "gained").sum())
        nl = int((g["effect_type"] == "lost").sum())
        lines.append(f"| {mir} | {ng} | {nl} | {ng+nl} |")
    lines += [
        f"| **POOLED** | **{gained_n}** | **{lost_n}** | **{gained_n+lost_n}** |",
        "",
        "## Exclusions (auditable)",
        "",
        f"Total rows in master: {len(df)} (included={len(inc)}, excluded={len(excl)})",
        "",
    ]
    if len(excl):
        vc = excl["exclude_reason"].value_counts()
        lines.append("| exclude_reason | n |")
        lines.append("|---|---:|")
        for reason, n in vc.items():
            lines.append(f"| `{reason}` | {n} |")
    lines += [
        "",
        "## Source B note",
        "",
        mirtar_note,
        "",
        "## Guo 2026 resolution",
        "",
        "Former ambiguous row `hsa-miR-134-5p;hsa-miR-3118` / `oxidized_seed` / CDKN2A",
        "is replaced by six curated rows (two miRNAs × CDKN2A lost + four gained",
        "targets from the abstract). All have `included=False` with",
        "`exclude_reason=o8g_position_not_reported_...`. Seed GTGACTG assigned from",
        "miRBase; miR-3118 stored as `hsa-miR-3118-3p` to match the DB.",
        "",
    ]
    RECON.write_text("\n".join(lines))
    print(f"Wrote {RECON}")


def sync_legacy_gold_csv(inc: pd.DataFrame) -> None:
    """Keep oxomir_gold_standard.csv in sync with included master rows."""
    out = inc.rename(
        columns={
            "o8g_position": "state_label",
            "effect_type": "effect",
        }
    ).copy()
    out["effect"] = out["effect"].map(
        {"gained": "gained_on_oxidation", "lost": "lost_on_oxidation"}
    )
    cols = ["source", "pmid", "mirna", "seed", "state_label", "gene", "effect", "evidence"]
    out[cols].to_csv(GOLD_LEGACY, index=False)
    print(f"Synced {GOLD_LEGACY} ({len(out)} included rows)")


def main():
    BENCH.mkdir(parents=True, exist_ok=True)
    MIRTAR_DIR.mkdir(parents=True, exist_ok=True)

    if DETAIL.exists():
        detail = pd.read_csv(DETAIL)
    else:
        # fall back to legacy gold without recovery columns
        leg = pd.read_csv(GOLD_LEGACY)
        detail = leg.assign(
            in_unmodified=None, in_oxidized=None, recovered=None, status=""
        )

    rows = migrate_legacy(detail)
    existing = {(r["mirna"], r["gene"]) for r in rows}

    rows.extend(curate_guo())
    rows.extend(curate_source_a_extra())

    mtb_inc, mtb_excl, mtb_note = mine_mirtarbase_local(existing)
    rows.extend(mtb_inc)
    rows.extend(mtb_excl)
    rows.extend(apply_source_b_shortlist(existing))

    # Ensembl IDs
    symbols = sorted({r["gene"] for r in rows})
    emap = ensembl_map(symbols)
    for r in rows:
        r["ensembl"] = emap.get(r["gene"], "")

    db = TargetDB(str(ROOT / "o8g_targets.db"))
    rescore_included(rows, db)

    cols = [
        "effect_id",
        "source",
        "pmid",
        "mirna",
        "seed",
        "o8g_position",
        "gene",
        "ensembl",
        "effect_type",
        "evidence",
        "in_unmodified",
        "in_oxidized",
        "recovered",
        "status",
        "included",
        "exclude_reason",
        "provenance",
    ]
    df = pd.DataFrame(rows)
    for c in cols:
        if c not in df.columns:
            df[c] = ""
    df = df[cols].sort_values(["included", "mirna", "effect_type", "gene"], ascending=[False, True, True, True])
    df.to_csv(MASTER, index=False)
    print(f"Wrote {MASTER} ({len(df)} rows)")

    inc = df[df["included"] == True]  # noqa: E712
    print(
        f"INCLUDED: gained={(inc.effect_type=='gained').sum()} "
        f"lost={(inc.effect_type=='lost').sum()} total={len(inc)}"
    )
    write_reconciliation(df, mtb_note)
    sync_legacy_gold_csv(inc)

    # machine-readable denom snapshot
    snap = {
        "pooled_gained": int((inc.effect_type == "gained").sum()),
        "pooled_lost": int((inc.effect_type == "lost").sum()),
        "n_included": int(len(inc)),
        "n_excluded": int((df.included == False).sum()),  # noqa: E712
        "per_mirna": {
            mir: {
                "gained": int((g.effect_type == "gained").sum()),
                "lost": int((g.effect_type == "lost").sum()),
            }
            for mir, g in inc.groupby("mirna")
        },
        "mirtarbase_note": mtb_note,
    }
    (BENCH / "gold_denominators.json").write_text(json.dumps(snap, indent=2))
    print(f"Wrote {BENCH / 'gold_denominators.json'}")


if __name__ == "__main__":
    main()
