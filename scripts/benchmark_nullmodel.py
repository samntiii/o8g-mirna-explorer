#!/usr/bin/env python3
"""
Null / permutation baselines for oxomiR gold-standard recovery (rigor upgrade).

Includes
--------
1. Per-miRNA composition-matched decoy-seed null (existing).
2. Pooled-across-miRNA decoy null (headline): sum per-miRNA decoy recoveries
   within each permutation so near-zero null means cannot invent 40× folds.
3. Orthogonal target-label permutation null: shuffle gold labels over the UTR
   universe (~19,159 genes), keep |gold| fixed, score recovery by the *real*
   predicted set. If decoy_p and labelperm_p agree at α=0.05 → agreement_flag.
   If they diverge, composition-matching is load-bearing (stated in output notes).
4. Wilson score 95% CIs on every recovery fraction (no Wald / normal approx).

Invariant: PrecisionMode filters apply to each seed state before gained/lost
setdiff (same as gold_recovery_by_mode / targets_filtered).

Outputs
-------
  paper/benchmarks/nullmodel_recovery.csv          (per-miRNA; Wilson CIs)
  paper/benchmarks/nullmodel_pooled.csv            (headline pooled)
  paper/benchmarks/gold_recovery_summary_by_mode.csv (Wilson CIs)
  paper/figures/figS_nullmodel.pdf
  paper/figures/figS_pooled_null.pdf

Env: O8G_NULL_N (default 10000), O8G_NULL_SEED (default 42)
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from o8g_db import TargetDB
from o8g_engine import SeedState, clean_seq, g_positions
from o8g_precision import PrecisionMode, PrecisionConfig, apply_precision_filter
from o8g_scanner import TargetScanner, SITE_RANK

BENCH = ROOT / "paper" / "benchmarks"
GOLD_MASTER = ROOT / "paper" / "benchmarks" / "gold_master.csv"
GOLD_LEGACY = ROOT / "paper" / "gold" / "oxomir_gold_standard.csv"
# Denominators MUST come from gold_master (included=True). Legacy CSV is a sync mirror.
FIGS = ROOT / "paper" / "figures"
UTR = ROOT / "utr3_human.parquet"


# ---------------------------------------------------------------------------
# Wilson score interval (Newcombe 1998; preferred for small n binomial)
# ---------------------------------------------------------------------------
def wilson_ci(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    """95% Wilson score CI for binomial proportion k/n.

    Cite: Wilson EB (1927) JASA; Newcombe RG (1998) Stat Med.
    Do NOT substitute a Wald interval here — n is small.
    """
    if n <= 0:
        return float("nan"), float("nan")
    phat = k / n
    z2 = z * z
    denom = 1.0 + z2 / n
    centre = (phat + z2 / (2.0 * n)) / denom
    half = (z / denom) * np.sqrt(phat * (1.0 - phat) / n + z2 / (4.0 * n * n))
    return float(max(0.0, centre - half)), float(min(1.0, centre + half))


def parse_ox_positions(label: str) -> tuple[int, ...]:
    if not label or label == "none":
        return ()
    return tuple(int(x) for x in label.replace("o8G@", "").split(",") if x)


def composition_shuffle(seed: str, rng: np.random.Generator) -> str:
    """Exact mononucleotide-composition shuffle (dinucleotide NOT Markov-matched)."""
    chars = list(clean_seq(seed))
    rng.shuffle(chars)
    return "".join(chars)


def decoy_ox_positions(
    decoy_seed: str, n_ox: int, rng: np.random.Generator
) -> tuple[int, ...]:
    gpos = g_positions(decoy_seed)
    if n_ox <= 0:
        return ()
    if len(gpos) < n_ox:
        raise ValueError(f"decoy has {len(gpos)} Gs but need {n_ox} oxidized")
    chosen = sorted(rng.choice(np.array(gpos), size=n_ox, replace=False).tolist())
    return tuple(int(x) for x in chosen)


def best_rank_in_utr(utr: str, state: SeedState) -> int:
    utr = utr.upper().replace("U", "T")
    best = 0
    for stype, motif in state.motifs.items():
        if motif and motif in utr:
            best = max(best, SITE_RANK[stype])
    return best


def gold_utr_table(scanner: TargetScanner, symbols: set[str]) -> pd.DataFrame:
    sym_to_i = {s: i for i, s in enumerate(scanner.symbols)}
    rows = []
    for s in sorted(symbols):
        if s not in sym_to_i:
            continue
        i = sym_to_i[s]
        rows.append(dict(symbol=s, gene_id=scanner.genes[i], utr=scanner.utrs[i]))
    return pd.DataFrame(rows)


def filtered_symbols_on_gold(
    gold_utrs: pd.DataFrame,
    state: SeedState,
    mode: PrecisionMode,
    *,
    conserved: set[str] | None,
    is_unmodified: bool,
) -> set[str]:
    if gold_utrs.empty:
        return set()
    ranks = [best_rank_in_utr(u, state) for u in gold_utrs["utr"]]
    df = pd.DataFrame(
        {
            "symbol": gold_utrs["symbol"].values,
            "site_rank": ranks,
            "is_conserved": [s in (conserved or set()) for s in gold_utrs["symbol"]],
        }
    )
    df = df[df["site_rank"] > 0]
    cfg = PrecisionConfig.from_mode(mode)
    filt = apply_precision_filter(
        df, cfg, conserved_symbols=conserved, is_unmodified_state=is_unmodified
    )
    return set(filt["symbol"]) if len(filt) else set()


def partition_from_states(
    gold_utrs: pd.DataFrame,
    seed: str,
    ox_pos: tuple[int, ...],
    mode: PrecisionMode,
    conserved: set[str],
) -> dict[str, set[str]]:
    """Per-state mode filter then setdiff — matches gold_recovery_by_mode."""
    unmod_state = SeedState(seed, ())
    ox_state = SeedState(seed, ox_pos)
    u = filtered_symbols_on_gold(
        gold_utrs, unmod_state, mode, conserved=conserved, is_unmodified=True
    )
    o = filtered_symbols_on_gold(
        gold_utrs, ox_state, mode, conserved=conserved, is_unmodified=False
    )
    return {"gained": o - u, "lost": u - o, "shared": u & o}


def recovery_count(parts: dict[str, set[str]], gold_genes: set[str], effect: str) -> int:
    return len(parts[effect] & gold_genes)


def hypergeom_p(k: int, n_draw: int, K: int, N: int) -> float:
    if N <= 0 or n_draw < 0 or K <= 0:
        return float("nan")
    if n_draw == 0:
        return 1.0 if k == 0 else 0.0
    return float(hypergeom.sf(k - 1, N, K, n_draw)) if k > 0 else 1.0


def load_gold_by_mirna() -> dict[str, dict]:
    """Load included effects from gold_master.csv (canonical denominators).

    Falls back to legacy oxomir_gold_standard.csv only if master is missing
    (should not happen in production — run scripts/build_gold_master.py first).
    """
    if GOLD_MASTER.exists():
        gold = pd.read_csv(GOLD_MASTER)
        gold = gold[gold["included"].astype(str).str.lower().isin(["true", "1"])].copy()
        # normalize column names
        if "effect_type" in gold.columns:
            gold["effect"] = gold["effect_type"].map(
                {"gained": "gained_on_oxidation", "lost": "lost_on_oxidation"}
            ).fillna(gold["effect_type"])
        if "o8g_position" in gold.columns:
            gold["state_label"] = gold["o8g_position"]
    else:
        gold = pd.read_csv(GOLD_LEGACY)
        gold = gold[gold["state_label"].notna() & (gold["state_label"] != "oxidized_seed")]
    gold = gold[gold["state_label"].notna() & (gold["state_label"].astype(str) != "") & (gold["state_label"] != "oxidized_seed")]
    out: dict[str, dict] = {}
    for mir, g in gold.groupby("mirna"):
        out[mir] = {
            "gained": set(g.loc[g["effect"] == "gained_on_oxidation", "gene"]),
            "lost": set(g.loc[g["effect"] == "lost_on_oxidation", "gene"]),
            "seed": g["seed"].iloc[0],
            "ox_label": g["state_label"].iloc[0],
        }
    return out


def conserved_for(db: TargetDB, mirna: str, seed: str) -> set[str]:
    try:
        from conservation import get_conserved_index, build_seed_family_map

        idx = get_conserved_index()
        fam = build_seed_family_map(db.path).get(seed)
        return set(idx.conserved_symbols_for_mirna(mirna, fam))
    except Exception:
        return set()


def observed_partition(db: TargetDB, mirna: str, seed: str, ox_label: str, mode: PrecisionMode):
    info = db.mirna_info(mirna)
    return db.retarget_partition(
        seed,
        ox_label,
        mode,
        mature_dna=info["seq_dna"] if info else None,
        mirna=mirna,
    )


def label_perm_null(
    pred_set: set[str],
    K: int,
    universe: np.ndarray,
    n_perm: int,
    rng: np.random.Generator,
    observed: int,
) -> tuple[np.ndarray, float]:
    """Shuffle gold labels over UTR universe; score recovery by fixed pred_set.

    Orthogonal to decoy-seed null: prediction chemistry stays fixed; only the
    identity of which K genes are called 'gold' is randomized. If this null and
    the decoy-seed null disagree on α=0.05 significance, composition-matching
    of the seed is load-bearing for the claim.
    """
    N = len(universe)
    if K <= 0 or N <= 0:
        return np.zeros(n_perm, dtype=np.int32), 1.0
    # index predicted genes in universe for fast membership
    pred_mask = np.isin(universe, np.array(list(pred_set), dtype=object))
    # sampling without replacement of K indices each perm
    counts = np.empty(n_perm, dtype=np.int32)
    for i in range(n_perm):
        idx = rng.choice(N, size=K, replace=False)
        counts[i] = int(pred_mask[idx].sum())
    emp_p = float(np.mean(counts >= observed))
    return counts, emp_p


def run_all_nulls(
    *,
    db: TargetDB,
    scanner: TargetScanner,
    gold: dict[str, dict],
    n_decoy: int,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame, dict, dict]:
    """Return (per_mirna_df, pooled_df, decoy_hist, pooled_hist)."""
    N_pop = int(scanner.genes.shape[0])
    universe = np.array(list(scanner.symbols), dtype=object)

    # Storage: (mirna, mode, effect) -> decoy counts / labelperm counts / obs
    decoy_hist: dict[tuple[str, str, str], np.ndarray] = {}
    label_hist: dict[tuple[str, str, str], np.ndarray] = {}
    obs_map: dict[tuple[str, str, str], int] = {}
    kgold_map: dict[tuple[str, str, str], int] = {}
    npred_map: dict[tuple[str, str, str], int] = {}
    pred_sets: dict[tuple[str, str, str], set[str]] = {}

    per_rows = []

    for mirna, meta in gold.items():
        seed = clean_seq(meta["seed"])
        ox_label = meta["ox_label"]
        real_ox = parse_ox_positions(ox_label)
        n_ox = len(real_ox)
        conserved = conserved_for(db, mirna, seed)
        all_gold = meta["gained"] | meta["lost"]
        gold_utrs = gold_utr_table(scanner, all_gold)
        present = set(gold_utrs["symbol"])

        print(
            f"\n=== {mirna}  seed={seed}  {ox_label}  "
            f"gained={len(meta['gained'])} lost={len(meta['lost'])} ===",
            flush=True,
        )

        for mode in PrecisionMode:
            parts_full = observed_partition(db, mirna, seed, ox_label, mode)
            parts_gold = partition_from_states(gold_utrs, seed, real_ox, mode, conserved)

            for effect in ("gained", "lost"):
                gold_genes = meta[effect] & present
                K = len(gold_genes)
                obs = recovery_count(parts_gold, gold_genes, effect)
                n_draw = len(parts_full[effect])
                key = (mirna, mode.value, effect)
                obs_map[key] = obs
                kgold_map[key] = K
                npred_map[key] = n_draw
                pred_sets[key] = set(parts_full[effect])

                # --- decoy-seed null ---
                null_counts = np.empty(n_decoy, dtype=np.int32)
                for i in range(n_decoy):
                    dseed = composition_shuffle(seed, rng)
                    if dseed == seed:
                        dseed = composition_shuffle(seed, rng)
                    ox_pos = decoy_ox_positions(dseed, n_ox, rng)
                    parts = partition_from_states(
                        gold_utrs, dseed, ox_pos, mode, conserved
                    )
                    null_counts[i] = recovery_count(parts, gold_genes, effect)
                decoy_hist[key] = null_counts

                # --- label-permutation null ---
                lp_counts, lp_p = label_perm_null(
                    pred_sets[key], K, universe, n_decoy, rng, obs
                )
                label_hist[key] = lp_counts

                null_mean = float(null_counts.mean())
                null_sd = float(null_counts.std(ddof=1)) if n_decoy > 1 else 0.0
                decoy_p = float(np.mean(null_counts >= obs))
                fold = (
                    (obs / null_mean)
                    if null_mean > 0
                    else (float("inf") if obs > 0 else float("nan"))
                )
                z = (obs - null_mean) / null_sd if null_sd > 0 else float("nan")
                hg_p = hypergeom_p(obs, n_draw, K, N_pop)
                rate = obs / K if K else float("nan")
                ci_lo, ci_hi = wilson_ci(obs, K)
                agree = (decoy_p < 0.05) == (lp_p < 0.05)

                print(
                    f"  {mode.value:10s} {effect:6s}: obs={obs}/{K}  "
                    f"decoy_p={decoy_p:.4g}  labelperm_p={lp_p:.4g}  "
                    f"agree={agree}  Wilson=[{ci_lo:.2f},{ci_hi:.2f}]",
                    flush=True,
                )

                per_rows.append(
                    dict(
                        mirna=mirna,
                        mode=mode.value,
                        effect_type=effect,
                        n_gold=K,
                        observed=obs,
                        recovery_rate=rate,
                        wilson_ci_low=ci_lo,
                        wilson_ci_high=ci_hi,
                        n_pred_set=n_draw,
                        null_mean=null_mean,
                        null_sd=null_sd,
                        decoy_p=decoy_p,
                        labelperm_p=lp_p,
                        agreement_flag=agree,
                        empirical_p=decoy_p,  # back-compat alias
                        hypergeom_p=hg_p,
                        z_score=z,
                        fold_enrichment=fold,
                        significant_emp_p05=decoy_p < 0.05,
                        N_population=N_pop,
                        n_decoy=n_decoy,
                    )
                )

    per_df = pd.DataFrame(per_rows)

    # --- pooled across miRNAs ---
    pooled_rows = []
    pooled_hist: dict[tuple[str, str], np.ndarray] = {}
    mirnas = list(gold.keys())
    print("\n=== POOLED (headline) ===", flush=True)
    for mode in PrecisionMode:
        for effect in ("gained", "lost"):
            keys = [(m, mode.value, effect) for m in mirnas]
            obs = int(sum(obs_map[k] for k in keys))
            K = int(sum(kgold_map[k] for k in keys))
            # sum per-miRNA decoy recoveries within each permutation index
            decoy_pool = sum(decoy_hist[k] for k in keys)
            label_pool = sum(label_hist[k] for k in keys)
            pooled_hist[(mode.value, effect)] = decoy_pool

            null_mean = float(decoy_pool.mean())
            null_sd = float(decoy_pool.std(ddof=1)) if n_decoy > 1 else 0.0
            decoy_p = float(np.mean(decoy_pool >= obs))
            lp_p = float(np.mean(label_pool >= obs))
            agree = (decoy_p < 0.05) == (lp_p < 0.05)
            diffs = obs - decoy_pool
            ci_low, ci_high = float(np.percentile(diffs, 2.5)), float(
                np.percentile(diffs, 97.5)
            )
            obs_minus_null = obs - null_mean
            rate = obs / K if K else float("nan")
            w_lo, w_hi = wilson_ci(obs, K)

            print(
                f"  POOLED {mode.value:10s} {effect:6s}: obs={obs}/{K}  "
                f"null={null_mean:.2f}±{null_sd:.2f}  "
                f"Δ={obs_minus_null:.2f} [{ci_low:.2f},{ci_high:.2f}]  "
                f"decoy_p={decoy_p:.4g}  labelperm_p={lp_p:.4g}  agree={agree}",
                flush=True,
            )
            pooled_rows.append(
                dict(
                    mode=mode.value,
                    effect_type=effect,
                    n_gold=K,
                    observed=obs,
                    recovery_rate=rate,
                    wilson_ci_low=w_lo,
                    wilson_ci_high=w_hi,
                    null_mean=null_mean,
                    null_sd=null_sd,
                    obs_minus_null=obs_minus_null,
                    ci_low=ci_low,
                    ci_high=ci_high,
                    decoy_p=decoy_p,
                    labelperm_p=lp_p,
                    agreement_flag=agree,
                    empirical_p=decoy_p,
                    n_decoy=n_decoy,
                    N_population=N_pop,
                )
            )

    # sanity: pooled observed == sum of per-miRNA observed
    for mode in PrecisionMode:
        for effect in ("gained", "lost"):
            s = int(
                per_df[
                    (per_df["mode"] == mode.value) & (per_df["effect_type"] == effect)
                ]["observed"].sum()
            )
            p = int(
                next(
                    r["observed"]
                    for r in pooled_rows
                    if r["mode"] == mode.value and r["effect_type"] == effect
                )
            )
            assert s == p, (mode, effect, s, p)

    return per_df, pd.DataFrame(pooled_rows), decoy_hist, label_hist, pooled_hist


def update_gold_summary_with_wilson():
    """Rebuild recovery summary (+ Wilson CIs) from gold_master denominators."""
    master_path = BENCH / "gold_master.csv"
    path = BENCH / "gold_recovery_by_mode.csv"
    if master_path.exists():
        master = pd.read_csv(master_path)
        master = master[master["included"].astype(str).str.lower().isin(["true", "1"])]
        sens = master.copy()
        sens["effect"] = sens["effect_type"].map(
            {"gained": "gained_on_oxidation", "lost": "lost_on_oxidation"}
        )
        sens["state_label"] = sens["o8g_position"]
        sens.to_csv(BENCH / "gold_recovery_detail.csv", index=False)

    if path.exists():
        detail = pd.read_csv(path)
        usable = detail[detail["status"] != "skipped"].copy()
    elif master_path.exists():
        usable = master.copy()
        usable["mode"] = "Sensitive"
        usable["effect"] = usable["effect_type"].map(
            {"gained": "gained_on_oxidation", "lost": "lost_on_oxidation"}
        )
        if "source" not in usable.columns:
            usable["source"] = "ALL"
    else:
        return

    if "source" not in usable.columns:
        usable = usable.assign(source="ALL")
    rows = []
    for (mode, source, effect), g in usable.groupby(["mode", "source", "effect"]):
        n = len(g)
        n_ok = int(g["recovered"].fillna(False).astype(bool).sum())
        lo, hi = wilson_ci(n_ok, n)
        rows.append(
            dict(
                mode=mode,
                source=source,
                effect=effect,
                n_gold=n,
                n_recovered=n_ok,
                recovery_rate=n_ok / n if n else float("nan"),
                wilson_ci_low=lo,
                wilson_ci_high=hi,
            )
        )
    for (mode, effect), g in usable.groupby(["mode", "effect"]):
        n = len(g)
        n_ok = int(g["recovered"].fillna(False).astype(bool).sum())
        lo, hi = wilson_ci(n_ok, n)
        rows.append(
            dict(
                mode=mode,
                source="ALL",
                effect=effect,
                n_gold=n,
                n_recovered=n_ok,
                recovery_rate=n_ok / n if n else float("nan"),
                wilson_ci_low=lo,
                wilson_ci_high=hi,
            )
        )
    out = pd.DataFrame(rows).drop_duplicates().sort_values(["mode", "source", "effect"])
    out.to_csv(BENCH / "gold_recovery_summary_by_mode.csv", index=False)
    print(f"Wrote {BENCH / 'gold_recovery_summary_by_mode.csv'} (Wilson CIs)", flush=True)


def plot_null_per_mirna(df: pd.DataFrame, decoy_hist: dict, out_pdf: Path):
    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 150,
            "pdf.fonttype": 42,
        }
    )
    modes = ["Sensitive", "Stringent", "Consensus"]
    mirnas = sorted(df["mirna"].unique())
    fig, axes = plt.subplots(
        len(mirnas), len(modes), figsize=(10.5, 2.35 * len(mirnas)), squeeze=False
    )
    for i, mir in enumerate(mirnas):
        for j, mode in enumerate(modes):
            ax = axes[i, j]
            for effect, color in (("gained", "#C44E52"), ("lost", "#4C72B0")):
                counts = decoy_hist.get((mir, mode, effect), [])
                if counts is None or (hasattr(counts, "__len__") and len(counts) == 0):
                    continue
                arr = np.asarray(counts)
                mx = int(arr.max()) if len(arr) else 0
                bins = np.arange(0, mx + 2) - 0.5
                ax.hist(
                    arr,
                    bins=bins,
                    alpha=0.45,
                    color=color,
                    label=f"null {effect}",
                    density=True,
                )
                row = df[
                    (df.mirna == mir) & (df.mode == mode) & (df.effect_type == effect)
                ]
                if len(row):
                    ax.axvline(int(row.iloc[0]["observed"]), color=color, ls="--", lw=1.4)
            ax.set_title(f"{mir.replace('hsa-', '')} · {mode}", fontsize=8)
            if i == len(mirnas) - 1:
                ax.set_xlabel("Gold genes recovered")
            if j == 0:
                ax.set_ylabel("Density")
            if i == 0 and j == len(modes) - 1:
                ax.legend(frameon=False, fontsize=6)
    fig.suptitle(
        "Fig. S  Per-miRNA decoy-seed null (composition shuffle)",
        fontsize=10,
        fontweight="bold",
        y=1.01,
        x=0.01,
        ha="left",
    )
    fig.tight_layout()
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_pdf.with_suffix(".png"), bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_pdf}", flush=True)


def plot_pooled_null(pooled: pd.DataFrame, pooled_hist: dict, out_pdf: Path):
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 150,
            "pdf.fonttype": 42,
        }
    )
    modes = ["Sensitive", "Stringent", "Consensus"]
    fig, axes = plt.subplots(2, 3, figsize=(11.0, 5.8), squeeze=False)
    for j, mode in enumerate(modes):
        for i, effect in enumerate(("gained", "lost")):
            ax = axes[i, j]
            counts = pooled_hist.get((mode, effect), np.array([]))
            row = pooled[
                (pooled["mode"] == mode) & (pooled["effect_type"] == effect)
            ]
            if not len(row):
                continue
            r = row.iloc[0]
            obs = int(r["observed"])
            color = "#C44E52" if effect == "gained" else "#4C72B0"
            if len(counts):
                mx = int(np.max(counts))
                bins = np.arange(0, max(mx, obs) + 2) - 0.5
                ax.hist(counts, bins=bins, color=color, alpha=0.55, density=True)
            ax.axvline(obs, color="k", ls="--", lw=1.6, label=f"obs={obs}")
            ax.axvline(
                float(r["null_mean"]),
                color="#666",
                ls=":",
                lw=1.2,
                label=f"null={r['null_mean']:.1f}",
            )
            sig = "*" if r["decoy_p"] < 0.05 else "ns"
            agree = "agree" if r["agreement_flag"] else "diverge"
            ax.set_title(
                f"{mode} · {effect}\n"
                f"decoy_p={r['decoy_p']:.3g}  label_p={r['labelperm_p']:.3g} ({agree}) {sig}",
                fontsize=8,
            )
            if i == 1:
                ax.set_xlabel("Pooled gold recovered")
            if j == 0:
                ax.set_ylabel("Density")
            ax.legend(frameon=False, fontsize=6)
    fig.suptitle(
        "Fig. S  Pooled null (sum of per-miRNA decoy recoveries; headline)",
        fontsize=10,
        fontweight="bold",
        y=1.01,
        x=0.01,
        ha="left",
    )
    fig.tight_layout()
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_pdf.with_suffix(".png"), bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_pdf}", flush=True)


def real_seed_gold_recovery(
    scanner: TargetScanner,
    db: TargetDB,
    mirna: str,
    mode: PrecisionMode,
) -> tuple[int, int]:
    meta = load_gold_by_mirna()[mirna]
    seed = clean_seq(meta["seed"])
    ox_pos = parse_ox_positions(meta["ox_label"])
    conserved = conserved_for(db, mirna, seed)
    gold_utrs = gold_utr_table(scanner, meta["gained"] | meta["lost"])
    parts = partition_from_states(gold_utrs, seed, ox_pos, mode, conserved)
    present = set(gold_utrs["symbol"])
    g = recovery_count(parts, meta["gained"] & present, "gained")
    l = recovery_count(parts, meta["lost"] & present, "lost")
    return g, l


def main():
    n_decoy = int(os.environ.get("O8G_NULL_N", "10000"))
    rng = np.random.default_rng(int(os.environ.get("O8G_NULL_SEED", "42")))
    BENCH.mkdir(parents=True, exist_ok=True)

    print(f"Loading scanner from {UTR} …", flush=True)
    scanner = TargetScanner.from_parquet(str(UTR))
    print(f"UTR genes N={scanner.genes.shape[0]}", flush=True)
    db = TargetDB(str(ROOT / "o8g_targets.db"))
    gold = load_gold_by_mirna()

    per_df, pooled_df, decoy_hist, label_hist, pooled_hist = run_all_nulls(
        db=db, scanner=scanner, gold=gold, n_decoy=n_decoy, rng=rng
    )

    per_path = BENCH / "nullmodel_recovery.csv"
    pool_path = BENCH / "nullmodel_pooled.csv"
    per_df.to_csv(per_path, index=False)
    pooled_df.to_csv(pool_path, index=False)
    print(f"\nWrote {per_path}", flush=True)
    print(f"Wrote {pool_path}", flush=True)

    # Persist all permutation draws (NAR: no illustrative histograms)
    draws_dir = BENCH / "null_draws"
    draws_dir.mkdir(parents=True, exist_ok=True)
    decoy_rows = []
    label_rows = []
    for (mir, mode, effect), arr in decoy_hist.items():
        for i, v in enumerate(arr):
            decoy_rows.append(
                dict(mirna=mir, mode=mode, effect_type=effect, draw_i=i, count=int(v))
            )
    for (mir, mode, effect), arr in label_hist.items():
        for i, v in enumerate(arr):
            label_rows.append(
                dict(mirna=mir, mode=mode, effect_type=effect, draw_i=i, count=int(v))
            )
    decoy_df = pd.DataFrame(decoy_rows)
    label_df = pd.DataFrame(label_rows)
    decoy_pq = draws_dir / "null_draws_decoy_per_mirna.parquet"
    label_pq = draws_dir / "null_draws_labelperm_per_mirna.parquet"
    decoy_df.to_parquet(decoy_pq, index=False)
    label_df.to_parquet(label_pq, index=False)
    print(f"Wrote {decoy_pq} ({len(decoy_df):,} rows)", flush=True)
    print(f"Wrote {label_pq} ({len(label_df):,} rows)", flush=True)

    # Pooled draws (sum across miRNAs within each draw_i)
    pooled_draw_rows = []
    for mode in PrecisionMode:
        for effect in ("gained", "lost"):
            keys = [(m, mode.value, effect) for m in gold.keys()]
            decoy_pool = sum(decoy_hist[k] for k in keys)
            label_pool = sum(label_hist[k] for k in keys)
            for i, (d, l) in enumerate(zip(decoy_pool, label_pool)):
                pooled_draw_rows.append(
                    dict(
                        mode=mode.value,
                        effect_type=effect,
                        draw_i=i,
                        decoy_count=int(d),
                        labelperm_count=int(l),
                    )
                )
    pooled_draws = pd.DataFrame(pooled_draw_rows)
    pooled_pq = draws_dir / "null_draws_pooled.parquet"
    pooled_draws.to_parquet(pooled_pq, index=False)
    print(f"Wrote {pooled_pq} ({len(pooled_draws):,} rows)", flush=True)

    # Sanity: means from parquet match pooled CSV
    for _, r in pooled_df.iterrows():
        sub = pooled_draws[
            (pooled_draws["mode"] == r["mode"])
            & (pooled_draws["effect_type"] == r["effect_type"])
        ]
        m = float(sub["decoy_count"].mean())
        assert abs(m - float(r["null_mean"])) < 1e-9, (r["mode"], r["effect_type"], m, r["null_mean"])
    print("Persisted-draw means match nullmodel_pooled.csv", flush=True)

    update_gold_summary_with_wilson()
    plot_null_per_mirna(per_df, decoy_hist, FIGS / "figS_nullmodel.pdf")
    plot_pooled_null(pooled_df, pooled_hist, FIGS / "figS_pooled_null.pdf")

    # headline printout
    print("\n=== HEADLINE POOLED VERDICTS ===", flush=True)
    for _, r in pooled_df.iterrows():
        sig = "SIGNIFICANT" if r["decoy_p"] < 0.05 else "ns"
        print(
            f"POOLED {r['mode']} {r['effect_type']}: {sig}  "
            f"obs={int(r['observed'])}/{int(r['n_gold'])}  "
            f"Δ={r['obs_minus_null']:.1f} [{r['ci_low']:.1f},{r['ci_high']:.1f}]  "
            f"decoy_p={r['decoy_p']:.4g}  labelperm_p={r['labelperm_p']:.4g}  "
            f"agree={r['agreement_flag']}",
            flush=True,
        )


if __name__ == "__main__":
    main()
