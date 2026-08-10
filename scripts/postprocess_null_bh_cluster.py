#!/usr/bin/env python3
"""
Post-process null draws: BH FDR, clustered bootstrap CIs, update pooled CSV.
Requires paper/benchmarks/null_draws/null_draws_*.parquet from benchmark_nullmodel.py.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "paper" / "benchmarks"
DRAWS = BENCH / "null_draws"


def bh_adjust(pvals: list[float]) -> list[float]:
    """Benjamini–Hochberg FDR (independent/positive regression)."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adj = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked[i] * n / rank
        prev = min(prev, val)
        adj[order[i]] = min(prev, 1.0)
    return adj.tolist()


def cluster_bootstrap_rate(
    per_mirna: pd.DataFrame,
    mode: str,
    effect: str,
    n_boot: int = 10000,
    seed: int = 42,
) -> tuple[float, float, float]:
    """Cluster bootstrap by miRNA for pooled recovery rate.

    Resample miRNAs with replacement; sum observed and n_gold within bootstrap.
    """
    sub = per_mirna[
        (per_mirna["mode"] == mode) & (per_mirna["effect_type"] == effect)
    ].copy()
    if sub.empty:
        return float("nan"), float("nan"), float("nan")
    obs = sub["observed"].to_numpy(dtype=float)
    kg = sub["n_gold"].to_numpy(dtype=float)
    rate_point = float(obs.sum() / kg.sum()) if kg.sum() else float("nan")
    rng = np.random.default_rng(seed)
    n = len(sub)
    rates = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        idx = rng.integers(0, n, size=n)
        o = obs[idx].sum()
        k = kg[idx].sum()
        rates[i] = o / k if k > 0 else float("nan")
    lo, hi = float(np.nanpercentile(rates, 2.5)), float(np.nanpercentile(rates, 97.5))
    return rate_point, lo, hi


def main():
    pooled = pd.read_csv(BENCH / "nullmodel_pooled.csv")
    per = pd.read_csv(BENCH / "nullmodel_recovery.csv")
    pooled_draws = pd.read_parquet(DRAWS / "null_draws_pooled.parquet")

    # Verify means
    for _, r in pooled.iterrows():
        sub = pooled_draws[
            (pooled_draws["mode"] == r["mode"])
            & (pooled_draws["effect_type"] == r["effect_type"])
        ]
        m = float(sub["decoy_count"].mean())
        if abs(m - float(r["null_mean"])) > 1e-6:
            raise SystemExit(
                f"MEAN MISMATCH {r['mode']} {r['effect_type']}: "
                f"draws={m} csv={r['null_mean']}"
            )

    # BH across 6 pooled cells × 2 nulls = 12 tests
    keys = []
    raw = []
    for _, r in pooled.iterrows():
        keys.append((r["mode"], r["effect_type"], "decoy"))
        raw.append(float(r["decoy_p"]))
        keys.append((r["mode"], r["effect_type"], "labelperm"))
        raw.append(float(r["labelperm_p"]))
    adj = bh_adjust(raw)
    bh_map = {k: a for k, a in zip(keys, adj)}

    # Per-miRNA BH (all decoy_p and labelperm_p)
    per_keys = []
    per_raw = []
    for _, r in per.iterrows():
        per_keys.append((r["mirna"], r["mode"], r["effect_type"], "decoy"))
        per_raw.append(float(r["decoy_p"]))
        per_keys.append((r["mirna"], r["mode"], r["effect_type"], "labelperm"))
        per_raw.append(float(r["labelperm_p"]))
    per_adj = bh_adjust(per_raw)
    per_bh = {k: a for k, a in zip(per_keys, per_adj)}

    rows = []
    for _, r in pooled.iterrows():
        d = r.to_dict()
        d["decoy_p_bh"] = bh_map[(r["mode"], r["effect_type"], "decoy")]
        d["labelperm_p_bh"] = bh_map[(r["mode"], r["effect_type"], "labelperm")]
        rate, clo, chi = cluster_bootstrap_rate(per, r["mode"], r["effect_type"])
        d["cluster_boot_rate"] = rate
        d["cluster_boot_ci_low"] = clo
        d["cluster_boot_ci_high"] = chi
        # empirical p from draws (should match decoy_p)
        sub = pooled_draws[
            (pooled_draws["mode"] == r["mode"])
            & (pooled_draws["effect_type"] == r["effect_type"])
        ]
        obs = int(r["observed"])
        d["decoy_p_from_draws"] = float(np.mean(sub["decoy_count"].to_numpy() >= obs))
        d["labelperm_p_from_draws"] = float(
            np.mean(sub["labelperm_count"].to_numpy() >= obs)
        )
        rows.append(d)

    out = pd.DataFrame(rows)
    out_path = BENCH / "nullmodel_pooled.csv"
    out.to_csv(out_path, index=False)
    print(f"Wrote {out_path} with BH + cluster bootstrap columns")

    # per-miRNA BH columns
    per2 = per.copy()
    per2["decoy_p_bh"] = [
        per_bh[(r.mirna, r.mode, r.effect_type, "decoy")] for r in per2.itertuples()
    ]
    per2["labelperm_p_bh"] = [
        per_bh[(r.mirna, r.mode, r.effect_type, "labelperm")] for r in per2.itertuples()
    ]
    per2.to_csv(BENCH / "nullmodel_recovery.csv", index=False)
    print(f"Updated nullmodel_recovery.csv with BH columns")

    print("\nPooled BH summary:")
    for _, r in out.iterrows():
        print(
            f"  {r['mode']:10s} {r['effect_type']:6s}  "
            f"decoy={r['decoy_p']:.4g} (BH={r['decoy_p_bh']:.4g})  "
            f"label={r['labelperm_p']:.4g} (BH={r['labelperm_p_bh']:.4g})  "
            f"cluster_rate={r['cluster_boot_rate']:.3f} "
            f"[{r['cluster_boot_ci_low']:.3f},{r['cluster_boot_ci_high']:.3f}]"
        )


if __name__ == "__main__":
    main()
