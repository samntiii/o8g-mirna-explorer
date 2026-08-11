"""Score oxidation states against user DEG UP/DOWN sets (lost/gained concordance)."""
from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from o8g_engine import enumerate_states, g_positions
from o8g_precision import PrecisionConfig, PrecisionMode


def _jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0


def _fisher_overlap(query: set[str], foreground: set[str], universe: set[str]) -> tuple[int, float, float]:
    """Return (k, odds_ratio, p_value) for enrichment of query in foreground within universe."""
    U = {g.upper() for g in universe}
    Q = {g.upper() for g in query} & U
    F = {g.upper() for g in foreground} & U
    k = len(Q & F)
    # 2x2: in F / not in F  ×  in Q / not in Q  — test enrichment of F among Q
    # Standard: rows = in query vs not; cols = in foreground vs not
    a = k
    b = len(Q) - k
    c = len(F) - k
    d = len(U) - len(Q) - len(F) + k
    if min(len(Q), len(F), len(U)) == 0:
        return k, float("nan"), 1.0
    # guard negative cells from set inconsistencies
    a, b, c, d = max(a, 0), max(b, 0), max(c, 0), max(d, 0)
    oddsratio, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    return k, float(oddsratio), float(p)


def _bh(pvals: list[float]) -> list[float]:
    m = len(pvals)
    if m == 0:
        return []
    order = np.argsort(pvals)
    ranked = np.asarray(pvals, dtype=float)[order]
    q = ranked * m / (np.arange(1, m + 1))
    q = np.minimum.accumulate(q[::-1])[::-1].clip(0, 1)
    out = np.empty(m, dtype=float)
    out[order] = q
    return out.tolist()


def list_ox_labels(seed: str, *, max_states: int = 32, single_g_only: bool = False) -> list[str]:
    gpos = g_positions(seed)
    n = 2 ** len(gpos)
    if n > max_states and not single_g_only:
        # prefer single-G when combinatorial explosion
        single_g_only = True
    states = enumerate_states(seed)
    labels = []
    for s in states:
        if s.label == "none":
            continue
        if single_g_only and len(s.oxidized_positions) != 1:
            continue
        labels.append(s.label)
    return labels


def score_mirna_states(
    db,
    *,
    mirna: str,
    up: set[str],
    down: set[str],
    universe: set[str],
    precision_cfg: PrecisionConfig | PrecisionMode | str,
    single_g_only: bool = False,
    max_states: int = 32,
    down_weight: float = 0.0,
) -> pd.DataFrame:
    """Rank oxidation states for one miRNA by UP∩lost (and optional DOWN∩gained) concordance."""
    info = db.mirna_info(mirna)
    seed = info["seed"]
    cfg = PrecisionConfig.from_mode(precision_cfg)
    labels = list_ox_labels(seed, max_states=max_states, single_g_only=single_g_only)
    rows = []
    p_up = []
    p_down = []
    for lab in labels:
        parts = db.retarget_partition(seed, lab, cfg, mirna=mirna)
        lost = {g.upper() for g in parts["lost"]}
        gained = {g.upper() for g in parts["gained"]}
        k_up, or_up, p_u = _fisher_overlap(up, lost, universe)
        k_dn, or_dn, p_d = _fisher_overlap(down, gained, universe)
        jac_up = _jaccard(up, lost)
        jac_dn = _jaccard(down, gained)
        rows.append(
            dict(
                mirna=mirna,
                seed=seed,
                ox_label=lab,
                n_unmod=len(parts["unmod"]),
                n_oxid=len(parts["oxid"]),
                n_lost=len(lost),
                n_gained=len(gained),
                n_up=len(up),
                n_down=len(down),
                n_up_lost=k_up,
                n_down_gained=k_dn,
                odds_up_lost=or_up,
                p_up_lost=p_u,
                jaccard_up_lost=jac_up,
                odds_down_gained=or_dn,
                p_down_gained=p_d,
                jaccard_down_gained=jac_dn,
                up_lost_genes=";".join(sorted(up & lost)[:40]),
                down_gained_genes=";".join(sorted(down & gained)[:40]),
            )
        )
        p_up.append(p_u)
        p_down.append(p_d)

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    df["q_up_lost"] = _bh(p_up)
    df["q_down_gained"] = _bh(p_down)
    # concordance: primary = UP~lost; optional DOWN~gained
    neglog_q = -np.log10(df["q_up_lost"].clip(lower=1e-300))
    score = neglog_q + df["jaccard_up_lost"]
    if down_weight:
        score = score + float(down_weight) * (
            -np.log10(df["q_down_gained"].clip(lower=1e-300)) + df["jaccard_down_gained"]
        )
    df["concordance_score"] = score
    return df.sort_values("concordance_score", ascending=False).reset_index(drop=True)


def score_panel(
    db,
    mirnas: Iterable[str],
    *,
    up: set[str],
    down: set[str],
    universe: set[str],
    precision_cfg: PrecisionConfig | PrecisionMode | str,
    single_g_only: bool = False,
    max_states: int = 32,
    down_weight: float = 0.0,
) -> pd.DataFrame:
    frames = []
    for mir in mirnas:
        try:
            frames.append(
                score_mirna_states(
                    db,
                    mirna=mir,
                    up=up,
                    down=down,
                    universe=universe,
                    precision_cfg=precision_cfg,
                    single_g_only=single_g_only,
                    max_states=max_states,
                    down_weight=down_weight,
                )
            )
        except Exception as e:
            frames.append(
                pd.DataFrame(
                    [
                        dict(
                            mirna=mir,
                            seed="",
                            ox_label="ERROR",
                            concordance_score=np.nan,
                            up_lost_genes=str(e)[:120],
                        )
                    ]
                )
            )
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    # Re-BH across the full panel for primary p
    if "p_up_lost" in out.columns and out["p_up_lost"].notna().any():
        mask = out["p_up_lost"].notna()
        qs = _bh(out.loc[mask, "p_up_lost"].tolist())
        out.loc[mask, "q_up_lost"] = qs
        neglog_q = -np.log10(out["q_up_lost"].clip(lower=1e-300))
        out["concordance_score"] = neglog_q + out["jaccard_up_lost"].fillna(0)
        if down_weight and "q_down_gained" in out.columns:
            out["concordance_score"] = out["concordance_score"] + float(down_weight) * (
                -np.log10(out["q_down_gained"].clip(lower=1e-300))
                + out["jaccard_down_gained"].fillna(0)
            )
    return out.sort_values("concordance_score", ascending=False).reset_index(drop=True)
