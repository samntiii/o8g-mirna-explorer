"""
o8g_binding.py
==============
Binding-efficiency (BE) score used to rank predicted targets in the explorer.

Formula (transparent, not a black-box ML score)
-----------------------------------------------
For each gene's best seed match under the current oxidation state:

  BE = W(site_type)
     + λ_mult · log2(1 + n_weighted)
     + λ_ctx  · σ(context_score)
     + λ_cons · 1[is_conserved]

where
  W(8mer)=1.00, W(7mer-m8)=0.70, W(7mer-A1)=0.40, W(6mer)=0.15
      (Grimson et al., Mol Cell 2007 site-type hierarchy; matches SITE_WEIGHT)
  n_weighted = n_8mer·1.0 + n_7mer_m8·0.7 + n_7mer_A1·0.4 + n_6mer·0.15
  σ(x) = 1 / (1 + exp(−x)) maps context_score onto (0,1)
  λ_mult=0.35, λ_ctx=0.25, λ_cons=0.15

Context features (when live-scanned) follow the lightweight context++ analog in
``o8g_scanner.py`` (AU flanks, 3′ supplementary pairing, end proximity, UTR
position — signs from Garcia NSMB 2011 / Agarwal eLife 2015).

This is an **engine-internal ranking metric**, not an official TargetScan
context++ percentile. Cite the coefficients above when reporting.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

# Site-type weights (identical spirit to o8g_scanner.SITE_WEIGHT)
W_SITE = {
    "8mer": 1.00,
    "7mer-m8": 0.70,
    "7mer-A1": 0.40,
    "6mer": 0.15,
}
W_BY_RANK = {4: 1.00, 3: 0.70, 2: 0.40, 1: 0.15}

LAMBDA_MULT = 0.35
LAMBDA_CTX = 0.25
LAMBDA_CONS = 0.15


def _sigmoid(x: np.ndarray | float) -> np.ndarray | float:
    x = np.asarray(x, dtype=float)
    return 1.0 / (1.0 + np.exp(-np.clip(x, -20, 20)))


def binding_efficiency_series(df: pd.DataFrame) -> pd.Series:
    """Compute BE for each row of a targets DataFrame."""
    n = len(df)
    if n == 0:
        return pd.Series(dtype=float)

    if "site_type" in df.columns:
        w = df["site_type"].map(W_SITE).fillna(0.15).astype(float).to_numpy()
    elif "site_rank" in df.columns:
        w = df["site_rank"].map(W_BY_RANK).fillna(0.15).astype(float).to_numpy()
    else:
        w = np.full(n, 0.15)

    n8 = df["n_8mer"].to_numpy(dtype=float) if "n_8mer" in df.columns else np.zeros(n)
    n7 = df["n_7mer_m8"].to_numpy(dtype=float) if "n_7mer_m8" in df.columns else np.zeros(n)
    n7a = df["n_7mer_A1"].to_numpy(dtype=float) if "n_7mer_A1" in df.columns else np.zeros(n)
    n6 = df["n_6mer"].to_numpy(dtype=float) if "n_6mer" in df.columns else np.zeros(n)
    n_w = n8 * 1.0 + n7 * 0.7 + n7a * 0.4 + n6 * 0.15
    # fall back to multiplicity score column if site counts missing
    if "score" in df.columns and n_w.sum() == 0:
        n_w = df["score"].to_numpy(dtype=float)

    mult = LAMBDA_MULT * np.log2(1.0 + np.maximum(n_w, 0.0))

    if "context_score" in df.columns:
        ctx = LAMBDA_CTX * _sigmoid(df["context_score"].fillna(0.0).to_numpy(dtype=float))
    else:
        ctx = np.zeros(n)

    if "is_conserved" in df.columns:
        cons = LAMBDA_CONS * df["is_conserved"].fillna(False).astype(float).to_numpy()
    else:
        cons = np.zeros(n)

    return pd.Series(w + mult + ctx + cons, index=df.index, name="binding_efficiency")


def add_binding_efficiency(
    df: pd.DataFrame,
    *,
    sort: bool = True,
    scanner=None,
    mature_dna: str | None = None,
    oxidized_positions: tuple[int, ...] | list[int] = (),
    mirna: str | None = None,
    is_unmodified: bool = True,
    with_thermo: bool = True,
) -> pd.DataFrame:
    """Return a copy with ``binding_efficiency`` (+ optional Vienna/TS metrics).

    When ``with_thermo`` and a live scanner + mature sequence are provided, also
    attaches ``dG_RNAduplex``, ``dG_RNAup``, and ``contextpp_TargetScan``
    (see ``o8g_thermo``). Default sort remains by binding_efficiency.
    """
    out = df.copy()
    if out.empty:
        out["binding_efficiency"] = pd.Series(dtype=float)
        return out
    out["binding_efficiency"] = binding_efficiency_series(out).to_numpy()

    if with_thermo and scanner is not None and mature_dna:
        try:
            from o8g_thermo import score_targets_thermo

            out = score_targets_thermo(
                out,
                scanner=scanner,
                mature_dna=mature_dna,
                oxidized_positions=oxidized_positions,
                mirna=mirna,
                is_unmodified=is_unmodified,
            )
        except Exception:
            pass

    if sort:
        cols = ["binding_efficiency"]
        asc = [False]
        if "site_rank" in out.columns:
            cols.append("site_rank")
            asc.append(False)
        if "symbol" in out.columns:
            cols.append("symbol")
            asc.append(True)
        out = out.sort_values(cols, ascending=asc).reset_index(drop=True)
    return out


FORMULA_CAPTION = (
    "Ranking uses **binding efficiency** "
    "BE = W(site) + 0.35·log₂(1+n_w) + 0.25·σ(context) + 0.15·conserved "
    "(Grimson 2007 site types; Garcia/Agarwal context++ signs). "
    "Also reported: **dG_RNAduplex** (ViennaRNA), **dG_RNAup** (duplex+opening), "
    "**contextpp_TargetScan** (official TS8 weighted context++, unmodified only)."
)
