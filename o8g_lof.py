"""Loss-of-function view: miRDB-anchored baseline attrition (optional file)."""
from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
import streamlit as st

from o8g_enrich import enrich, enrich_within_pool
from o8g_sections import SectionContext

ROOT = Path(__file__).resolve().parent


def enrich_within_baseline(query_genes, baseline_genes, library: str = "GO_BP", **kw):
    """Alias kept for callers expecting the LoF-shaped API."""
    return enrich_within_pool(query_genes, baseline_genes, library=library, **kw)


def _load_mirdb(mirna: str, score_min: float) -> set[str] | None:
    """Return miRDB symbols for mirna at score >= cutoff, or None if file missing."""
    candidates = [
        ROOT / "mirdb_ref.parquet",
        ROOT / "paper" / "data" / "mirdb_ref.parquet",
        ROOT / "mirdb_custom_cache.db",
    ]
    present = [p for p in candidates if p.exists()]
    if not present:
        return None
    # Prefer parquet extract if present
    for p in present:
        if p.suffix == ".parquet":
            try:
                df = pd.read_parquet(p)
                cols = {c.lower(): c for c in df.columns}
                mir_col = cols.get("mirna") or cols.get("mir") or list(df.columns)[0]
                score_col = cols.get("score") or cols.get("mirdb_score")
                sym_col = cols.get("symbol") or cols.get("gene") or cols.get("gene_symbol")
                sub = df[df[mir_col].astype(str) == mirna]
                if score_col is not None:
                    sub = sub[sub[score_col] >= score_min]
                return set(sub[sym_col].astype(str).str.upper())
            except Exception:
                continue
    return set()  # file present but unreadable shape → empty with warning upstream


def render(ctx: SectionContext) -> None:
    st.markdown(f"#### Loss of function (miRDB-anchored attrition) · `{ctx.precision_mode}`")
    st.warning(
        "miRDB catalogs **unmodified** matures only. The anchored oxidized set is "
        "`baseline ∩ oxidized_strong`, so **gained is identically 0 by construction** — "
        "structural, not measured. Use this view for loss/attrition only; do not treat "
        "gained≡0 as a biological finding."
    )
    score = st.slider("miRDB score cutoff", min_value=50, max_value=100, value=80, step=5)
    baseline = _load_mirdb(ctx.mirna, float(score))
    if baseline is None:
        st.error(
            "miRDB reference missing. Expected `mirdb_ref.parquet` (or "
            "`mirdb_custom_cache.db`) next to the app. Attrition plot unavailable."
        )
        # Still show Explorer-only loss fractions so the section does not crash.
        baseline = ctx.strong_set("none")
        st.caption("Falling back to Explorer unmodified strong set as a temporary baseline.")
        anchored = True
    else:
        anchored = True
        if len(baseline) == 0:
            st.warning("miRDB file present but returned 0 genes for this miRNA/cutoff.")

    ox_labels = [s for s in ctx.state_labels if s != "none"]
    rows = []
    for lab in ox_labels:
        ox = ctx.strong_set(lab)
        anchored_ox = baseline & ox
        lost = baseline - ox
        rows.append(
            {
                "state": lab,
                "baseline": len(baseline),
                "anchored_oxidized": len(anchored_ox),
                "lost": len(lost),
                "retained_frac": (len(anchored_ox) / len(baseline)) if baseline else 0.0,
                "gained_structural_zero": 0,
            }
        )
    tab = pd.DataFrame(rows)
    st.dataframe(tab, width="stretch", hide_index=True)
    if len(tab):
        st.bar_chart(tab.set_index("state")["lost"])
    st.caption(
        f"Vulnerability summary for {ctx.mirna}: median lost = "
        f"{tab['lost'].median() if len(tab) else 0:.0f} of {len(baseline)} baseline genes."
    )

    # Dual enrichment on none→first ox state
    if ox_labels and len(baseline) >= 5:
        lab = "o8G@7" if "o8G@7" in ox_labels else ox_labels[0]
        lost = baseline - ctx.strong_set(lab)
        pool = baseline  # within-baseline control
        eg = enrich(list(lost), ctx.library, background=ctx.universe, top=15)
        ep = enrich_within_baseline(list(lost), pool, library=ctx.library, top=15)
        c1, c2 = st.columns(2)
        with c1:
            st.markdown(f"**Lost vs UTR universe** ({lab})")
            st.metric("q<0.05 terms", int((eg["q_value"] < 0.05).sum()) if len(eg) else 0)
            st.dataframe(eg[["term", "overlap", "q_value"]] if len(eg) else eg, hide_index=True, height=200)
        with c2:
            st.markdown(f"**Lost within miRDB baseline** ({lab})")
            st.metric("q<0.05 terms", int((ep["q_value"] < 0.05).sum()) if len(ep) else 0)
            st.dataframe(ep[["term", "overlap", "q_value"]] if len(ep) else ep, hide_index=True, height=200)
