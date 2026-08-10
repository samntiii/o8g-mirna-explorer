"""Statistics view: length-matched ORA or preranked GSEA (opt-in; slow)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import streamlit as st

from o8g_enrich import enrich, enrich_within_pool
from o8g_sections import SectionContext


def render(ctx: SectionContext) -> None:
    st.markdown(f"#### Statistics · `{ctx.precision_mode}`")
    st.caption(
        "Differential enrichment with an explicit within-pool control. "
        "GSEA is optional and requires `gseapy`."
    )
    s_un = st.selectbox(
        "State A (baseline)",
        ctx.state_labels,
        index=ctx.state_labels.index("none") if "none" in ctx.state_labels else 0,
        key="stats_a",
    )
    ox_opts = [s for s in ctx.state_labels if s != s_un] or ctx.state_labels
    s_ox = st.selectbox(
        "State B",
        ox_opts,
        index=ox_opts.index("o8G@7") if "o8G@7" in ox_opts else 0,
        key="stats_b",
    )
    a = ctx.strong_set(s_un)
    b = ctx.strong_set(s_ox)
    lost, gained, pool = a - b, b - a, a | b
    m1, m2, m3 = st.columns(3)
    m1.metric("Lost", len(lost))
    m2.metric("Gained", len(gained))
    m3.metric("Pool |A∪B|", len(pool))

    run = st.checkbox("Run (resampling / enrichment can be slow)", value=False, key="stats_run")
    if not run:
        st.info("Check **Run** to compute ORA (and optional GSEA).")
        return

    method = st.radio("Method", ["ORA (hypergeometric)", "GSEA (preranked)"], horizontal=True)
    if method.startswith("ORA"):
        eg_lost = enrich(list(lost), ctx.library, background=ctx.universe, top=25)
        ep_lost = enrich_within_pool(list(lost), pool, library=ctx.library, top=25)
        eg_gain = enrich(list(gained), ctx.library, background=ctx.universe, top=25)
        ep_gain = enrich_within_pool(list(gained), pool, library=ctx.library, top=25)
        st.info(
            "Left = vs 3′UTR universe; right = within target pool. "
            "A genome-only number must never stand alone for differential claims."
        )
        for title, eg, ep in (
            ("Lost", eg_lost, ep_lost),
            ("Gained", eg_gain, ep_gain),
        ):
            st.markdown(f"##### {title}")
            c1, c2 = st.columns(2)
            with c1:
                st.caption(f"vs UTR universe · sig q<0.05: {int((eg['q_value']<0.05).sum()) if len(eg) else 0}")
                st.dataframe(eg[["term", "overlap", "odds_ratio", "q_value"]] if len(eg) else eg,
                             hide_index=True, height=220)
            with c2:
                st.caption(f"within pool · sig q<0.05: {int((ep['q_value']<0.05).sum()) if len(ep) else 0}")
                st.dataframe(ep[["term", "overlap", "odds_ratio", "q_value"]] if len(ep) else ep,
                             hide_index=True, height=220)
    else:
        try:
            import gseapy as gp  # noqa: F401
        except ImportError:
            st.error("`gseapy` is not installed. `pip install gseapy>=1.1` then rerun.")
            return
        # Simple rank: +1 gained, -1 lost, 0 else — illustrative prerank
        ranks = {}
        for g in pool:
            if g in gained:
                ranks[g] = 1.0
            elif g in lost:
                ranks[g] = -1.0
            else:
                ranks[g] = 0.0
        rnk = pd.Series(ranks).sort_values(ascending=False)
        st.caption(f"Preranked {len(rnk)} genes in pool (gained=+1, lost=-1).")
        try:
            from o8g_enrich import LIBRARY_FILES, GENESET_DIR
            import os

            gmt = os.path.join(GENESET_DIR, LIBRARY_FILES.get(ctx.library, ""))
            if not os.path.exists(gmt):
                st.error(f"GMT for {ctx.library} not found at {gmt}")
                return
            res = gp.prerank(
                rnk=rnk,
                gene_sets=gmt,
                threads=1,
                min_size=5,
                max_size=500,
                permutation_num=100,
                outdir=None,
                seed=42,
                verbose=False,
            )
            out = res.res2d if hasattr(res, "res2d") else pd.DataFrame(res.results).T
            st.dataframe(out.head(30), hide_index=True)
        except Exception as e:
            st.error(f"GSEA failed ({type(e).__name__}: {e})")
