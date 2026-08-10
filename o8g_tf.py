"""Transcription-factor targets among gained / lost / retained sets."""
from __future__ import annotations

import os

import pandas as pd
import streamlit as st

from o8g_enrich import TF_LIBRARIES, LIBRARY_FILES, GENESET_DIR, enrich, enrich_within_pool
from o8g_sections import SectionContext


def _tf_symbols_from_libraries() -> set[str]:
    """Union of TF names appearing as term prefixes / gene members in TF GMTs."""
    out: set[str] = set()
    for lib in TF_LIBRARIES:
        path = os.path.join(GENESET_DIR, LIBRARY_FILES[lib])
        if not os.path.exists(path):
            continue
        with open(path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                # ChEA/ENCODE term often starts with TF name
                head = parts[0].split("_")[0].upper()
                if head.isalpha() and 2 <= len(head) <= 10:
                    out.add(head)
                for g in parts[2:]:
                    g = g.split(",")[0].strip().upper()
                    if g:
                        out.add(g)
    return out


def render(ctx: SectionContext) -> None:
    st.markdown(f"#### Transcription factors as direct targets · `{ctx.precision_mode}`")
    missing = [
        LIBRARY_FILES[k]
        for k in TF_LIBRARIES
        if not os.path.exists(os.path.join(GENESET_DIR, LIBRARY_FILES[k]))
    ]
    if missing:
        st.warning(
            "TF gene-set libraries not installed (expected under `genesets/`): "
            + ", ".join(missing)
            + ". Direct-target TF table still runs on a heuristic name list when any GMT is present; "
            "otherwise only the state partition is shown."
        )

    s_un = st.selectbox("Unmodified / baseline state", ctx.state_labels,
                        index=ctx.state_labels.index("none") if "none" in ctx.state_labels else 0,
                        key="tf_unmod")
    ox_opts = [s for s in ctx.state_labels if s != s_un] or ctx.state_labels
    s_ox = st.selectbox("Oxidized state", ox_opts,
                        index=ox_opts.index("o8G@7") if "o8G@7" in ox_opts else 0,
                        key="tf_ox")
    un = ctx.strong_set(s_un)
    ox = ctx.strong_set(s_ox)
    lost, gained, retained = un - ox, ox - un, un & ox

    c1, c2, c3 = st.columns(3)
    c1.metric("Lost", len(lost))
    c2.metric("Gained", len(gained))
    c3.metric("Retained", len(retained))

    tfs = _tf_symbols_from_libraries()
    amplify = st.checkbox(
        "One-hop regulon amplification (targets of TFs that are themselves gained/lost)",
        value=False,
        help="Optional; uses TF→target edges from installed TF GMTs when available.",
    )

    def _tf_table(genes: set[str], label: str) -> pd.DataFrame:
        hits = sorted(g for g in genes if g.upper() in tfs)
        return pd.DataFrame({"symbol": hits, "category": label})

    tab = pd.concat(
        [_tf_table(lost, "lost"), _tf_table(gained, "gained"), _tf_table(retained, "retained")],
        ignore_index=True,
    )
    st.dataframe(tab, width="stretch", height=280, hide_index=True)
    st.download_button(
        "⬇ Download TF targets (CSV)",
        tab.to_csv(index=False),
        file_name=f"{ctx.mirna}_{s_un}_vs_{s_ox}_TF_targets.csv",
        mime="text/csv",
    )

    if amplify and tfs:
        st.caption("Amplification lists genes co-annotated with selected gained/lost TFs in TF GMTs.")
        seed_tfs = {g.upper() for g in (lost | gained) if g.upper() in tfs}
        amplified: set[str] = set()
        for lib in TF_LIBRARIES:
            path = os.path.join(GENESET_DIR, LIBRARY_FILES[lib])
            if not os.path.exists(path):
                continue
            with open(path) as fh:
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 3:
                        continue
                    head = parts[0].split("_")[0].upper()
                    members = {g.split(",")[0].strip().upper() for g in parts[2:] if g.strip()}
                    if head in seed_tfs or (seed_tfs & members):
                        amplified |= members
        amplified -= un | ox
        st.write(f"Extra one-hop genes (not in either strong set): {len(amplified)}")
        if amplified:
            st.dataframe(pd.DataFrame({"symbol": sorted(amplified)}).head(200), hide_index=True)

    pool = un | ox
    if len(lost) >= 5:
        e = enrich_within_pool(list(lost), pool, library=ctx.library, top=15)
        st.markdown("Within-pool enrichment of **lost** genes (control for differential TF/pathway claims)")
        st.dataframe(e[["term", "overlap", "odds_ratio", "q_value"]] if len(e) else e, hide_index=True)
