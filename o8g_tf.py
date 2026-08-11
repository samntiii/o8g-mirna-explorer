"""Transcription-factor enrichment among gained / lost / retained target sets."""
from __future__ import annotations

import os
import re

import pandas as pd
import streamlit as st

from o8g_enrich import (
    TF_LIBRARIES,
    TF_MOTIF_LIBRARIES,
    LIBRARY_FILES,
    GENESET_DIR,
    available_tf_libraries,
    enrich,
    enrich_within_pool,
)
from o8g_sections import SectionContext


def _parse_tf_name(term: str) -> str:
    """Best-effort TF symbol from Enrichr term labels."""
    t = str(term).strip()
    # e.g. "HINFP (human)", "Gata1 (mouse)", "SP1_HUMAN", "FOXO1_..."
    t = re.sub(r"\s*\(.*\)\s*$", "", t)
    head = t.split("_")[0].split(" ")[0]
    return head.upper()


def _tf_symbols_from_libraries(libs: list[str]) -> set[str]:
    out: set[str] = set()
    for lib in libs:
        path = os.path.join(GENESET_DIR, LIBRARY_FILES[lib])
        if not os.path.exists(path):
            continue
        with open(path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                out.add(_parse_tf_name(parts[0]))
    return {g for g in out if g.isalpha() and 2 <= len(g) <= 12}


def render(ctx: SectionContext) -> None:
    st.markdown(f"#### Transcription-factor enrichment · `{ctx.precision_mode}`")
    st.caption(
        "Primary output is **enriched TFs** (motif / ChIP regulons) among lost or gained "
        "targets — not just whether a TF gene is itself a direct miRNA target. "
        "Motif libraries: TRANSFAC+JASPAR, JASPAR 2025, Genome Browser PWMs "
        "(UCSC/ENCODE-style; Ensembl Regulatory Build–adjacent). "
        "Rebuild missing GMTs with `python scripts/fetch_tf_genesets.py`."
    )

    present = available_tf_libraries()
    missing = [LIBRARY_FILES[k] for k in TF_LIBRARIES if k not in present]
    if missing:
        st.warning("Missing TF gene-set files under `genesets/`: " + ", ".join(missing[:6])
                   + ("…" if len(missing) > 6 else ""))
    if not present:
        st.error("No TF libraries installed — cannot run enrichment.")
        return

    default_libs = [k for k in TF_MOTIF_LIBRARIES if k in present] or present[:3]
    libs = st.multiselect(
        "TF / motif libraries",
        options=present,
        default=default_libs,
        key="tf_libs",
        help="Motif PWMs first; add ChEA/ENCODE/TRRUST for ChIP / curated edges.",
    )
    if not libs:
        st.info("Select at least one TF library.")
        return

    s_un = st.selectbox(
        "Unmodified / baseline state",
        ctx.state_labels,
        index=ctx.state_labels.index("none") if "none" in ctx.state_labels else 0,
        key="tf_unmod",
    )
    ox_opts = [s for s in ctx.state_labels if s != s_un] or ctx.state_labels
    s_ox = st.selectbox(
        "Oxidized state",
        ox_opts,
        index=ox_opts.index("o8G@7") if "o8G@7" in ox_opts else 0,
        key="tf_ox",
    )
    un = ctx.strong_set(s_un)
    ox = ctx.strong_set(s_ox)
    lost, gained, retained = un - ox, ox - un, un & ox
    pool = un | ox

    c1, c2, c3 = st.columns(3)
    c1.metric("Lost", len(lost))
    c2.metric("Gained", len(gained))
    c3.metric("Retained", len(retained))

    bg_mode = st.radio(
        "Enrichment background",
        ["UTR universe", "Within target pool (A∪B)"],
        horizontal=True,
        key="tf_bg",
        help="Within-pool is the correct control for differential (oxidation-change) claims.",
    )
    query_which = st.multiselect(
        "Gene sets to enrich",
        ["lost", "gained", "retained"],
        default=["lost", "gained"],
        key="tf_query_sets",
    )
    top_n = st.slider("Top terms per library", 5, 50, 15, key="tf_top")
    run = st.checkbox("Run TF enrichment", value=False, key="tf_run")
    if not run:
        st.info("Check **Run TF enrichment** to score selected libraries (can be slow).")
    else:
        gene_sets = {"lost": lost, "gained": gained, "retained": retained}
        frames = []
        for qname in query_which:
            genes = gene_sets.get(qname, set())
            if len(genes) < 5:
                st.caption(f"Skip `{qname}` — need ≥5 genes (have {len(genes)}).")
                continue
            for lib in libs:
                if bg_mode.startswith("Within"):
                    e = enrich_within_pool(list(genes), pool, library=lib, top=top_n)
                else:
                    e = enrich(list(genes), library=lib, background=ctx.universe, top=top_n)
                if e is None or e.empty:
                    continue
                e = e.copy()
                # enrich() already sets library; assign query/tf columns without insert()
                e = e.assign(query=qname, library=lib, tf_symbol=e["term"].map(_parse_tf_name))
                frames.append(e)
        if not frames:
            st.warning("No enrichment rows returned.")
        else:
            out = pd.concat(frames, ignore_index=True)
            # compact display
            show = out[
                [c for c in ["query", "library", "tf_symbol", "term", "overlap",
                             "odds_ratio", "p_value", "q_value"] if c in out.columns]
            ].sort_values(["query", "q_value", "p_value"])
            st.markdown("##### Enriched TF regulons")
            st.dataframe(show, hide_index=True, height=380, use_container_width=True)
            sig = show[show["q_value"] < 0.05] if "q_value" in show.columns else show.iloc[0:0]
            st.caption(f"Terms with q < 0.05: {len(sig)} / {len(show)}")
            st.download_button(
                "⬇ Download TF enrichment (CSV)",
                show.to_csv(index=False),
                file_name=f"{ctx.mirna}_{s_un}_vs_{s_ox}_TF_enrichment.csv",
                mime="text/csv",
                key="tf_enrich_dl",
            )

    with st.expander("Direct TF genes among targets (secondary)"):
        st.caption(
            "Intersection of target symbols with TF names parsed from selected libraries. "
            "This is **not** motif enrichment — use the table above for that."
        )
        tfs = _tf_symbols_from_libraries(libs)
        rows = []
        for label, genes in (("lost", lost), ("gained", gained), ("retained", retained)):
            for g in sorted(x for x in genes if x.upper() in tfs):
                rows.append({"symbol": g, "category": label})
        tab = pd.DataFrame(rows)
        st.dataframe(tab if len(tab) else pd.DataFrame(columns=["symbol", "category"]),
                     hide_index=True, height=220, use_container_width=True)
        if len(tab):
            st.download_button(
                "⬇ Download direct TF targets (CSV)",
                tab.to_csv(index=False),
                file_name=f"{ctx.mirna}_{s_un}_vs_{s_ox}_TF_direct.csv",
                mime="text/csv",
                key="tf_direct_dl",
            )
