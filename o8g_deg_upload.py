"""Streamlit UI: RNA-seq / DEG upload → oxomiR concordance ranking."""
from __future__ import annotations

import io

import pandas as pd
import plotly.express as px
import streamlit as st

from o8g_deg import (
    detect_columns,
    read_upload,
    build_deg_sets,
    is_excel_filename,
    list_excel_sheets,
    guess_deg_sheet,
)
from o8g_deg_score import score_panel
from o8g_engine import g_positions
from o8g_sections import SectionContext

DEFAULT_PANEL = [
    "hsa-miR-1-3p",
    "hsa-miR-124-3p",
    "hsa-let-7a-5p",
    "hsa-miR-122-5p",
]


def render(ctx: SectionContext) -> None:
    st.markdown(f"#### RNA-seq / DEG upload · `{ctx.precision_mode}`")
    st.caption(
        "Exploratory concordance: UP DEGs ↔ **lost** strong targets; optional DOWN ↔ gained. "
        "Uploads stay in session memory only (not written to disk). Correlative — not causal. "
        "External catalogs (miRDB etc.) are unmodified-only; Consensus still refuses when "
        "conservation tables are missing."
    )
    st.warning(
        "Not a DE pipeline. Upload differential-expression results (symbol + log2FC + padj), "
        "not raw counts. Search space is a **selected miRNA panel**, not the full catalog."
    )

    up_file = st.file_uploader(
        "DEG table (CSV / TSV / Excel)",
        type=["csv", "tsv", "txt", "xlsx", "xls"],
        key="deg_upload_file",
    )
    if up_file is None:
        st.info("Upload a DEG table to begin. Required columns: gene symbol, log2FC, adjusted p-value.")
        return

    raw = up_file.getvalue()
    fname = up_file.name or "upload.csv"
    sheet_name: str | None = None
    if is_excel_filename(fname):
        try:
            sheets = list_excel_sheets(raw, fname)
        except Exception as e:
            st.error(str(e))
            return
        if len(sheets) > 1:
            guessed = guess_deg_sheet(sheets) or sheets[0]
            g_i = sheets.index(guessed) if guessed in sheets else 0
            sheet_name = st.selectbox(
                "Excel sheet",
                sheets,
                index=g_i,
                key="deg_excel_sheet",
                help="Workbook has multiple sheets — pick the DEG results table.",
            )
            st.caption(f"{len(sheets)} sheets detected · default guess: `{guessed}`")
        elif len(sheets) == 1:
            sheet_name = sheets[0]
            st.caption(f"Excel sheet: `{sheet_name}`")

    try:
        df = read_upload(raw, fname, sheet_name=sheet_name)
    except Exception as e:
        st.error(str(e))
        return

    detected = detect_columns(list(df.columns))
    cols = list(df.columns)
    c1, c2, c3 = st.columns(3)
    with c1:
        sym_i = cols.index(detected["symbol"]) if detected["symbol"] in cols else 0
        symbol_col = st.selectbox("Gene symbol column", cols, index=sym_i, key="deg_sym_col")
    with c2:
        lfc_i = cols.index(detected["lfc"]) if detected["lfc"] in cols else min(1, len(cols) - 1)
        lfc_col = st.selectbox("log2FC column", cols, index=lfc_i, key="deg_lfc_col")
    with c3:
        padj_i = cols.index(detected["padj"]) if detected["padj"] in cols else min(2, len(cols) - 1)
        padj_col = st.selectbox("padj / FDR column", cols, index=padj_i, key="deg_padj_col")

    t1, t2 = st.columns(2)
    lfc_thr = t1.number_input("|log2FC| threshold", min_value=0.0, max_value=5.0, value=0.5, step=0.1)
    padj_thr = t2.number_input("padj threshold", min_value=0.0, max_value=1.0, value=0.05, step=0.01)

    try:
        deg = build_deg_sets(
            df,
            symbol_col=symbol_col,
            lfc_col=lfc_col,
            padj_col=padj_col,
            lfc_thr=float(lfc_thr),
            padj_thr=float(padj_thr),
            universe={u.upper() for u in ctx.universe},
        )
    except Exception as e:
        st.error(str(e))
        return

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("|UP|", len(deg.up))
    m2.metric("|DOWN|", len(deg.down))
    m3.metric("Rows kept", len(deg.table))
    m4.metric("Unmatched (UTR)", deg.unmatched)
    for note in deg.notes:
        st.caption(note)
    if not deg.up and not deg.down:
        st.warning("No UP/DOWN genes after thresholds ∩ UTR universe. Relax thresholds or check columns.")
        return

    # --- miRNA panel ---
    all_mir = []
    try:
        all_mir = ctx.db.mirnas()["mirna"].astype(str).tolist()
    except Exception:
        all_mir = [ctx.mirna]
    starter = [m for m in DEFAULT_PANEL if m in all_mir]
    if ctx.mirna and ctx.mirna not in starter:
        starter = [ctx.mirna] + starter

    pick = st.multiselect(
        "miRNA panel",
        options=all_mir,
        default=starter,
        key="deg_mirna_panel",
    )
    if st.button("Add current sidebar miRNA", key="deg_add_sidebar"):
        if ctx.mirna and ctx.mirna not in pick:
            pick = list(pick) + [ctx.mirna]
            st.session_state["deg_mirna_panel"] = pick
            st.rerun()

    if not pick:
        st.info("Select at least one miRNA.")
        return

    # Cap warning for combinatorial seeds
    warn_seeds = []
    for mir in pick:
        info = ctx.db.mirna_info(mir)
        if not info:
            continue
        n = 2 ** len(g_positions(info["seed"]))
        if n > 32:
            warn_seeds.append(f"{mir} (2^{len(g_positions(info['seed']))}={n})")
    single_g = st.checkbox(
        "Scan single-G oxidation states only (recommended when 2^k > 32)",
        value=bool(warn_seeds),
        key="deg_single_g",
    )
    if warn_seeds and not single_g:
        st.warning(
            "Large state spaces: " + ", ".join(warn_seeds) + ". "
            "Enable single-G only or the enumerator will auto-cap."
        )

    use_oboe = st.checkbox(
        "Re-rank by OBOE oxo-G prior (local RNABERT / GC fallback)",
        value=False,
        key="deg_oboe",
        help="Multiplies concordance by mean P(o8G) of oxidized seed positions "
        "(fine-tuned OBOE RNABERT when available).",
    )

    down_w = st.slider(
        "DOWN∩gained weight in concordance score (0 = UP∩lost only)",
        min_value=0.0,
        max_value=1.0,
        value=0.0,
        step=0.1,
        key="deg_down_w",
    )

    run = st.checkbox("Run concordance scoring (can be slow for large panels)", value=False, key="deg_run")
    if not run:
        st.info("Check **Run** to score selected miRNAs × oxidation states against your DEGs.")
        return

    cfg = getattr(ctx, "precision_cfg", None) or ctx.precision_mode
    with st.spinner(f"Scoring {len(pick)} miRNA(s)…"):
        ranked = score_panel(
            ctx.db,
            pick,
            up=deg.up,
            down=deg.down,
            universe={u.upper() for u in ctx.universe},
            precision_cfg=cfg,
            single_g_only=single_g,
            max_states=32,
            down_weight=float(down_w),
            scanner=getattr(ctx, "scanner", None),
        )

    if ranked.empty:
        st.warning("No states scored.")
        return

    if use_oboe and "ox_label" in ranked.columns:
        try:
            import importlib
            import o8g_oboe as _oboe
            import o8g_oboe_model as _oboe_model

            _oboe_model = importlib.reload(_oboe_model)
            _oboe = importlib.reload(_oboe)

            priors = []
            for mir in ranked["mirna"].unique():
                info = ctx.db.mirna_info(mir)
                if not info:
                    continue
                rt = _oboe.rank_oxidation_states(info["seq_dna"])
                if rt.empty:
                    continue
                rt = rt.assign(mirna=mir)
                priors.append(rt[["mirna", "ox_label", "mean_oboe_prior"]])
            if priors:
                pr = pd.concat(priors, ignore_index=True)
                ranked = ranked.merge(pr, on=["mirna", "ox_label"], how="left")
                ranked["concordance_score"] = ranked["concordance_score"] * (
                    0.25 + ranked["mean_oboe_prior"].fillna(0.5)
                )
                ranked = ranked.sort_values("concordance_score", ascending=False).reset_index(drop=True)
        except Exception as e:
            st.caption(f"OBOE re-rank skipped: {e}")

    show_cols = [
        c
        for c in [
            "mirna",
            "ox_label",
            "n_lost",
            "n_up_lost",
            "odds_up_lost",
            "q_up_lost",
            "jaccard_up_lost",
            "n_down_gained",
            "q_down_gained",
            "mean_oboe_prior",
            "concordance_score",
        ]
        if c in ranked.columns
    ]
    st.dataframe(ranked[show_cols], hide_index=True, height=360, use_container_width=True)

    top = ranked.head(15).copy()
    if len(top) and "concordance_score" in top.columns:
        top["label"] = top["mirna"].astype(str) + " · " + top["ox_label"].astype(str)
        fig = px.bar(
            top.sort_values("concordance_score"),
            x="concordance_score",
            y="label",
            orientation="h",
            title="Top concordance scores (UP∩lost primary)",
            labels={"concordance_score": "Concordance score", "label": ""},
        )
        fig.update_layout(height=max(280, 28 * len(top)), margin=dict(l=10, r=10, t=40, b=10))
        st.plotly_chart(fig, use_container_width=True)

    with st.expander("Intersection gene lists / download"):
        row_lab = st.selectbox(
            "Detail row",
            options=[f"{r.mirna} | {r.ox_label}" for r in ranked.itertuples()],
            key="deg_detail_row",
        )
        mir_sel, ox_sel = [x.strip() for x in row_lab.split("|", 1)]
        detail = ranked[(ranked["mirna"] == mir_sel) & (ranked["ox_label"] == ox_sel)].iloc[0]
        st.markdown(
            f"**UP ∩ lost** (n={int(detail.get('n_up_lost', 0))}): "
            f"`{detail.get('up_lost_genes', '')}`"
        )
        st.markdown(
            f"**DOWN ∩ gained** (n={int(detail.get('n_down_gained', 0) or 0)}): "
            f"`{detail.get('down_gained_genes', '')}`"
        )
        buf = io.StringIO()
        ranked.to_csv(buf, index=False)
        st.download_button(
            "Download full ranked table (CSV)",
            data=buf.getvalue(),
            file_name="deg_oxomir_concordance.csv",
            mime="text/csv",
            key="deg_dl",
        )
