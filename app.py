"""
o8G-miRNA Retargeting Explorer  —  Streamlit app
=================================================
Browse human miRNA seeds, enumerate 8-oxoguanine (o8G) seed-oxidation states,
inspect the predicted mRNA target list for each state, and compare target
lists by pathway enrichment (volcano + heatmap).

Run:  streamlit run app.py
Data: o8g_targets.db, o8g_states.parquet, genesets/  (same folder)
Engine modules: o8g_engine.py, o8g_scanner.py, o8g_enrich.py, o8g_plots.py, o8g_db.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import streamlit as st

from o8g_engine import extract_seed, enumerate_states, g_positions, SeedState
from o8g_db import TargetDB
from o8g_enrich import enrich, compare_states, available_libraries
from o8g_genes import ID_TYPES, GeneResolver
from o8g_precision import PrecisionMode, PrecisionConfig
import o8g_plots as plots
import importlib as _importlib
plots = _importlib.reload(plots)   # pick up edits to the plotting module on rerun

st.set_page_config(page_title="o8G-miRNA Retargeting Explorer", layout="wide",
                   initial_sidebar_state="expanded")

DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "o8g_targets.db")
RANK_LABEL = {1: "6mer", 2: "7mer-A1", 3: "7mer-m8", 4: "8mer"}


# Bump when TargetDB / GeneResolver constructor API changes so Streamlit
# drops stale @cache_resource instances across hot-reloads.
_CACHE_VER = 6


@st.cache_resource
def get_db(_ver: int = _CACHE_VER):
    return TargetDB(DB_PATH)

@st.cache_resource
def get_scanner(_ver: int = _CACHE_VER):
    """Lazy scanner for on-the-fly 6mer-level scans (only built if needed)."""
    from o8g_scanner import TargetScanner
    utr = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utr3_human.parquet")
    if not os.path.exists(utr):
        return None
    return TargetScanner.from_parquet(utr)

@st.cache_resource
def get_resolver(_ver: int = _CACHE_VER):
    try:
        return GeneResolver()
    except FileNotFoundError:
        return None

db = get_db()
# Heal stale cached TargetDB instances from before reverse-index / precision APIs.
if not hasattr(db, "targets_filtered") or not hasattr(db, "retarget_partition"):
    get_db.clear()
    db = get_db()
resolver = get_resolver()


# ---------- helpers ----------
def seed_html(seed: str, oxidized: set) -> str:
    """Render seed with G colored; oxidized Gs boxed red (o8G), normal G blue."""
    out = []
    for i, b in enumerate(seed):
        pos = i + 2
        if b == "G" and pos in oxidized:
            out.append(f"<span style='background:#c44e52;color:white;padding:2px 5px;"
                       f"border-radius:3px;font-weight:700' title='o8G at position {pos}'>{b}<sub>o8</sub></span>")
        elif b == "G":
            out.append(f"<span style='background:#dbe4f0;color:#4c72b0;padding:2px 5px;"
                       f"border-radius:3px;font-weight:700' title='G at position {pos}'>{b}</span>")
        else:
            out.append(f"<span style='padding:2px 4px;color:#333'>{b}</span>")
    return ("<div style='font-family:monospace;font-size:26px;letter-spacing:2px'>"
            "<span style='color:#999;font-size:13px'>5'-</span>" + "".join(out) +
            "<span style='color:#999;font-size:13px'>-3'</span></div>"
            "<div style='font-family:monospace;font-size:10px;color:#999;letter-spacing:2px'>"
            "&nbsp;&nbsp;&nbsp;&nbsp;" + "".join(f"&nbsp;{p}&nbsp;&nbsp;" for p in range(2,9)) + "</div>")


# ---------- sidebar: miRNA selection ----------
st.sidebar.title("🧬 o8G-miRNA Explorer")
st.sidebar.caption("Predicting miRNA target retargeting by 8-oxoguanine seed oxidation")

mir_df = db.mirnas()
query = st.sidebar.text_input("Search miRNA (optional filter)", value="",
                              placeholder="e.g. miR-1, miR-124, let-7 — leave blank for all",
                              help="Type part of a miRNA name to narrow the list. "
                                   "Leave blank to browse all miRNAs in the dropdown below.")
matches = mir_df[mir_df["mirna"].str.contains(query, case=False, na=False)] if query else mir_df
if len(matches) == 0:
    st.sidebar.warning("No miRNA matches that search.")
    st.stop()
# preselect hsa-miR-1-3p when present so the app opens on a worked example
opts = matches["mirna"].tolist()
default_idx = 0
for _pref in ("hsa-miR-1-3p", "hsa-miR-1-2-3p", "hsa-miR-1"):
    if _pref in opts:
        default_idx = opts.index(_pref)
        break
mirna = st.sidebar.selectbox(
    f"Select miRNA ({len(matches)} of {len(mir_df)} shown)",
    opts, index=default_idx)
info = db.mirna_info(mirna)
seed = info["seed"]
gpos = g_positions(seed)

st.sidebar.markdown("---")
st.sidebar.markdown(f"**Accession:** `{info['accession']}`")
st.sidebar.markdown(f"**Mature (DNA):** `{info['seq_dna']}`")
st.sidebar.markdown(f"**Seed (pos 2–8):** `{seed}`")
st.sidebar.markdown(f"**Guanines in seed:** {len(gpos)} at position(s) {gpos or '—'}")
st.sidebar.markdown(f"**Oxidation states:** {2**len(gpos)}")

lib = st.sidebar.selectbox("Pathway library", available_libraries(), index=0)

st.sidebar.markdown("### Precision mode")
_MODE_OPTIONS = [m.value for m in PrecisionMode]  # Sensitive / Stringent / Consensus
mode_name = st.sidebar.radio(
    "Target-list filter",
    _MODE_OPTIONS,
    index=0,  # Sensitive — discovery default
    key=f"precision_mode_v{_CACHE_VER}",
    help=(
        "Sensitive: 7mer-m8+8mer. Stringent: 8mer only. "
        "Consensus: + TargetScan conserved (unmodified baseline). "
        "Filters apply to each oxidation state before gained/lost."
    ),
)
precision_cfg = PrecisionConfig.from_mode(mode_name)
st.sidebar.caption(
    "Paper / claims: use **Stringent** or **Consensus**. "
    "Gained/lost are always computed after filtering both states."
)
scanner = None
# Legacy slider kept as an additional floor within the mode (rarely needed)
min_rank = st.sidebar.select_slider(
    "Extra min site quality floor",
    options=[1, 2, 3, 4],
    value=3,
    format_func=lambda r: {1: "6mer", 2: "7mer-A1", 3: "7mer-m8", 4: "8mer"}[r],
)
st.sidebar.caption("Applied on top of the precision mode (usually leave at 7mer-m8).")


# ---------- main ----------
st.title(f"{mirna}")
if len(gpos) == 0:
    st.info("This miRNA's seed contains no guanine, so 8-oxoG oxidation cannot alter its "
            "seed pairing — only one target list exists. Choose a G-containing seed to explore retargeting.")

st.markdown("#### Seed sequence")
st.caption("Guanines are oxidation-eligible. Blue = normal G (pairs C); red = o8G (pairs A).")

states = enumerate_states(seed)
state_labels = [s.label for s in states]

tab_browse, tab_compare, tab_all, tab_gene = st.tabs(
    ["Single state — targets", "Compare two states", "All states",
     "Gene → miRNA/oxomiR"])

# ===== TAB 1: single state target list =====
with tab_browse:
    c1, c2 = st.columns([3, 2])
    with c1:
        ox_choice = st.multiselect("Oxidize which guanine position(s)?", gpos, default=[],
            help="Select seed G positions to model as 8-oxoG (o8G). None = fully normal seed.")
        cur = SeedState(seed, tuple(sorted(ox_choice)))
        st.markdown(seed_html(seed, set(ox_choice)), unsafe_allow_html=True)
        st.markdown(f"**State:** `{cur.label}`  •  target motifs (5'→3' on mRNA):")
        st.code(f"8mer   : {cur.motifs['8mer']}\n7mer-m8: {cur.motifs['7mer-m8']}\n"
                f"7mer-A1: {cur.motifs['7mer-A1']}\n6mer   : {cur.motifs['6mer']}")
    with c2:
        meta = db.states_for_seed(seed)
        row = meta[meta["label"] == cur.label]
        if len(row):
            r = row.iloc[0]
            st.metric("Strong-site target genes (7mer-m8 + 8mer)", int(r["n_strong"]))
            st.caption(f"Site inventory: {int(r['n_8mer'])} 8mer · {int(r['n_7mer_m8'])} 7mer-m8 · "
                       f"{int(r['n_7mer_A1'])} 7mer-A1 · {int(r['n_6mer'])} 6mer")

    tdf = db.targets_filtered(
        seed,
        cur.label,
        precision_cfg,
        scanner=scanner,
        mature_dna=info["seq_dna"],
        mirna=mirna,
    )
    if "site_rank" in tdf.columns:
        tdf = tdf[tdf["site_rank"] >= min_rank]
    st.markdown(
        f"**{len(tdf)} predicted target genes** · mode `{precision_cfg.mode.value}` "
        f"(extra floor ≥ {RANK_LABEL[min_rank]})"
    )
    show_cols = [c for c in [
        "symbol", "gene_id", "site_type", "site_rank", "n_8mer", "n_7mer_m8",
        "score", "context_score", "is_conserved",
    ] if c in tdf.columns]
    st.dataframe(tdf[show_cols] if show_cols else tdf, width='stretch', height=320)
    st.download_button("⬇ Download target list (CSV)", tdf.to_csv(index=False),
                       file_name=f"{mirna}_{cur.label}_{precision_cfg.mode.value}_targets.csv", mime="text/csv")

    if len(tdf) >= 5:
        with st.spinner("Running pathway enrichment…"):
            e = enrich(tdf["symbol"].tolist(), lib, top=20)
        st.markdown(f"**Top {lib} pathways** (hypergeometric over-representation)")
        st.dataframe(e[["term","overlap","term_size","odds_ratio","p_value","q_value"]],
                     width='stretch', height=280)

# ===== TAB 2: compare two states → volcano =====
with tab_compare:
    if len(gpos) == 0:
        st.info("No oxidation states to compare for a G-free seed.")
    else:
        cc1, cc2 = st.columns(2)
        sA = cc1.selectbox("State A", state_labels, index=0)
        sB = cc2.selectbox("State B", state_labels, index=min(3, len(state_labels)-1))
        stA = next(s for s in states if s.label == sA)
        stB = next(s for s in states if s.label == sB)
        cc1.markdown(seed_html(seed, set(stA.oxidized_positions)), unsafe_allow_html=True)
        cc2.markdown(seed_html(seed, set(stB.oxidized_positions)), unsafe_allow_html=True)

        tA = db.targets_filtered(
            seed, sA, precision_cfg, scanner=scanner, mature_dna=info["seq_dna"], mirna=mirna
        )
        tB = db.targets_filtered(
            seed, sB, precision_cfg, scanner=scanner, mature_dna=info["seq_dna"], mirna=mirna
        )
        if "site_rank" in tA.columns:
            tA = tA[tA["site_rank"] >= max(min_rank, 3 if precision_cfg.mode != PrecisionMode.STRINGENT else 4)]
            tB = tB[tB["site_rank"] >= max(min_rank, 3 if precision_cfg.mode != PrecisionMode.STRINGENT else 4)]
        gA = tA["symbol"].tolist()
        gB = tB["symbol"].tolist()
        sa, sb = set(gA), set(gB)
        m1,m2,m3 = st.columns(3)
        m1.metric(f"Targets in {sA}", len(sa))
        m2.metric(f"Targets in {sB}", len(sb))
        m3.metric("Shared", len(sa & sb), delta=f"{len(sa^sb)} state-specific")

        # ---- gene-level differential (parallel to the pathway volcano) ----
        st.markdown("#### Differential target genes")
        st.caption(
            f"**Lost** = targeted in `{sA}` but not `{sB}` · "
            f"**Gained** = new targets in `{sB}` · **Shared** = targeted in both. "
            f"Precision mode `{precision_cfg.mode.value}` applied to each state before the set difference."
        )        # merge on symbol, carrying each state's best site type
        mA = tA.drop_duplicates("symbol").set_index("symbol")
        mB = tB.drop_duplicates("symbol").set_index("symbol")
        merged = pd.DataFrame(index=sorted(sa | sb))
        merged[f"site_{sA}"] = mA["site_type"].reindex(merged.index)
        merged[f"site_{sB}"] = mB["site_type"].reindex(merged.index)
        merged["gene_id"] = mA["gene_id"].reindex(merged.index).fillna(
            mB["gene_id"].reindex(merged.index))
        merged["category"] = np.where(
            merged.index.isin(sa & sb), "Shared",
            np.where(merged.index.isin(sa - sb), "Lost", "Gained"))
        merged = merged.reset_index().rename(columns={"index": "symbol"})

        lost   = merged[merged["category"] == "Lost"][["symbol", "gene_id", f"site_{sA}"]]
        gained = merged[merged["category"] == "Gained"][["symbol", "gene_id", f"site_{sB}"]]
        shared = merged[merged["category"] == "Shared"][["symbol", "gene_id", f"site_{sA}", f"site_{sB}"]]

        gc1, gc2, gc3 = st.columns(3)
        with gc1:
            st.markdown(f"**🔴 Lost on {sB}** ({len(lost)})")
            st.dataframe(lost, width='stretch', height=300, hide_index=True)
        with gc2:
            st.markdown(f"**🔵 Gained on {sB}** ({len(gained)})")
            st.dataframe(gained, width='stretch', height=300, hide_index=True)
        with gc3:
            st.markdown(f"**⚪ Shared** ({len(shared)})")
            st.dataframe(shared, width='stretch', height=300, hide_index=True)

        st.download_button(
            "⬇ Download differential gene table (CSV)",
            merged[["symbol", "gene_id", f"site_{sA}", f"site_{sB}", "category"]].to_csv(index=False),
            file_name=f"{mirna}_{sA}_vs_{sB}_differential_genes.csv", mime="text/csv")

        st.markdown("#### Differential pathway enrichment")
        chart = st.radio("Display", ["Diverging bar", "Volcano"], horizontal=True,
            help="Diverging bar names each pathway and shows direction directly — clearest "
                 "when few terms are significant. Volcano is better with many terms.")
        if len(sa) >= 5 and len(sb) >= 5:
            with st.spinner("Computing differential pathway enrichment…"):
                cmp = compare_states(gA, gB, library=lib, label_A=sA, label_B=sB, min_overlap=3)
                cmp["log2_or"] = (-cmp["log2_or_ratio"]).clip(-8, 8)   # x>0 → toward A
                cmp["q_value"] = cmp[[f"q_{sA}", f"q_{sB}"]].min(axis=1)
                cmp["neglog10_q"] = cmp["neglog10_q_best"]
            if chart == "Diverging bar":
                fig = plots.diverging_bar_plotly(cmp, label_A=sA, label_B=sB, x_col="log2_or")
                st.plotly_chart(fig, width='stretch')
                st.caption(f"Each bar is a pathway; length = −log10 q. Right (red) = enriched in "
                           f"**{sA}**; left (blue) = enriched in **{sB}**. Dashed lines = q<0.05.")
            else:
                fig = plots.volcano_plotly(cmp, label_A=sA, label_B=sB, x_col="log2_or")
                st.plotly_chart(fig, width='stretch')
                st.caption(f"Each point is a pathway. Right (red) = enriched in **{sA}**; "
                           f"left (blue) = enriched in **{sB}**. Dashed line = q<0.05.")
            st.dataframe(cmp.sort_values("neglog10_q_best", ascending=False)
                         [["term","library",f"overlap_{sA}",f"overlap_{sB}",
                           f"odds_{sA}",f"odds_{sB}","q_value"]].head(40),
                         width='stretch', height=300)
        else:
            st.warning("Need ≥5 strong-site targets in each state for enrichment.")

# ===== TAB 3: all states heatmap =====
with tab_all:
    if len(gpos) == 0:
        st.info("Only one state exists for a G-free seed.")
    else:
        st.markdown(f"Pathway enrichment (−log10 q) across all **{len(states)}** seed-oxidation "
                    f"states of {mirna}, library **{lib}**.")
        view = st.radio("Display", ["Heatmap", "Dot plot"], horizontal=True,
            help="Heatmap = enrichment intensity per state. Dot plot adds gene-overlap as dot "
                 "size, so you also see how many genes support each term in each state.")
        top_terms = st.slider("Pathways to show", 8, 40, 22)
        sg = {s.label: db.target_symbols(seed, s.label, min_rank=max(min_rank,3)) for s in states}
        states_order = [s.label for s in states]
        if view == "Heatmap":
            with st.spinner("Enriching all states…"):
                mat = plots.enrichment_matrix(sg, enrich, library=lib, top_terms=top_terms)
            if mat.shape[0] >= 2:
                fig = plots.heatmap_plotly(mat)
                st.plotly_chart(fig, width='stretch')
                st.download_button("⬇ Download enrichment matrix (CSV)", mat.to_csv(),
                                   file_name=f"{mirna}_state_enrichment_{lib}.csv", mime="text/csv")
            else:
                st.warning("Not enough enriched pathways to build a heatmap at this setting.")
        else:
            with st.spinner("Enriching all states…"):
                long = plots.state_dotplot_data(sg, enrich, library=lib, top_terms=top_terms)
            if long["term"].nunique() >= 2:
                fig = plots.dotplot_plotly(long, states_order=states_order)
                st.plotly_chart(fig, width='stretch')
                st.caption("Dot size = number of overlapping genes; color = −log10 q. "
                           "Compare a row across columns to see how oxidation shifts each pathway.")
                st.download_button("⬇ Download enrichment table (CSV)", long.to_csv(index=False),
                                   file_name=f"{mirna}_state_dotplot_{lib}.csv", mime="text/csv")
            else:
                st.warning("Not enough enriched pathways to build a dot plot at this setting.")

# ===== TAB 4: reverse gene → miRNA/oxomiR =====
with tab_gene:
    st.markdown("#### Find miRNA / oxomiR states that target a gene")
    st.caption(
        "Offline reverse lookup against strong sites (7mer-m8 + 8mer) in the precomputed "
        "database. Gene IDs resolve via a local NCBI + HGNC map (no live API)."
    )
    if resolver is None:
        st.error(
            "Gene alias map missing. Build it with:\n\n"
            "`python scripts/fetch_gene_aliases.py`\n\n"
            "Then rebuild the reverse index if needed:\n\n"
            "`python scripts/build_reverse_index.py`"
        )
    elif getattr(db, "reverse_path", None) is None:
        st.error(
            "Reverse index missing. Build it with:\n\n"
            "`python scripts/build_reverse_index.py`"
        )
    else:
        gc1, gc2 = st.columns([3, 1])
        gene_q = gc1.text_input(
            "Gene",
            value="HDAC4",
            placeholder="HDAC4 · ENSG00000068024 · 9759",
            help="Symbol, Ensembl gene ID, or Entrez Gene ID.",
        )
        id_type = gc2.selectbox("ID type", ID_TYPES, index=0)

        if gene_q.strip():
            res = resolver.resolve(gene_q, id_type)
            if res.error:
                st.warning(res.error)
            else:
                if res.ambiguous:
                    st.info(
                        f"Ambiguous {res.detected_type} match ({len(res.hits)} genes). "
                        "Showing the first; refine the query or pick Ensembl/Entrez."
                    )
                    pick_labels = [
                        f"{h.symbol} · {h.ensembl}"
                        + (f" · Entrez {h.entrez}" if h.entrez else "")
                        for h in res.hits
                    ]
                    pick = st.selectbox("Choose gene", pick_labels, index=0)
                    hit = res.hits[pick_labels.index(pick)]
                else:
                    hit = res.hits[0]

                meta_bits = [
                    f"**{hit.symbol}**",
                    f"Ensembl `{hit.ensembl}`",
                ]
                if hit.entrez is not None:
                    meta_bits.append(f"Entrez `{hit.entrez}`")
                meta_bits.append(f"matched as {res.detected_type if id_type == 'Auto' else id_type}")
                st.markdown(" · ".join(meta_bits))
                built = resolver.meta.get("built_at", "?")
                st.caption(
                    f"ID map built {built} · sources: {resolver.meta.get('sources', 'local')}"
                )

                with st.spinner("Querying reverse index…"):
                    hits = db.states_targeting_gene(hit.gene_idx)

                if hits.empty:
                    st.info("No strong-site miRNA/oxomiR states target this gene in the database.")
                else:
                    filt = st.radio(
                        "Show states",
                        ["All", "Unmodified only", "Oxidized only",
                         "Gained on oxidation only", "Also in unmodified"],
                        horizontal=True,
                    )
                    view = hits
                    if filt == "Unmodified only":
                        view = hits[hits["state_label"] == "none"]
                    elif filt == "Oxidized only":
                        view = hits[hits["state_label"] != "none"]
                    elif filt == "Gained on oxidation only":
                        view = hits[hits["vs_unmodified"] == "gained on oxidation"]
                    elif filt == "Also in unmodified":
                        view = hits[hits["vs_unmodified"] == "also in unmodified"]

                    n_seeds = view["seed"].nunique()
                    n_ox = int((view["state_label"] != "none").sum())
                    m1, m2, m3 = st.columns(3)
                    m1.metric("Matching states", len(view))
                    m2.metric("Unique seeds", int(n_seeds))
                    m3.metric("Oxidized-state rows", n_ox)

                    show = view[
                        [
                            "mirna",
                            "seed",
                            "state_label",
                            "oxidized_positions",
                            "site_type",
                            "site_rank",
                            "motif_7mer_m8",
                            "motif_8mer",
                            "vs_unmodified",
                            "all_mirnas",
                        ]
                    ].rename(
                        columns={
                            "state_label": "oxidation_state",
                            "vs_unmodified": "vs_unmodified_seed",
                        }
                    )
                    st.dataframe(show, width="stretch", height=420, hide_index=True)
                    st.caption(
                        "**vs_unmodified_seed:** *unmodified* = targeted by the normal seed; "
                        "*also in unmodified* = still targeted when oxidized; "
                        "*gained on oxidation* = only appears under an o8G state for that seed "
                        "(retargeting)."
                    )
                    st.download_button(
                        "⬇ Download miRNA/oxomiR hits (CSV)",
                        show.to_csv(index=False),
                        file_name=f"{hit.symbol}_targeting_oxomirs.csv",
                        mime="text/csv",
                    )

st.markdown("---")
st.caption("Prediction: seed positions 2–8; unoxidized G pairs C, 8-oxoG (o8G) pairs A "
           "(Hoogsteen). Targets = TargetScan-style 6mer/7mer/8mer sites in human 3′UTRs "
           "(longest per gene, Ensembl). Enrichment = hypergeometric ORA, BH-corrected. "
           "Large seed-match target lists yield broad enrichment — interpret pathway shifts, "
           "not absolute significance, and validate candidates experimentally. "
           "Gene ID resolution uses a locally cached NCBI/HGNC map (no runtime API).")
