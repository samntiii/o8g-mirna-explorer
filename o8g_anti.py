"""Antagomir design view — state discrimination, not AGO-loaded seed targeting."""
from __future__ import annotations

import pandas as pd
import streamlit as st

from o8g_energy import design_antagomir, position2_warning
from o8g_sections import SectionContext


def render(ctx: SectionContext) -> None:
    st.markdown(f"#### Antagomir design · `{ctx.precision_mode}`")
    st.warning(
        "Oxidized-state ViennaRNA ΔG uses a **G→U proxy** (ViennaRNA has no o8G:A "
        "Hoogsteen parameter). `dG_normal` and collateral tables are trustworthy; "
        "`dG_oxo` / `ddG` / fold-preference are directionally right with placeholder magnitudes."
    )
    st.info(
        "An antagomir **is not AGO-loaded and has no seed**. This view scores "
        "**state discrimination** only — it must never present an mRNA target list "
        "derived from the oligo."
    )

    ox_lab = st.selectbox(
        "Design against oxidation state",
        [s for s in ctx.state_labels if s != "none"] or ctx.state_labels,
        index=0,
        key="anti_state",
    )
    # parse o8G@2,7 → (2,7)
    positions: tuple[int, ...] = ()
    if ox_lab and ox_lab != "none":
        positions = tuple(
            int(x) for x in ox_lab.replace("o8G@", "").split(",") if x.strip().isdigit()
        )

    if position2_warning(positions):
        st.error(
            "Position-2-only oxidation states are **poor antagomir targets**: "
            "median fold preference ≈8×, ~64% of designs below 10×, versus ≈150× "
            "better for interior positions. Prefer an interior o8G state."
        )

    design = design_antagomir(ctx.info.get("seq_dna", ""), positions)
    st.code(design.sequence)
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("dG_normal", f"{design.dG_normal:.2f}" if design.dG_normal is not None else "—")
    m2.metric("dG_oxo (proxy)", f"{design.dG_oxo:.2f}" if design.dG_oxo is not None else "—")
    m3.metric("ddG", f"{design.ddG:.2f}" if design.ddG is not None else "—")
    m4.metric(
        "Fold preference",
        f"{design.fold_preference:.1f}×" if design.fold_preference is not None else "—",
    )
    st.caption(design.note)

    # Feasibility verdict
    fp = design.fold_preference or 0.0
    if fp >= 50:
        st.success("Feasibility: promising discrimination (fold ≥ 50×, directional).")
    elif fp >= 10:
        st.warning("Feasibility: modest discrimination (10–50×).")
    else:
        st.error("Feasibility: weak discrimination (<10×) — reconsider state or chemistry.")

    # Collateral: other miRNAs sharing the reverse-complement 7mer seed window
    st.markdown("##### Collateral mature scan (sequence identity, not target prediction)")
    try:
        mir_df = ctx.db.list_mirnas() if hasattr(ctx.db, "list_mirnas") else None
        if mir_df is None and hasattr(ctx.db, "mirnas"):
            mir_df = ctx.db.mirnas()
        if mir_df is None and hasattr(ctx.db, "_con"):
            mir_df = pd.read_sql("SELECT mirna, seq_dna, seed FROM mirnas", ctx.db._con)
        if mir_df is not None and len(mir_df):
            oligo = design.sequence
            hits = []
            for _, row in mir_df.iterrows():
                seq = str(row.get("seq_dna", "")).upper()
                if not seq or row["mirna"] == ctx.mirna:
                    continue
                # simple 7-nt window identity vs oligo
                best = 0
                for i in range(max(1, len(oligo) - 6)):
                    win = oligo[i : i + 7]
                    if win and win in seq:
                        best = max(best, 7)
                    elif any(win[j:] in seq or win[: 7 - j] in seq for j in range(1, 3)):
                        best = max(best, 5)
                if best >= 6:
                    hits.append({"mirna": row["mirna"], "seed": row.get("seed", ""), "match_nt": best})
            hit_df = pd.DataFrame(hits).sort_values("match_nt", ascending=False) if hits else pd.DataFrame()
            st.dataframe(hit_df.head(30) if len(hit_df) else hit_df, hide_index=True)
            st.caption(f"{len(hit_df)} other matures with ≥6 nt oligo window match (collateral risk).")
        else:
            st.caption("miRNA table unavailable for collateral scan.")
    except Exception as e:
        st.caption(f"Collateral scan skipped ({type(e).__name__}: {e}).")
