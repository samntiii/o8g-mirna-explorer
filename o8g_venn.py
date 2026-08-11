"""Overlap view: Venn (≤3 sets) or UpSet (>3); region gene table + CSV."""
from __future__ import annotations

import importlib

import pandas as pd
import streamlit as st

from o8g_sections import SectionContext


def render(ctx: SectionContext) -> None:
    st.markdown(f"#### Overlap across oxidation states · `{ctx.precision_mode}`")
    st.caption(
        "Sets are filtered with the sidebar precision mode before overlap. "
        "2–3 states → Venn; 4+ states → UpSet."
    )
    default = [s for s in ("none", "o8G@7") if s in ctx.state_labels] or ctx.state_labels[:2]
    chosen = st.multiselect(
        "States to compare",
        ctx.state_labels,
        default=default,
        key="overlap_states",
    )
    if len(chosen) < 2:
        st.info("Select at least two states.")
        return

    sets = {lab: set(ctx.strong_set(lab)) for lab in chosen}
    sizes = {k: len(v) for k, v in sets.items()}
    st.write(sizes)

    import o8g_plots as plots

    plots = importlib.reload(plots)

    # ----- plot first -----
    try:
        if len(chosen) == 2:
            a, b = chosen
            fig = plots.venn2_plotly(sets[a], sets[b], label_a=a, label_b=b)
            st.plotly_chart(fig, width="stretch")
        elif len(chosen) == 3:
            a, b, c = chosen
            fig = plots.venn3_plotly(sets[a], sets[b], sets[c], label_a=a, label_b=b, label_c=c)
            st.plotly_chart(fig, width="stretch")
        else:
            fig = plots.upset_plotly(sets)
            st.plotly_chart(fig, width="stretch")
    except Exception as e:
        st.warning(f"Could not render overlap plot ({type(e).__name__}: {e}).")

    # ----- region gene table (2–3) or size table (4+) -----
    labels = list(sets.keys())
    if len(labels) <= 3:
        # Build membership without materializing every gene row twice
        if len(labels) == 2:
            a, b = labels
            regions = {
                f"only {a}": sorted(sets[a] - sets[b]),
                f"only {b}": sorted(sets[b] - sets[a]),
                "shared": sorted(sets[a] & sets[b]),
            }
        else:
            a, b, c = labels
            regions = {
                f"only {a}": sorted(sets[a] - sets[b] - sets[c]),
                f"only {b}": sorted(sets[b] - sets[a] - sets[c]),
                f"only {c}": sorted(sets[c] - sets[a] - sets[b]),
                f"{a}∩{b} only": sorted((sets[a] & sets[b]) - sets[c]),
                f"{a}∩{c} only": sorted((sets[a] & sets[c]) - sets[b]),
                f"{b}∩{c} only": sorted((sets[b] & sets[c]) - sets[a]),
                "all three": sorted(sets[a] & sets[b] & sets[c]),
            }
        cols = st.columns(min(4, len(regions)))
        for i, (name, genes) in enumerate(regions.items()):
            cols[i % len(cols)].metric(name, len(genes))

        region_name = st.selectbox("Show genes in region", list(regions.keys()), key="overlap_region")
        genes = regions[region_name]
        show = pd.DataFrame({"symbol": genes})
        for lab in labels:
            show[lab] = [g in sets[lab] for g in genes]
        st.dataframe(show, width="stretch", height=280, hide_index=True)
        st.download_button(
            "⬇ Download region genes (CSV)",
            show.to_csv(index=False),
            file_name=f"{ctx.mirna}_overlap_{region_name.replace(' ', '_')}.csv",
            mime="text/csv",
            key="overlap_region_dl",
        )
    else:
        st.dataframe(
            pd.DataFrame({"state": list(sizes), "n_genes": list(sizes.values())}),
            hide_index=True,
        )
        st.caption("Select ≤3 states above for a per-region gene table; 4+ uses UpSet.")
