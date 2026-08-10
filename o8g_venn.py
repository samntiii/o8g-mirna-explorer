"""Overlap view: Venn (≤3 sets) or UpSet (>3); region gene table + CSV."""
from __future__ import annotations

import pandas as pd
import streamlit as st

from o8g_sections import SectionContext


def render(ctx: SectionContext) -> None:
    st.markdown(f"#### Overlap across oxidation states · `{ctx.precision_mode}`")
    st.caption(
        "Sets are filtered with the sidebar precision mode before overlap. "
        "≤3 states → Venn-style region counts; more → UpSet via the plot helper."
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
    sets = {lab: ctx.strong_set(lab) for lab in chosen}
    sizes = {k: len(v) for k, v in sets.items()}
    st.write({f"{k}": v for k, v in sizes.items()})

    if len(chosen) <= 3:
        # region table for 2–3 sets
        labels = list(sets.keys())
        universe = set().union(*sets.values())
        rows = []
        for g in sorted(universe):
            membership = tuple(g in sets[lab] for lab in labels)
            rows.append({"symbol": g, **{lab: membership[i] for i, lab in enumerate(labels)}})
        tab = pd.DataFrame(rows)
        if len(labels) == 2:
            a, b = labels
            regions = {
                f"only {a}": tab[tab[a] & ~tab[b]],
                f"only {b}": tab[~tab[a] & tab[b]],
                "shared": tab[tab[a] & tab[b]],
            }
        else:
            a, b, c = labels
            regions = {
                f"only {a}": tab[tab[a] & ~tab[b] & ~tab[c]],
                f"only {b}": tab[~tab[a] & tab[b] & ~tab[c]],
                f"only {c}": tab[~tab[a] & ~tab[b] & tab[c]],
                f"{a}∩{b} only": tab[tab[a] & tab[b] & ~tab[c]],
                f"{a}∩{c} only": tab[tab[a] & ~tab[b] & tab[c]],
                f"{b}∩{c} only": tab[~tab[a] & tab[b] & tab[c]],
                "all three": tab[tab[a] & tab[b] & tab[c]],
            }
        cols = st.columns(min(3, len(regions)))
        for i, (name, rdf) in enumerate(regions.items()):
            cols[i % len(cols)].metric(name, len(rdf))
        region_name = st.selectbox("Show genes in region", list(regions.keys()))
        show = regions[region_name][["symbol"] + labels]
        st.dataframe(show, width="stretch", height=280, hide_index=True)
        st.download_button(
            "⬇ Download region genes (CSV)",
            show.to_csv(index=False),
            file_name=f"{ctx.mirna}_overlap_{region_name.replace(' ', '_')}.csv",
            mime="text/csv",
        )
    else:
        try:
            import o8g_plots as plots

            fig = plots.upset_plotly(sets) if hasattr(plots, "upset_plotly") else None
            if fig is not None:
                st.plotly_chart(fig, width="stretch")
            else:
                st.dataframe(
                    pd.DataFrame({"state": list(sizes), "n_genes": list(sizes.values())}),
                    hide_index=True,
                )
                st.caption("UpSet helper unavailable — showing set sizes only.")
        except Exception as e:
            st.warning(f"Could not render UpSet ({type(e).__name__}: {e}).")
            st.dataframe(
                pd.DataFrame({"state": list(sizes), "n_genes": list(sizes.values())}),
                hide_index=True,
            )
