"""
o8g_plots.py
============
Plotting for the o8G retargeting analysis.

volcano_static / volcano_plotly:
    Differential pathway volcano (one seed-state vs another).
    x = log2 odds-ratio (Fisher differential),  y = -log10 BH q.
    Points right of 0 skew toward state A (lost on oxidation); left toward B.

diverging_bar_plotly:
    Two-state differential enrichment as a diverging horizontal bar — one named
    row per pathway, length = -log10 q, direction by color. Clearer than a
    volcano when few terms are significant.

state_dotplot_data / dotplot_plotly:
    Dot matrix across all 2^k seed-oxidation states: dot size = gene overlap,
    color = -log10 q. Extends the heatmap with gene-count support.

heatmap_static / heatmap_plotly:
    Pathways (rows) x seed-states (columns) enrichment matrix
    (-log10 q of per-state over-representation), hierarchically ordered.
"""
from __future__ import annotations
import numpy as np
import pandas as pd


# ----------------------------------------------------------------- volcano
def volcano_static(diff: pd.DataFrame, label_A="normal", label_B="o8G",
                   q_thresh=0.05, top_n=8, title=None, ax=None):
    import matplotlib.pyplot as plt
    d = diff.copy()
    d = d[np.isfinite(d["log2_or"]) & np.isfinite(d["neglog10_q"])]
    x = d["log2_or"].to_numpy()
    y = d["neglog10_q"].to_numpy()
    sig = d["q_value"].to_numpy() < q_thresh
    toward_A = x > 0
    if ax is None:
        fig, ax = plt.subplots(figsize=(7.0, 5.2))
    else:
        fig = ax.figure
    # non-sig grey
    ax.scatter(x[~sig], y[~sig], s=14, c="#c8ccd0", alpha=0.6, linewidths=0, zorder=1)
    # sig, colored by direction
    cA, cB = "#c44e52", "#4c72b0"
    mA = sig & toward_A
    mB = sig & ~toward_A
    ax.scatter(x[mA], y[mA], s=26, c=cA, alpha=0.9, linewidths=0, zorder=3,
               label=f"toward {label_A} (lost on oxidation)")
    ax.scatter(x[mB], y[mB], s=26, c=cB, alpha=0.9, linewidths=0, zorder=3,
               label=f"toward {label_B} (gained on oxidation)")
    ax.axhline(-np.log10(q_thresh), ls="--", lw=0.9, c="#888", zorder=2)
    ax.axvline(0, ls="-", lw=0.7, c="#bbb", zorder=1)
    ax.text(ax.get_xlim()[0], -np.log10(q_thresh), f"q={q_thresh} ", va="bottom",
            ha="left", fontsize=7, color="#888")
    # label the strongest significant terms on each side, staggered to avoid overlap
    sd = d[sig].copy()
    if len(sd):
        for side, sub in (("A", sd[sd["log2_or"] > 0]), ("B", sd[sd["log2_or"] <= 0])):
            sub = sub.nlargest(top_n, "neglog10_q")
            # stagger vertically from the top down
            ymax = d["neglog10_q"].max()
            for i, (_, r) in enumerate(sub.iterrows()):
                ha = "left" if side == "A" else "right"
                dx = 6 if side == "A" else -6
                ax.annotate(_short(r["term"], 34), (r["log2_or"], r["neglog10_q"]),
                            fontsize=6.2, ha=ha, va="center",
                            xytext=(dx, -i * 11 - 4), textcoords="offset points",
                            arrowprops=dict(arrowstyle="-", lw=0.4, color="#aaa"))
    ax.set_xlabel(f"log2 odds-ratio   (\u2190 toward {label_B}   |   toward {label_A} \u2192)")
    ax.set_ylabel("-log10 BH q-value")
    ax.set_title(title or f"Differential pathway shift: {label_A} vs {label_B}",
                 fontsize=9, pad=26)
    ax.legend(frameon=False, fontsize=7, loc="lower center",
              bbox_to_anchor=(0.5, 1.005), ncol=2, handletextpad=0.3, columnspacing=1.2)
    ax.margins(0.06)
    return fig, ax


def volcano_plotly(diff: pd.DataFrame, label_A="normal", label_B="o8G", q_thresh=0.05,
                   x_col="log2_or"):
    import plotly.graph_objects as go
    d = diff.copy()
    d = d[np.isfinite(d[x_col]) & np.isfinite(d["neglog10_q"])].rename(columns={x_col: "log2_or"})
    d["sig"] = d["q_value"] < q_thresh
    d["grp"] = np.where(~d["sig"], "ns",
                np.where(d["log2_or"] > 0, f"toward {label_A}", f"toward {label_B}"))
    cmap = {"ns": "#c8ccd0", f"toward {label_A}": "#c44e52", f"toward {label_B}": "#4c72b0"}
    fig = go.Figure()
    for grp, sub in d.groupby("grp"):
        fig.add_trace(go.Scatter(
            x=sub["log2_or"], y=sub["neglog10_q"], mode="markers", name=grp,
            marker=dict(size=7, color=cmap.get(grp, "#999")),
            text=sub["term"] + "<br>lib=" + sub["library"].astype(str),
            hovertemplate="%{text}<br>log2 OR=%{x:.2f}<br>-log10 q=%{y:.2f}<extra></extra>"))
    fig.add_hline(y=-np.log10(q_thresh), line_dash="dash", line_color="#888")
    fig.add_vline(x=0, line_color="#ccc")
    fig.update_layout(template="simple_white",
                      xaxis_title=f"log2 odds-ratio (\u2190 {label_B} | {label_A} \u2192)",
                      yaxis_title="-log10 BH q", title=f"Pathway shift: {label_A} vs {label_B}",
                      legend=dict(orientation="h", y=1.06), height=520)
    return fig


# ------------------------------------------------------------- diverging bar
def diverging_bar_plotly(diff: pd.DataFrame, label_A="normal", label_B="o8G",
                         q_thresh=0.05, top_n=20, x_col="log2_or"):
    """Two-state differential enrichment as a diverging horizontal bar.

    Bar length = -log10 q; bars extend right (red) for terms enriched toward
    state A, left (blue) for state B. One named row per pathway — the readable
    alternative to a volcano when the number of significant terms is small.
    """
    import plotly.graph_objects as go
    d = diff.copy()
    d = d[np.isfinite(d[x_col]) & np.isfinite(d["neglog10_q"])].rename(columns={x_col: "log2_or"})
    sig = d[d["q_value"] < q_thresh]
    show = sig if len(sig) >= 3 else d.nlargest(min(top_n, len(d)), "neglog10_q")
    if len(show) > top_n:
        show = show.nlargest(top_n, "neglog10_q")
    show = show.copy()
    show["signed"] = np.where(show["log2_or"] > 0, show["neglog10_q"], -show["neglog10_q"])
    show = show.sort_values("signed")
    colors = np.where(show["signed"] > 0, "#c44e52", "#4c72b0")
    labels = [_short(t, 48) for t in show["term"]]
    oA = show.get(f"overlap_{label_A}", pd.Series([np.nan]*len(show)))
    oB = show.get(f"overlap_{label_B}", pd.Series([np.nan]*len(show)))
    hover = [f"{t}<br>{label_A}: {a:g} genes, {label_B}: {b:g} genes<br>-log10 q={q:.2f}"
             for t, a, b, q in zip(show["term"], oA, oB, show["neglog10_q"])]
    fig = go.Figure(go.Bar(
        x=show["signed"], y=labels, orientation="h",
        marker=dict(color=colors, line=dict(width=0.4, color="white")),
        text=hover, hovertemplate="%{text}<extra></extra>"))
    thr = -np.log10(q_thresh)
    for xv in (thr, -thr):
        fig.add_vline(x=xv, line_dash="dash", line_color="#888", line_width=1)
    fig.add_vline(x=0, line_color="#333", line_width=1)
    fig.update_layout(
        template="simple_white", height=max(320, 26*len(show)+120),
        xaxis_title=f"\u2212log10 q   (\u2190 toward {label_B}   |   toward {label_A} \u2192)",
        title=f"Differential pathway enrichment: {label_A} vs {label_B}",
        margin=dict(l=10, r=10, t=50, b=40))
    return fig


# ------------------------------------------------------------- dot matrix
def state_dotplot_data(state_gene_lists: dict[str, list[str]], enrich_fn,
                       library="GO_BP", top_terms=20):
    """Long-form table (term, state, neglog10_q, overlap) for a dot matrix
    across every seed-oxidation state. Keeps the top terms by max enrichment."""
    frames = []
    for lab, genes in state_gene_lists.items():
        e = enrich_fn(genes, library).copy()
        e["state"] = lab
        e["neglog10_q"] = -np.log10(e["q_value"].clip(lower=1e-300))
        frames.append(e[["term", "state", "neglog10_q", "overlap"]])
    long = pd.concat(frames, ignore_index=True)
    rowmax = long.groupby("term")["neglog10_q"].max()
    sig = rowmax[rowmax > 1.30]
    keep = (sig if len(sig) >= 3 else rowmax).nlargest(top_terms).index
    return long[long["term"].isin(keep)].reset_index(drop=True)


def dotplot_plotly(long_df: pd.DataFrame, states_order=None, cbar_label="-log10 q"):
    """Dot matrix: states (x) x pathways (y); dot size = gene overlap,
    color = -log10 q. Scales to all 2^k oxidation states at once."""
    import plotly.graph_objects as go
    d = long_df.copy()
    d["term_s"] = [_short(t, 48) for t in d["term"]]
    term_order = d.groupby("term_s")["neglog10_q"].max().sort_values().index.tolist()
    if states_order is None:
        states_order = list(dict.fromkeys(d["state"]))
    smax = max(d["overlap"].max(), 1)
    fig = go.Figure(go.Scatter(
        x=d["state"], y=d["term_s"], mode="markers",
        marker=dict(size=d["overlap"], sizemode="area",
                    sizeref=2.0*smax/(28.0**2), sizemin=3,
                    color=d["neglog10_q"], colorscale="Magma", showscale=True,
                    colorbar=dict(title=cbar_label), line=dict(width=0.5, color="#ddd")),
        text=[f"{t}<br>state={s}<br>{o:g} genes<br>-log10 q={q:.2f}"
              for t, s, o, q in zip(d["term"], d["state"], d["overlap"], d["neglog10_q"])],
        hovertemplate="%{text}<extra></extra>"))
    fig.update_yaxes(categoryorder="array", categoryarray=term_order)
    fig.update_xaxes(categoryorder="array", categoryarray=states_order, tickangle=45)
    fig.update_layout(template="simple_white", height=max(420, 26*len(term_order)+140),
                      xaxis_title="seed-oxidation state",
                      title="Pathway enrichment across seed-oxidation states (dot size = gene overlap)",
                      margin=dict(l=10, r=10, t=50, b=60))
    return fig


# ----------------------------------------------------------------- heatmap
def enrichment_matrix(state_gene_lists: dict[str, list[str]], enrich_fn,
                      library="GO_BP", top_terms=25, value="neglog10_q"):
    """Build pathways x states matrix of -log10 q.

    state_gene_lists: {state_label: [gene symbols]}
    enrich_fn: callable(genes, library) -> DataFrame with term,q_value
    """
    cols = {}
    for lab, genes in state_gene_lists.items():
        e = enrich_fn(genes, library)
        e = e.set_index("term")
        cols[lab] = -np.log10(e["q_value"].clip(lower=1e-300))
    mat = pd.DataFrame(cols).fillna(0.0)
    # prefer terms that reach significance (q<0.05 -> -log10 q > 1.30) in >=1 state;
    # rank those by max enrichment. Fall back to top-by-max if too few are significant.
    rowmax = mat.max(axis=1)
    sig_rows = rowmax[rowmax > 1.30]
    if len(sig_rows) >= 3:
        keep = sig_rows.nlargest(top_terms).index
    else:
        keep = rowmax.nlargest(top_terms).index
    return mat.loc[keep]


def heatmap_static(mat: pd.DataFrame, title=None, cbar_label="-log10 q", cluster=True):
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import linkage, leaves_list
    m = mat.copy()
    if cluster and m.shape[0] > 2:
        try:
            order = leaves_list(linkage(m.values, method="average"))
            m = m.iloc[order]
        except Exception:
            pass
    fig, ax = plt.subplots(figsize=(max(5.0, 0.7*m.shape[1]+3.5), max(4.0, 0.26*m.shape[0]+1.2)))
    im = ax.imshow(m.values, aspect="auto", cmap="magma")
    ax.set_xticks(range(m.shape[1])); ax.set_xticklabels(m.columns, rotation=45, ha="right", fontsize=7.5)
    ax.set_yticks(range(m.shape[0])); ax.set_yticklabels([_short(t, 46) for t in m.index], fontsize=6.8)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02); cb.set_label(cbar_label, fontsize=8)
    ax.set_title(title or "Pathway enrichment across seed-oxidation states", fontsize=9)
    ax.set_xlabel("seed-oxidation state", fontsize=8)
    return fig, ax


def heatmap_plotly(mat: pd.DataFrame, cbar_label="-log10 q"):
    import plotly.graph_objects as go
    fig = go.Figure(go.Heatmap(
        z=mat.values, x=list(mat.columns), y=[_short(t, 46) for t in mat.index],
        colorscale="Magma", colorbar=dict(title=cbar_label)))
    fig.update_layout(template="simple_white", height=max(420, 22*mat.shape[0]+120),
                      xaxis_title="seed-oxidation state",
                      title="Pathway enrichment across seed-oxidation states")
    return fig


def _short(s, n=40):
    s = str(s)
    # strip trailing GO/Reactome ids
    import re
    s = re.sub(r"\s*\((GO:\d+)\)$", "", s)
    s = re.sub(r"\s+R-HSA-\d+$", "", s)
    return s if len(s) <= n else s[:n-1] + "\u2026"


def upset_plotly(sets: dict[str, set], max_intersections: int = 40):
    """UpSet-style plot of gene-set intersections (plotly; no extra dependency).

    Top panel: intersection sizes (sorted desc).
    Bottom panel: which sets participate (dot matrix).
    """
    import itertools
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    names = list(sets.keys())
    if len(names) < 2:
        fig = go.Figure()
        fig.update_layout(title="Need ≥2 databases for an UpSet plot")
        return fig

    # exclusive intersections for every non-empty subset
    rows = []
    for r in range(1, len(names) + 1):
        for combo in itertools.combinations(names, r):
            combo = list(combo)
            inter = set.intersection(*(sets[n] for n in combo))
            # exclusive: remove genes also in databases outside the combo
            outside = [n for n in names if n not in combo]
            if outside:
                union_out = set.union(*(sets[n] for n in outside)) if outside else set()
                inter = inter - union_out
            if not inter:
                continue
            rows.append({"combo": combo, "size": len(inter), "degree": len(combo)})

    if not rows:
        fig = go.Figure()
        fig.update_layout(title="No non-empty intersections")
        return fig

    rows.sort(key=lambda x: (-x["size"], -x["degree"]))
    rows = rows[:max_intersections]
    labels = [" ∩ ".join(r["combo"]) for r in rows]
    sizes = [r["size"] for r in rows]

    fig = make_subplots(
        rows=2,
        cols=1,
        row_heights=[0.58, 0.42],
        shared_xaxes=True,
        vertical_spacing=0.10,
    )
    fig.add_trace(
        go.Bar(x=list(range(len(rows))), y=sizes, marker_color="#4C72B0", name="size",
               text=sizes, textposition="outside"),
        row=1, col=1,
    )
    # membership matrix as scatter
    xs, ys, colors = [], [], []
    for i, r in enumerate(rows):
        for j, name in enumerate(names):
            xs.append(i)
            ys.append(j)
            colors.append("#4C72B0" if name in r["combo"] else "#E8E8E8")
    fig.add_trace(
        go.Scatter(
            x=xs, y=ys, mode="markers",
            marker=dict(
                size=11,
                color=colors,
                line=dict(width=1, color="#888888"),
            ),
            showlegend=False,
            hoverinfo="text",
            text=[names[j] for j in ys],
        ),
        row=2, col=1,
    )
    fig.update_xaxes(showticklabels=False, row=1, col=1)
    fig.update_xaxes(
        tickmode="array",
        tickvals=list(range(len(rows))),
        ticktext=[_short(lb, 24) for lb in labels],
        tickangle=-50,
        row=2, col=1,
    )
    fig.update_yaxes(title_text="n genes", row=1, col=1)
    # Drop y tick labels on the membership panel — they collide with the dots;
    # set identity is in the hover text on each marker.
    fig.update_yaxes(
        tickmode="array",
        tickvals=list(range(len(names))),
        showticklabels=False,
        title_text="",
        row=2,
        col=1,
    )
    fig.update_layout(
        template="simple_white",
        height=400 + 22 * len(names),
        margin=dict(l=60, r=20, t=50, b=130),
        title="UpSet — exclusive intersections (hover dots for set names)",
    )
    return fig


def venn2_plotly(set_a, set_b, label_a: str = "A", label_b: str = "B"):
    """Two-set Venn with exclusive / shared counts (plotly shapes; no extra dep)."""
    import plotly.graph_objects as go

    a, b = set(set_a), set(set_b)
    only_a, only_b, both = len(a - b), len(b - a), len(a & b)
    fig = go.Figure()
    fig.add_shape(
        type="circle", xref="x", yref="y",
        x0=0, y0=0, x1=2, y1=2,
        line_color="#c44e52", fillcolor="rgba(196,78,82,0.25)",
    )
    fig.add_shape(
        type="circle", xref="x", yref="y",
        x0=1.2, y0=0, x1=3.2, y1=2,
        line_color="#4c72b0", fillcolor="rgba(76,114,176,0.25)",
    )
    fig.add_annotation(x=0.7, y=1.0, text=f"<b>{only_a:,}</b>", showarrow=False, font=dict(size=18))
    fig.add_annotation(x=2.5, y=1.0, text=f"<b>{only_b:,}</b>", showarrow=False, font=dict(size=18))
    fig.add_annotation(x=1.6, y=1.0, text=f"<b>{both:,}</b>", showarrow=False, font=dict(size=18))
    fig.add_annotation(x=0.7, y=2.25, text=_short(label_a, 24), showarrow=False, font=dict(size=13, color="#c44e52"))
    fig.add_annotation(x=2.5, y=2.25, text=_short(label_b, 24), showarrow=False, font=dict(size=13, color="#4c72b0"))
    fig.add_annotation(x=0.7, y=0.35, text="only A", showarrow=False, font=dict(size=10, color="#666"))
    fig.add_annotation(x=2.5, y=0.35, text="only B", showarrow=False, font=dict(size=10, color="#666"))
    fig.add_annotation(x=1.6, y=0.35, text="shared", showarrow=False, font=dict(size=10, color="#666"))
    fig.update_xaxes(visible=False, range=[-0.2, 3.4])
    fig.update_yaxes(visible=False, range=[-0.2, 2.6], scaleanchor="x", scaleratio=1)
    fig.update_layout(
        template="simple_white",
        height=340,
        margin=dict(l=20, r=20, t=40, b=20),
        title=f"|A|={len(a):,} · |B|={len(b):,} · shared={both:,}",
    )
    return fig


def venn3_plotly(set_a, set_b, set_c, label_a="A", label_b="B", label_c="C"):
    """Three-set Venn with region counts (plotly shapes; schematic, not area-proportional)."""
    import plotly.graph_objects as go

    a, b, c = set(set_a), set(set_b), set(set_c)
    only_a = len(a - b - c)
    only_b = len(b - a - c)
    only_c = len(c - a - b)
    ab = len((a & b) - c)
    ac = len((a & c) - b)
    bc = len((b & c) - a)
    abc = len(a & b & c)

    fig = go.Figure()
    # equilateral-ish triangle of circles
    circles = [
        (0.0, 0.55, "#c44e52", label_a),
        (1.1, 0.55, "#4c72b0", label_b),
        (0.55, -0.35, "#55a868", label_c),
    ]
    r = 1.15
    for x0, y0, color, _lab in circles:
        fig.add_shape(
            type="circle",
            xref="x",
            yref="y",
            x0=x0,
            y0=y0,
            x1=x0 + 2 * r,
            y1=y0 + 2 * r,
            line_color=color,
            fillcolor=color.replace(")", ",0.22)").replace("rgb", "rgba")
            if color.startswith("rgb")
            else (
                "rgba(196,78,82,0.22)"
                if color == "#c44e52"
                else "rgba(76,114,176,0.22)"
                if color == "#4c72b0"
                else "rgba(85,168,104,0.22)"
            ),
        )
    # region count annotations (approximate centers)
    anns = [
        (0.55, 2.35, only_a, "#c44e52"),
        (2.75, 2.35, only_b, "#4c72b0"),
        (1.65, -0.05, only_c, "#55a868"),
        (1.65, 2.05, ab, "#333"),
        (0.95, 1.15, ac, "#333"),
        (2.35, 1.15, bc, "#333"),
        (1.65, 1.45, abc, "#111"),
    ]
    for x, y, n, col in anns:
        fig.add_annotation(
            x=x, y=y, text=f"<b>{n:,}</b>", showarrow=False, font=dict(size=14, color=col)
        )
    fig.add_annotation(x=0.4, y=3.0, text=_short(label_a, 18), showarrow=False, font=dict(size=12, color="#c44e52"))
    fig.add_annotation(x=2.9, y=3.0, text=_short(label_b, 18), showarrow=False, font=dict(size=12, color="#4c72b0"))
    fig.add_annotation(x=1.65, y=-0.55, text=_short(label_c, 18), showarrow=False, font=dict(size=12, color="#55a868"))
    fig.update_xaxes(visible=False, range=[-0.3, 3.6])
    fig.update_yaxes(visible=False, range=[-0.9, 3.3], scaleanchor="x", scaleratio=1)
    fig.update_layout(
        template="simple_white",
        height=420,
        margin=dict(l=20, r=20, t=40, b=20),
        title=(
            f"|A|={len(a):,} · |B|={len(b):,} · |C|={len(c):,} · "
            f"A∩B∩C={abc:,}"
        ),
    )
    return fig
