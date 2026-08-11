"""Section router helpers for optional explorer views.

Design constraint: sections never query the DB directly. Every target set
arrives via ``ctx.strong_set(state_label)``, which routes through the
precision ladder in ``app.filtered_targets``.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import streamlit as st


@dataclass
class SectionContext:
    db: object
    mirna: str
    info: dict
    seed: str
    gpos: list
    state_labels: list[str]
    library: str
    precision_mode: str
    strong_set: Callable[[str], set[str]]
    universe: set[str]
    matched_background: set[str]
    external_refs: object
    precision_cfg: object | None = None
    scanner: object | None = None


def render_overlap(ctx: SectionContext) -> None:
    import importlib
    import o8g_venn as _mod

    _mod = importlib.reload(_mod)
    _mod.render(ctx)


def render_tf(ctx: SectionContext) -> None:
    import importlib
    import o8g_tf as _mod

    _mod = importlib.reload(_mod)
    _mod.render(ctx)


def render_lof(ctx: SectionContext) -> None:
    from o8g_lof import render as _render

    _render(ctx)


def render_stats(ctx: SectionContext) -> None:
    from o8g_stats import render as _render

    _render(ctx)


def render_antagomir(ctx: SectionContext) -> None:
    import importlib
    import o8g_energy as _energy
    import o8g_anti as _mod

    importlib.reload(_energy)
    _mod = importlib.reload(_mod)
    _mod.render(ctx)


def render_deg_upload(ctx: SectionContext) -> None:
    from o8g_deg_upload import render as _render

    _render(ctx)
