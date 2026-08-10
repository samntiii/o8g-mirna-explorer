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


def render_overlap(ctx: SectionContext) -> None:
    from o8g_venn import render as _render

    _render(ctx)


def render_tf(ctx: SectionContext) -> None:
    from o8g_tf import render as _render

    _render(ctx)


def render_lof(ctx: SectionContext) -> None:
    from o8g_lof import render as _render

    _render(ctx)


def render_stats(ctx: SectionContext) -> None:
    from o8g_stats import render as _render

    _render(ctx)


def render_antagomir(ctx: SectionContext) -> None:
    from o8g_anti import render as _render

    _render(ctx)
