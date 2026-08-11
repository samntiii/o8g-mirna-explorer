"""
o8g_precision.py
================
Precision ladder for oxomiR target lists.

CORE CONSTRAINT
---------------
Every filter is applied independently to the unmodified and oxidized target
tables; gained / lost / shared are set differences *after* filtering.
Never filter the delta itself.

Modes (PrecisionMode)
---------------------
Sensitive  — rank >= 3 (7mer-m8 + 8mer). Discovery default in the UI.
Stringent  — rank == 4 (8mer only). Default for paper / Fig.4 exports.
Consensus  — rank >= 3 AND TargetScan-conserved on the *unmodified* baseline;
             oxidized state uses the same rank gate (conservation is undefined
             for novel o8G motifs in TargetScan — see note below). Deltas =
             setdiff after filters. Use for main-text quantitative claims.
TargetScan — rank >= 3 AND TargetScan *predicted* strong sites (8mer+7mer-m8
             from Predicted_Targets_Info) on the unmodified baseline. Same
             partition logic as Consensus, but does **not** require conserved-
             family tables — use when you want catalog predictions alone.

Deprecated alias
----------------
"HighConf" previously meant SITE_WEIGHT score ≥ 1.0, but without stored
multiplicity blobs that gate collapsed to 8mer-only (= Stringent) and recovered
the fewest gold effects. It is merged into Stringent; ``from_mode("HighConf")``
still resolves for backward compatibility.

Conservation note
-----------------
TargetScan conserved tables annotate *canonical* (unoxidized) seed matches.
Applying is_conserved==True identically to oxidized motifs would zero oxomiR
gains (those complementary words are absent from TargetScan). Consensus
therefore applies conservation to the unmodified baseline set and the same
*rank* filter to both states.

Feature flags (env or constructors)
-----------------------------------
O8G_USE_CONSERVATION=1|0
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from enum import Enum

import pandas as pd


class PrecisionMode(str, Enum):
    SENSITIVE = "Sensitive"
    STRINGENT = "Stringent"
    TARGETSCAN = "TargetScan"
    CONSENSUS = "Consensus"


# Old UI / CSV label → current mode (HighConf was a Stringent duplicate).
_MODE_ALIASES = {
    "HighConf": PrecisionMode.STRINGENT,
    "highconf": PrecisionMode.STRINGENT,
}


def _coerce_mode(mode: PrecisionMode | str | object) -> PrecisionMode:
    """Normalize to PrecisionMode; survive Streamlit module reloads (enum identity)."""
    if isinstance(mode, PrecisionMode):
        return mode
    # Enum-like from a reloaded module (same name/value, different class object)
    if hasattr(mode, "value") and not hasattr(mode, "mode"):
        try:
            return PrecisionMode(str(getattr(mode, "value")))
        except Exception:
            pass
    if isinstance(mode, str):
        mode = _MODE_ALIASES.get(mode, mode)
        return mode if isinstance(mode, PrecisionMode) else PrecisionMode(mode)
    # Nested / mistaken PrecisionConfig passed as mode
    if hasattr(mode, "mode"):
        return _coerce_mode(getattr(mode, "mode"))
    raise TypeError(f"Cannot coerce precision mode from {type(mode)!r}: {mode!r}")


@dataclass(frozen=True)
class PrecisionConfig:
    mode: PrecisionMode = PrecisionMode.SENSITIVE
    use_conservation: bool = True

    @staticmethod
    def from_mode(mode: PrecisionMode | str | "PrecisionConfig", **kwargs) -> "PrecisionConfig":
        # Streamlit hot-reload can produce a *different* PrecisionConfig class object;
        # treat duck-typed configs as already built.
        if isinstance(mode, PrecisionConfig) or (
            type(mode).__name__ == "PrecisionConfig" and hasattr(mode, "mode")
        ):
            if kwargs:
                return PrecisionConfig(
                    mode=_coerce_mode(getattr(mode, "mode")),
                    use_conservation=kwargs.get(
                        "use_conservation", getattr(mode, "use_conservation", True)
                    ),
                )
            return PrecisionConfig(
                mode=_coerce_mode(getattr(mode, "mode")),
                use_conservation=bool(getattr(mode, "use_conservation", True)),
            )
        return PrecisionConfig(mode=_coerce_mode(mode), **kwargs)

    @staticmethod
    def ui_default() -> "PrecisionConfig":
        return PrecisionConfig(mode=PrecisionMode.SENSITIVE)

    @staticmethod
    def paper_default() -> "PrecisionConfig":
        return PrecisionConfig(mode=PrecisionMode.STRINGENT)

    @staticmethod
    def claim_default() -> "PrecisionConfig":
        return PrecisionConfig(mode=PrecisionMode.CONSENSUS)


def env_flag(name: str, default: bool = True) -> bool:
    v = os.environ.get(name)
    if v is None:
        return default
    return v not in ("0", "false", "False", "no")


def apply_precision_filter(
    df: pd.DataFrame,
    cfg: PrecisionConfig,
    *,
    conserved_symbols: set[str] | None = None,
    is_unmodified_state: bool = True,
) -> pd.DataFrame:
    """Filter one state's target table. Does not touch deltas.

    Expected columns (subset OK): symbol, site_rank, is_conserved.
    """
    if df is None or df.empty:
        return df.iloc[0:0].copy() if df is not None else pd.DataFrame()

    cfg = PrecisionConfig.from_mode(cfg)
    out = df
    use_cons = cfg.use_conservation and env_flag("O8G_USE_CONSERVATION", True)

    # Compare by *value* string — never by enum identity (Streamlit reload-safe).
    mode_s = str(_coerce_mode(cfg.mode).value)
    if mode_s == PrecisionMode.SENSITIVE.value:
        out = out[out["site_rank"] >= 3]
    elif mode_s == PrecisionMode.STRINGENT.value:
        out = out[out["site_rank"] >= 4]
    elif mode_s in (PrecisionMode.CONSENSUS.value, PrecisionMode.TARGETSCAN.value):
        out = out[out["site_rank"] >= 3]
        if use_cons and is_unmodified_state:
            if mode_s == PrecisionMode.CONSENSUS.value and "is_conserved" in out.columns:
                out = out[out["is_conserved"].astype(bool)]
            elif conserved_symbols is not None:
                # TargetScan predictions, or Consensus when is_conserved column absent
                syms = {str(s).upper() for s in conserved_symbols}
                out = out[out["symbol"].astype(str).str.upper().isin(syms)]
        # oxidized: rank gate only (TargetScan has no o8G motifs)
    else:
        raise ValueError(f"Unknown precision mode: {mode_s!r}")

    return out.reset_index(drop=True)


def partition_after_filter(
    unmod: pd.DataFrame,
    oxid: pd.DataFrame,
    cfg: PrecisionConfig,
    *,
    conserved_symbols: set[str] | None = None,
) -> dict[str, set[str]]:
    """Apply filters to each state; deltas = setdiff after filtering.

    Consensus special-case (documented): conservation shrinks the *unmodified
    baseline* for claims, but **gained** uses the same rank gate on both states
    without conservation on the oxidized side *and* subtracts the full
    rank-filtered unmodified set — so gains remain true o8G retargeting events,
    not artifacts of dropping non-conserved unmodified targets.
    """
    mode_s = str(getattr(cfg.mode, "value", cfg.mode))
    if mode_s not in (PrecisionMode.CONSENSUS.value, PrecisionMode.TARGETSCAN.value):
        u = apply_precision_filter(
            unmod, cfg, conserved_symbols=conserved_symbols, is_unmodified_state=True
        )
        o = apply_precision_filter(
            oxid, cfg, conserved_symbols=conserved_symbols, is_unmodified_state=False
        )
        su = set(u["symbol"]) if len(u) else set()
        so = set(o["symbol"]) if len(o) else set()
        return {
            "unmod": su,
            "oxid": so,
            "shared": su & so,
            "lost": su - so,
            "gained": so - su,
        }

    # Consensus / TargetScan: anchor shrinks unmodified baseline; gains use rank-only
    u_rank = apply_precision_filter(
        unmod,
        PrecisionConfig(mode=PrecisionMode.SENSITIVE),
        is_unmodified_state=True,
    )
    o_rank = apply_precision_filter(
        oxid,
        PrecisionConfig(mode=PrecisionMode.SENSITIVE),
        is_unmodified_state=False,
    )
    su_rank = set(u_rank["symbol"]) if len(u_rank) else set()
    so_rank = set(o_rank["symbol"]) if len(o_rank) else set()
    u_anch = apply_precision_filter(
        unmod, cfg, conserved_symbols=conserved_symbols, is_unmodified_state=True
    )
    su_anch = set(u_anch["symbol"]) if len(u_anch) else set()
    return {
        "unmod": su_anch,
        "oxid": so_rank,
        "shared": su_anch & so_rank,
        "lost": su_anch - so_rank,
        "gained": so_rank - su_rank,
    }


def assert_retargeting_signal(
    parts: dict[str, set[str]],
    *,
    min_gained: int = 50,
    min_lost: int = 50,
    min_shared: int = 50,
    label: str = "miR-1 o8G@7",
) -> None:
    """Regression: oxidation still yields substantial partitions after filtering."""
    g, l, s = len(parts["gained"]), len(parts["lost"]), len(parts["shared"])
    if g < min_gained or l < min_lost or s < min_shared:
        raise AssertionError(
            f"{label} retargeting collapsed under filter: "
            f"gained={g}, lost={l}, shared={s} "
            f"(need ≥{min_gained}/{min_lost}/{min_shared})"
        )
