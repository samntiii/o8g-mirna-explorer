#!/usr/bin/env python3
"""Stub: benchmark_refsets driver.

The real driver lives at ``paper/scripts/benchmark_refsets.py`` (gitignored
``paper/`` tree used for the NAR manuscript). Reusable external-catalog loaders
are in ``o8g_refsets``.

This path previously was a dangling symlink into ``../paper/scripts/``, which
broke fresh clones, ``ls -L``, tarball builds, and CI checkouts.
"""
from __future__ import annotations

import sys

if __name__ == "__main__":
    print(
        "benchmark_refsets: driver is paper/scripts/benchmark_refsets.py "
        "(not shipped in this clone). Use o8g_refsets for loaders.",
        file=sys.stderr,
    )
    sys.exit(2)
