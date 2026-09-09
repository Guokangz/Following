#!/usr/bin/env python3
"""Generate all current linear-paper reproductions (Figs. 3-6)."""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.linear import reproduce_fig3, reproduce_fig4, reproduce_fig5, reproduce_fig6


def main() -> None:
    for label, workflow in (
        ("Fig. 3", reproduce_fig3.main),
        ("Fig. 4", reproduce_fig4.main),
        ("Fig. 5", reproduce_fig5.main),
        ("Fig. 6", reproduce_fig6.main),
    ):
        print(f"\n=== {label} ===")
        workflow()


if __name__ == "__main__":
    main()

