#!/usr/bin/env python3
"""Generate all reproduced numerical and conceptual nonlinear figures."""

from reproduce_conceptual_figures import main as reproduce_conceptual_figures
from reproduce_fig3 import main as reproduce_fig3
from reproduce_fig4 import main as reproduce_fig4


def main() -> None:
    reproduce_fig3()
    reproduce_fig4()
    reproduce_conceptual_figures()


if __name__ == "__main__":
    main()

