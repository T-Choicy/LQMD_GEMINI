#!/usr/bin/env python3
"""Command-line entry point for A/Z/N cross-section distributions."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


# Allow direct execution from the repository root without requiring installation.
FIGURE_ROOT = Path(__file__).resolve().parents[1]
if str(FIGURE_ROOT) not in sys.path:
    sys.path.insert(0, str(FIGURE_ROOT))

from figtools.config import apply_style  # noqa: E402
from figtools.io import derive_title, load_result  # noqa: E402
from figtools.plots import plot_cross_section_distributions  # noqa: E402


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="multi-block .dat input file")
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        help="output image path (default: figure/output/<stem>_cross_section.png)",
    )
    parser.add_argument("--yscale", choices=("linear", "log"), default="linear")
    parser.add_argument("--quantities", nargs="+", choices=("Z", "N", "A"), default=["Z", "N", "A"])
    parser.add_argument("--title", default=None, help="figure title; defaults to a title derived from input")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    output = args.output or Path("figure") / "output" / f"{args.input.stem}_cross_section.png"
    font_name = apply_style()
    result = load_result(args.input)
    plot_cross_section_distributions(
        result,
        output=output,
        quantities=args.quantities,
        yscale=args.yscale,
        title=args.title if args.title is not None else derive_title(args.input),
    )
    print(f"Saved {output}")
    print(f"Font: {font_name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
