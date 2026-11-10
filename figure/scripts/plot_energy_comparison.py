#!/usr/bin/env python3
"""Create independent PPT-ready Z/N/A energy comparison plots."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


# Allow direct execution from the repository root without requiring installation.
FIGURE_ROOT = Path(__file__).resolve().parents[1]
if str(FIGURE_ROOT) not in sys.path:
    sys.path.insert(0, str(FIGURE_ROOT))

from figtools.config import apply_style  # noqa: E402
from figtools.io import discover_energy_series  # noqa: E402
from figtools.plots import plot_energy_comparison  # noqa: E402


COMPARISON_ENERGIES = (80, 100, 200, 300, 600, 1000)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("data/p+232Th"))
    parser.add_argument("--mode", choices=("COLD", "HOT"), default="COLD")
    parser.add_argument("--quantities", nargs="+", choices=("Z", "N", "A"), default=["Z", "N", "A"])
    parser.add_argument("--output-dir", type=Path, default=Path("figure/output"))
    parser.add_argument("--yscale", choices=("linear", "log"), default="log")
    parser.add_argument("--title", default=None, help="optional title applied to each figure")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    discovered = discover_energy_series(args.root, mode=args.mode)
    # The confirmed PPT comparison consists of these six energies.  The
    # scanner remains generic, so extra directories (for example 800MeV in
    # the repository) do not silently change the figure specification.
    series = [
        (label, result)
        for label, result in discovered
        if int(label.split()[0]) in COMPARISON_ENERGIES
    ]
    if not series:
        raise SystemExit(f"No CNPCS_{args.mode}.dat energy files found below {args.root}")

    font_name = apply_style()
    reaction = args.root.name or "energy_series"
    for quantity in args.quantities:
        output = args.output_dir / f"{reaction}_{quantity}_energy.png"
        plot_energy_comparison(
            series,
            quantity,
            output=output,
            title=args.title,
            yscale=args.yscale,
        )
        print(f"Saved {output}")
    print(f"Font: {font_name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
