"""Reusable plotting utilities for spallation analysis results."""

from . import config, io, plots
from .config import apply_style
from .io import Result, derive_title, discover_energy_series, load_result, parse_multiblock
from .plots import plot_cross_section_distributions, plot_energy_comparison

__version__ = "0.1.0"

__all__ = [
    "Result",
    "apply_style",
    "config",
    "derive_title",
    "discover_energy_series",
    "io",
    "load_result",
    "parse_multiblock",
    "plot_cross_section_distributions",
    "plot_energy_comparison",
    "plots",
]
