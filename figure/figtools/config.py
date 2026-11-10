"""Plot style and constants shared by the figure tools."""

from __future__ import annotations

from typing import Iterable

from matplotlib import font_manager, rcParams

FIG_DPI = 300
PANEL_SIZE = (2.8, 3.0)
COLORS = {
    "Z": "#0072B2",
    "N": "#D55E00",
    "A": "#009E73",
}
LINE_WIDTH = 1.1
MARKER_SIZE = 2.5
AXIS_LABEL_SIZE = 10
TICK_LABEL_SIZE = 9
PPT_AXIS_LABEL_SIZE = 22
PPT_TICK_LABEL_SIZE = 22
PPT_LEGEND_SIZE = 18
PPT_TITLE_SIZE = 22
PPT_LINE_WIDTH = 1.8
ENERGY_COLORS = (
    "#CFE3F5",
    "#A5C8EB",
    "#78ABDE",
    "#4B88C6",
    "#25639E",
    "#123F6B",
)


def _installed_font_names(fonts: Iterable[object]) -> set[str]:
    """Return normalized names from Matplotlib's installed-font records."""

    names: set[str] = set()
    for font in fonts:
        name = getattr(font, "name", None)
        if name:
            names.add(str(name).casefold())
    return names


def detect_font() -> str:
    """Select the first available font from the project's sans-serif fallback chain."""

    installed = _installed_font_names(font_manager.fontManager.ttflist)
    for candidate in ("Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"):
        if candidate.casefold() in installed:
            return candidate
    # Matplotlib ships DejaVu Sans in normal installations.  Returning it keeps
    # rcParams deterministic even in unusually minimal environments.
    return "DejaVu Sans"


def apply_style() -> str:
    """Apply the shared style and return the font selected on this machine."""

    actual_font = detect_font()
    rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [actual_font, "DejaVu Sans"],
            "axes.linewidth": 0.8,
            "axes.grid": True,
            "grid.color": "#B0B0B0",
            "grid.linewidth": 0.5,
            "grid.alpha": 0.45,
            "xtick.labelsize": TICK_LABEL_SIZE,
            "ytick.labelsize": TICK_LABEL_SIZE,
            "axes.labelsize": AXIS_LABEL_SIZE,
            "figure.dpi": FIG_DPI,
        }
    )
    return actual_font
