"""Matplotlib Figure-based plots for spallation results."""

from __future__ import annotations

from os import PathLike
from pathlib import Path
from typing import Any, Sequence

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator

from .config import (
    AXIS_LABEL_SIZE,
    COLORS,
    LINE_WIDTH,
    MARKER_SIZE,
    PANEL_SIZE,
    ENERGY_COLORS,
    PPT_AXIS_LABEL_SIZE,
    PPT_LEGEND_SIZE,
    PPT_LINE_WIDTH,
    PPT_TICK_LABEL_SIZE,
    PPT_TITLE_SIZE,
    apply_style,
)


def plot_cross_section_distributions(
    result: Any,
    output: str | PathLike[str] | None = None,
    quantities: Sequence[str] = ("Z", "N", "A"),
    yscale: str = "linear",
    title: str | None = None,
    dpi: int = 300,
) -> Figure:
    """Create a one-row cross-section distribution figure and optionally save it."""

    if yscale not in {"linear", "log"}:
        raise ValueError("yscale must be 'linear' or 'log'")
    normalized_quantities = tuple(str(quantity).upper() for quantity in quantities)
    if not normalized_quantities:
        raise ValueError("quantities must contain at least one quantity")
    invalid = set(normalized_quantities) - set(COLORS)
    if invalid:
        raise ValueError(f"unsupported quantities: {sorted(invalid)}")

    apply_style()
    fig = Figure(figsize=(PANEL_SIZE[0] * len(normalized_quantities), PANEL_SIZE[1]))
    FigureCanvasAgg(fig)
    axes = fig.subplots(1, len(normalized_quantities), squeeze=False)[0]

    for ax, quantity in zip(axes, normalized_quantities):
        block = result.get(quantity) if hasattr(result, "get") else None
        if block is None:
            ax.set_visible(False)
            continue

        x, y = block
        ax.plot(
            x,
            y,
            color=COLORS[quantity],
            linewidth=LINE_WIDTH,
            marker="o",
            markersize=MARKER_SIZE,
        )
        ax.set_xlabel(quantity, fontsize=AXIS_LABEL_SIZE)
        ax.set_yscale(yscale)
        ax.grid(True, color="#B0B0B0", linewidth=0.5, alpha=0.45)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.supylabel("Cross Section (mb)", fontsize=AXIS_LABEL_SIZE)
    if title:
        fig.suptitle(title)
    fig.tight_layout(rect=(0.04, 0.0, 1.0, 0.94 if title else 1.0))

    if output is not None:
        output_path = Path(output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    return fig


def plot_energy_comparison(
    series: Sequence[tuple[str, Any]],
    quantity: str,
    output: str | PathLike[str] | None = None,
    title: str | None = None,
    yscale: str = "log",
    dpi: int = 300,
) -> Figure:
    """Plot one quantity across an energy series in a PPT-friendly Figure."""

    quantity_name = str(quantity).upper()
    if quantity_name not in {"Z", "N", "A"}:
        raise ValueError("quantity must be one of 'Z', 'N', or 'A'")
    if yscale not in {"linear", "log"}:
        raise ValueError("yscale must be 'linear' or 'log'")

    apply_style()
    entries = list(series)
    fig = Figure(figsize=(7.5, 6.0))
    FigureCanvasAgg(fig)
    ax = fig.subplots()
    plotted = 0

    for index, (label, result) in enumerate(entries):
        block = result.get(quantity_name) if hasattr(result, "get") else None
        if block is None:
            continue
        x, y = block
        ax.plot(
            x,
            y,
            color=ENERGY_COLORS[index] if index < len(ENERGY_COLORS) else ENERGY_COLORS[-1],
            linewidth=PPT_LINE_WIDTH,
            drawstyle="steps-mid",
            label=str(label),
        )
        plotted += 1

    ax.set_xlabel(quantity_name, fontsize=PPT_AXIS_LABEL_SIZE, fontweight="bold")
    ax.set_ylabel("Cross Section (mb)", fontsize=PPT_AXIS_LABEL_SIZE, fontweight="bold")
    ax.tick_params(
        axis="both",
        labelsize=PPT_TICK_LABEL_SIZE,
        width=1.8,
        length=5.0,
        top=False,
        right=False,
        labeltop=False,
        labelright=False,
    )
    ax.set_yscale(yscale)
    ax.grid(False)
    ax.set_xlim(left=0)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(
        axis="both",
        which="minor",
        width=1.8,
        length=3.0,
        top=False,
        right=False,
        labeltop=False,
        labelright=False,
    )
    for spine in ax.spines.values():
        spine.set_linewidth(1.8)
    for tick_label in ax.get_xticklabels() + ax.get_yticklabels():
        tick_label.set_fontweight("bold")
    if plotted:
        ax.legend(
            loc="best",
            frameon=False,
            prop={"size": PPT_LEGEND_SIZE, "weight": "bold"},
        )
    if title:
        ax.set_title(title, fontsize=PPT_TITLE_SIZE, fontweight="bold")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94 if title else 1.0))

    if output is not None:
        output_path = Path(output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    return fig
