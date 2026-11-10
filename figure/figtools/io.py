"""Readers for the multi-block ASCII output produced by anahyp."""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path
import re
import numpy as np


ArrayPair = tuple[np.ndarray, np.ndarray]


@dataclass
class Result(Mapping[str, ArrayPair]):
    """Parsed result blocks keyed by normalized quantity name.

    Each value is an ``(x, y)`` NumPy-array pair.  One-column blocks use the
    zero-based row number as ``x``.  ``Result`` also behaves like a read-only
    mapping for convenient access such as ``result["Z"]``.
    """

    blocks: dict[str, ArrayPair]
    source: Path | None = None

    def __getitem__(self, key: str) -> ArrayPair:
        return self.blocks[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self.blocks)

    def __len__(self) -> int:
        return len(self.blocks)


def _numeric_tokens(line: str) -> list[float] | None:
    """Parse a line as numeric data, preserving scientific notation support."""

    tokens = line.split()
    if not tokens:
        return None
    try:
        return [float(token) for token in tokens]
    except ValueError:
        return None


def _block_key(header: str) -> str:
    """Map a potentially multi-line header to the public normalized key."""

    normalized = " ".join(header.split()).casefold()
    if "special" in normalized:
        return "special"
    if "excitation" in normalized:
        if "mass" in normalized:
            return "excitation_mass"
        return "excitation_energy"
    for quantity in ("z", "n", "a"):
        if re.search(rf"function\s+of\s+{quantity}(?:\b|\s|:)", normalized):
            return quantity.upper()
    return normalized or "unknown"


def parse_multiblock(source: str | Path | Iterator[str]) -> Result:
    """Parse a multi-block result file or text iterator into a :class:`Result`.

    A line is data only when every whitespace-separated token converts with
    ``float()``.  This is intentionally the discriminator so values such as
    ``2.1382461E-03`` are never mistaken for headers.
    """

    source_path: Path | None = None
    if isinstance(source, Path):
        source_path = source
        lines = source.read_text(encoding="utf-8").splitlines()
    elif isinstance(source, str):
        possible_path = Path(source)
        if "\n" not in source and "\r" not in source and possible_path.is_file():
            source_path = possible_path
            lines = possible_path.read_text(encoding="utf-8").splitlines()
        else:
            lines = source.splitlines()
    else:
        lines = source

    blocks: dict[str, list[tuple[float, float]]] = {}
    current_header: list[str] = []
    current_rows: list[tuple[float, float]] = []

    def flush() -> None:
        nonlocal current_header, current_rows
        if not current_header:
            current_rows = []
            return
        key = _block_key("\n".join(current_header))
        if current_rows:
            blocks.setdefault(key, []).extend(current_rows)
        current_header = []
        current_rows = []

    for raw_line in lines:
        line = raw_line.strip()
        if not line:
            # Blank lines separate blocks only after data has started.  A
            # title followed by a column header is intentionally kept joined.
            if current_rows:
                flush()
            continue

        values = _numeric_tokens(line)
        if values is not None and len(values) in (1, 2):
            if not current_header:
                # Be forgiving of a data-only fragment while keeping the
                # mapping useful to callers.
                current_header = ["unknown"]
            if len(values) == 2:
                current_rows.append((values[0], values[1]))
            else:
                current_rows.append((float(len(current_rows)), values[0]))
            continue

        # A header after a header belongs to the same logical block (for
        # example, "EXCITATION ..." followed by "MeV mb/MeV").
        if current_rows:
            flush()
        current_header.append(line)

    flush()

    parsed: dict[str, ArrayPair] = {}
    for key, rows in blocks.items():
        data = np.asarray(rows, dtype=float)
        parsed[key] = (data[:, 0], data[:, 1])
    return Result(parsed, source=source_path)


def load_result(path: str | Path) -> Result:
    """Load and parse one anahyp ``.dat`` file."""

    return parse_multiblock(Path(path))


def discover_energy_series(root: str | Path, mode: str = "COLD") -> list[tuple[str, Result]]:
    """Discover energy-specific CNPCS files and return them in energy order.

    The expected layout is ``root/<energy>MeV/ana_results/CNPCS_<mode>.dat``.
    Missing or invalid roots simply produce an empty list.
    """

    mode_name = str(mode).upper()
    if mode_name not in {"COLD", "HOT"}:
        raise ValueError("mode must be 'COLD' or 'HOT'")

    root_path = Path(root)
    matches: list[tuple[int, Path]] = []
    expected_name = f"CNPCS_{mode_name}.dat".casefold()
    for candidate in root_path.glob("*/ana_results/CNPCS_*.dat"):
        if candidate.name.casefold() != expected_name:
            continue
        energy_match = re.search(r"(\d+)\s*MeV", candidate.parent.parent.name, flags=re.IGNORECASE)
        if energy_match:
            matches.append((int(energy_match.group(1)), candidate))

    matches.sort(key=lambda item: item[0])
    return [(f"{energy} MeV", load_result(path)) for energy, path in matches]


def derive_title(path: str | Path) -> str:
    """Derive a concise reaction/energy/state title from a result path."""

    path_obj = Path(path)
    parts = path_obj.parts
    reaction = next((part for part in reversed(parts) if "+" in part), None)
    energy = next(
        (
            match.group(0)
            for part in reversed(parts)
            if (match := re.search(r"\d+\s*MeV", part, flags=re.IGNORECASE))
        ),
        None,
    )
    if energy:
        energy = re.sub(r"\s+", "", energy)
        energy = energy[:-3] + "MeV"
    state_match = re.search(r"(?:COLD|HOT)", path_obj.stem, flags=re.IGNORECASE)
    state = state_match.group(0).upper() if state_match else None

    fields = [field for field in (reaction, energy, state) if field]
    return ", ".join(fields) if fields else path_obj.stem
