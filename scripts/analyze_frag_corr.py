#!/usr/bin/env python3
"""Summarize primary-excitation correlations from FRAG_CORR.dat."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


FIELDS = (
    "Z0 A0 J0 E0 Zf Af Ef Jf Id WGT".split()
)


def parse_target(value):
    try:
        z_text, n_text = value.split(",", 1)
        return int(z_text), int(n_text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("target must be Z,N (for example 90,137)") from exc


def read_records(path):
    records = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            text = raw_line.strip()
            if not text or text.startswith("#"):
                continue
            parts = text.split()
            if len(parts) != len(FIELDS):
                raise ValueError(
                    f"{path}:{line_number}: expected {len(FIELDS)} fields, got {len(parts)}"
                )
            try:
                values = {
                    "Z0": int(parts[0]), "A0": int(parts[1]),
                    "J0": float(parts[2]), "E0": float(parts[3]),
                    "Zf": int(parts[4]), "Af": int(parts[5]),
                    "Ef": float(parts[6]), "Jf": float(parts[7]),
                    "Id": int(parts[8]), "WGT": float(parts[9]),
                }
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: invalid numeric field") from exc
            if values["Zf"] < 0 or values["Af"] <= 0:
                raise ValueError(
                    f"{path}:{line_number}: invalid selected final fragment"
                )
            records.append(values)
    return records


def energy_bin(e0, edges):
    if e0 < edges[0]:
        return "low"
    if e0 < edges[1]:
        return "mid"
    return "high"


def summarize(records, selected, edges):
    groups = defaultdict(lambda: {"records": 0, "weight": 0.0,
                                  "bins": defaultdict(float)})
    for item in records:
        isotope = item["Zf"], item["Af"] - item["Zf"]
        if selected and isotope not in selected:
            continue
        group = groups[isotope]
        group["records"] += 1
        group["weight"] += item["WGT"]
        group["bins"][energy_bin(item["E0"], edges)] += item["WGT"]
    return groups


def print_table(title, groups, edges):
    print(title)
    print("Zf Nf Af records weight p_low p_mid p_high")
    for (zf, nf), group in sorted(groups.items()):
        total = group["weight"]
        probabilities = [group["bins"][name] / total if total else 0.0
                         for name in ("low", "mid", "high")]
        print(f"{zf:2d} {nf:3d} {zf + nf:3d} {group['records']:7d} "
              f"{total:10.5g} {probabilities[0]:.6f} {probabilities[1]:.6f} "
              f"{probabilities[2]:.6f}")
    print(f"bins: low E0 < {edges[0]:g}; mid {edges[0]:g} <= E0 < "
          f"{edges[1]:g}; high E0 >= {edges[1]:g}")


def write_csv(path, groups, edges):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Zf", "Nf", "Af", "records", "weight", "bin",
                         "p_e0_given_isotope", "lower_mev", "upper_mev"])
        limits = (("low", "", edges[0]), ("mid", edges[0], edges[1]),
                  ("high", edges[1], ""))
        for (zf, nf), group in sorted(groups.items()):
            total = group["weight"]
            for name, lower, upper in limits:
                probability = group["bins"][name] / total if total else 0.0
                writer.writerow([zf, nf, zf + nf, group["records"], group["weight"],
                                 name, probability, lower, upper])


def print_background(records, edges):
    weight = sum(item["WGT"] for item in records)
    bins = defaultdict(float)
    for item in records:
        bins[energy_bin(item["E0"], edges)] += item["WGT"]
    probabilities = [bins[name] / weight if weight else 0.0
                     for name in ("low", "mid", "high")]
    print("Overall final-leaf background")
    print("records weight p_low p_mid p_high")
    print(f"{len(records):7d} {weight:10.5g} {probabilities[0]:.6f} "
          f"{probabilities[1]:.6f} {probabilities[2]:.6f}")
    print(f"bins: low E0 < {edges[0]:g}; mid {edges[0]:g} <= E0 < "
          f"{edges[1]:g}; high E0 >= {edges[1]:g}")


def main():
    parser = argparse.ArgumentParser(
        description="Report weighted p(E0 | final isotope) from FRAG_CORR.dat."
    )
    parser.add_argument("input", nargs="?", default="FRAG_CORR.dat", type=Path)
    parser.add_argument("--target", action="append", type=parse_target, default=[],
                        help="final isotope Z,N; repeat to select several")
    parser.add_argument("--edges", nargs=2, type=float, metavar=("LOW_MID", "MID_HIGH"),
                        default=(20.0, 100.0), help="two E0 bin edges in MeV")
    parser.add_argument("--csv", type=Path, help="write selected-isotope table as CSV")
    parser.add_argument("--expected-rows", type=int,
                        help="fail unless the data-row count has this value")
    args = parser.parse_args()
    if args.edges[0] >= args.edges[1]:
        parser.error("--edges must be strictly increasing")

    records = read_records(args.input)
    if args.expected_rows is not None and len(records) != args.expected_rows:
        parser.error(
            f"expected {args.expected_rows} data rows, found {len(records)}"
        )
    all_groups = summarize(records, set(), args.edges)
    selected = set(args.target)
    selected_groups = summarize(records, selected, args.edges)

    print(f"FRAG_CORR records: {len(records)}")
    if args.expected_rows is not None:
        print(f"Rows / expected rows: {len(records)} / {args.expected_rows}")
    print("Coverage: each record is the largest valid final leaf of the MST "
          "maxloc(ICL) source that reached GEMINI. It is not a CNPCS-COLD "
          "denominator. With kktot=1, row order follows event L1.")
    print_table("All final isotopes", all_groups, args.edges)
    if selected:
        print_table("Requested final isotopes", selected_groups, args.edges)
    print_background(records, args.edges)
    if args.csv:
        write_csv(args.csv, selected_groups if selected else all_groups, args.edges)


if __name__ == "__main__":
    main()
