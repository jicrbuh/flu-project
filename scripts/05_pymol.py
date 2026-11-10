#!/usr/bin/env python3
"""Create a PyMOL session script from the drift summary.

This is intentionally a scaffold: it reads a drift file and emits a .pml
showing residue selections for each season. Fill in the exact antigenic site
residue lists from a validated HA epitope map before final analysis.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create a PyMOL session script from the drift summary.")
    parser.add_argument("--drift-csv", type=Path, default=ROOT / "results" / "drift.csv")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results")
    return parser.parse_args()


def build_pml(drift: pd.DataFrame) -> str:
    seasons = [col for col in drift.columns if col != "mature"]
    lines = [
        "load 7KNV.pdb, prot",
        "select chainA, chain A",
        "hide everything, all",
        "show cartoon, chainA",
        "show sticks, chainA",
        "color white, chainA",
        "",
        "# NOTE: populate exact antigenic site residue ranges from a validated HA epitope map.",
        "# For now this script keeps the region-based color logic and records the trend.",
        "",
    ]

    for season in seasons:
        residues = sorted(int(resid) for resid in drift.index if pd.notna(drift.loc[resid, season]) and drift.loc[resid, season] > 0.05)
        if not residues:
            continue
        residue_list = ", ".join(str(res) for res in residues)
        lines.append(f"select {season.replace('-', '_')}, chainA and resi {residue_list}")
        lines.append(f"color {season_to_color(season)}, {season.replace('-', '_')}")
        lines.append("")

    # Example HA1/HA2 split: edit once a precise HA reference map is validated.
    lines.extend([
        "select head, chainA and resi 1-300",
        "select stem, chainA and resi 301-550",
        "color lightblue, head",
        "color orange, stem",
        "",
        "set cartoon_transparency, 0.2",
        "set sphere_scale, 0.25",
        "dss",
    ])
    return "\n".join(lines) + "\n"


def season_to_color(season: str) -> str:
    palette = [
        "red",
        "blue",
        "green",
        "yellow",
        "magenta",
        "cyan",
        "orange",
        "purple",
        "salmon",
        "gray",
    ]
    idx = abs(hash(season)) % len(palette)
    return palette[idx]


def main() -> int:
    args = parse_args()
    drift_path = args.drift_csv
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    if not drift_path.exists():
        raise FileNotFoundError(f"Drift file not found: {drift_path}. Run make variability first.")

    drift = pd.read_csv(drift_path, index_col=0)
    out_path = output_dir / "drift.pml"
    pml = build_pml(drift)
    out_path.write_text(pml, encoding="utf-8")
    print(f"Wrote PyMOL file to {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
