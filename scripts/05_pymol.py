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
    parser.add_argument("--numbering-csv", type=Path, default=None, help="Numbering map to map mature positions to PDB residue IDs")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results")
    return parser.parse_args()


def build_pml(drift: pd.DataFrame, numbering_map: pd.DataFrame | None = None) -> str:
    seasons = [col for col in drift.columns if col != "mature"]
    residue_map = {}
    if numbering_map is not None and not numbering_map.empty:
        if "mature_index" in numbering_map.columns and "pdb_resi" in numbering_map.columns:
            residue_map = dict(zip(numbering_map["mature_index"].astype(int), numbering_map["pdb_resi"].astype(int)))

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
        residues = []
        for mature in drift.index:
            mature_int = int(mature)
            score = drift.loc[mature, season]
            if pd.notna(score) and score > 0.05:
                pdb_resi = residue_map.get(mature_int, mature_int)
                residues.append(int(pdb_resi))
        residues = sorted(set(residues))
        if not residues:
            continue
        residue_list = ", ".join(str(res) for res in residues)
        lines.append(f"select {season.replace('-', '_')}, chainA and resi {residue_list}")
        lines.append(f"color {season_to_color(season)}, {season.replace('-', '_')}")
        lines.append("")

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
    numbering_map = None
    if args.numbering_csv is not None and args.numbering_csv.exists():
        numbering_map = pd.read_csv(args.numbering_csv)

    out_path = output_dir / "drift.pml"
    pml = build_pml(drift, numbering_map)
    out_path.write_text(pml, encoding="utf-8")
    print(f"Wrote PyMOL file to {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
