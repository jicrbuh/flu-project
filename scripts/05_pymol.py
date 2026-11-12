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
    parser.add_argument("--pdb-id", default="7KNA", help="PDB identifier to load in PyMOL")
    parser.add_argument("--chain", default="A", help="PDB chain to select and color")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results")
    return parser.parse_args()


def build_pml(drift: pd.DataFrame, numbering_map: pd.DataFrame | None = None, pdb_id: str = "7KNA", chain: str = "A") -> str:
    seasons = [col for col in drift.columns if col != "mature"]
    residue_map = {}
    if numbering_map is not None and not numbering_map.empty:
        if "mature_index" in numbering_map.columns and "pdb_resi" in numbering_map.columns:
            residue_map = dict(zip(numbering_map["mature_index"].astype(int), numbering_map["pdb_resi"].astype(int)))

    chain_id = chain.upper()
    lines = [
        f"fetch {pdb_id}, async=0",
        "hide everything, all",
        f"show cartoon, {pdb_id}",
        f"color white, {pdb_id}",
        "",
        "# The residue list below is derived from the mature HA numbering in the drift table and",
        "# then translated to PDB numbering using numbering_map.csv where available.",
        "# Use the count_atoms check below to confirm the selection contains modeled residues.",
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
        sel_name = to_pymol_selection_name(season)
        residue_list = "+".join(str(res) for res in residues)
        lines.append(f"select {sel_name}, {pdb_id} and chain {chain_id} and name CA and resi {residue_list}")
        lines.append(f"show spheres, {sel_name}")
        lines.append(f"color {season_to_color(season)}, {sel_name}")
        lines.append(f"print('count_atoms {sel_name}: %d' % cmd.count_atoms('{sel_name}'))")
        lines.append("")

    lines.extend([
        "",
        "# Domain coloring is intentionally kept minimal and should be validated against a published",
        "# HA domain definition before publication. The season-specific residue selections are the primary output.",
        "",
        "set cartoon_transparency, 0.2",
        "set sphere_scale, 0.6",
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


def to_pymol_selection_name(season: str) -> str:
    name = season.replace("-", "").replace("_", "")
    name = ''.join(ch for ch in name if ch.isalnum())
    if not name:
        return "s"
    return f"s{name}"


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
    pml = build_pml(drift, numbering_map, pdb_id=args.pdb_id, chain=args.chain)
    out_path.write_text(pml, encoding="utf-8")
    print(f"Wrote PyMOL file to {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
