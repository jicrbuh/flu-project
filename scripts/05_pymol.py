#!/usr/bin/env python3
"""Create a PyMOL script from the structure-aware seasonal drift report."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create a PyMOL script from validated structural drift positions.")
    parser.add_argument("--drift-csv", type=Path, default=ROOT / "results" / "structural_drift.csv")
    parser.add_argument("--pdb-id", default="7KNA", help="PDB identifier to load in PyMOL")
    parser.add_argument("--chain", default="A", help="PDB chain to select and color")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results")
    return parser.parse_args()


def build_pml(drift: pd.DataFrame, pdb_id: str = "7KNA", chain: str = "A") -> str:
    required = {"season", "fraction_different", "pdb_chain", "pdb_resi", "structure_status"}
    missing = required - set(drift.columns)
    if missing:
        raise ValueError(f"Structure-aware drift table is missing columns: {sorted(missing)}")

    chain_id = chain.upper()
    lines = [
        f"fetch {pdb_id}, async=0",
        "hide everything, all",
        f"show cartoon, {pdb_id}",
        f"color white, {pdb_id}",
        "",
        "# Selections include only above-threshold residues modeled in the validated chain mapping.",
        "# Positions absent from the model remain listed in structural_drift.csv but are not selected.",
        "",
    ]

    seasons = sorted(drift["season"].dropna().unique())
    for season_index, season in enumerate(seasons):
        season_rows = drift[
            (drift["season"] == season)
            & (drift["fraction_different"] > 0.05)
        ]
        modeled = season_rows[
            (season_rows["structure_status"] == "modeled")
            & (season_rows["pdb_chain"].astype(str).str.upper() == chain_id)
            & season_rows["pdb_resi"].notna()
            & (season_rows["pdb_resi"].astype(str) != "")
        ]
        residues = sorted(set(modeled["pdb_resi"].astype(str)), key=pdb_residue_sort_key)
        missing_count = len(season_rows) - len(modeled)
        if not residues:
            lines.append(f"# {season}: no above-threshold modeled residues; {missing_count} absent from the structure")
            lines.append("")
            continue

        sel_name = to_pymol_selection_name(str(season))
        residue_list = "+".join(residues)
        lines.append(
            f"select {sel_name}, {pdb_id} and chain {chain_id} and name CA and resi {residue_list}"
        )
        lines.append(f"show spheres, {sel_name}")
        lines.append(f"color {season_color(season_index)}, {sel_name}")
        lines.append(f"print('count_atoms {sel_name}: %d' % cmd.count_atoms('{sel_name}'))")
        lines.append(f"# {season}: {len(residues)} mapped positions; {missing_count} above-threshold positions absent")
        lines.append("")

    lines.extend(
        [
            "set cartoon_transparency, 0.2",
            "set sphere_scale, 0.6",
            "dss",
        ]
    )
    return "\n".join(lines) + "\n"


def season_color(index: int) -> str:
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
    return palette[index % len(palette)]


def pdb_residue_sort_key(value: str) -> tuple[int, str]:
    numeric = "".join(character for character in value if character in "-0123456789")
    return int(numeric or 0), value


def to_pymol_selection_name(season: str) -> str:
    name = "".join(ch for ch in season.replace("-", "").replace("_", "") if ch.isalnum())
    return f"s{name}" if name else "s"


def main() -> int:
    args = parse_args()
    if not args.drift_csv.exists():
        raise FileNotFoundError(f"Structure-aware drift report not found: {args.drift_csv}")

    drift = pd.read_csv(args.drift_csv, dtype={"pdb_resi": str})
    args.output_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.output_dir / "drift.pml"
    out_path.write_text(build_pml(drift, pdb_id=args.pdb_id, chain=args.chain), encoding="utf-8")
    print(f"Wrote PyMOL file to {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
