#!/usr/bin/env python3
"""Generate a PyMOL script that renders seasonal drift positions from four angles."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
REQUIRED_COLUMNS = {
    "mature",
    "season",
    "fraction_different",
    "pdb_chain",
    "pdb_resi",
    "structure_status",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate four PyMOL views of mapped seasonal HA drift positions."
    )
    parser.add_argument(
        "--drift-csv",
        type=Path,
        default=ROOT / "results" / "structural_drift.csv",
    )
    parser.add_argument(
        "--structure",
        type=Path,
        default=ROOT / "data" / "structure" / "structure.cif",
    )
    parser.add_argument("--season", default="2024_25")
    parser.add_argument("--chain", default="A")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "figures")
    parser.add_argument("--width", type=int, default=1600)
    parser.add_argument("--height", type=int, default=1200)
    parser.add_argument("--label-size", type=int, default=9)
    return parser.parse_args()


def read_drift_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
        return list(reader)


def selected_positions(
    rows: Iterable[dict[str, str]], season: str, chain: str
) -> list[tuple[int, str]]:
    positions: dict[int, str] = {}
    for row in rows:
        if row["season"] != season:
            continue
        if row["structure_status"] != "modeled":
            continue
        if row["pdb_chain"].upper() != chain.upper():
            continue
        try:
            fraction = float(row["fraction_different"])
            mature = int(row["mature"])
        except (TypeError, ValueError) as error:
            raise ValueError(f"Invalid drift row: {row}") from error
        if fraction <= 0:
            continue

        pdb_resi = row["pdb_resi"].strip()
        if not re.fullmatch(r"-?\d+[A-Za-z]?", pdb_resi):
            raise ValueError(f"Invalid PDB residue identifier: {pdb_resi!r}")
        existing_resi = positions.get(mature)
        if existing_resi is not None and existing_resi != pdb_resi:
            raise ValueError(
                f"Mature position {mature} maps to multiple PDB residues "
                f"for chain {chain}: {existing_resi}, {pdb_resi}"
            )
        positions[mature] = pdb_resi

    return sorted(positions.items(), key=lambda item: item[0])


def quote_pml_path(path: Path) -> str:
    return '"' + path.as_posix().replace("\\", "\\\\").replace('"', '\\"') + '"'


def build_pml(
    positions: list[tuple[int, str]],
    structure: Path,
    output_dir: Path,
    season: str = "2024_25",
    chain: str = "A",
    width: int = 1600,
    height: int = 1200,
    label_size: int = 9,
) -> str:
    if not positions:
        raise ValueError(f"No mapped, non-zero drift positions for {season} chain {chain}")
    if width <= 0 or height <= 0 or label_size <= 0:
        raise ValueError("Image dimensions and label size must be positive")

    safe_season = re.sub(r"[^A-Za-z0-9_-]", "", season)
    residues = sorted({resi for _, resi in positions}, key=pdb_residue_sort_key)
    residue_list = "+".join(residues)
    trimer_selection = "ha and chain A+B+C"
    lines = [
        f"load {quote_pml_path(structure.resolve())}, ha",
        "hide everything, all",
        "bg_color white",
        "set ray_opaque_background, on",
        "set antialias, 2",
        f"show cartoon, {trimer_selection}",
        "color red, ha and chain A",
        "color green, ha and chain B",
        "color blue, ha and chain C",
        f"select drift_sites, ha and chain {chain.upper()} and name CA and resi {residue_list}",
        "show spheres, drift_sites",
        "color yellow, drift_sites",
        "set sphere_scale, 0.65",
        "set label_color, black",
        "set label_outline_color, white",
        f"set label_size, {label_size}",
        "set label_connector, on",
    ]

    for mature, pdb_resi in positions:
        lines.append(
            f'label (ha and chain {chain.upper()} and name CA and resi {pdb_resi}), "{mature}"'
        )

    lines.extend(
        [
            f"orient {trimer_selection}",
            f"zoom {trimer_selection}, 5",
        ]
    )
    for angle in range(1, 5):
        image_path = output_dir.resolve() / f"{safe_season}_selected_residues_angle{angle}.png"
        lines.append(f"png {image_path.as_posix()}, {width}, {height}, ray=1")
        if angle < 4:
            lines.append("turn y, 90")

    lines.append(
        f"print('Rendered {len(positions)} mapped {season} drift positions on chain {chain.upper()}')"
    )
    return "\n".join(lines) + "\n"


def pdb_residue_sort_key(value: str) -> tuple[int, str]:
    match = re.fullmatch(r"(-?\d+)([A-Za-z]?)", value)
    if match is None:
        raise ValueError(f"Invalid PDB residue identifier: {value!r}")
    return int(match.group(1)), match.group(2)


def main() -> int:
    args = parse_args()
    if not args.drift_csv.is_file():
        raise FileNotFoundError(f"Structure-aware drift report not found: {args.drift_csv}")
    if not args.structure.is_file():
        raise FileNotFoundError(f"Structure file not found: {args.structure}")

    positions = selected_positions(
        read_drift_rows(args.drift_csv),
        season=args.season,
        chain=args.chain,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pml_path = args.output_dir / f"{re.sub(r'[^A-Za-z0-9_-]', '', args.season)}_structure_views.pml"
    pml_path.write_text(
        build_pml(
            positions,
            structure=args.structure,
            output_dir=args.output_dir,
            season=args.season,
            chain=args.chain,
            width=args.width,
            height=args.height,
            label_size=args.label_size,
        ),
        encoding="utf-8",
    )
    print(f"Wrote {pml_path}")
    print(f"Mapped positions: {len(positions)}")
    print(f"Open in PyMOL with: pymol {pml_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
