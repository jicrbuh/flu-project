#!/usr/bin/env python3
"""Select representative influenza HA sequences per season.

- deduplicate exact duplicates by sequence
- remove partial/short sequences
- subsample N per season with a fixed random seed
- write selected FASTA files to data/selected/
"""

from __future__ import annotations

import argparse
import csv
import random
from pathlib import Path

from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw"
SELECTED_DIR = ROOT / "data" / "selected"
META_PATH = ROOT / "data" / "meta.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Subsample sequences per season.")
    parser.add_argument("--n-per-season", type=int, default=10)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--min-length", type=int, default=200)
    parser.add_argument("--input-dir", type=Path, default=RAW_DIR, help="Directory with per-season FASTA files")
    parser.add_argument("--output-dir", type=Path, default=SELECTED_DIR, help="Directory for selected sequences")
    return parser.parse_args()


def iter_raw_files(input_dir: Path) -> list[Path]:
    return sorted(input_dir.glob("*.faa"))


def dedupe_and_filter(path: Path, min_length: int = 200) -> list[tuple[str, str]]:
    seen: set[str] = set()
    selected: list[tuple[str, str]] = []
    for rec in SeqIO.parse(path, "fasta"):
        seq = str(rec.seq).replace("\n", "")
        if len(seq) < min_length:
            continue
        if seq in seen:
            continue
        seen.add(seq)
        selected.append((rec.id, seq))
    return selected


def subsample(records: list[tuple[str, str]], n_per_season: int = 10, seed: int = 42) -> list[tuple[str, str]]:
    rng = random.Random(seed)
    if len(records) <= n_per_season:
        return records
    return rng.sample(records, n_per_season)


def main() -> int:
    args = parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for raw_path in iter_raw_files(args.input_dir):
        season = raw_path.stem
        records = dedupe_and_filter(raw_path, min_length=args.min_length)
        selected = subsample(records, n_per_season=args.n_per_season, seed=args.seed)

        out_path = args.output_dir / f"selected_{season}.faa"
        with out_path.open("w") as fh:
            for rec_id, sequence in selected:
                fh.write(f">{rec_id}\n{sequence}\n")

        print(f"Selected {len(selected)} sequences for season {season} -> {out_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
