#!/usr/bin/env python3
"""Select representative influenza HA sequences per season.

- deduplicate exact duplicates by sequence
- remove partial/short sequences
- subsample N per season with a fixed random seed
- write selected FASTA files and a merged metadata table
"""

from __future__ import annotations

import argparse
import csv
import random
from pathlib import Path

import pandas as pd
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


def read_meta_for_season(meta_path: Path) -> pd.DataFrame:
    if not meta_path.exists():
        return pd.DataFrame(columns=["id", "strain", "season", "passage", "country"])
    return pd.read_csv(meta_path, dtype=str).fillna("unknown")


def main() -> int:
    args = parse_args()

    all_meta: list[pd.DataFrame] = []
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for raw_path in iter_raw_files(args.input_dir):
        if raw_path.name.startswith("selected_"):
            continue
        season = raw_path.stem
        records = dedupe_and_filter(raw_path, min_length=args.min_length)
        selected = subsample(records, n_per_season=args.n_per_season, seed=args.seed)
        selected_ids = {rec_id for rec_id, _ in selected}

        out_path = args.output_dir / f"selected_{season}.faa"
        with out_path.open("w") as fh:
            for rec_id, sequence in selected:
                fh.write(f">{rec_id}\n{sequence}\n")

        meta_path = args.input_dir / f"{season}_meta.csv"
        season_meta = read_meta_for_season(meta_path)
        season_meta = season_meta[season_meta["id"].isin(selected_ids)] if not season_meta.empty else season_meta
        if not season_meta.empty:
            all_meta.append(season_meta)

        print(f"Selected {len(selected)} sequences for season {season} -> {out_path}")

    if all_meta:
        merged = pd.concat(all_meta, ignore_index=True, sort=False)
        merged = merged[["id", "strain", "season", "passage", "country"]]
        merged_path = args.output_dir / "selected_meta.csv"
        merged.to_csv(merged_path, index=False)
        print(f"Wrote merged metadata to {merged_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
