#!/usr/bin/env python3
"""Fetch influenza HA sequences by season and write FASTA + metadata rows.

This is a lightweight Entrez-based scaffold for one season at a time. It is
intended to be rerun as needed without modifying the source logic.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Iterable

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
META_PATH = DATA_DIR / "meta.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fetch influenza HA sequences by season.")
    parser.add_argument("--season", required=True, help="Season label, e.g. 2025-26")
    parser.add_argument("--term", default="Influenza A H1N1 hemagglutinin", help="Entrez search term")
    parser.add_argument("--email", default="", help="Email for Entrez requests")
    parser.add_argument("--retmax", type=int, default=200)
    parser.add_argument("--output-dir", type=Path, default=RAW_DIR, help="Directory for downloaded FASTA files")
    parser.add_argument("--meta-output", type=Path, default=None, help="Optional metadata CSV path to write for this season")
    return parser.parse_args()


def season_from_collection_date(collection_date: str) -> str:
    """Return season label in the form YYYY-YY, with Oct-Sep split."""
    try:
        year, month, day = collection_date[:10].split("-")
    except ValueError as exc:
        raise ValueError(f"Bad collection_date: {collection_date!r}") from exc
    month_int = int(month)
    year_int = int(year)
    season_start = year_int if month_int >= 10 else year_int - 1
    season_end = season_start + 1
    return f"{season_start}-{season_end % 100:02d}"


def ensure_meta_headers() -> None:
    META_PATH.parent.mkdir(parents=True, exist_ok=True)
    if not META_PATH.exists():
        with META_PATH.open("w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["id", "strain", "season", "passage", "country"])


def append_meta_rows(rows: Iterable[dict]) -> None:
    ensure_meta_headers()
    with META_PATH.open("a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "strain", "season", "passage", "country"])
        for row in rows:
            writer.writerow(row)


def fetch_fasta_records(term: str, email: str, retmax: int = 200) -> list[SeqRecord]:
    if not email:
        raise ValueError("Entrez email is required. Pass --email or set ENTREZ_EMAIL.")
    Entrez.email = email
    Entrez.api_key = None

    search = Entrez.esearch(db="protein", term=term, retmax=retmax, sort="relevance")
    search_handle = Entrez.read(search)
    ids = search_handle.get("IdList", [])
    if not ids:
        return []

    fetch_handle = Entrez.efetch(db="protein", id=",".join(ids), rettype="fasta", retmode="text")
    records = list(SeqIO.parse(fetch_handle, "fasta"))
    return records


def write_records(records: list[SeqRecord], season: str, output_dir: Path, meta_output: Path | None = None) -> list[dict]:
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{season}.faa"
    rows: list[dict] = []

    with output_path.open("w") as fh:
        for rec in records:
            fh.write(f">{rec.id}\n{str(rec.seq)}\n")
            desc = rec.description
            strain = desc.split("[")[-1].split("]")[0] if "[" in desc else desc
            row = {
                "id": rec.id,
                "strain": strain,
                "season": season,
                "passage": "unknown",
                "country": "unknown",
            }
            rows.append(row)

    meta_output = meta_output or output_dir / f"{season}_meta.csv"
    meta_output.parent.mkdir(parents=True, exist_ok=True)
    with meta_output.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "strain", "season", "passage", "country"])
        writer.writeheader()
        writer.writerows(rows)

    return rows


def main() -> int:
    args = parse_args()

    records = fetch_fasta_records(args.term, args.email, retmax=args.retmax)
    meta_output = args.meta_output or args.output_dir / f"{args.season}_meta.csv"
    rows = write_records(records, args.season, args.output_dir, meta_output=meta_output)
    append_meta_rows(rows)

    print(f"Fetched {len(records)} records for season {args.season} into {args.output_dir / f'{args.season}.faa'}")
    print(f"Wrote season metadata to {meta_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
