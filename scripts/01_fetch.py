#!/usr/bin/env python3
"""Fetch seasonal influenza HA records and write per-season FASTA + metadata.

This script is designed to fetch the actual HA sequences for a season window,
store the metadata from the GenBank source qualifiers, and write one output set
per season instead of appending to a global project-level metadata file.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd
from Bio import Entrez, SeqIO


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fetch influenza HA sequences for one season.")
    parser.add_argument("--season", required=True, help="Season label, e.g. 2024_25")
    parser.add_argument("--email", required=True, help="Email for Entrez requests")
    parser.add_argument("--output-dir", default=".", type=Path, help="Directory for the season FASTA and metadata")
    parser.add_argument("--retmax", type=int, default=3000, help="Maximum number of records to keep")
    parser.add_argument("--term", default=None, help="Optional override for the Entrez search term")
    parser.add_argument("--meta-output", default=None, type=Path, help="Optional explicit metadata path to write for this season")
    return parser.parse_args()


def parse_season(season: str) -> tuple[int, pd.Timestamp, pd.Timestamp]:
    season_clean = season.replace("-", "_")
    year_part = season_clean.split("_")[0]
    year0 = int(year_part)
    start = pd.Timestamp(f"{year0}-10-01")
    end = pd.Timestamp(f"{year0 + 1}-09-30")
    return year0, start, end


def is_human_host(host_value: str) -> bool:
    normalized = (host_value or "").strip().lower()
    return normalized in {"homo sapiens", "human"}


def build_default_term(season: str) -> str:
    year0, start, end = parse_season(season)
    return (
        '"Influenza A virus"[Organism] AND H1N1[All Fields] AND hemagglutinin[All Fields] '
        f'AND {start.strftime("%Y/%m/%d")}:{end.strftime("%Y/%m/%d")}[PDAT]'
    )


def fetch_records(term: str, email: str, retmax: int = 3000) -> tuple[list[str], str]:
    Entrez.email = email
    Entrez.api_key = None

    search = Entrez.esearch(db="nuccore", term=term, usehistory="y", retmax=0)
    result = Entrez.read(search)
    count = min(int(result.get("Count", 0)), retmax)
    if count == 0:
        return [], result.get("WebEnv", "")

    webenv = result.get("WebEnv", "")
    query_key = result.get("QueryKey", "")
    ids: list[str] = []
    for start in range(0, count, 200):
        handle = Entrez.efetch(
            db="nuccore",
            rettype="gb",
            retmode="text",
            retstart=start,
            retmax=200,
            webenv=webenv,
            query_key=query_key,
        )
        records = list(SeqIO.parse(handle, "gb"))
        ids.extend(rec.id for rec in records)
    return ids, webenv


def fetch_by_id(record_ids: list[str], webenv: str, email: str) -> list:
    Entrez.email = email
    if not record_ids:
        return []

    chunks: list[list[str]] = [record_ids[i : i + 200] for i in range(0, len(record_ids), 200)]
    records: list = []
    for chunk in chunks:
        handle = Entrez.efetch(
            db="nuccore",
            id=",".join(chunk),
            rettype="gb",
            retmode="text",
            webenv=webenv,
        )
        records.extend(list(SeqIO.parse(handle, "gb")))
    return records


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    term = args.term or build_default_term(args.season)
    year0, start_date, end_date = parse_season(args.season)

    Entrez.email = args.email
    search = Entrez.esearch(db="nuccore", term=term, usehistory="y", retmax=0)
    result = Entrez.read(search)
    count = min(int(result.get("Count", 0)), args.retmax)
    print(f"{count} candidate records")

    fasta_path = args.output_dir / f"{args.season}.faa"
    meta_path = args.meta_output or (args.output_dir / f"{args.season}_meta.csv")

    kept = 0
    with fasta_path.open("w", encoding="utf-8") as fasta_handle, meta_path.open("w", newline="", encoding="utf-8") as meta_handle:
        writer = csv.writer(meta_handle)
        writer.writerow(["id", "strain", "season", "passage", "country", "collection_date"])

        webenv = result.get("WebEnv", "")
        query_key = result.get("QueryKey", "")
        for start in range(0, count, 200):
            fetch_handle = Entrez.efetch(
                db="nuccore",
                rettype="gb",
                retmode="text",
                retstart=start,
                retmax=200,
                webenv=webenv,
                query_key=query_key,
            )
            for rec in SeqIO.parse(fetch_handle, "gb"):
                source = next((feat.qualifiers for feat in rec.features if feat.type == "source"), {})
                cds = next((feat for feat in rec.features if feat.type == "CDS"), None)
                if cds is None or "translation" not in cds.qualifiers:
                    continue

                host = (source.get("host", [""])[0] or "").strip()
                if host and not is_human_host(host):
                    continue

                collection_date = source.get("collection_date", [""])[0]
                d = pd.to_datetime(collection_date, errors="coerce")
                if pd.isna(d) or len(collection_date) <= 4 or not (start_date <= d <= end_date):
                    continue

                protein = cds.qualifiers["translation"][0]
                if len(protein) < 550:
                    continue

                strain = source.get("strain", [""])[0] or rec.annotations.get("source", "") or "unknown"
                serotype = source.get("serotype", [""])[0] or ""
                if "H1N1" not in serotype and "H1N1" not in rec.description and "H1N1" not in strain:
                    continue

                passage = source.get("lab_host", [""])[0] or source.get("note", ["unknown"])[0] or "unknown"
                country = source.get("geo_loc_name", ["unknown"])[0] or "unknown"

                fasta_handle.write(f">{rec.id}\n{protein}\n")
                writer.writerow([rec.id, strain, args.season, passage, country, collection_date])
                kept += 1

    print(f"kept {kept}")
    print(f"Wrote {fasta_path}")
    print(f"Wrote {meta_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
