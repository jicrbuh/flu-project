#!/usr/bin/env python3
"""Compute reference-anchored variability across seasons.

The output is a table with columns roughly matching the project design:
  mature, season, aa, ref_aa, id, differs, etc.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio import AlignIO

ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute reference-anchored variability across seasons.")
    parser.add_argument("--alignment", type=Path, default=ROOT / "results" / "aln.faa")
    parser.add_argument("--meta", type=Path, default=ROOT / "data" / "meta.csv")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    align_path = args.alignment
    meta_path = args.meta
    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    if not align_path.exists():
        raise FileNotFoundError(f"Alignment file not found: {align_path}. Run make align first.")
    if not meta_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {meta_path}. Build it with the fetch step first.")

    aln = AlignIO.read(str(align_path), "fasta")
    if len(aln) < 2:
        raise ValueError("Alignment needs a reference plus at least one sequence.")

    ref = aln[0]
    meta = pd.read_csv(meta_path).fillna("unknown").set_index("id")

    rows = []
    ref_position = 0

    for col in range(aln.get_alignment_length()):
        if ref.seq[col] == "-":
            continue
        ref_position += 1

        # compute mature peptide numbering with the same signal-peptide offset used elsewhere
        mature = ref_position - 17
        if mature < 1:
            continue

        for rec in aln[1:]:
            aa = rec.seq[col]
            if aa == "-":
                continue
            rows.append(
                {
                    "idx": ref_position,
                    "mature": mature,
                    "ref_aa": ref.seq[col],
                    "aa": aa,
                    "id": rec.id,
                }
            )

    df = pd.DataFrame(rows)
    if df.empty:
        raise ValueError("No non-gap reference-anchored variants were found.")

    df = df.join(meta, on="id")
    df["differs"] = df["aa"] != df["ref_aa"]
    df = df[(df["aa"] != "-") & (df["mature"] >= 1)]

    summary = (
        df.groupby(["mature", "season"])["differs"]
        .mean()
        .unstack(fill_value=0.0)
        .sort_index()
    )

    variability_path = out_dir / "variability.csv"
    summary.to_csv(variability_path)

    drift = summary[(summary > 0.05).any(axis=1)]
    drift_path = out_dir / "drift.csv"
    drift.to_csv(drift_path)

    print(f"Wrote variability summary to {variability_path}")
    print(f"Wrote drift subset to {drift_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
