#!/usr/bin/env python3
"""Compute reference-anchored variability and per-season drift details."""

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
    parser.add_argument("--signal-peptide-length", type=int, default=17)
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
        mature = ref_position - args.signal_peptide_length
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

    drift_positions = set(drift.index.astype(int))
    seasonal_counts = df.groupby("season")["id"].nunique()
    detail_rows = []
    drift_df = df[df["mature"].isin(drift_positions)]
    for (mature, season), group in drift_df.groupby(["mature", "season"], sort=True):
        n_compared = int(group["id"].nunique())
        n_different = int(group["differs"].sum())
        ref_aa = str(group["ref_aa"].iloc[0])
        substitutions = group.loc[group["differs"], "aa"].value_counts().sort_index()
        substitution_text = ";".join(
            f"{ref_aa}>{aa}:{count}" for aa, count in substitutions.items()
        )
        detail_rows.append(
            {
                "mature": int(mature),
                "ref_index": int(mature) + args.signal_peptide_length,
                "ref_aa": ref_aa,
                "season": season,
                "sequences_in_season": int(seasonal_counts.get(season, 0)),
                "sequences_compared": n_compared,
                "sequences_different": n_different,
                "fraction_different": n_different / n_compared if n_compared else 0.0,
                "substitutions": substitution_text or "none",
            }
        )

    detail_path = out_dir / "drift_detail.csv"
    pd.DataFrame(
        detail_rows,
        columns=[
            "mature",
            "ref_index",
            "ref_aa",
            "season",
            "sequences_in_season",
            "sequences_compared",
            "sequences_different",
            "fraction_different",
            "substitutions",
        ],
    ).to_csv(detail_path, index=False)

    print(f"Wrote variability summary to {variability_path}")
    print(f"Wrote drift subset to {drift_path}")
    print(f"Wrote detailed drift table to {detail_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
