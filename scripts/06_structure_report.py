#!/usr/bin/env python3
"""Join seasonal drift details to the coordinate-validated PDB mapping."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create a structure-aware seasonal drift report.")
    parser.add_argument("--drift-detail", type=Path, required=True)
    parser.add_argument("--numbering-map", type=Path, required=True)
    parser.add_argument("--pdb-id", default="7KNA")
    parser.add_argument("--chain", default="A")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results")
    return parser.parse_args()


def build_structure_report(
    drift_detail: pd.DataFrame,
    numbering_map: pd.DataFrame,
) -> tuple[pd.DataFrame, dict]:
    required_detail = {
        "mature",
        "ref_index",
        "ref_aa",
        "season",
        "sequences_in_season",
        "sequences_compared",
        "sequences_different",
        "fraction_different",
        "substitutions",
    }
    required_map = {
        "ref_index",
        "mature_index",
        "pdb_chain",
        "pdb_resi",
        "pdb_aa",
        "identity_match",
        "mapping_status",
    }
    missing_detail = required_detail - set(drift_detail.columns)
    missing_map = required_map - set(numbering_map.columns)
    if missing_detail:
        raise ValueError(f"Drift detail is missing columns: {sorted(missing_detail)}")
    if missing_map:
        raise ValueError(f"Numbering map is missing columns: {sorted(missing_map)}")
    if numbering_map["mature_index"].duplicated().any():
        raise ValueError("Numbering map has duplicate mature_index values")

    mapping = numbering_map[
        [
            "ref_index",
            "mature_index",
            "pdb_chain",
            "pdb_resi",
            "pdb_aa",
            "identity_match",
            "mapping_status",
        ]
    ].rename(columns={"ref_index": "mapped_ref_index"})
    report = drift_detail.merge(
        mapping,
        left_on="mature",
        right_on="mature_index",
        how="left",
        validate="many_to_one",
    )
    if report["mapping_status"].isna().any():
        missing_positions = sorted(report.loc[report["mapping_status"].isna(), "mature"].unique())
        raise ValueError(f"Drift positions absent from numbering map: {missing_positions}")
    if (report["ref_index"] != report["mapped_ref_index"]).any():
        raise ValueError("Drift and numbering tables use inconsistent reference coordinates")

    report["structure_status"] = report["mapping_status"]
    report = report.drop(columns=["mapped_ref_index", "mature_index", "mapping_status"])

    drift_positions = set(report["mature"].astype(int))
    modeled_positions = set(
        report.loc[report["structure_status"] == "modeled", "mature"].astype(int)
    )
    unmapped_positions = sorted(drift_positions - modeled_positions)
    identity_mismatches = int(
        report.loc[
            (report["structure_status"] == "modeled")
            & (report["identity_match"] == "no"),
            "mature",
        ].nunique()
    )
    stats = {
        "drift_positions": len(drift_positions),
        "modeled_positions": len(modeled_positions),
        "unmodeled_positions": unmapped_positions,
        "season_position_rows": len(report),
        "mapped_season_position_rows": int((report["structure_status"] == "modeled").sum()),
        "reference_structure_mismatches": identity_mismatches,
    }
    return report, stats


def main() -> int:
    args = parse_args()
    drift_detail = pd.read_csv(args.drift_detail)
    numbering_map = pd.read_csv(args.numbering_map, dtype={"pdb_resi": str})
    report, stats = build_structure_report(drift_detail, numbering_map)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    report_path = args.output_dir / "structural_drift.csv"
    validation_path = args.output_dir / "structure_validation.txt"
    report.to_csv(report_path, index=False)

    unmapped = ",".join(str(position) for position in stats["unmodeled_positions"]) or "none"
    validation_path.write_text(
        "STRUCTURAL DRIFT REPORT GENERATED\n"
        f"PDB: {args.pdb_id.upper()}\n"
        f"Chain: {args.chain}\n"
        f"Distinct drift positions: {stats['drift_positions']}\n"
        f"Positions with modeled residues: {stats['modeled_positions']}\n"
        f"Positions absent from the model: {unmapped}\n"
        f"Season-position rows mapped to modeled residues: {stats['mapped_season_position_rows']}/{stats['season_position_rows']}\n"
        f"Modeled reference/PDB amino-acid mismatches: {stats['reference_structure_mismatches']}\n"
        "Unmodeled positions are retained in structural_drift.csv and must not be interpreted from the structure.\n"
        "This report does not assign antigenic sites or establish antigenic effect.\n",
        encoding="utf-8",
    )

    print(f"Wrote structure-aware drift table to {report_path}")
    print(f"Wrote structural report summary to {validation_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
