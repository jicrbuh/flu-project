#!/usr/bin/env python3
"""Map reference HA numbering to PDB residue numbering.

Workflow:
1. Load a reference FASTA.
2. Download the matching PDB chain FASTA from RCSB.
3. Pairwise-align the two sequences.
4. Compute a mature-peptide numbering map: idx - 17 + offset.
5. Save a CSV ready for downstream structure work.
"""

from __future__ import annotations

import argparse
import csv
from io import StringIO
from pathlib import Path

import requests
from Bio import SeqIO
from Bio.Align import PairwiseAligner

ROOT = Path(__file__).resolve().parents[1]
REF_PATH = ROOT / "data" / "ref" / "ref_A_Missouri_11_2025.faa"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute HA reference-to-PDB numbering.")
    parser.add_argument("--reference", type=Path, default=REF_PATH, help="Reference FASTA path")
    parser.add_argument("--pdb-id", default="7KNA", help="PDB identifier to fetch from RCSB")
    parser.add_argument("--chain", default="A", help="PDB chain to use")
    parser.add_argument("--expected-protein", default="Hemagglutinin", help="Expected protein name in the PDB FASTA header")
    parser.add_argument("--expected-organism", default="Influenza A virus", help="Expected organism in the PDB FASTA header")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results", help="Directory to store numbering outputs")
    return parser.parse_args()


def load_reference(path: Path) -> str:
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        raise FileNotFoundError(f"No FASTA records found in {path}")
    return str(records[0].seq)


def validate_pdb_record(record, pdb_id: str, chain_id: str, expected_protein: str, expected_organism: str) -> None:
    label = f"{record.id} {record.description}".upper()
    description = record.description.upper()
    protein_ok = expected_protein.upper() in description
    organism_ok = expected_organism.upper() in description
    chain_ok = (
        f"{pdb_id}_{chain_id}".upper() in label
        or f"CHAIN {chain_id}" in label
        or f"CHAINS" in label and str(chain_id).upper() in label
    )
    if not chain_ok:
        raise ValueError(f"Fetched record {record.id} did not match PDB {pdb_id} chain {chain_id}")
    if not protein_ok:
        raise ValueError(
            f"PDB {pdb_id} fetched a record whose description does not match expected protein '{expected_protein}': {record.description}"
        )
    if not organism_ok:
        raise ValueError(
            f"PDB {pdb_id} fetched a record whose description does not match expected organism '{expected_organism}': {record.description}"
        )


def fetch_pdb_fasta(
    pdb_id: str,
    chain_id: str = "A",
    expected_protein: str = "Hemagglutinin",
    expected_organism: str = "Influenza A virus",
) -> str:
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}/display"
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    records = list(SeqIO.parse(StringIO(response.text), "fasta"))
    if not records:
        raise ValueError(f"No FASTA records returned for PDB {pdb_id} from RCSB")

    candidate = None
    for record in records:
        if record.id.startswith(f"{pdb_id}_{chain_id}") or record.id.endswith(f"_{chain_id}"):
            candidate = record
            break
    if candidate is None:
        candidate = records[0]

    validate_pdb_record(candidate, pdb_id, chain_id, expected_protein, expected_organism)
    return str(candidate.seq)


def build_numbering_map(reference_seq: str, pdb_seq: str) -> list[dict]:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(reference_seq, pdb_seq)
    if not alignments:
        raise RuntimeError("Pairwise alignment failed.")
    alignment = alignments[0]
    ref_aligned, pdb_aligned = str(alignment[0]), str(alignment[1])

    ref_index = 0
    pdb_index = 0
    offset = None
    rows: list[dict] = []

    for ref_char, pdb_char in zip(ref_aligned, pdb_aligned):
        if ref_char != "-":
            ref_index += 1
        if pdb_char != "-":
            pdb_index += 1

        if ref_char != "-" and pdb_char != "-":
            mature_index = ref_index - 17
            if mature_index >= 1 and offset is None:
                offset = pdb_index - mature_index
            if mature_index >= 1:
                pdb_resi = mature_index + offset
                rows.append(
                    {
                        "ref_index": ref_index,
                        "mature_index": mature_index,
                        "pdb_resi": pdb_resi,
                        "ref_aa": ref_char,
                        "pdb_aa": pdb_char,
                    }
                )

    if offset is None:
        raise RuntimeError("No mature-peptide residues were aligned to the PDB chain.")

    return rows


def main() -> int:
    args = parse_args()

    reference_seq = load_reference(args.reference)
    pdb_seq = fetch_pdb_fasta(
        args.pdb_id,
        args.chain,
        expected_protein=args.expected_protein,
        expected_organism=args.expected_organism,
    )
    rows = build_numbering_map(reference_seq, pdb_seq)

    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "numbering_map.csv"
    check_path = out_dir / "numbering_check.txt"

    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["ref_index", "mature_index", "pdb_resi", "ref_aa", "pdb_aa"])
        writer.writeheader()
        writer.writerows(rows)

    first = rows[0]
    later = rows[min(len(rows) - 1, 50)]
    check_path.write_text(
        "Reference-to-PDB numbering diagnostic\n"
        f"First mature residue: ref_index={first['ref_index']}, mature_index={first['mature_index']}, pdb_resi={first['pdb_resi']}\n"
        f"Later check: ref_index={later['ref_index']}, mature_index={later['mature_index']}, pdb_resi={later['pdb_resi']}\n"
        "Verify one landmark directly in PyMOL, e.g. select a residue near the receptor binding site and compare numbering.\n",
        encoding="utf-8",
    )

    print(f"Saved numbering map to {out_path}")
    print(f"Saved numbering check summary to {check_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
