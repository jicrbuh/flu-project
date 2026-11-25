#!/usr/bin/env python3
"""Map reference HA positions to modeled PDB residues."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import requests
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1

ROOT = Path(__file__).resolve().parents[1]
REF_PATH = ROOT / "data" / "ref" / "ref_A_Missouri_11_2025.faa"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Map reference HA positions to modeled PDB residues.")
    parser.add_argument("--reference", type=Path, default=REF_PATH, help="Reference FASTA path")
    parser.add_argument("--structure", type=Path, default=None, help="PDB mmCIF file; fetched from RCSB if omitted")
    parser.add_argument("--pdb-id", default="7KNA", help="PDB identifier")
    parser.add_argument("--chain", default="A", help="PDB author chain ID")
    parser.add_argument("--signal-peptide-length", type=int, default=17)
    parser.add_argument("--min-aligned", type=int, default=250, help="Minimum aligned HA residues required")
    parser.add_argument("--min-identity", type=float, default=0.80, help="Minimum aligned sequence identity")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results")
    return parser.parse_args()


def load_reference(path: Path) -> str:
    records = list(SeqIO.parse(path, "fasta"))
    if len(records) != 1:
        raise ValueError(f"Expected exactly one reference FASTA record in {path}, found {len(records)}")
    return str(records[0].seq).upper()


def download_structure(pdb_id: str, output_path: Path) -> Path:
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_bytes(response.content)
    return output_path


def get_chain_residues(structure_path: Path, chain_id: str):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", str(structure_path))
    model = next(structure.get_models())
    if chain_id not in model:
        available = ", ".join(chain.id for chain in model)
        raise ValueError(f"Chain {chain_id!r} not found in {structure_path}; available chains: {available}")

    residues = [
        residue
        for residue in model[chain_id]
        if residue.id[0] == " " and "CA" in residue
    ]
    if not residues:
        raise ValueError(f"No modeled protein residues with CA atoms found in chain {chain_id!r}")
    return residues


def build_numbering_map(
    reference_seq: str,
    structure_residues,
    chain_id: str,
    signal_peptide_length: int = 17,
    min_aligned: int = 250,
    min_identity: float = 0.80,
) -> tuple[list[dict], dict]:
    structure_seq = "".join(
        seq1(residue.resname, custom_map={"MSE": "M"}, undef_code="X")
        for residue in structure_residues
    )
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(reference_seq, structure_seq)
    if not alignments:
        raise RuntimeError("Reference-to-structure sequence alignment failed.")
    alignment = alignments[0]

    ref_to_residue = {}
    identical = 0
    aligned = 0
    ref_blocks, structure_blocks = alignment.aligned
    for (ref_start, ref_end), (structure_start, structure_end) in zip(ref_blocks, structure_blocks):
        if ref_end - ref_start != structure_end - structure_start:
            raise RuntimeError("Alignment produced incompatible paired blocks.")
        for offset in range(ref_end - ref_start):
            ref_index = int(ref_start + offset) + 1
            residue = structure_residues[int(structure_start + offset)]
            ref_to_residue[ref_index] = residue
            aligned += 1
            ref_aa = reference_seq[ref_index - 1]
            structure_aa = seq1(residue.resname, custom_map={"MSE": "M"}, undef_code="X")
            identical += ref_aa == structure_aa

    identity = identical / aligned if aligned else 0.0
    if aligned < min_aligned:
        raise ValueError(
            f"Only {aligned} reference residues aligned to {chain_id}; expected at least {min_aligned}"
        )
    if identity < min_identity:
        raise ValueError(
            f"Reference-to-structure identity is {identity:.3f}; expected at least {min_identity:.3f}"
        )

    rows = []
    for ref_index, ref_aa in enumerate(reference_seq, start=1):
        mature_index = ref_index - signal_peptide_length
        if mature_index < 1:
            continue

        residue = ref_to_residue.get(ref_index)
        if residue is None:
            rows.append(
                {
                    "ref_index": ref_index,
                    "mature_index": mature_index,
                    "ref_aa": ref_aa,
                    "pdb_chain": chain_id,
                    "pdb_resi": "",
                    "pdb_aa": "",
                    "identity_match": "",
                    "mapping_status": "not_modeled",
                }
            )
            continue

        pdb_aa = seq1(residue.resname, custom_map={"MSE": "M"}, undef_code="X")
        pdb_resi = f"{residue.id[1]}{residue.id[2].strip()}"
        rows.append(
            {
                "ref_index": ref_index,
                "mature_index": mature_index,
                "ref_aa": ref_aa,
                "pdb_chain": chain_id,
                "pdb_resi": pdb_resi,
                "pdb_aa": pdb_aa,
                "identity_match": "yes" if ref_aa == pdb_aa else "no",
                "mapping_status": "modeled",
            }
        )

    stats = {
        "aligned_residues": aligned,
        "identical_residues": identical,
        "identity": identity,
        "reference_coverage": len(ref_to_residue) / len(reference_seq),
        "mapped_reference_positions": len(ref_to_residue),
    }
    return rows, stats


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    structure_path = args.structure
    if structure_path is None:
        structure_path = download_structure(args.pdb_id, args.output_dir / f"{args.pdb_id.upper()}.cif")
    if not structure_path.exists():
        raise FileNotFoundError(f"Structure file not found: {structure_path}")

    reference_seq = load_reference(args.reference)
    residues = get_chain_residues(structure_path, args.chain)
    rows, stats = build_numbering_map(
        reference_seq,
        residues,
        args.chain,
        signal_peptide_length=args.signal_peptide_length,
        min_aligned=args.min_aligned,
        min_identity=args.min_identity,
    )

    out_dir = args.output_dir
    map_path = out_dir / "numbering_map.csv"
    check_path = out_dir / "numbering_check.txt"

    with map_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "ref_index",
                "mature_index",
                "ref_aa",
                "pdb_chain",
                "pdb_resi",
                "pdb_aa",
                "identity_match",
                "mapping_status",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    modeled = sum(row["mapping_status"] == "modeled" for row in rows)
    check_path.write_text(
        "STRUCTURE CHAIN MAPPING PASSED\n"
        f"PDB: {args.pdb_id.upper()}\n"
        f"Chain: {args.chain}\n"
        f"Reference residues aligned: {stats['aligned_residues']}\n"
        f"Exact amino-acid matches: {stats['identical_residues']}\n"
        f"Aligned sequence identity: {stats['identity']:.4f}\n"
        f"Reference coverage: {stats['reference_coverage']:.4f}\n"
        f"Mature positions mapped to modeled residues: {modeled}\n"
        f"Signal-peptide offset used for mature numbering: {args.signal_peptide_length}\n"
        "This validates chain sequence mapping, not the biological signal-peptide cleavage site.\n",
        encoding="utf-8",
    )

    print(f"Saved coordinate-based numbering map to {map_path}")
    print(f"Saved chain validation summary to {check_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
