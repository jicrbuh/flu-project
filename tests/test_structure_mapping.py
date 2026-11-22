import sys
import unittest
from importlib import import_module
from pathlib import Path

import pandas as pd
from Bio.PDB import Residue

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

build_numbering_map = import_module("03_numbering").build_numbering_map
build_pml = import_module("05_pymol").build_pml
build_structure_report = import_module("06_structure_report").build_structure_report


THREE_LETTER = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
}


def make_residue(amino_acid: str, residue_number: int) -> Residue.Residue:
    return Residue.Residue((" ", residue_number, " "), THREE_LETTER[amino_acid], " ")


class StructureMappingTests(unittest.TestCase):
    def test_reference_positions_use_actual_structure_residue_ids(self):
        residues = [
            make_residue(aa, residue_number)
            for aa, residue_number in zip("CDEFG", [101, 102, 105, 106, 107])
        ]

        rows, stats = build_numbering_map(
            "AACDEFGH",
            residues,
            chain_id="A",
            signal_peptide_length=2,
            min_aligned=5,
            min_identity=1.0,
        )

        by_reference = {row["ref_index"]: row for row in rows}
        self.assertEqual(by_reference[3]["pdb_resi"], "101")
        self.assertEqual(by_reference[7]["pdb_resi"], "107")
        self.assertEqual(by_reference[8]["mapping_status"], "not_modeled")
        self.assertEqual(stats["aligned_residues"], 5)

    def test_structure_report_keeps_unmodeled_positions_explicit(self):
        drift = pd.DataFrame(
            [
                {
                    "mature": 1,
                    "ref_index": 18,
                    "ref_aa": "A",
                    "season": "2021_22",
                    "sequences_in_season": 10,
                    "sequences_compared": 10,
                    "sequences_different": 2,
                    "fraction_different": 0.2,
                    "substitutions": "A>V:2",
                },
                {
                    "mature": 2,
                    "ref_index": 19,
                    "ref_aa": "C",
                    "season": "2021_22",
                    "sequences_in_season": 10,
                    "sequences_compared": 10,
                    "sequences_different": 1,
                    "fraction_different": 0.1,
                    "substitutions": "C>Y:1",
                },
            ]
        )
        numbering = pd.DataFrame(
            [
                {
                    "ref_index": 18,
                    "mature_index": 1,
                    "pdb_chain": "A",
                    "pdb_resi": "11",
                    "pdb_aa": "A",
                    "identity_match": "yes",
                    "mapping_status": "modeled",
                },
                {
                    "ref_index": 19,
                    "mature_index": 2,
                    "pdb_chain": "A",
                    "pdb_resi": "",
                    "pdb_aa": "",
                    "identity_match": "",
                    "mapping_status": "not_modeled",
                },
            ]
        )

        report, stats = build_structure_report(drift, numbering)

        self.assertEqual(report["structure_status"].tolist(), ["modeled", "not_modeled"])
        self.assertEqual(stats["unmodeled_positions"], [2])

    def test_pymol_selection_does_not_guess_numbers_for_unmodeled_sites(self):
        report = pd.DataFrame(
            [
                {
                    "season": "2021_22",
                    "fraction_different": 0.2,
                    "pdb_chain": "A",
                    "pdb_resi": "11",
                    "structure_status": "modeled",
                },
                {
                    "season": "2021_22",
                    "fraction_different": 0.1,
                    "pdb_chain": "A",
                    "pdb_resi": "",
                    "structure_status": "not_modeled",
                },
            ]
        )

        script = build_pml(report, pdb_id="7KNA", chain="A")

        self.assertIn("resi 11", script)
        self.assertNotIn("resi 2", script)
        self.assertIn("1 above-threshold positions absent", script)


if __name__ == "__main__":
    unittest.main()
