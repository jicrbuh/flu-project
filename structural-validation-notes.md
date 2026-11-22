# Structural validation notes for the HA drift view

## Current mapping workflow

The workflow aligns the reference HA sequence to residues with modeled CA atoms in the configured PDB chain. It records the chain's actual author residue IDs, amino-acid identities, and whether each reference position is modeled. It does not infer a constant offset from the first aligned residue.

For the current `7KNA`, chain `A` run:

- 483 reference residues align to modeled chain residues at 94.2% amino-acid identity.
- The first mature reference residue (`ref_index=18`) maps to PDB residue `A:11`.
- Of the 68 drift positions, 60 have modeled residues and 8 do not.
- The 8 unmodeled drift positions are retained in `structural_drift.csv`, but omitted from PyMOL selections.
- The model and the Missouri reference differ at 20 of the mapped drift positions; `pdb_aa` and `identity_match` make those differences visible in the table.

These are sequence-to-coordinate mapping checks, not evidence that a drift substitution changes antibody recognition.

## Generated outputs

- `results/numbering_map.csv` maps each mature reference position to its PDB chain/residue ID, or marks it `not_modeled`.
- `results/numbering_check.txt` reports aligned-residue count, sequence identity, and reference coverage. The configured minimums are 250 aligned residues and 80% identity.
- `results/drift_detail.csv` reports per-season substitutions, counts, fractions, and denominators for drift positions.
- `results/structural_drift.csv` joins the detailed drift table to the coordinate map and retains unmodeled positions explicitly.
- `results/drift.pml` selects only above-threshold residues that are modeled in the configured chain and restricts selections to CA atoms.
- `results/structure_validation.txt` summarizes mapped and unmodeled drift positions.

## Interpretation limits

- The project uses a 17-residue signal-peptide offset for mature indexing. The mapping workflow applies that convention consistently but does not independently validate the biological cleavage site.
- The PDB model is not the same isolate as the Missouri reference. A `pdb_aa` mismatch is a difference between the structure and reference, not a drift substitution in the seasonal samples.
- No published antigenic-site map is applied. Structural proximity or sequence variability alone does not establish an antigenic effect.
- The seasonal sample is capped at 10 sequences per season, so reported fractions are based on a small sample and should be read with their denominators.

Before making antigenic claims, select an experimentally supported H1 epitope map, translate its numbering to the project reference, and support functional claims with relevant antibody or neutralization evidence.
