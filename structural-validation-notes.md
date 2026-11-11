# Structural validation notes for the HA drift view

This note records the main issues we found while validating the PyMOL output against the 7KNA structure and the drift table.

## 1. The drift data is real, but the visual is not yet a validated biological result

The drift file is generated and non-empty. We confirmed that `results/drift.csv` contains real seasonal drift values and multiple nonzero per-season positions.

That means the underlying per-residue drift analysis is operationally valid. However, the downstream structural view is a rendering layer and must still be validated against the actual modeled structure.

## 2. PyMOL selection syntax was wrong

The generated PyMOL script originally used comma-separated residue syntax such as:

```
resi 19, 54, 96
```

PyMOL expects residue lists written as:

```
resi 19+54+96
```

This was a real syntax bug. It can silently drop many residues or fail to select them as intended.

## 3. Selections were not restricted to CA atoms

The script originally selected whole residues rather than alpha-carbons. For a clear atom-level visualization, the selection should be restricted to `name CA`:

```
7KNA and chain A and name CA and resi 19+54+96
```

This matters because the expected number of visible atoms should match the actual modeled atoms in the structure, and not all drift residues are necessarily modeled.

## 4. The seasonal colors were overwritten by the final domain coloring

The broad domain coloring step was applied after the drift selections, so the HA head and stem color assignments could cover the same atoms and hide the season-specific drift colors.

This is a rendering issue, not a drift computation issue. The season-specific drift colors should be the primary viewpoint and any broad domain coloring should be used carefully and clearly as a background annotation.

## 5. The head/stem split was not biologically validated

The original split:

```
head = resi 1-300
stem = resi 301-550
```

is too crude for a real HA structural interpretation. HA domain boundaries and the PDB numbering layout do not necessarily align with a simple 1–300 / 301–550 split. This was only a rough visual helper, not a scientifically supported domain assignment.

## 6. The PDB structure is not a simple HA-only ectodomain in a single clean chain

The live PyMOL check showed warnings such as:

```
AssignSS-Warning: Ignoring incomplete residue /7KNA/C/C/795 ...
```

This indicates that the structure includes incomplete residues and that chain C continues beyond the expected HA mature residue range. That means the model is not a clean single-chain HA ectodomain with a trivial 1–550 numbering map.

This is important because the script assumes the residue numbering corresponds cleanly to the mature HA residues used in the drift analysis.

## 7. The numbering map must be validated against the actual model before trusting any residue position

The drift table is in mature HA numbering, and the structural rendering only makes sense if the numbers are correctly translated to PDB numbering. The mapping in `results/numbering_map.csv` must be checked against the atomic model.

This should be done by comparing:

- the residues in `7KNA and chain A and name CA`
- the translated mature-to-PDB mapping in `numbering_map.csv`
- the expected Missouri residue identities at the relevant drift positions

Until that identity check is done, the structure rendering should be treated as a visual scaffold, not as a confirmed biological mapping.

## 8. Many residue positions are not actually present in the modeled structure

The live PyMOL CA counts for the seasonal selections were lower than the number of positions in the drift table. For example, the season selection counts were roughly 25, 18, 19, 17, and 11 CA atoms instead of the full drift counts in the CSV.

This tells us that many of the drift positions are either:

- not in the PDB model,
- not present in the selected chain,
- not present as modeled CA atoms,
- or not part of the HA ectodomain represented by this PDB entry.

A structural visualization must account for that before it is interpreted as a biological distribution.

## 9. The view is useful as a check, but not as a result

The current PyMOL view is useful for:

- confirming that the script is executing,
- confirming that residue selections are being created,
- checking whether the selected positions fall on the modeled HA body,
- identifying that some positions are in the head and others in the stem-like region.

It is not yet enough to support a publication claim about drift distribution without:

- an explicit chain/domain validation,
- a residue-by-residue match to the numbered model,
- a clear exclusion of non-modeled positions,
- and a proper domain definition.

## 10. Recommended next step

The right next move is not to call the PyMOL result “done.” Instead, the next step is to create a strict structural validation: 

1. print the CA range of the relevant chain,
2. compare it against the mature numbering map,
3. count atoms in each seasonal selection,
4. exclude or mark residues that are not modeled,
5. then generate the final visual only after the model identity is confirmed.

This is the needed epistemic check before the draft should present the structural picture as a real biological result.
