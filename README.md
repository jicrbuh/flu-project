# Flu HA Drift Project

This repository contains a small workflow for exploring seasonal influenza HA sequence drift, aligning sequences to a reference, and generating summary outputs for downstream analysis and visualization.

The project combines:

- a Python-based fetch/filter pipeline for seasonal influenza HA records
- a Nextflow workflow orchestrating the steps
- Quarto-based project pages for overview and results
- output tables and a PyMOL script for structural visualization

## Repository layout

- `main.nf` — main Nextflow workflow
- `scripts/` — Python processing scripts
  - `01_fetch.py` — fetch seasonal GenBank records
  - `02_select.py` — select and filter candidate sequences
  - `03_numbering.py` — align the reference to modeled PDB-chain residues and map coordinates
  - `04_variability.py` — compute variability, drift, and per-season substitutions
  - `05_pymol.py` — generate a PyMOL script from structure-aware drift positions
  - `06_structure_report.py` — combine drift counts with the validated structure mapping
- `data/` — input/reference and fetched records
- `results/` — workflow outputs such as alignment, drift tables, and validation files
- `entries/` and `project-overview.qmd` — Quarto content for the website
- `_quarto.yml` — Quarto site configuration

## Requirements

- Python 3
- Nextflow
- Quarto (for the website)
- MAFFT (for sequence alignment)
- Biopython, pandas, and other Python dependencies in `requirements.txt`
- Optional: PyMOL for viewing generated `.pml` output

On macOS, install MAFFT with Homebrew:

```bash
brew install mafft
```

Install Python dependencies:

```bash
python3 -m pip install -r requirements.txt
```

## Configuration

Create a local environment file from the example if needed:

```bash
cp .env.example .env
```

Set your email for NCBI Entrez access in `.env`:

```bash
ENTREZ_EMAIL=you@example.com
```

If you have an NCBI API key, you can also export it for the fetch step:

```bash
export NCBI_API_KEY="your_api_key_here"
```

## Run the workflow

From the project root:

```bash
nextflow run main.nf
```

This runs the seasonal fetch/select/alignment/variability/coordinate-mapping/structure-report/PyMOL pipeline and publishes outputs into the project folders.

## Outputs

Typical outputs include:

- `data/raw/<season>.faa`
- `data/raw/<season>_meta.csv`
- `data/selected/selected_*.faa`
- `results/aln.faa`
- `results/variability.csv`
- `results/drift.csv`
- `results/drift_detail.csv` — substitution counts, sample denominators, and drift fractions by residue and season
- `results/numbering_map.csv` — reference positions mapped to actual modeled PDB residue IDs, including unmodeled positions
- `results/numbering_check.txt` — chain-alignment coverage and identity checks
- `results/structural_drift.csv` — detailed drift table joined to PDB coordinates, with absent residues marked
- `results/structure_validation.txt` — mapping/report summary and unmapped drift positions
- `results/drift.pml`
- `results/validation_summary.txt`

## Quarto site

To preview the website locally:

```bash
quarto preview
```

To render the site:

```bash
quarto render
```

## Notes

This project is meant for reproducible analysis and exploratory structural interpretation. The drift outputs are generated from sequence data, but structural interpretation should be treated with caution and validated against the actual model and residue mapping before drawing strong biological conclusions.

The structure mapping aligns the reference sequence to modeled residues (residues with CA atoms) in the configured PDB chain. It reports actual PDB residue IDs rather than assuming a constant numbering offset. Drift positions absent from the model remain in `structural_drift.csv` and are excluded from the PyMOL selections. The configured 17-residue signal-peptide offset is still a numbering convention; this structural check does not independently establish its biological cleavage site, assign antigenic sites, or prove an antigenic effect.
