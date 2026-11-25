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

## Render structure figures locally

After running the workflow, create a PyMOL script for the current mapped drift data and open it:

```bash
python3 scripts/07_make_pymol_figures.py --season 2024_25
pymol figures/2024_25_structure_views.pml
```

The script writes four PNG views under `figures/`. They show chains A/B/C in red/green/blue, mapped non-zero drift positions for the requested season as yellow spheres on the selected chain, and mature HA residue numbers as labels. These figures are rendered from the current `results/structural_drift.csv`, so they may differ from the older images on the published site.

To inspect or adjust the generated PyMOL commands before rendering, open `figures/2024_25_structure_views.pml` in a text editor. If labels still overlap, generate smaller labels:

```bash
python3 scripts/07_make_pymol_figures.py --season 2024_25 --label-size 8
```

Run the generated `.pml` with PyMOL again to regenerate the PNGs.

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

### Next structural step: consider another PDB candidate

The current `7KNA` scaffold should be replaced or supplemented if a more appropriate H1 structure is available. [PDB 9GSP](https://www.rcsb.org/structure/9GSP) is a stronger current candidate than 7KNA: it is a full-length H1 hemagglutinin with 96.6% identity to the Missouri reference and 488 modeled CA-bearing residues in chain A.

9GSP is A/Victoria/2570/2019, not an exact A/Missouri/11/2025 structure. It has 566 residues and differs from Missouri at 19 positions, so the replacement should preserve `pdb_aa`, `identity_match`, and unmodeled-position reporting rather than treating the PDB sequence as the Missouri sequence.

Before changing `pdb_id` in `main.nf`, verify:

- the candidate matches H1N1 and the intended HA chain,
- sequence identity and coverage are at least as good as or better than 7KNA,
- the modeled chain and residue numbering are complete enough for the drift positions,
- the structure is experimental and not merely a computed model, unless a validated homology model is the explicit goal,
- 9GSP or another candidate improves on 7KNA rather than introducing worse gaps or mismatches.

If no exact Missouri structure is found, document that the selected PDB is a related H1 structural scaffold and keep all biological claims exploratory.
