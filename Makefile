PYTHON ?= python3
ENTREZ_EMAIL ?= $(shell [ -f .env ] && grep -E '^ENTREZ_EMAIL=' .env 2>/dev/null | head -n 1 | cut -d= -f2- || true)

.PHONY: fetch select align variability numbering structure-report pymol all clean

all: results/drift.pml

fetch: data/raw/2025_26.faa

select: data/selected/selected_2025_26.faa

align: results/aln.faa

variability: results/variability.csv

numbering: results/numbering_map.csv

structure-report: results/structural_drift.csv

pymol: results/drift.pml

data/raw/2025_26.faa:
	@test -n "$(ENTREZ_EMAIL)" || { echo "Set ENTREZ_EMAIL in .env before running make fetch"; exit 1; }
	mkdir -p data/raw
	$(PYTHON) scripts/01_fetch.py --season 2025_26 --email "$(ENTREZ_EMAIL)" --output-dir data/raw --meta-output data/raw/2025_26_meta.csv

data/raw/2025_26_meta.csv: data/raw/2025_26.faa
	@test -s $@

data/selected/selected_2025_26.faa: data/raw/2025_26.faa data/raw/2025_26_meta.csv scripts/02_select.py
	mkdir -p data/selected
	$(PYTHON) scripts/02_select.py --n-per-season 10 --seed 42 --min-length 566 --input-dir data/raw --output-dir data/selected

data/selected/selected_meta.csv: data/selected/selected_2025_26.faa
	@test -s $@

results/all.faa: data/ref/ref_A_Missouri_11_2025.faa data/selected/selected_2025_26.faa
	mkdir -p results
	cat data/ref/ref_A_Missouri_11_2025.faa data/selected/selected_*.faa > results/all.faa

results/aln.faa: results/all.faa
	mafft --auto results/all.faa > results/aln.faa

results/variability.csv: results/aln.faa data/selected/selected_meta.csv scripts/04_variability.py
	$(PYTHON) scripts/04_variability.py --alignment results/aln.faa --meta data/selected/selected_meta.csv --signal-peptide-length 17 --output-dir results

results/drift_detail.csv: results/variability.csv
	@test -s results/drift_detail.csv

results/numbering_map.csv: data/ref/ref_A_Missouri_11_2025.faa scripts/03_numbering.py
	$(PYTHON) scripts/03_numbering.py --reference data/ref/ref_A_Missouri_11_2025.faa --pdb-id 7KNA --chain A --signal-peptide-length 17 --output-dir results

results/structural_drift.csv: results/drift_detail.csv results/numbering_map.csv scripts/06_structure_report.py
	$(PYTHON) scripts/06_structure_report.py --drift-detail results/drift_detail.csv --numbering-map results/numbering_map.csv --pdb-id 7KNA --chain A --output-dir results

results/drift.pml: results/structural_drift.csv scripts/05_pymol.py
	$(PYTHON) scripts/05_pymol.py --drift-csv results/structural_drift.csv --pdb-id 7KNA --chain A --output-dir results

clean:
	rm -rf data/raw data/selected results
