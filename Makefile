PYTHON ?= python3
ENTREZ_EMAIL ?= $(shell [ -f .env ] && grep -E '^ENTREZ_EMAIL=' .env 2>/dev/null | head -n 1 | cut -d= -f2- || true)

.PHONY: fetch select align variability pymol all numbering clean

all: variability pymol

# Convenience target for the numbering workstream that starts the project logic.
numbering:
	$(PYTHON) scripts/03_numbering.py --reference data/ref/ref_A_Missouri_11_2025.faa --pdb-id 7KNV --chain A

fetch: data/raw/2025-26.faa

select: data/selected/2025-26.faa

align: results/aln.faa

variability: results/variability.csv

pymol: results/drift.pml

# This is the actual fetch step, built as a file target so make only reruns when needed.
data/raw/2025-26.faa: scripts/01_fetch.py
	@test -n "$(ENTREZ_EMAIL)" || { echo "Set ENTREZ_EMAIL in .env before running make fetch"; exit 1; }
	$(PYTHON) scripts/01_fetch.py --season 2025-26 --term "Influenza A H1N1 hemagglutinin" --email "$(ENTREZ_EMAIL)"
	@mkdir -p data/raw

# Select representative sequences after fetch completes.
data/selected/2025-26.faa: data/raw/2025-26.faa scripts/02_select.py
	$(PYTHON) scripts/02_select.py --n-per-season 10 --seed 42 --min-length 200

# Build reference + selected alignment.
results/all.faa: data/ref/*.faa data/selected/2025-26.faa
	mkdir -p results
	cat data/ref/*.faa data/selected/*.faa > results/all.faa

results/aln.faa: results/all.faa
	mafft --auto results/all.faa > results/aln.faa

results/variability.csv: results/aln.faa scripts/04_variability.py
	$(PYTHON) scripts/04_variability.py

results/drift.pml: results/variability.csv scripts/05_pymol.py
	$(PYTHON) scripts/05_pymol.py

clean:
	rm -rf data/raw data/selected results
