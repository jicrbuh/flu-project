# Flu HA Drift Project Checklist

Goal: map hemagglutinin (HA) sequence changes across seasons onto a 3D structure, check whether they cluster at antibody epitopes (head) vs conserved regions (stem). Same skills as the pembrolizumab project, antigen side.

Related: [Compugen role](https://cgen.com/careers/career-detail/default.aspx?uid=F3.27F)

## Background
- [x] Read WHO 2026-27 NH recommendation (27 Feb 2026) — **done** — strains: A/Missouri/11/2025 (H1N1)pdm09, A/Darwin/1454/2025 [egg] or A/Darwin/1415/2025 [non-egg] (H3N2), B/Tokyo/EIS13-175/2025 [egg] or B/Pennsylvania/14/2025 [non-egg] (B/Victoria)
- [x] Understand egg vs cell vs recombinant vs mRNA platforms (egg adaptation mutations in HA explain why strains differ by platform)
- [ ] Understand HA basics: head vs stem, receptor-binding site, antigenic sites, why HA drifts
  - [ ] Why it's called hemagglutinin (red cell clumping assay) vs what it does (attachment + membrane fusion)
  - [ ] HA0 → HA1 + HA2 cleavage; low-pH conformational change; fusion peptide
  - [ ] Look at an HA trimer in PyMOL: color head (HA1) vs stem (HA2), find the receptor-binding site
- [x] Understand why NA is secondary (vaccine standardization is HA-based; verify NA content per vaccine)
- [ ] Read about broadly neutralizing antibodies targeting the stem

## Data
- [ ] GISAID EpiFlu: optional. Registration needs identity verification + Database Access Agreement; webmail addresses may be slowed. DAA forbids redistributing GISAID sequences (no raw sequences in a public GitHub repo — scripts + EPI_ISL ID lists + acknowledgment table only). Approval for an individual without an institution: unknown
- [ ] Pull HA sequences for the 2026-27 strains (H1N1 first)
- [ ] Pull HA sequences from past seasons from open sources first: NCBI Virus / Influenza Virus Resource (GenBank), BV-BRC. NOTE: Nextstrain seasonal flu builds use GISAID data and are restricted in re-sharing, so they are NOT a fully open source (earlier note here was wrong). Check what the tree page's 'Download Data' offers
- [ ] Check GenBank/BV-BRC coverage of the 2026-27 strains and subclade K (recent strains may be missing)
- [ ] Choose one subtype to start (H1N1pdm09 or H3N2)

## Sequence analysis
- [ ] Install MAFFT and Biopython
- [ ] Align HA sequences across seasons
- [ ] Identify variable positions; note which changed in which season
- [ ] Separate egg-adaptation changes from genuine antigenic drift where possible

## Structure mapping
- [ ] Pick an HA structure from the PDB for the chosen subtype
- [ ] Check numbering (H3 numbering convention vs PDB resi; verify with a landmark residue, same lesson as Chothia/PyMOL mismatch)
- [ ] Color variable positions on the structure in PyMOL
- [ ] Mark known antigenic sites / receptor-binding site (look up published epitope maps)
- [ ] Check: do drifting positions cluster on the head near epitopes? Is the stem conserved?
- [ ] Overlay an antibody-bound HA structure (head-binding vs stem-binding) and see which variable positions fall in the footprint

## Writeup
- [ ] Put pipeline in a GitHub repo (alignment script, mapping script, PyMOL session)
- [ ] Short README: question, data sources, method, figure, findings, caveats
- [ ] Note limitations (sampling bias in sequence databases, egg vs non-egg, one subtype only)

## Stretch
- [ ] Repeat for second subtype or B/Victoria
- [ ] Link to antibody engineering: which epitope would a therapeutic or broadly neutralizing antibody need to avoid drifting positions?
- [ ] Compare mutation density head vs stem quantitatively

## Concepts covered so far (discussion; verify details against sources)
- [x] Strain selection process (GISRS surveillance, sequencing, ferret HI, human serology, manufacturing lead time)
- [x] Platform lists: egg vs cell/recombinant/nucleic acid; why picks differ (egg adaptation)
- [x] Notation: A(H1N1)pdm09, strain naming A/Place/#/Year, B/Victoria lineage, "-like"
- [x] Trivalent vs quadrivalent (Yamagata dropped; verify)
- [x] Zoonotic H5; receptor specificity (α2,3 vs α2,6); why few human-to-human cases
- [x] Clade nomenclature: H1 6B.1A.5a.2a.1 / D.3.1; H3 2a.3a.1 / J / K (= J.2.4.1)
- [x] Why vaccines aren't 100% effective; mucosal (FluMist) vs injected; pre-existing and anti-vector immunity
- [x] Antiviral resistance surveillance (H275Y, I223V + S247N in NA; baloxavir/PA)
- [x] Passage level and SIAT/MDCK cell culture; reading an HI table

## Findings to verify / follow up
- [ ] Authoritative subclade of A/Missouri/11/2025 (one source says D.3.1, another D.3.1.1) — check WHO report or Nextstrain
- [ ] Authoritative defining HA substitutions for H3N2 subclade K (T135K, S144N vs K2N, Q173R — sources disagree)
- [ ] Does S144N create an N-glycosylation motif (depends on residue 146)? Links to liability scan idea
- [ ] H1pdm09 numbering (signal peptide removed) vs PDB residue numbers; H3 numbering convention
- [ ] Check GISAID metadata for passage history; prefer original-specimen sequences to avoid culture artifacts

## Reading to add
- [ ] Canadian SPSN paper (Eurosurveillance 2024, 5a.2a/5a.2a.1 HA substitutions modeled on PDB 7KNV in PyMOL) — closest published analogue to this project
- [ ] England early characterisation + VE paper on subclade K (autumn 2025)

## Project ideas branching from discussion
- [ ] H3N2 subclade K: map T135K / S144N (H3 numbering) on an H3 HA structure; check antigenic sites
- [ ] N1 neuraminidase: map H275Y, I223V, S247N relative to bound oseltamivir (find PDB entry)
- [ ] H5 stretch: compare avian H5 HA vs seasonal H1 head and stem; receptor-binding positions
