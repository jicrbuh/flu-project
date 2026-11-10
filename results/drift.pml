load 7KNV.pdb, prot
select chainA, chain A
hide everything, all
show cartoon, chainA
show sticks, chainA
color white, chainA

# NOTE: populate exact antigenic site residue ranges from a validated HA epitope map.
# For now this script keeps the region-based color logic and records the trend.

select 2025_26, chainA and resi 1, 2, 5, 6, 8, 9, 10, 11, 12, 13, 15, 18, 19, 20, 21, 22, 23, 24, 25, 26, 28, 29, 30, 31, 32, 33, 34, 38, 39, 40, 42, 44, 45
color purple, 2025_26

select head, chainA and resi 1-300
select stem, chainA and resi 301-550
color lightblue, head
color orange, stem

set cartoon_transparency, 0.2
set sphere_scale, 0.25
dss
