load 7KNV.pdb, prot
select chainA, chain A
hide everything, all
show cartoon, chainA
show sticks, chainA
color white, chainA

# NOTE: populate exact antigenic site residue ranges from a validated HA epitope map.
# For now this script keeps the region-based color logic and records the trend.

select 2025_26, chainA and resi 63, 64, 67, 68, 70, 71, 72, 73, 74, 75, 77, 80, 81, 82, 83, 84, 85, 86, 87, 88, 90, 91, 92, 93, 94, 95, 96, 100, 101, 102, 104, 106, 107
color magenta, 2025_26

select head, chainA and resi 1-300
select stem, chainA and resi 301-550
color lightblue, head
color orange, stem

set cartoon_transparency, 0.2
set sphere_scale, 0.25
dss
