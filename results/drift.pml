load 7KNV.pdb, prot
select chainA, chain A
hide everything, all
show cartoon, chainA
show sticks, chainA
color white, chainA

# NOTE: populate exact antigenic site residue ranges from a validated HA epitope map.
# For now this script keeps the region-based color logic and records the trend.

select 2025_26, chainA and resi 217, 218, 221, 222, 224, 225, 226, 227, 228, 229, 231, 234, 235, 236, 237, 238, 239, 240, 241, 242, 244, 245, 246, 247, 248, 249, 250, 254, 255, 256, 258, 260, 261
color blue, 2025_26

select head, chainA and resi 1-300
select stem, chainA and resi 301-550
color lightblue, head
color orange, stem

set cartoon_transparency, 0.2
set sphere_scale, 0.25
dss
