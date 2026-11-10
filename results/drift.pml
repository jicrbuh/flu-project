load 7KNV.pdb, prot
select chainA, chain A
hide everything, all
show cartoon, chainA
show sticks, chainA
color white, chainA

# NOTE: populate exact antigenic site residue ranges from a validated HA epitope map.
# For now this script keeps the region-based color logic and records the trend.

select 2021_22, chainA and resi 19, 54, 96, 113, 120, 130, 137, 142, 156, 161, 173, 186, 187, 207, 216, 224, 250, 259, 260, 270, 277, 308, 356, 372, 399, 418, 451, 454, 460, 506, 510, 520
color yellow, 2021_22

select 2022_23, chainA and resi 19, 39, 61, 73, 97, 120, 125, 137, 138, 141, 142, 152, 173, 184, 186, 190, 195, 216, 260, 277, 278, 308, 321, 356, 372, 411, 418, 451, 455, 460, 500, 513, 520
color red, 2022_23

select 2023_24, chainA and resi 38, 44, 45, 47, 85, 94, 96, 113, 120, 137, 142, 146, 169, 216, 260, 277, 286, 302, 326, 356, 372, 418, 427, 451, 460, 474, 520, 527, 533
color gray, 2023_24

select 2024_25, chainA and resi 36, 57, 66, 82, 83, 113, 119, 137, 141, 142, 166, 169, 186, 216, 260, 262, 277, 323, 324, 356, 372, 394, 450, 451, 455, 460, 468, 503, 510, 520
color orange, 2024_25

select 2025_26, chainA and resi 30, 36, 57, 110, 113, 139, 141, 146, 205, 260, 277, 283, 302, 333, 451, 458, 496, 500
color cyan, 2025_26

select head, chainA and resi 1-300
select stem, chainA and resi 301-550
color lightblue, head
color orange, stem

set cartoon_transparency, 0.2
set sphere_scale, 0.25
dss
