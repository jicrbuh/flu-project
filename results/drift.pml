fetch 7KNA, async=0
hide everything, all
show cartoon, 7KNA
color white, 7KNA

# Selections include only above-threshold residues modeled in the validated chain mapping.
# Positions absent from the model remain listed in structural_drift.csv but are not selected.

select s202122, 7KNA and chain A and name CA and resi 29+64+106+123+130+139A+146+151+165+170+182+195+196+216+225+233+259+268+269+278+285+316+522+538+565+584+617+620+626
show spheres, s202122
color red, s202122
print('count_atoms s202122: %d' % cmd.count_atoms('s202122'))
# 2021_22: 29 mapped positions; 3 above-threshold positions absent

select s202223, 7KNA and chain A and name CA and resi 29+58+104+118+130+144+146+151+161+166+203+225+269+278+282+285+316+329+522+538+584+593+617+626+646+663
show spheres, s202223
color blue, s202223
print('count_atoms s202223: %d' % cmd.count_atoms('s202223'))
# 2022_23: 26 mapped positions; 1 above-threshold positions absent

select s202324, 7KNA and chain A and name CA and resi 48+57+80+96+130+131+146+151+178+206+218+225+269+285+522+538+611+616+617+626+638+646
show spheres, s202324
color green, s202324
print('count_atoms s202324: %d' % cmd.count_atoms('s202324'))
# 2023_24: 22 mapped positions; 3 above-threshold positions absent

select s202425, 7KNA and chain A and name CA and resi 45+55+64+78+93+123+130+146+148+151+175+178+199+225+269+270+285+291+310+522+538+617+626
show spheres, s202425
color yellow, s202425
print('count_atoms s202425: %d' % cmd.count_atoms('s202425'))
# 2024_25: 23 mapped positions; 3 above-threshold positions absent

select s202526, 7KNA and chain A and name CA and resi 57+80+123+138+148+161+164+166+178+214+232+268+269+290+291+310+578
show spheres, s202526
color magenta, s202526
print('count_atoms s202526: %d' % cmd.count_atoms('s202526'))
# 2025_26: 17 mapped positions; 0 above-threshold positions absent

set cartoon_transparency, 0.2
set sphere_scale, 0.6
dss
