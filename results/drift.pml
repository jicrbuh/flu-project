fetch 7KNA, async=0
hide everything, all
show cartoon, 7KNA
color white, 7KNA

# Rough background domain coloring for a quick structural overview. This is a visual aid only,
# not a validated HA head/stem domain annotation.
select head, 7KNA and chain A and resi 1-300
select stem, 7KNA and chain A and resi 301-550
show surface, head
show surface, stem
color teal, head
color orange, stem
set surface_transparency, 0.75
set surface_quality, 1

# The residue list below is derived from the mature HA numbering in the drift table and
# then translated to PDB numbering using numbering_map.csv where available.
# Use the count_atoms check below to confirm the selection contains modeled residues.

select s202122, 7KNA and chain A and name CA and resi 19+54+96+113+120+130+137+142+156+161+173+186+187+207+216+224+250+259+260+270+277+308+356+372+399+418+451+454+460+506+510+520
show spheres, s202122
color cyan, s202122
print('count_atoms s202122: %d' % cmd.count_atoms('s202122'))

select s202223, 7KNA and chain A and name CA and resi 39+48+96+108+113+120+135+137+142+152+216+241+256+260+274+277+308+356+372+418+427+451+460+485+520
show spheres, s202223
color orange, s202223
print('count_atoms s202223: %d' % cmd.count_atoms('s202223'))

select s202324, 7KNA and chain A and name CA and resi 45+47+68+86+120+137+142+169+216+260+277+286+329+356+372+451+460+474+501+510+511+520+527+549
show spheres, s202324
color red, s202324
print('count_atoms s202324: %d' % cmd.count_atoms('s202324'))

select s202425, 7KNA and chain A and name CA and resi 45+83+113+137+138+142+152+166+169+216+260+262+273+277+356+372+451+454+460+472+496+510+520
show spheres, s202425
color red, s202425
print('count_atoms s202425: %d' % cmd.count_atoms('s202425'))

select s202526, 7KNA and chain A and name CA and resi 45+48+68+113+139+141+152+155+169+283+302+383+462
show spheres, s202526
color orange, s202526
print('count_atoms s202526: %d' % cmd.count_atoms('s202526'))


# Domain coloring is intentionally kept minimal and should be validated against a published
# HA domain definition before publication. The season-specific residue selections are the primary output.

set cartoon_transparency, 0.2
set sphere_scale, 0.6
dss
