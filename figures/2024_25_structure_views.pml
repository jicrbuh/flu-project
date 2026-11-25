load "/Users/Chen/workspace/flu-project/data/structure/structure.cif", ha
hide everything, all
bg_color white
set ray_opaque_background, on
set antialias, 2
show cartoon, ha and chain A+B+C
color red, ha and chain A
color green, ha and chain B
color blue, ha and chain C
select drift_sites, ha and chain A and name CA and resi 45+55+64+78+93+123+130+146+148+151+175+178+199+225+269+270+285+291+310+522+538+617+626
show spheres, drift_sites
color yellow, drift_sites
set sphere_scale, 0.65
set label_color, black
set label_outline_color, white
set label_size, 9
set label_connector, on
label (ha and chain A and name CA and resi 45), "35"
label (ha and chain A and name CA and resi 55), "45"
label (ha and chain A and name CA and resi 64), "54"
label (ha and chain A and name CA and resi 78), "68"
label (ha and chain A and name CA and resi 93), "83"
label (ha and chain A and name CA and resi 123), "113"
label (ha and chain A and name CA and resi 130), "120"
label (ha and chain A and name CA and resi 146), "137"
label (ha and chain A and name CA and resi 148), "139"
label (ha and chain A and name CA and resi 151), "142"
label (ha and chain A and name CA and resi 175), "166"
label (ha and chain A and name CA and resi 178), "169"
label (ha and chain A and name CA and resi 199), "190"
label (ha and chain A and name CA and resi 225), "216"
label (ha and chain A and name CA and resi 269), "260"
label (ha and chain A and name CA and resi 270), "261"
label (ha and chain A and name CA and resi 285), "277"
label (ha and chain A and name CA and resi 291), "283"
label (ha and chain A and name CA and resi 310), "302"
label (ha and chain A and name CA and resi 522), "356"
label (ha and chain A and name CA and resi 538), "372"
label (ha and chain A and name CA and resi 617), "451"
label (ha and chain A and name CA and resi 626), "460"
orient ha and chain A+B+C
zoom ha and chain A+B+C, 5
png /Users/Chen/workspace/flu-project/figures/2024_25_selected_residues_angle1.png, 1600, 1200, ray=1
turn y, 90
png /Users/Chen/workspace/flu-project/figures/2024_25_selected_residues_angle2.png, 1600, 1200, ray=1
turn y, 90
png /Users/Chen/workspace/flu-project/figures/2024_25_selected_residues_angle3.png, 1600, 1200, ray=1
turn y, 90
png /Users/Chen/workspace/flu-project/figures/2024_25_selected_residues_angle4.png, 1600, 1200, ray=1
print('Rendered 23 mapped 2024_25 drift positions on chain A')
