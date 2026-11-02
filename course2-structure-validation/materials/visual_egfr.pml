# visual_egfr.pml - teacher demo script (course 2, task card 3)
# Run in PyMOL:   @visual_egfr.pml
# Requires: 1M17.pdb loaded OR run from the materials directory.
# Teacher demo only. Students reproduce the same views step by step.

# ---------------------------------------------------------------
# View 1: domain coloring - kinase domain highlighted
# ---------------------------------------------------------------
load 1M17.pdb
hide everything
show cartoon, 1M17
color grey80, 1M17
select kinase, 1M17 and polymer.protein
select erlotinib, 1M17 and resn AQ4
color slate, kinase
show sticks, erlotinib
color orange, erlotinib
zoom kinase
set ray_opaque_background, 1
png view1_domains.png, width=1600, height=1200, dpi=300, ray=1

# ---------------------------------------------------------------
# View 2: pocket surface half-wrapping the ligand
# ---------------------------------------------------------------
hide everything
show cartoon, 1M17
select pocket_res, byres ((1M17 and polymer.protein) within 6 of erlotinib)
show surface, pocket_res
color white, pocket_res
set transparency, 0.5
show sticks, erlotinib
color orange, erlotinib
zoom erlotinib
viewport 1200, 900
png view2_pocket_surface.png, width=1600, height=1200, dpi=300, ray=1

# ---------------------------------------------------------------
# View 3: key residue T790 labelled near the ligand
# ---------------------------------------------------------------
hide everything
show cartoon, 1M17
select t790, 1M17 and chain A and resi 766 and resn THR
show sticks, t790
color magenta, t790
label t790 and name CA, "T790 (1M17: T766)"
set label_size, 26
set label_color, black
show sticks, erlotinib
color orange, erlotinib
zoom (t790 or erlotinib), 8
distance dist_t790_lig, t790, erlotinib
hide labels, dist_t790_lig
png view3_residue_t790.png, width=1600, height=1200, dpi=300, ray=1

print "Demo script finished: view1_domains.png / view2_pocket_surface.png / view3_residue_t790.png written."