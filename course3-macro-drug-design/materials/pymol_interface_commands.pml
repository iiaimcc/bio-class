# pymol_interface_commands.pml - course 3 interface analysis on 1YY9
# Teacher: verify chain ids first (1YY9: A=EGFR ectodomain, C=Fab light, D=Fab heavy;
# B/E/F/G are N-glycan chains, not protein).
# Usage: run in PyMOL after `load 1YY9.pdb` (or @this file inside a session).
# CDR ranges below are IMGT, converted to PDB resi (see CDR标注表.md).

# ---- setup ----
# Self-contained: load the structure if it is not already in the session.
load 1YY9.pdb, 1YY9
hide everything
show cartoon, 1YY9
color marine, chain A
color palecyan, chain C
color salmon, chain D

# ---- CDR highlight (IMGT -> PDB resi) ----
select cdrH, chain D and resi 26-33+51-57+96-108
select cdrL, chain C and resi 27-32+50-52+89-97
color red, cdrH
color orange, cdrL
zoom (cdrH or cdrL)
png cdr_view.png, width=1600, height=1200, dpi=300, ray=1

# ---- interface contacts (4.5 A) ----
select egfr_chain, chain A and polymer.protein
select antibody, (chain C or chain D) and polymer.protein
select absurf, antibody within 4.5 of egfr_chain
select agsurf, egfr_chain within 4.5 of antibody
# Whole contacting residues: CA atoms of many contacting residues are >4.5 A
# from the partner, so expand with byres before listing.
select absurf_res, byres absurf
select agsurf_res, byres agsurf
show sticks, (absurf_res or agsurf_res) and (name CB or name CA)
# count + list contacting residues (21 antibody residues, 22 EGFR residues):
count_atoms absurf_res
count_atoms agsurf_res
iterate absurf_res and name CA, print("%s%s-%s" % (chain, resn, resi))
iterate agsurf_res and name CA, print("%s%s-%s" % (chain, resn, resi))

# ---- H-bonds: measure a few candidate polar pairs after visual check ----
# distance hb1, chain D and resi 100 and name OD1, chain A and resi 467 and name NE2
# (update atom/residue pairs to what you observe; then:)
# show labels, hb1

# ---- sequence positions export (for mutation mapping) ----
# select vh, chain D and resi 1-119
# iterate vh and name CA, print("VH %3d -> %s" % (resi, res1))
