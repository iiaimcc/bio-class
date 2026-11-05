load /home/cluster1/caoji/zhesheng1/research/bio-class/course3-macro-drug-design/materials/fallback/inputs/1YY9.pdb, ref
load /home/cluster1/caoji/zhesheng1/work/boltz_ref/out/boltz_results_boltz_egfr_cetuximab_msa/predictions/boltz_egfr_cetuximab_msa/boltz_egfr_cetuximab_msa_model_0.pdb, pred
alter pred and chain A, resi=str(int(resi)+309)
sort
align pred and chain A and polymer.protein and name CA, ref and chain A and polymer.protein and name CA
bg_color white
set ray_opaque_background, 0
set cartoon_transparency, 0.55, ref
color grey60, ref
color skyblue, pred and chain A
color orange, pred and chain C
color marine, pred and chain D
orient pred
ray 1400, 1000
png /home/cluster1/caoji/zhesheng1/work/boltz_ref/align_vs_1yy9_global.png, dpi=150
hide ref
orient pred and chain A within 8 of (pred and chain C+D)
ray 1400, 1000
png /home/cluster1/caoji/zhesheng1/work/boltz_ref/align_vs_1yy9_interface.png, dpi=150
