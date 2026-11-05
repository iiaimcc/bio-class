"""Compare the Boltz-predicted EGFR DIII + cetuximab Fab complex with the experimental 1YY9."""
from pymol import cmd
import json

PRED = "/home/cluster1/caoji/zhesheng1/work/boltz_ref/out/boltz_results_boltz_egfr_cetuximab_msa/predictions/boltz_egfr_cetuximab_msa/boltz_egfr_cetuximab_msa_model_0.pdb"
REF = "/home/cluster1/caoji/zhesheng1/research/bio-class/course3-macro-drug-design/materials/fallback/inputs/1YY9.pdb"

cmd.reinitialize()
cmd.load(REF, "ref")
cmd.load(PRED, "pred")

# Prediction chains: A = DIII (1..171), C = VL (1..107), D = VH (1..119).
# 1YY9 numbering: A 310-480 = the same DIII, C/D Fab chains.
cmd.alter("pred and chain A", "resi=str(int(resi)+309)")
cmd.sort()

rms = cmd.align("pred and chain A and polymer.protein and name CA",
                "ref and chain A and polymer.protein and name CA")
print(f"align chain A (DIII): RMSD={rms[0]:.3f} A over {rms[1]} CA atoms")

for label, sel_ref, sel_pred in [
    ("VL", "ref and chain C and polymer.protein", "pred and chain C and polymer.protein"),
    ("VH", "ref and chain D and polymer.protein", "pred and chain D and polymer.protein"),
]:
    r = cmd.align(f"{sel_pred} and name CA and resi 1-107", f"{sel_ref} and name CA and resi 1-107")
    print(f"align {label}: RMSD={r[0]:.3f} A over {r[1]} CA atoms")


def interface(prefix, fab_chains, target_chain):
    sel = f"({prefix} and chain {target_chain} and polymer.protein) within 4.5 of ({prefix} and chain {fab_chains} and polymer.protein)"
    atoms = cmd.get_model(sel).atom
    residues = sorted({int(a.resi) for a in atoms})
    return residues, atoms


ref_res, ref_atoms = interface("ref", "C+D", "A")
pred_res, pred_atoms = interface("pred", "C+D", "A")
print(f"interface residues on EGFR: ref={len(ref_res)} pred={len(pred_res)}")
print(f"shared: {len(set(ref_res) & set(pred_res))}")
print(f"ref only: {sorted(set(ref_res) - set(pred_res))}")
print(f"pred only: {sorted(set(pred_res) - set(ref_res))}")
print(f"interface atom pairs: ref={len(ref_atoms)} pred={len(pred_atoms)}")

summary = {
    "prediction": "Boltz-2 (boltz 2.2.1, checkpoint default)",
    "input": "EGFR domain III (P00533 334-504) + cetuximab VH (1YY9 D 1-119) + VL (1YY9 C 1-107)",
    "msa": "ColabFold MMseqs2 server, pairs greedy; A 3090 / C 10728 / D 8718 sequences after merge",
    "chains": {"A": "EGFR DIII", "C": "cetuximab VL", "D": "cetuximab VH"},
    "confidence": json.load(open("/home/cluster1/caoji/zhesheng1/work/boltz_ref/out/boltz_results_boltz_egfr_cetuximab_msa/predictions/boltz_egfr_cetuximab_msa/confidence_boltz_egfr_cetuximab_msa_model_0.json")),
    "align_vs_1YY9": {
        "EGFR_DIII_CA_rmsd_A": round(rms[0], 3),
        "EGFR_DIII_CA_atoms": rms[1],
    },
    "interface_4p5A": {
        "reference_residues_on_EGFR": ref_res,
        "predicted_residues_on_EGFR": pred_res,
        "shared": sorted(set(ref_res) & set(pred_res)),
        "reference_atom_count": len(ref_atoms),
        "predicted_atom_count": len(pred_atoms),
    },
}
with open("/home/cluster1/caoji/zhesheng1/work/boltz_ref/compare_summary.json", "w") as fh:
    json.dump(summary, fh, indent=1)

cmd.bg_color("white")
cmd.set("ray_opaque_background", 0)
cmd.show("cartoon", "ref")
cmd.color("grey70", "ref")
cmd.show("cartoon", "pred")
cmd.color("skyblue", "pred and chain A")
cmd.color("orange", "pred and chain C")
cmd.color("marine", "pred and chain D")
cmd.orient("pred")
cmd.ray(1200, 900)
cmd.png("/home/cluster1/caoji/zhesheng1/work/boltz_ref/align_vs_1YY9_global.png", dpi=150)
cmd.orient("pred and chain A within 8 of (pred and chain C+D)")
cmd.ray(1200, 900)
cmd.png("/home/cluster1/caoji/zhesheng1/work/boltz_ref/align_vs_1YY9_interface.png", dpi=150)
print("wrote PNGs and compare_summary.json")
