"""Interface analysis for the 60-90 aa family (main run)."""
import json
from pymol import cmd

EPITOPE = [349,353,382,384,408,409,411,412,415,417,418,438,440,441,443,465,466,467,468,469,471,473]
BASE = "out_main/intermediate_designs_inverse_folded/refold_cif/egfr_diii_binder_%s.cif"
DESIGNS = {"main_082": "082", "main_078": "078", "main_020": "020", "main_101": "101"}
out = {}
for name, idx in DESIGNS.items():
    cmd.reinitialize()
    cmd.load(BASE % idx, "cx")
    target = binder = None
    for ch in cmd.get_chains("cx"):
        seq = cmd.get_fastastr(f"cx and chain {ch}").splitlines()[1]
        if seq.startswith("RKVCNG"):
            target = ch
        else:
            binder = ch
    sel = f"(cx and chain {target}) within 4.5 of (cx and chain {binder})"
    atoms = cmd.get_model(sel).atom
    resis = sorted({int(a.resi) for a in atoms})
    pdb_num = [r + 309 for r in resis]
    shared = sorted(set(pdb_num) & set(EPITOPE))
    out[name] = {"chains": {"binder": binder, "target": target},
                 "n_interface_residues_target_side": len(resis),
                 "n_interface_atoms_target_side": len(atoms),
                 "interface_residues_pdb_numbering": pdb_num,
                 "epitope_overlap": shared,
                 "epitope_overlap_fraction": round(len(shared)/len(EPITOPE), 3)}
    print(name, out[name]["n_interface_residues_target_side"], "res,", out[name]["n_interface_atoms_target_side"], "atoms, overlap", len(shared), "/22")
    if name == "main_082":
        cmd.bg_color("white"); cmd.set("ray_opaque_background", 0)
        cmd.hide("everything"); cmd.show("cartoon", "cx")
        cmd.color("grey70", f"cx and chain {target}"); cmd.color("orange", f"cx and chain {binder}")
        cmd.color("red", f"cx and chain {target} and resi " + "+".join(str(r) for r in shared))
        cmd.orient("cx"); cmd.ray(1400, 1000)
        cmd.png("/home/cluster1/caoji/zhesheng1/work/boltzgen_egfr/design_main082_interface.png", dpi=150)
json.dump(out, open("design_interface_summary_main.json","w"), indent=1)
print("ok")
