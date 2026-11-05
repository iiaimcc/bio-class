"""Interface analysis for selected BoltzGen designs against the EGFR DIII epitope."""
import json
from pymol import cmd

REF = "/home/cluster1/caoji/zhesheng1/research/bio-class/course3-macro-drug-design/materials/fallback/inputs/1YY9.pdb"
# cetuximab epitope on EGFR (4.5 A, computed earlier; PDB numbering)
EPITOPE = [349,353,382,384,408,409,411,412,415,417,418,438,440,441,443,465,466,467,468,469,471,473]
BASE = "out_long/intermediate_designs_inverse_folded/refold_cif/egfr_diii_binder_long_%s.cif"
DESIGNS = {"long_14": "14", "long_29": "29", "long_86": "86", "long_99": "99"}

out = {}
for name, idx in DESIGNS.items():
    cmd.reinitialize()
    cmd.load(BASE % idx, "cx")
    chains = cmd.get_chains("cx")
    # target = the chain whose sequence starts with RKVCNG (EGFR DIII); binder = the other
    target, binder = None, None
    for ch in chains:
        seq = cmd.get_fastastr(f"cx and chain {ch}").splitlines()[1]
        if seq.startswith("RKVCNG"):
            target = ch
        else:
            binder = ch
    sel = f"(cx and chain {target}) within 4.5 of (cx and chain {binder})"
    atoms = cmd.get_model(sel).atom
    resis = sorted({int(a.resi) for a in atoms})
    pdb_num = [r + 309 for r in resis]           # input 1..171 = PDB 310..480
    shared = sorted(set(pdb_num) & set(EPITOPE))
    out[name] = {
        "chains": {"binder": binder, "target": target},
        "n_interface_residues_target_side": len(resis),
        "n_interface_atoms_target_side": len(atoms),
        "interface_residues_pdb_numbering": pdb_num,
        "epitope_overlap": shared,
        "epitope_overlap_fraction": round(len(shared) / len(EPITOPE), 3),
    }
    print(name, out[name]["n_interface_residues_target_side"], "residues,",
          out[name]["n_interface_atoms_target_side"], "atoms,",
          "overlap", len(shared), "/", len(EPITOPE))
    if name == "long_14":
        cmd.bg_color("white")
        cmd.set("ray_opaque_background", 0)
        cmd.hide("everything")
        cmd.show("cartoon", "cx")
        cmd.color("grey70", f"cx and chain {target}")
        cmd.color("orange", f"cx and chain {binder}")
        cmd.color("red", f"cx and chain {target} and resi " + "+".join(str(r) for r in shared))
        cmd.orient("cx")
        cmd.ray(1400, 1000)
        cmd.png("/home/cluster1/caoji/zhesheng1/work/boltzgen_egfr/design_long14_interface.png", dpi=150)

json.dump(out, open("/home/cluster1/caoji/zhesheng1/work/boltzgen_egfr/design_interface_summary.json", "w"), indent=1)
print("wrote summary")
