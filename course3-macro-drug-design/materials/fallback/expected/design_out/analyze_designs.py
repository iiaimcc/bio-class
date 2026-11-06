"""Interface analysis for the selected 80-110 aa BoltzGen designs (long family).

Usage:
    python analyze_designs.py [design_out_dir]

Defaults to this file's own directory, so the script runs from a fresh clone.
Needs PyMOL (pymol-open-source) on the PATH or importable as a module.
"""
import csv
import json
import os
import sys

from pymol import cmd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = sys.argv[1] if len(sys.argv) > 1 else HERE

# Cetuximab epitope on EGFR (4.5 A contacts, PDB numbering), computed from 1YY9.
EPITOPE = [349, 353, 382, 384, 408, 409, 411, 412, 415, 417, 418,
           438, 440, 441, 443, 465, 466, 467, 468, 469, 471, 473]
# Target numbering inside the BoltzGen runs: 1-171 == 1YY9 PDB 310-480.
PDB_OFFSET = 309
DESIGNS = ["14", "29", "86", "99"]


def analyse(path):
    cmd.reinitialize()
    cmd.load(path, "cx")
    target = binder = None
    for ch in cmd.get_chains("cx"):
        seq = cmd.get_fastastr(f"cx and chain {ch}").splitlines()[1]
        if seq.startswith("RKVCNG"):
            target = ch
        else:
            binder = ch
    if target is None or binder is None:
        raise RuntimeError(f"could not identify chains in {path}: {cmd.get_chains('cx')}")
    atoms = cmd.get_model(f"(cx and chain {target}) within 4.5 of (cx and chain {binder})").atom
    resis = sorted({int(a.resi) for a in atoms})
    pdb_num = [r + PDB_OFFSET for r in resis]
    shared = sorted(set(pdb_num) & set(EPITOPE))
    return {
        "chains": {"binder": binder, "target": target},
        "n_interface_residues_target_side": len(resis),
        "n_interface_atoms_target_side": len(atoms),
        "interface_residues_pdb_numbering": pdb_num,
        "epitope_overlap": shared,
        "epitope_overlap_fraction": round(len(shared) / len(EPITOPE), 3),
    }


def main():
    summary = {}
    for idx in DESIGNS:
        path = os.path.join(OUT_DIR, "refolded", f"refold_long_{idx}.cif")
        summary[f"long_{idx}"] = analyse(path)
        s = summary[f"long_{idx}"]
        print(f"long_{idx}: {s['n_interface_residues_target_side']} residues, "
              f"{s['n_interface_atoms_target_side']} atoms, "
              f"epitope overlap {len(s['epitope_overlap'])}/{len(EPITOPE)}")

    with open(os.path.join(OUT_DIR, "interface_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)

    table = os.path.join(OUT_DIR, "sequence_candidates.csv")
    if os.path.exists(table):
        with open(table) as fh:
            rows = {r["design_id"]: r for r in csv.DictReader(fh)}
        for key, s in summary.items():
            row = rows.get(key)
            if row:
                row["interface_residues_target_side"] = s["n_interface_residues_target_side"]
                row["epitope_overlap_of_22"] = len(s["epitope_overlap"])
        with open(table, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(next(iter(rows.values())).keys()))
            writer.writeheader()
            writer.writerows(rows.values())
    print("wrote interface_summary.json")


if __name__ == "__main__":
    main()
