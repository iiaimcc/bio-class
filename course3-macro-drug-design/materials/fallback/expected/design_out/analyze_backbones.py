"""Analyze the packaged BoltzGen backbones: python analyze_backbones.py [CIF ...].

These files use chain A for the binder and B for EGFR DIII (residues 1-171).
Contacts with the aligned crystal reference and clashes inside the design
complex are reported separately. Requires a Python environment with PyMOL.
"""
import json
import math
import sys
from pathlib import Path

from pymol import cmd

HERE = Path(__file__).resolve().parent
TARGET = HERE.parents[1] / "inputs" / "EGFR_DIII_1YY9.pdb"
EPITOPE = {349, 353, 382, 384, 408, 409, 411, 412, 415, 417, 418,
           438, 440, 441, 443, 465, 466, 467, 468, 469, 471, 473}


def atom_label(atom):
    return {"chain": atom.chain, "resi": atom.resi, "resn": atom.resn, "atom": atom.name}


def clashes(target, binder):
    target_atoms = cmd.get_model(f"({target}) within 2.0 of ({binder})").atom
    binder_atoms = cmd.get_model(f"({binder}) within 2.0 of ({target})").atom
    pairs = []
    for b in binder_atoms:
        for t in target_atoms:
            distance = math.dist(b.coord, t.coord)
            if distance < 2.0:
                pairs.append({"binder": atom_label(b), "target": atom_label(t),
                              "distance_angstrom": round(distance, 3)})
    return len(binder_atoms), pairs


def analyze(path):
    cmd.reinitialize()
    cmd.load(str(TARGET), "egfr")
    cmd.load(str(path), "bb")
    rmsd = cmd.align("bb and chain B and polymer.protein",
                     "egfr and chain A and polymer.protein")[0]
    target = "egfr and chain A and polymer.protein"
    binder = "bb and chain A and polymer.protein"
    epitope = target + " and resi " + "+".join(map(str, sorted(EPITOPE)))
    contacts = cmd.get_model(f"({target}) within 4.5 of ({binder})").atom
    residues = {int(atom.resi) for atom in contacts}
    near = cmd.count_atoms(f"({binder}) within 8 of ({epitope})")
    clash_count, pairs = clashes(target, binder)
    complex_clash_count, complex_pairs = clashes("bb and chain B and polymer.protein", binder)
    return {"target_alignment_rmsd_angstrom": round(rmsd, 4),
            "interface_target_atoms": len(contacts),
            "interface_target_residues": sorted(residues),
            "epitope_overlap": sorted(residues & EPITOPE),
            "near_epitope_binder_atoms": near,
            "clash_binder_atoms": clash_count, "clash_pairs": pairs,
            "complex_clash_binder_atoms": complex_clash_count,
            "complex_clash_pairs": complex_pairs}


if __name__ == "__main__":
    paths = [Path(p) for p in sys.argv[1:]] or sorted((HERE / "backbones").glob("design_*.cif"))
    results = {path.stem: analyze(path) for path in paths}
    print(json.dumps(results, indent=2))
