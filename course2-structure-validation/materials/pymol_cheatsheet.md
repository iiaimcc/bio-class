# PyMOL cheatsheet (course 2)

One-page quick reference for the in-class tasks. Everything can also be done
via GUI buttons, but typing these in the PyMOL command line is faster.

## Basics

| Action | Command |
|---|---|
| Load structure | `load 1M17.pdb` |
| Show as cartoon | `show cartoon` |
| Show sticks for a selection | `show sticks, <sel>` |
| Select 1M17 erlotinib | `select erlotinib, 1M17 and resn AQ4` |
| Select canonical T790 in 1M17 | `select t790, 1M17 and chain A and resi 766 and resn THR` |
| Select 1M17 protein construct | `select kinase, 1M17 and chain A and polymer.protein` |
| Color | `color slate, kinase` / `color grey80, 1M17` |
| Surface | `show surface, <sel>` ; `set transparency, 0.5` |
| Label mapped residue | `label t790 and name CA, "T790 (1M17: T766)"`; `set label_size, 26` |
| Align two structures | `align target_sel, reference_sel` (RMSD printed) |
| pLDDT coloring (AlphaFold-style, 4 bins) | `color orange, <sel> and b<50` / `color yellow, <sel> and (b<70 and not b<50)` / `color skyblue, <sel> and (b<90 and not b<70)` / `color marine, <sel> and not b<90` |
| B-factor ramp (注意方向) | `spectrum b, blue_white_red, <sel>, minimum=0, maximum=100`：**低值=蓝、高值=红**，与 AlphaFold 官网配色相反，用它时务必看色标 |
| Save image | `png shot.png, dpi=300` |
| Run a script | `@visual_egfr.pml` |
| Reset view | `reset` |

## Keyboard / mouse (GUI)

- Left drag: rotate; middle or right drag: zoom/pan; Shift+middle: slab.
- `Space` toggles selection highlight; `Ctrl+click` on object: center view.

## Quick diagnosis

| Symptom | Fix |
|---|---|
| ligand not visible | `show sticks, organic`; then `color byelement, organic` |
| `resn ERT` selects 0 atoms | official 1M17 uses `resn AQ4`; confirm in RCSB/PDB `HETNAM` |
| `resi 790` is LEU | canonical T790 maps to 1M17 A:T766; verify `resn THR` |
| labels invisible | `set label_color, black` + `set label_size, 40` |
| surface covers everything | `set transparency, 0.6` |
| RMSD huge in align | check chains: align kinase domain only (`and resi 709-968`) |
| png file not written | use absolute path: `png /home/.../shot.png, dpi=300` |