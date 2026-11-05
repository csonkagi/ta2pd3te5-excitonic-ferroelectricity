"""Assign each Wannier function to the nearest atom of its element and report the centre-atom distance.

Shows that the Te-px-seeded functions of the maximally localized model sit 1.7-2.4 A from their seed
atoms (Pd-Te bond hybrids), so no parent-atom assignment restores Pmmn (Results, Wannier-symmetry section).
"""
from pathlib import Path
import importlib.util
import re
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
W = ROOT / "data" / "wannier"

spec = importlib.util.spec_from_file_location("v6", ROOT / "code" / "ta2pd3te5_channel_mf_v6.py")
v6 = importlib.util.module_from_spec(spec); spec.loader.exec_module(v6)

Rvec, HR, degen, nw = v6.load_hr(str(W / "wannier90_hr.dat"))
L, C = v6.load_wout(str(W / "wannier90_wout_final_state.txt"), nw)
orbs = ["dz2"] * 4 + ["px"] * 10 + ["dxz"] * 2
elem = ["Ta"] * 4 + ["Te"] * 10 + ["Pd"] * 2

txt = (W / "wannier90.win").read_text()
blk = re.search(r"begin atoms_cart(.*?)end atoms_cart", txt, re.S).group(1)
atoms = [(l.split()[0], np.array(l.split()[1:4], float)) for l in blk.strip().splitlines()]

used = set()
for i in range(nw):
    best = None
    for j, (e, r) in enumerate(atoms):
        if e != elem[i] or j in used:
            continue
        d = r - C[i]; n = np.round(d / L); dd = np.linalg.norm(d - n * L)
        if best is None or dd < best[0]:
            best = (dd, j)
    used.add(best[1])
    print(f"WF{i:2d} {orbs[i]:3s} -> {elem[i]} atom {best[1]:2d}  |centre-atom| = {best[0]:.3f} A")
