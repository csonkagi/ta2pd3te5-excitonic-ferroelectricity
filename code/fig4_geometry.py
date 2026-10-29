"""Shared geometry for Fig. 4: atoms, projection seeds and maximally localized Wannier centres of monolayer Ta2Pd3Te5.

Coordinates are folded into one cell: x in [-a/4, 3a/4) (chain axis), y in [-b/2, b/2), z in [-c/2, c/2).
Each Wannier centre is placed in the periodic image nearest to its projection seed; seeds follow the Wannier90
projection order (Ta d_z2 on the Ta atoms, Te p_x on the Te atoms, Pd d_xz on the two listed Pd sites).
"""
import re, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "code"))
import ta2pd3te5_channel_mf_v6 as v6

HR = ROOT/"data/wannier/wannier90_hr.dat"; WOUT = ROOT/"data/wannier/wannier90_wout_final_state.txt"; WIN = ROOT/"data/wannier/wannier90.win"
Rvec, HRm, degen, nw = v6.load_hr(str(HR)); L, C = v6.load_wout(str(WOUT), nw)
_blk = re.search(r"begin atoms_cart(.*?)end atoms_cart", WIN.read_text(), re.S | re.I).group(1)
ATOMS_RAW = [(l.split()[0], np.array(l.split()[1:4], float)) for l in _blk.strip().splitlines()]
SPREADS = np.array([float(x) for x in re.findall(r"WF centre and spread\s+\d+\s+\([^)]*\)\s+([\d.]+)", WOUT.read_text())][-nw:])
COL = {"Ta": "#2f6db5", "Pd": "#9a9a9a", "Te": "#1b9e5a"}
ELW = ["Ta"]*4 + ["Te"]*10 + ["Pd"]*2

def fold(r):
    r = r.copy()
    r[0] = ((r[0] + L[0]/4) % L[0]) - L[0]/4
    r[1] = ((r[1] + L[1]/2) % L[1]) - L[1]/2
    r[2] = ((r[2] + L[2]/2) % L[2]) - L[2]/2
    return r

ATOMS = [(e, fold(r)) for e, r in ATOMS_RAW]
_seeds = [r for e, r in ATOMS_RAW if e == "Ta"] + [r for e, r in ATOMS_RAW if e == "Te"] + \
         [np.array([0.5, 0.5, 0.0108789126880766])*L, np.array([0, 0, 0.9891210873119234])*L]
SEEDS = np.array([fold(s) for s in _seeds])
CENTRES = np.array([SEEDS[i] + (lambda d: d - np.round(d/L)*L)(C[i] - _seeds[i]) for i in range(nw)])
OFFSET = np.linalg.norm(CENTRES - SEEDS, axis=1)

def nearest_sigma_x(x):
    xm = min((0.0, L[0]/2, -L[0]/2, L[0]), key=lambda q: abs(x - q)); return xm, abs(x - xm)

def bonds(cut=3.05):
    out = []; inside = lambda r: -L[0]/4 - 0.3 <= r[0] <= 0.75*L[0] + 0.3
    for i, (ei, ri) in enumerate(ATOMS):
        for j, (ej, rj) in enumerate(ATOMS):
            if j <= i: continue
            for nx in (-1, 0, 1):
                rr = rj + np.array([nx*L[0], 0, 0])
                if np.linalg.norm(rr - ri) < cut and {ei, ej} != {"Te"} and {ei, ej} != {"Pd"} and inside(ri) and inside(rr):
                    out.append((ri, rr))
    return out
