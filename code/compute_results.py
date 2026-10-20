#!/usr/bin/env python3
"""Reproduce every number quoted in the paper from the Wannier90 files in data/wannier.

Writes to results/:
  symmetry_residuals.csv      Table III  (r_g for the seven nontrivial Pmmn operations)
  gamma_eigenvalues.csv       Table I, right-hand columns (v1, c1, c2 under the exact operations)
  wannier_centres.csv         Fig. 4 and Sec. III.E (centre, seed, offset, spread, distance to sigma_x planes)
  v0c_merged_channels.csv     Fig. 5  (critical couplings for r_c = 3.4, 4.0, 4.6 A)
Usage: python code/compute_results.py   (from the repository root; about 5 min on a laptop)
"""
import csv, re, sys, argparse
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "code"))
import ta2pd3te5_channel_mf_v6 as v6

HR, WOUT, WIN = ROOT/"data/wannier/wannier90_hr.dat", ROOT/"data/wannier/wannier90_wout_final_state.txt", ROOT/"data/wannier/wannier90.win"
OUT = ROOT/"results"; OUT.mkdir(exist_ok=True)
ORBS = ["dz2"]*4 + ["px"]*10 + ["dxz"]*2
NAMES = {"C2x":"2_1x screw","C2y":"2_1y screw","C2z":"C2z rotation","I":"inversion","sx":"sigma_x mirror","sy":"sigma_y mirror","sz":"n_z glide"}

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--skip-v0c", action="store_true"); a = ap.parse_args()
    Rvec, HRm, degen, nw = v6.load_hr(str(HR)); L, C = v6.load_wout(str(WOUT), nw)
    ops = v6.pmmn_ops(L, [0.5, 0.5, 0.0])
    # --- Table III
    rows, maps = [], {}
    for n in v6.OPNAMES[1:]:
        om = v6.orbital_map(ops[n], C, L, ORBS); maps[n] = om
        r = v6.op_residual(ops[n], om, Rvec, HRm, degen) * 1e3
        rows.append(dict(operation=n, type=NAMES[n], residual_meV=round(r, 1), exact=r < 20))
    with open(OUT/"symmetry_residuals.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    # --- Table I, right-hand columns
    HG = np.einsum("r,rmn->mn", 1/degen, HRm); e, V = np.linalg.eigh(HG)
    def U(n):
        p, T, s, _ = maps[n]; M = np.zeros((nw, nw))
        for i in range(nw): M[p[i], i] = s[i]
        return M
    with open(OUT/"gamma_eigenvalues.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["state", "E_minus_EVBM_meV", "2_1x", "sigma_y", "n_z"])
        for k, lab in [(11, "v1"), (12, "c1"), (13, "c2")]:
            w.writerow([lab, round((e[k]-e[11])*1e3, 1)] + [round(float(np.vdot(V[:, k], U(n) @ V[:, k]).real), 6) for n in ("C2x", "sy", "sz")])
    # --- Wannier centres vs projection seeds
    txt = WIN.read_text(); blk = re.search(r"begin atoms_cart(.*?)end atoms_cart", txt, re.S | re.I).group(1)
    atoms = [(l.split()[0], np.array(l.split()[1:4], float)) for l in blk.strip().splitlines()]
    spreads = [float(x) for x in re.findall(r"WF centre and spread\s+\d+\s+\([^)]*\)\s+([\d.]+)", WOUT.read_text())][-nw:]
    seeds = [r for e_, r in atoms if e_ == "Ta"] + [r for e_, r in atoms if e_ == "Te"] + \
            [np.array([0.5, 0.5, 0.0108789126880766])*L, np.array([0, 0, 0.9891210873119234])*L]
    with open(OUT/"wannier_centres.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["wf", "seed_orbital", "centre_x", "centre_y", "centre_z", "spread_A2", "rms_radius_A", "offset_from_seed_A", "distance_to_sigma_x_plane_A"])
        for i in range(nw):
            d = C[i]-seeds[i]; d -= np.round(d/L)*L
            xm = np.mod(C[i][0], L[0]); dpl = min(abs(xm), abs(xm-L[0]/2), abs(xm-L[0]))
            w.writerow([i+1, ORBS[i]] + list(np.round(C[i], 4)) + [spreads[i], round(spreads[i]**0.5, 2), round(float(np.linalg.norm(d)), 3), round(dpl, 3)])
    print("wrote symmetry_residuals.csv, gamma_eigenvalues.csv, wannier_centres.csv")
    if a.skip_v0c: return
    # --- Fig. 5: critical couplings of the merged channels
    out = []
    for rc in (3.4, 4.0, 4.6):
        args = argparse.Namespace(hr=str(HR), wout=str(WOUT), orbitals="dz2*4,px*10,dxz*2", cb=[0, 1, 2, 3], nocc=12,
                                  mesh=[36, 8], tglide=[0.5, 0.5, 0.0], symtol=0.02, rcut=rc)
        m = v6.Model(args)
        exact = []
        for n in v6.OPNAMES:
            om = v6.orbital_map(ops[n], m.C, m.L, m.orbs)
            if om[3] < 0.05 and v6.op_residual(ops[n], om, m.Rvec, m.HR, m.degen) < args.symtol:
                D = v6.bond_rep(m.bonds, om, ops[n][0])
                if D is not None: exact.append((n, D))
        names = [n for n, _ in exact]; chi = m.chi_prim(); classes = {}
        for ir, ch in v6.CHAR.items():
            classes.setdefault(tuple(ch[v6.OPNAMES.index(n)] for n in names), []).append(ir)
        for key, irs in classes.items():
            if not any(x in irs for x in ("Au", "B3u")): continue
            P = sum(c*D for c, (_, D) in zip(key, exact))/len(exact); Uu, sv, _ = np.linalg.svd(P); Q = Uu[:, sv > 0.5]
            lam = np.linalg.eigvalsh(Q.T @ chi @ Q).max()
            out.append(dict(r_c_A=rc, n_bonds=len(m.bonds), channel="+".join(irs), lambda_max_per_eV=round(lam, 4), V0c_eV=round(1/lam, 3)))
            print(out[-1])
    with open(OUT/"v0c_merged_channels.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
    print("wrote v0c_merged_channels.csv")

if __name__ == "__main__":
    main()
