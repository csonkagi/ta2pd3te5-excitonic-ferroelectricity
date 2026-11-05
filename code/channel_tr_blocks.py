#!/usr/bin/env python3
"""
channel_tr_blocks.py -- time-reversal-resolved critical couplings for every excitonic channel.

Uses the model of ta2pd3te5_channel_mf_v6.py (imported, not modified). Each Ta-X bond carries two
real order-parameter components, S_b = (c^+ v + h.c.)/sqrt2 (time-reversal even) and
A_b = i(c^+ v - h.c.)/sqrt2 (time-reversal odd); v6 orders them [S_1..S_nb, A_1..A_nb]. For a
time-reversal-symmetric H0 the static susceptibility does not couple the two sets, so every channel
splits into a T-even and a T-odd block. This script reports V0c = 1/lambda_max for both blocks.

Spin interpretation (no spin-orbit coupling, interaction acting equally in all four spin channels):
  T-even block of channel G -> spin-singlet order G (ordinary excitonic order) and, degenerate with
                               it, the magnetic (T-odd) spin triplets G x S
  T-odd  block of channel G -> spin-singlet bond-current order G and, degenerate with it, the
                               T-even spin-current triplets G x S
For G = Au the T-odd block therefore gives the threshold of the polar triplet Au x S_x = B3u.

Usage (same model options as v6; several cutoffs in one run):
  python channel_tr_blocks.py --v6 ta2pd3te5_channel_mf_v6.py --hr wannier90_hr.dat \
         --wout wannier90.wout --rcut 3.4 4.0 4.6
Writes channel_tr_blocks.csv.
"""
import argparse
import csv
import importlib.util
import types

import numpy as np


def load_v6(path):
    spec = importlib.util.spec_from_file_location("v6", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    for name in ("Model", "pmmn_ops", "orbital_map", "op_residual", "bond_rep", "CHAR", "OPNAMES"):
        if not hasattr(mod, name):
            raise SystemExit(f"STOP: {path} has no '{name}'; this script needs the v6 model API")
    return mod


def residual(v6, op, om, m):
    try:
        return v6.op_residual(op, om, m.Rvec, m.HR, m.degen)
    except TypeError:   # newer v6 versions may take extra arguments
        return v6.op_residual(op, om, m)


def block_max(P, chi, sel):
    """lambda_max of chi on the image of P restricted to the components in sel."""
    Pr = P[:, sel]
    U, sv, _ = np.linalg.svd(Pr, full_matrices=False)
    Q = U[:, sv > 0.5]
    if Q.shape[1] == 0:
        return np.nan, 0
    return np.linalg.eigvalsh(Q.T @ chi @ Q).max(), Q.shape[1]


def run(v6, a, rcut):
    args = types.SimpleNamespace(**{**vars(a), "rcut": rcut})
    m = v6.Model(args)
    nb = len(m.bonds)
    ops = v6.pmmn_ops(m.L, a.tglide)
    exact = []
    for n in v6.OPNAMES:
        om = v6.orbital_map(ops[n], m.C, m.L, m.orbs)
        r = residual(v6, ops[n], om, m)
        D = v6.bond_rep(m.bonds, om, ops[n][0]) if om[3] < 0.05 else None
        if r < a.symtol and D is not None:
            exact.append((n, D))
    names = [n for n, _ in exact]
    classes = {}
    for ir, ch in v6.CHAR.items():
        key = tuple(ch[v6.OPNAMES.index(n)] for n in names)
        classes.setdefault(key, []).append(ir)
    chi = m.chi_prim()
    coupling = np.abs(chi[:nb, nb:]).max() / np.abs(chi).max()
    even, odd = np.arange(nb), np.arange(nb, 2 * nb)
    rows = []
    for key, irs in classes.items():
        P = sum(c * D for c, (_, D) in zip(key, exact)) / len(exact)
        le, de = block_max(P, chi, even)
        lo, do = block_max(P, chi, odd)
        lall = max(le, lo)
        rows.append(dict(rcut=rcut, bonds=nb, exact_group="{" + ",".join(names) + "}", channel="+".join(irs),
                         V0c_all=1 / lall, V0c_Teven=1 / le, V0c_Todd=1 / lo, dim_even=de, dim_odd=do,
                         leading="T-even" if le >= lo else "T-odd", chi_even_odd_coupling=coupling))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--v6", default="ta2pd3te5_channel_mf_v6.py")
    ap.add_argument("--hr", default="wannier90_hr.dat")
    ap.add_argument("--wout", default="wannier90.wout")
    ap.add_argument("--orbitals", default="dz2*4,px*10,dxz*2")
    ap.add_argument("--cb", type=int, nargs="+", default=[0, 1, 2, 3])
    ap.add_argument("--nocc", type=int, default=12)
    ap.add_argument("--mesh", type=int, nargs=2, default=[36, 8])
    ap.add_argument("--tglide", type=float, nargs=3, default=[0.5, 0.5, 0.0])
    ap.add_argument("--symtol", type=float, default=0.02)
    ap.add_argument("--rcut", type=float, nargs="+", default=[3.4, 4.0, 4.6])
    ap.add_argument("--out", default="channel_tr_blocks.csv")
    a = ap.parse_args()
    v6 = load_v6(a.v6)
    rows = []
    for rc in a.rcut:
        r = run(v6, a, rc)
        rows += r
        print(f"\nr_c = {rc} A: {r[0]['bonds']} bonds, exact group {r[0]['exact_group']}, "
              f"T-even/T-odd coupling in chi {r[0]['chi_even_odd_coupling']:.1e} (must be ~0)")
        print(f"  {'channel':<22}{'V0c all':>9}{'T-even':>9}{'T-odd':>9}   leading")
        for x in sorted(r, key=lambda x: x["V0c_all"]):
            print(f"  {x['channel']:<22}{x['V0c_all']:9.4f}{x['V0c_Teven']:9.4f}{x['V0c_Todd']:9.4f}   {x['leading']}")
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwritten {a.out}. For channel Au (or the merged channel containing Au): T-even = singlet order and "
          "magnetic triplets; T-odd = singlet bond currents and the polar (T-even) triplet Au x S_x = B3u.")


if __name__ == "__main__":
    main()
