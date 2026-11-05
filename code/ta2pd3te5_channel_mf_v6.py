#!/usr/bin/env python3
"""
ta2pd3te5_channel_mf_v6.py  --  symmetry-resolved excitonic mean field for monolayer Ta2Pd3Te5
==============================================================================================
Consolidates and replaces v3, v4_SDWF_path_B and v5cco.

What it computes
  For a Wannier90 tight-binding Hamiltonian H0 of the Pmmn monolayer and a bond-local
  electron-hole interaction of strength V0 between conduction (Ta) and valence (Te/Pd)
  Wannier functions, it determines, per irreducible representation of the symmetry group
  that H0 actually possesses, (i) the critical coupling V0c = 1/lambda_max(chi) from the
  static q = 0 pair susceptibility and (ii) the self-consistent condensation energy
  F(Delta) - E_normal above V0c.

Corrections relative to v5cco.py / v4
  1. Operators are generated from the Pmmn operations, the Wannier centres (wout) and the
     orbital types; each one is TESTED on H(k) including nonsymmorphic phases. Only exact
     operations (residual < --symtol) enter the projectors. v5cco used hand tables that do
     not commute with H0 (0.24-0.53 eV); v4 used a diagonal C2z.
  2. Delta is bond-indexed, Delta_cv(k) = sum_R z_(c,v,R) exp(2 pi i k.R), for Ta-X bonds
     within --rcut. A cell-local (R = 0) Delta is not closed under glide/screw operations.
  3. Sign and reference: Delta = -V0 * P[g(Delta) - g0], where g are bond expectation values
     and g0 their normal-state values. v5cco used +V0*rho and no reference, so the
     normal-state hybridisation (and H0's symmetry-breaking content) acted as a source term.
  4. Energy: F = E_band(Delta) - lambda.g0 + |lambda|^2/(2 V0), evaluated at the final Delta.
  5. Convergence and seed dependence are reported; the linear-response V0c is the primary,
     mixing-independent result.

Inputs
  --hr    wannier90_hr.dat  (monolayer: all R_z = 0)
  --wout  wannier90.wout or an excerpt containing the lattice vectors and the
          "Final State" WF centres
  --orbitals  orbital type per WF (default: 4 dz2, 10 px, 2 dxz), axes x = chain (a), z = normal
Requires numpy. Example:
  python ta2pd3te5_channel_mf_v6.py --hr wannier90_hr.dat --wout wannier90.wout --V0 0.5 1 1.5 2
"""
import argparse, csv, re, time
import numpy as np

# ------------------------------------------------------------------ input
def load_hr(path):
    lines = open(path).readlines()
    nw, nr = int(lines[1]), int(lines[2]); ndl = -(-nr // 15)
    degen = np.array([float(x) for l in lines[3:3 + ndl] for x in l.split()][:nr])
    raw = np.loadtxt(lines[3 + ndl:])
    keys = [tuple(r) for r in raw[:, :3].astype(int)]; idx, Rl = {}, []
    for k in keys:
        if k not in idx: idx[k] = len(Rl); Rl.append(k)
    HR = np.zeros((len(Rl), nw, nw), complex)
    for k, row in zip(keys, raw):
        HR[idx[k], int(row[3]) - 1, int(row[4]) - 1] = row[5] + 1j * row[6]
    return np.array(Rl), HR, degen, nw

def load_wout(path, nw):
    txt = open(path).read()
    L = np.array([float(m) for m in re.findall(r"a_[123]\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", txt)[:3] for m in m]).reshape(3, 3)
    cen = re.findall(r"WF centre and spread\s+\d+\s+\(\s*([-\d.]+),\s*([-\d.]+),\s*([-\d.]+)\s*\)", txt)
    if len(cen) < nw: raise ValueError("wout: fewer WF centres than Wannier functions")
    return np.diag(L), np.array(cen[-nw:], float)   # orthorhombic cell assumed

# ------------------------------------------------------------------ symmetry
PARITY = {"dz2": lambda R: 1, "s": lambda R: 1, "px": lambda R: R[0], "py": lambda R: R[1],
          "pz": lambda R: R[2], "dxz": lambda R: R[0] * R[2], "dyz": lambda R: R[1] * R[2],
          "dxy": lambda R: R[0] * R[1], "dx2-y2": lambda R: 1}
OPNAMES = ["E", "C2x", "C2y", "C2z", "I", "sx", "sy", "sz"]
ROT = {"E": (1, 1, 1), "C2x": (1, -1, -1), "C2y": (-1, 1, -1), "C2z": (-1, -1, 1),
       "I": (-1, -1, -1), "sx": (-1, 1, 1), "sy": (1, -1, 1), "sz": (1, 1, -1)}
CHAR = {"Ag": [1, 1, 1, 1, 1, 1, 1, 1], "B1g": [1, -1, -1, 1, 1, -1, -1, 1],
        "B2g": [1, -1, 1, -1, 1, -1, 1, -1], "B3g": [1, 1, -1, -1, 1, 1, -1, -1],
        "Au": [1, 1, 1, 1, -1, -1, -1, -1], "B1u": [1, -1, -1, 1, -1, 1, 1, -1],
        "B2u": [1, -1, 1, -1, -1, 1, -1, 1], "B3u": [1, 1, -1, -1, -1, -1, 1, 1]}

def pmmn_ops(L, tglide):
    """Pmmn (origin choice 1): mirrors through x=0, y=0; I, C2x, C2y, sz carry translation tglide."""
    t = np.array(tglide) * L
    return {n: (np.array(ROT[n]), t if n in ("I", "C2x", "C2y", "sz") else np.zeros(3)) for n in OPNAMES}

def orbital_map(op, C, L, orbs):
    Rd, t = op; img = C * Rd + t; p, T, err = [], [], 0.0
    for i in range(len(C)):
        d = img[i] - C; n = np.round(d / L); dist = np.linalg.norm(d - n * L, axis=1)
        same = np.array([o == orbs[i] for o in orbs]); dist[~same] = np.inf
        j = int(np.argmin(dist)); p.append(j); T.append(n[j].astype(int)); err = max(err, dist[j])
    s = np.array([PARITY[orbs[i]](Rd) for i in range(len(C))], float)
    return p, np.array(T), s, err

def residual_matrix(op, omap, Rvec, HR, degen, k):
    """Delta_g(k) = U_g(k) H0(k) U_g(k)^dagger - H0(gk) for one reduced k (2D)."""
    Rd = op[0]; p, T, s, _ = omap; nw = len(p); gk = k * Rd[:2]
    Hk = np.einsum("r,rmn->mn", np.exp(2j * np.pi * Rvec[:, :2] @ k) / degen, HR)
    Hg = np.einsum("r,rmn->mn", np.exp(2j * np.pi * Rvec[:, :2] @ gk) / degen, HR)
    U = np.zeros((nw, nw), complex)
    for i in range(nw): U[p[i], i] = s[i] * np.exp(-2j * np.pi * gk @ T[i][:2])
    return U @ Hk @ U.conj().T - Hg

def residual_kpoints(mesh=(36, 8), npath=101, nrandom=6):
    """k points of the residual test: the uniform mesh, the X-Gamma-Y path and six random points."""
    ks = [np.array([i / mesh[0], j / mesh[1]]) for i in range(mesh[0]) for j in range(mesh[1])]
    ks += [np.array([0.5 * (1 - t), 0.0]) for t in np.linspace(0, 1, npath)]
    ks += [np.array([0.0, 0.5 * t]) for t in np.linspace(0, 1, npath)]
    ks += list(np.random.default_rng(0).random((nrandom, 2)))
    return ks

def op_residual(op, omap, Rvec, HR, degen, ks=None):
    """r_g = max over the test k points of max |Delta_g(k)| (eV)."""
    ks = residual_kpoints() if ks is None else ks
    return max(np.abs(residual_matrix(op, omap, Rvec, HR, degen, k)).max() for k in ks)

def op_residual_bound(op, omap, Rvec, HR, degen):
    """k-independent bound sum_R max|Delta_g(R)| (eV); Delta_g(R) are the exact lattice Fourier
    coefficients of the residual matrix, from an FFT on a grid larger than their range."""
    n = 2 * (int(np.abs(Rvec[:, :2]).max()) + 3) + 4
    D = np.array([[residual_matrix(op, omap, Rvec, HR, degen, np.array([i / n, j / n])) for j in range(n)]
                  for i in range(n)])
    DR = np.fft.fft2(D, axes=(0, 1)) / n**2
    return np.abs(DR).max(axis=(2, 3)).sum()

# ------------------------------------------------------------------ bonds
def build_bonds(C, L, cb, vb, rcut, nmax=3):
    bonds = []
    for c in cb:
        for v in vb:
            for rx in range(-nmax, nmax + 1):
                for ry in range(-nmax, nmax + 1):
                    if np.linalg.norm(C[v] + np.array([rx, ry, 0]) * L - C[c]) < rcut:
                        bonds.append((c, v, rx, ry))
    return bonds

def bond_rep(bonds, omap, Rd):
    """Signed permutation of the primitive operators (S_b, A_b) under one operation."""
    p, T, s, _ = omap; index = {b: i for i, b in enumerate(bonds)}; nb = len(bonds)
    D = np.zeros((2 * nb, 2 * nb))
    for i, (c, v, rx, ry) in enumerate(bonds):
        R2 = np.array([rx, ry]) * Rd[:2] + T[v][:2] - T[c][:2]
        j = index.get((p[c], p[v], int(R2[0]), int(R2[1])))
        if j is None: return None
        D[j, i] = D[nb + j, nb + i] = s[c] * s[v]
    return D

# ------------------------------------------------------------------ model
class Model:
    def __init__(s, a):
        s.Rvec, HR, degen, s.nw = load_hr(a.hr)
        if np.any(s.Rvec[:, 2] != 0): raise ValueError("R_z != 0 in hr.dat: not a monolayer model")
        s.L, s.C = load_wout(a.wout, s.nw)
        s.orbs = sum(([o] * int(n) for o, n in (x.split("*") for x in a.orbitals.split(","))), [])
        if len(s.orbs) != s.nw: raise ValueError("--orbitals does not match num_wann")
        s.nocc, s.cb = a.nocc, a.cb; s.vb = [i for i in range(s.nw) if i not in a.cb]
        nkx, nky = a.mesh
        kx, ky = np.meshgrid(np.arange(nkx) / nkx, np.arange(nky) / nky, indexing="ij")
        s.k = np.column_stack([kx.ravel(), ky.ravel()]); s.wk = 1 / len(s.k)
        s.H0 = np.einsum("kr,rmn->kmn", np.exp(2j * np.pi * s.k @ s.Rvec[:, :2].T) / degen, HR)
        s.HR, s.degen = HR, degen
        s.e0, s.U0 = np.linalg.eigh(s.H0); s.En = s.wk * s.e0[:, :s.nocc].sum()
        s.bonds = build_bonds(s.C, s.L, s.cb, s.vb, a.rcut)
        s.ph = np.exp(2j * np.pi * s.k @ np.array([[b[2], b[3]] for b in s.bonds]).T)  # (nk, nb)
        s.g0 = s.grad(np.zeros(2 * len(s.bonds)))[1]

    def delta(s, lam):
        nb = len(s.bonds); z = (lam[:nb] + 1j * lam[nb:]) / np.sqrt(2)
        D = np.zeros((len(s.k), s.nw, s.nw), complex)
        for i, (c, v, *_ ) in enumerate(s.bonds):
            D[:, c, v] += z[i] * s.ph[:, i]
        return D + D.conj().transpose(0, 2, 1)

    def grad(s, lam):
        e, V = np.linalg.eigh(s.H0 + s.delta(lam)); occ = V[:, :, :s.nocc]
        Eb = s.wk * e[:, :s.nocc].sum()
        rho_vc = np.array([np.einsum("ki,ki->k", occ[:, v], occ[:, c].conj()) for c, v, *_ in s.bonds]).T
        w = s.wk * np.sum(s.ph * rho_vc, axis=0)
        return Eb, np.sqrt(2) * np.concatenate([w.real, -w.imag])

    def F(s, lam, V0):
        Eb, g = s.grad(lam)
        return Eb - lam @ s.g0 + lam @ lam / (2 * V0), g

    def chi_prim(s):
        occ, emp = s.U0[:, :, :s.nocc], s.U0[:, :, s.nocc:]
        de = s.e0[:, s.nocc:, None] - s.e0[:, None, :s.nocc]
        M = []
        for i, (c, v, *_ ) in enumerate(s.bonds):
            z = s.ph[:, i][:, None, None]
            Mcv = emp[:, c, :, None].conj() * occ[:, v, None, :] * z + emp[:, v, :, None].conj() * occ[:, c, None, :] * z.conj()
            Mcv2 = emp[:, c, :, None].conj() * occ[:, v, None, :] * z - emp[:, v, :, None].conj() * occ[:, c, None, :] * z.conj()
            M.append(Mcv / np.sqrt(2)); M.append(1j * Mcv2 / np.sqrt(2))
        nb = len(s.bonds); order = [2 * i for i in range(nb)] + [2 * i + 1 for i in range(nb)]
        M = np.array(M)[order] / np.sqrt(de)[None]
        Mf = M.reshape(len(M), -1)
        return 2 * s.wk * np.real(Mf @ Mf.conj().T)

    def scf(s, Q, V0, lam, alpha, tol, maxit):
        for it in range(maxit):
            _, g = s.grad(lam)
            new = -V0 * Q @ (Q.T @ (g - s.g0))
            d = np.abs(new - lam).max(); lam = alpha * new + (1 - alpha) * lam
            if d < tol: return lam, True, it + 1
        return lam, False, maxit

# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hr", default="wannier90_hr.dat"); ap.add_argument("--wout", default="wannier90.wout")
    ap.add_argument("--orbitals", default="dz2*4,px*10,dxz*2")
    ap.add_argument("--cb", type=int, nargs="+", default=[0, 1, 2, 3])
    ap.add_argument("--nocc", type=int, default=12)
    ap.add_argument("--mesh", type=int, nargs=2, default=[36, 8])
    ap.add_argument("--tglide", type=float, nargs=3, default=[0.5, 0.5, 0.0])
    ap.add_argument("--symtol", type=float, default=0.02, help="eV; max residual for an operation to count as exact")
    ap.add_argument("--rcut", type=float, default=3.4, help="Angstrom; Ta-X bond cutoff between WF centres")
    ap.add_argument("--V0", type=float, nargs="+", default=[0.5, 1.0, 1.5, 2.0])
    ap.add_argument("--seeds", type=int, default=2); ap.add_argument("--amp", type=float, default=0.05)
    ap.add_argument("--alpha", type=float, default=0.3); ap.add_argument("--tol", type=float, default=1e-6)
    ap.add_argument("--maxit", type=int, default=500); ap.add_argument("--out", default="channel_mf_v6.csv")
    ap.add_argument("--channels", choices=["all", "odd"], default="all",
                    help="all: every channel the exact group separates; odd: only channels containing Au or B3u")
    a = ap.parse_args(); t0 = time.time()

    m = Model(a)
    gap = m.e0[:, a.nocc].min() - m.e0[:, a.nocc - 1].max()
    print(f"H0: nw={m.nw}, nR={len(m.Rvec)}, mesh={a.mesh}, nocc={a.nocc}, gap at fixed filling {gap*1e3:.1f} meV")
    print(f"Bonds: {len(m.bonds)} Ta-X bonds within {a.rcut} A  ({2*len(m.bonds)} real order-parameter components)")

    ops = pmmn_ops(m.L, a.tglide); exact = []
    print("\nSymmetry test of H0 (geometric operators from WF centres):")
    for n in OPNAMES:
        om = orbital_map(ops[n], m.C, m.L, m.orbs); r = op_residual(ops[n], om, m.Rvec, m.HR, m.degen)
        D = bond_rep(m.bonds, om, ops[n][0]) if om[3] < 0.05 else None
        ok = r < a.symtol and D is not None
        print(f"  {n:>3}: centre mismatch {om[3]:6.3f} A   residual {r*1e3:8.1f} meV   {'EXACT' if ok else 'broken'}")
        if ok: exact.append((n, D))
    names = [n for n, _ in exact]
    # irreps of the exact subgroup = classes of D2h irreps with equal characters on it
    classes = {}
    for ir, ch in CHAR.items():
        key = tuple(ch[OPNAMES.index(n)] for n in names); classes.setdefault(key, []).append(ir)
    print(f"\nExact group of H0: {{{', '.join(names)}}} (order {len(names)})")
    if len(names) < 8:
        print("  D2h irreps are NOT separable in this basis. Merged channels:")
    for key, irs in classes.items(): print(f"    {' + '.join(irs)}")

    chi = m.chi_prim(); rows = []; rng = np.random.default_rng(7)
    for key, irs in classes.items():
        label = "+".join(irs)
        if a.channels == "odd" and not any(x in irs for x in ("Au", "B3u")): continue
        P = sum(c * D for c, (_, D) in zip(key, exact)) / len(exact)
        U, sv, _ = np.linalg.svd(P); Q = U[:, sv > 0.5]
        lam = np.linalg.eigvalsh(Q.T @ chi @ Q).max(); V0c = 1 / lam
        print(f"\n=== channel {label}: dim {Q.shape[1]}, lambda_max(chi) = {lam:.4f} 1/eV  ->  V0c = {V0c:.3f} eV")
        for V0 in a.V0:
            best = None
            for sd in range(a.seeds):
                l0 = Q @ rng.normal(size=Q.shape[1]); l0 *= a.amp / np.linalg.norm(l0)
                lam_f, conv, nit = m.scf(Q, V0, l0, a.alpha, a.tol, a.maxit)
                F, _ = m.F(lam_f, V0)
                r = dict(channel=label, V0=V0, seed=sd, converged=conv, iters=nit,
                         Delta_meV=np.linalg.norm(lam_f) * 1e3, dF_meV=(F - m.En) * 1e3, V0c=V0c)
                rows.append(r)
                if conv and (best is None or r["dF_meV"] < best["dF_meV"]): best = r
            print(f"    V0={V0:5.2f}  " + (f"|Delta|={best['Delta_meV']:9.3f} meV   F-E_n={best['dF_meV']:+10.4f} meV"
                                         if best else "no seed converged"))
    if rows:
        with open(a.out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"\nWrote {a.out}   elapsed {time.time()-t0:.1f} s")

if __name__ == "__main__":
    main()
