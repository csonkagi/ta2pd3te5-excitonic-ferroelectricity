"""Figure 5: critical couplings V0c of the four merged channels of the maximally localized model.

Reads results/v0c_merged_channels.csv (columns r_c_A, n_bonds, channel, lambda_max_per_eV, V0c_eV),
written by code/compute_results.py from data/wannier/wannier90_hr.dat for all four merged channels.
Markers are the computed values. The dotted curves are monotone piecewise-cubic Hermite interpolants
(PCHIP) through those values: smooth, exact at the computed points, and free of the spurious maxima
or minima an ordinary spline can introduce between three points. They guide the eye; they are not data.
"""
from pathlib import Path
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
plt.rcParams.update({"font.family": "Carlito", "font.size": 14, "mathtext.fontset": "dejavusans",
                     "axes.linewidth": 1.0})

# channel key in the CSV -> (label, colour, marker, where the value labels go)
CHANNELS = {
    "B3g+Au":  (r"$A_u\oplus B_{3g}$", "#1f5aa6", "o", "below"),
    "Ag+B3u":  (r"$A_g\oplus B_{3u}$", "#c0392b", "s", "above"),
    "B1g+B2u": (r"$B_{1g}\oplus B_{2u}$", "#6f7f5a", "D", "below"),
    "B2g+B1u": (r"$B_{2g}\oplus B_{1u}$", "#7d6b91", "^", "above"),
}

rows = list(csv.DictReader(open(ROOT / "results/v0c_merged_channels.csv")))
rc = sorted({float(r["r_c_A"]) for r in rows})
nbonds = [int(next(r["n_bonds"] for r in rows if float(r["r_c_A"]) == x)) for x in rc]
v0c = {ch: np.array([float(next(r["V0c_eV"] for r in rows if float(r["r_c_A"]) == x and r["channel"] == ch))
                     for x in rc]) for ch in CHANNELS}

fig, ax = plt.subplots(figsize=(7.6, 5.4))
REF = 1.1                       # eV, bulk Hartree-Fock interaction of Huang et al.
MIN_LABEL_GAP = 0.055           # eV, minimum vertical spacing of the channel names
ax.axhline(REF, ls="--", c="0.45", lw=1.1, zorder=0)
ax.text(3.27, 1.092, r"$V=1.1$ eV: bulk Hartree–Fock interaction (Huang et al.)",
        fontsize=12, color="0.35", va="top")

xs = np.linspace(rc[0], rc[-1], 200)
for ch, (label, col, mk, where) in CHANNELS.items():
    y = v0c[ch]
    ax.plot(xs, PchipInterpolator(rc, y)(xs), ls=(0, (1.2, 2.2)), lw=2.4, c=col, zorder=2,
            solid_capstyle="round", dash_capstyle="round")
    ax.plot(rc, y, ls="none", marker=mk, ms=8.5, mfc=col, mec="white", mew=1.0, c=col, zorder=3)
    for x, v in zip(rc, y):
        below = where == "below" and v - 0.06 > REF   # a label below the point must not reach the reference line
        dy, va = (-0.022, "top") if below else (0.022, "bottom")
        ax.text(x, v + dy, f"{v:.2f}", ha="center", va=va, fontsize=11.5, color=col, zorder=4)

# channel names at the right end, spread so that close curves keep readable labels
ends = sorted(((v0c[ch][-1], ch) for ch in CHANNELS))
ypos = []
for v, ch in ends:
    ypos.append(max(v, ypos[-1] + MIN_LABEL_GAP) if ypos else v)
shift = max(0.0, (ypos[-1] - ends[-1][0]) / 2)      # centre the spread group on the curves
for (v, ch), yl in zip(ends, ypos):
    label, col, _, _ = CHANNELS[ch]
    ax.text(rc[-1] + 0.07, yl - shift if len(ends) > 1 and yl - shift >= v - MIN_LABEL_GAP else yl,
            label, color=col, fontsize=13.5, ha="left", va="center")

ax.set_xticks(rc)
ax.set_xticklabels([f"{r:.1f} Å\n{n} bonds" for r, n in zip(rc, nbonds)], fontsize=13)
ax.set_xlabel(r"Ta–X bond cutoff $r_c$ and number of bond pairs", fontsize=14)
ax.set_ylabel(r"critical coupling $V_0^{c}$ (eV)", fontsize=15)
ax.tick_params(axis="y", labelsize=13)
ax.set_xlim(3.25, 5.15)
ax.set_ylim(1.0, 1.95)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
fig.tight_layout()
fig.savefig(ROOT / "figures/fig5_v0c.png", dpi=300)
plt.close(fig)
