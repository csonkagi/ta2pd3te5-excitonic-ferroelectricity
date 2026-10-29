"""Static Fig. 4: 3D perspective of the Wannier centres, matching the default view of figures/fig4_interactive.html."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from fig4_geometry import ROOT, L, nw, ATOMS, SEEDS, CENTRES, COL, ELW, nearest_sigma_x
plt.rcParams.update({"font.family": "Carlito", "font.size": 11, "mathtext.fontset": "dejavusans"})

# camera of the interactive default view: Plotly eye (1.4, -1.6, 0.9)
eye = np.array([1.4, -1.6, 0.9]); ELEV = np.degrees(np.arctan2(eye[2], np.hypot(eye[0], eye[1]))); AZIM = np.degrees(np.arctan2(eye[1], eye[0]))
fig = plt.figure(figsize=(7.4, 5.6)); ax = fig.add_axes([0.05, 0.17, 0.84, 0.86], projection="3d")
ex, ey, ez = L[0] + 1.0, L[1], 6.0
ax.set_box_aspect((ex, ey, ez), zoom=1.2); ax.view_init(elev=ELEV, azim=AZIM)
ax.set_xlim(-L[0]/4 - 0.5, 0.75*L[0] + 0.5); ax.set_ylim(-L[1]/2, L[1]/2); ax.set_zlim(-3, 3)
ax.set_xticks([0, L[0]/2]); ax.set_xticklabels(["0", "a/2"]); ax.set_yticks([-8, -4, 0, 4, 8]); ax.set_zticks([-2, 0, 2])
ax.set_xlabel("x (Å)", labelpad=0); ax.set_ylabel("y (Å)", labelpad=8); ax.set_zlabel("z (Å)", labelpad=4)
for a in (ax.xaxis, ax.yaxis, ax.zaxis): a.pane.fill = False
ax.grid(False)
yy, zz = np.meshgrid([-L[1]/2, L[1]/2], [-2.8, 2.8])
for xm in (0, L[0]/2): ax.plot_surface(np.full_like(yy, xm), yy, zz, color="#888888", alpha=0.12, linewidth=0, shade=False)
for e, r in ATOMS: ax.scatter(*r, s=95, c=COL[e], alpha=0.45, edgecolors=COL[e], linewidths=0.7, depthshade=False)
for i in range(nw):
    if ELW[i] == "Ta": continue
    d = CENTRES[i] - SEEDS[i]
    ax.quiver(*SEEDS[i], *d, color=COL["Te"] if ELW[i] == "Te" else "#555555", arrow_length_ratio=0.12, linewidth=1.7)
for i in range(nw):
    xm, dd = nearest_sigma_x(CENTRES[i, 0])
    if dd > 0.1: ax.plot([CENTRES[i, 0], xm], [CENTRES[i, 1]]*2, [CENTRES[i, 2]]*2, c="#d62728", lw=3.0)
MK = {"Ta": "o", "Te": "s", "Pd": "D"}
for i in range(nw):
    ax.scatter(*CENTRES[i], s=30, c=COL[ELW[i]], marker=MK[ELW[i]], edgecolors="k", linewidths=0.6, depthshade=False)
h = [Line2D([], [], marker="o", ls="", mfc=COL[e], mec=COL[e], alpha=0.5, ms=10, label=f"{e} atom") for e in ("Ta", "Te", "Pd")] + \
    [Line2D([], [], marker=MK["Ta"], ls="", mfc=COL["Ta"], mec="k", ms=6, label=r"Ta $d_{z^2}$ centre"),
     Line2D([], [], marker=MK["Te"], ls="", mfc=COL["Te"], mec="k", ms=6, label=r"Te $p_x$-seeded centre"),
     Line2D([], [], marker=MK["Pd"], ls="", mfc=COL["Pd"], mec="k", ms=6, label=r"Pd $d_{xz}$ centre"),
     Line2D([], [], color=COL["Te"], lw=1.7, label="seed atom → centre (true length)"),
     Patch(facecolor="#888888", alpha=0.3, label=r"$\sigma_x$ mirror planes"),
     Line2D([], [], color="#d62728", lw=3, label=r"offset from nearest $\sigma_x$ plane")]
fig.legend(handles=h, loc="upper center", ncol=3, fontsize=9.5, frameon=False, bbox_to_anchor=(0.5, 0.15))
out = ROOT/"figures/fig4_symmetry.png"; fig.savefig(out, dpi=300, bbox_inches="tight", pad_inches=0.2); print("wrote", out)
