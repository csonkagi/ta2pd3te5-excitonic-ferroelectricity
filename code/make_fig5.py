"""Figure 5: critical couplings of the merged channels, read from results/v0c_merged_channels.csv."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "code"))
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"Carlito","font.size":14,"mathtext.fontset":"dejavusans","axes.linewidth":1.0})
import csv
rows=list(csv.DictReader(open(ROOT/"results/v0c_merged_channels.csv")))
rc=sorted({float(r["r_c_A"]) for r in rows}); nb=[int(next(r["n_bonds"] for r in rows if float(r["r_c_A"])==x)) for x in rc]
Au=[float(next(r["V0c_eV"] for r in rows if float(r["r_c_A"])==x and "Au" in r["channel"])) for x in rc]
B3=[float(next(r["V0c_eV"] for r in rows if float(r["r_c_A"])==x and "B3u" in r["channel"])) for x in rc]
fig,ax=plt.subplots(figsize=(7.2,4.8))
ax.axhline(1.1,ls="--",c="0.4",lw=1.2,zorder=0)
ax.text(3.27,1.088,r"$V=1.1$ eV: bulk Hartree–Fock interaction (Huang et al.)",fontsize=12.5,color="0.3",va="top")
ax.plot(rc,B3,"s-",c="#c0392b",lw=2.2,ms=9,zorder=3); ax.plot(rc,Au,"o-",c="#1f5aa6",lw=2.2,ms=9,zorder=3)
for x,y in zip(rc,B3): ax.text(x,y+0.025,f"{y:.2f}",ha="center",va="bottom",fontsize=13,color="#c0392b")
for x,y in zip(rc[:2],Au[:2]): ax.text(x,y-0.03,f"{y:.2f}",ha="center",va="top",fontsize=13,color="#1f5aa6")
ax.text(rc[2]+0.04,Au[2]+0.012,f"{Au[2]:.2f}",ha="left",va="bottom",fontsize=13,color="#1f5aa6")
ax.text(4.42,1.47,"$A_g\\oplus B_{3u}$\n(contains polar $B_{3u}$)",color="#c0392b",fontsize=14,ha="center",va="center",linespacing=1.3)
ax.text(3.58,1.20,"$A_u\\oplus B_{3g}$\n(contains band-edge $A_u$)",color="#1f5aa6",fontsize=14,ha="center",va="center",linespacing=1.3)
ax.set_xticks(rc); ax.set_xticklabels([f"{r:.1f} Å\n{n} bonds" for r,n in zip(rc,nb)],fontsize=13)
ax.set_xlabel(r"Ta–X bond cutoff $r_c$ and number of bond pairs",fontsize=14)
ax.set_ylabel(r"critical coupling $V_0^{c}$ (eV)",fontsize=15); ax.tick_params(axis="y",labelsize=13)
ax.set_xlim(3.22,4.78); ax.set_ylim(1.0,1.66)
fig.tight_layout(); fig.savefig(ROOT/"figures/fig5_v0c.png",dpi=300); plt.close(fig)
