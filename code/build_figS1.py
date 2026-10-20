"""Figure S1: ARPES (Zhang et al., PRX 14, 011047, CC BY 4.0) with r2SCAN+rVV10 bands, panels rebuilt on calibrated axes."""
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image
from scipy import ndimage
plt.rcParams.update({"font.family":"Carlito","font.size":9,"mathtext.fontset":"dejavusans"})
src=np.array(Image.open(ROOT/'data/external/zhang2024_prx_fig2_with_r2scan_overlay.png').convert('RGB')).astype(float)
colx=[60,377,693]; rowy=[22,382]; W,H=251,287
kx=lambda x:(x-126.0)/417.5; ey=lambda y:(95.5-y)/637.5
L,R,T,B=17,3,3,13
panels=[]; boxes=[]
for r,y0 in enumerate(rowy):
    for c,x0 in enumerate(colx):
        p=src[y0+T:y0+H-B, x0+L:x0+W-R]
        if r==0:   # bounding box of the original temperature label glyphs (and halo)
            reg=p[:80,:]; x_off=int(0.45*p.shape[1])
            sub=reg[:,x_off:]
            glyph=(sub.mean(2)<130)&((sub.max(2)-sub.min(2))<30)
            lab,n=ndimage.label(glyph); sz=ndimage.sum(glyph,lab,range(1,n+1))
            glyph=np.isin(lab,1+np.where(sz>15)[0])
            ys,xs=np.where(glyph)
            boxes.append((x_off+xs.min()-6, ys.min()-6, x_off+xs.max()+6, ys.max()+6))
        panels.append(p.astype(np.uint8))
bx0=135-4-4; by0=0; bx1=p.shape[1]; by1=42
ext=[kx(L),kx(W-R),ey(H-B),ey(T)]
k0,k1=kx(L+bx0),ext[1]; e1,e0=ext[3],ey(T+by1)
fig,axs=plt.subplots(2,3,figsize=(7.0,4.9),sharex=True,sharey=True)
Tlab=["350 K","300 K","200 K"]; Elab=["30 eV","67 eV"]
for i,ax in enumerate(axs.flat):
    r,c=divmod(i,3)
    ax.imshow(panels[i],extent=ext,aspect="auto",interpolation="lanczos")
    ax.add_patch(Rectangle((k0,e0),k1-k0,e1-e0,facecolor="white",edgecolor="0.6",lw=0.5,zorder=3))
    ax.text((k0+k1)/2,(e0+e1)/2,f"{Tlab[c]}\n{Elab[r]}",ha="center",va="center",fontsize=8,zorder=4,linespacing=1.1)
    ax.text(0.03,0.97,f"({'abcdef'[i]})",transform=ax.transAxes,va="top",ha="left",fontsize=9,zorder=4,
            bbox=dict(facecolor="white",edgecolor="none",pad=1.2,alpha=0.85))
    ax.set_xlim(ext[0],ext[1]); ax.set_ylim(ext[2],ext[3])
    ax.set_xticks([-0.2,0,0.2]); ax.set_yticks([-0.2,-0.1,0,0.1])
    if c==0: ax.set_ylabel(r"$E-E_\mathrm{F}$ (eV)")
    if r==1: ax.set_xlabel(r"$k_b$ (Å$^{-1}$)")
fig.tight_layout(pad=0.4,w_pad=0.6,h_pad=0.6)
fig.savefig(ROOT/"figures/figS1_arpes_r2scan.png",dpi=400)
print("label box k",round(k0,3),round(k1,3),"E",round(e0,3),round(e1,3))
