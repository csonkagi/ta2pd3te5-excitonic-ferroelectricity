"""Figure 3 (bands, pairings at Gamma, interband velocity matrix elements); also writes the earlier two-panel 2D projections of Fig. 4 (figures/fig4_2d_projections.png). The final Fig. 4 is made by make_fig4.py."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "code"))
import numpy as np, re, importlib.util, matplotlib, os
FIG3_LABELPOS=tuple(float(v) for v in os.environ.get('FIG3_LABELPOS','0.655,0.035,0.17').split(','))
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import ta2pd3te5_channel_mf_v6 as v6
plt.rcParams.update({"font.family":"Carlito","font.size":10,"mathtext.fontset":"dejavusans","axes.linewidth":0.8})
Rvec,HR,degen,nw=v6.load_hr(str(ROOT/"data/wannier/wannier90_hr.dat")); L,C=v6.load_wout(str(ROOT/"data/wannier/wannier90_wout_final_state.txt"),nw)
orbs=["dz2"]*4+["px"]*10+["dxz"]*2
txt=open(ROOT/"data/wannier/wannier90.win").read()
blk=re.search(r"begin atoms_cart(.*?)end atoms_cart",txt,re.S).group(1)
atoms=[(l.split()[0],np.array(l.split()[1:4],float)) for l in blk.strip().splitlines()]
COL={"Ta":"#2f6db5","Pd":"#8c8c8c","Te":"#1b9e5a"}
def Hk(k): return np.einsum("r,rmn->mn",np.exp(2j*np.pi*Rvec[:,:2]@k)/degen,HR)
# ---------------- Fig. 3: bands, pairings at Gamma, interband dipole
Rc=Rvec[:,:2]*L[:2]
def Hv(k):
    ph=np.exp(2j*np.pi*Rvec[:,:2]@k)/degen
    return np.einsum("r,rmn->mn",ph,HR), np.einsum("r,rmn->mn",1j*Rc[:,0]*ph,HR)
X=np.array([0.5,0]); G=np.array([0,0]); Y=np.array([0,0.5]); b=np.array([2*np.pi/L[0],2*np.pi/L[1]])
ks=[];xs=[];x0=0
for A,B in [(X,G),(G,Y)]:
    for t in np.linspace(0,1,400): ks.append(A+(B-A)*t); xs.append(x0+np.linalg.norm((B-A)*b)*t)
    x0+=np.linalg.norm((B-A)*b)
xs=np.array(xs); ks=np.array(ks)
E=np.array([np.linalg.eigvalsh(Hv(k)[0]) for k in ks])
# follow v1, c1, c2 from Gamma outward by eigenvector overlap; near-degenerate pairs are first split into
# symmetry-pure states with the operations that are exact at that k (n_z everywhere; sigma_y, 2_1x on Gamma-X),
# so that true crossings between different symmetry sectors are not mistaken for mixing
OPS=v6.pmmn_ops(L,[0.5,0.5,0]); MAPS={n:v6.orbital_map(OPS[n],C,L,orbs) for n in ("sz","sy","C2x")}
def Uk(n,k):
    p,T,sg,_=MAPS[n]; Rd=OPS[n][0]; gk=k*Rd[:2]; M=np.zeros((nw,nw),complex)
    for i in range(nw): M[p[i],i]=sg[i]*np.exp(-2j*np.pi*gk@T[i][:2])
    return M
def pure_eigs(k):
    H,vx=Hv(k); e,U=np.linalg.eigh(H)
    keep=["sz"]+(["sy","C2x"] if abs(k[1])<1e-12 else [])
    coef={"sz":1.0,"sy":0.37,"C2x":0.19}
    S=sum(coef[n]*Uk(n,k) for n in keep)
    i=0
    while i<nw:
        j=i+1
        while j<nw and e[j]-e[j-1]<5e-3: j+=1
        if j-i>1:
            P=U[:,i:j]; w,Vs=np.linalg.eig(P.conj().T@S@P)
            Q,_=np.linalg.qr(P@Vs); U[:,i:j]=Q
        i=j
    return e,U,vx
iG=int(np.argmin(np.linalg.norm(ks,axis=1))); M1=np.zeros(len(ks)); M2=np.zeros(len(ks))
for direction in (range(iG,-1,-1), range(iG,len(ks))):
    prev=None
    for i in direction:
        e,U,vx=pure_eigs(ks[i])
        idx=[11,12,13] if prev is None else [int(np.argmax(abs(U.conj().T@q))) for q in prev]
        vec=[U[:,j] for j in idx]; prev=vec
        M1[i]=abs(vec[1].conj()@vx@vec[0])**2; M2[i]=abs(vec[2].conj()@vx@vec[0])**2
eG=np.linalg.eigvalsh(Hv(G)[0]); e0=eG[11]; xG=np.linalg.norm((G-X)*b)
cA,cB="#1f5aa6","#c0392b"
LX,LY1,LY2=FIG3_LABELPOS
fig,(ax,axd)=plt.subplots(2,1,figsize=(6.6,6.4),sharex=True,gridspec_kw=dict(height_ratios=[1,0.42],hspace=0.07))
for n in range(nw): ax.plot(xs,E[:,n]-e0,c="#3a6fb5" if n>=12 else "#555555",lw=1.2)
for x in [0,xG,x0]:
    ax.axvline(x,c="0.6",lw=0.6); axd.axvline(x,c="0.6",lw=0.6)
ax.axhline(0,c="0.6",lw=0.6,ls=":")
lab=[(11,r"$v_1\;(\Gamma_4^+)$"),(12,r"$c_1\;(\Gamma_4^-)$"),(13,r"$c_2\;(\Gamma_2^-)$")]
for idx,t in lab:
    ey=eG[idx]-e0
    ax.plot([xG],[ey],"o",ms=5,mfc="white",mec="k",mew=1.0,zorder=6)
    ax.plot([xG,x0+0.04],[ey,ey],ls=":",lw=0.9,color="0.3",clip_on=False,zorder=4)
    ax.text(x0+0.05,ey,t,fontsize=12,ha="left",va="center",clip_on=False)
g1=eG[12]-e0; g2=eG[13]-e0
for xa,top,col in [(xG-0.035,g1,cA),(xG-0.085,g2,cB)]:
    ax.annotate("",xy=(xa,top),xytext=(xa,0),arrowprops=dict(arrowstyle="-|>",lw=1.6,color=col,shrinkA=0,shrinkB=0,mutation_scale=11),zorder=5)
ax.text(LX,LY1,r"$c_1\!\leftarrow\!v_1$: $A_u$, 108 meV, nonpolar",color=cA,fontsize=10.5,ha="right",va="center")
ax.text(LX,LY2,r"$c_2\!\leftarrow\!v_1$: $B_{3u}$, 237 meV, polar",color=cB,fontsize=10.5,ha="right",va="center")
ax.set_xlim(0,x0); ax.set_ylim(-1.2,1.2)
ax.set_ylabel(r"$E-E_{\mathrm{VBM}}$ (eV)",fontsize=13); ax.tick_params(axis="y",labelsize=12)
ax.text(0.012,1.12,"(a)",fontsize=12,va="top")
keep=(xs>0.12)&(xs<x0-0.03)          # omit the ends: near X (three near-degenerate valence bands) and at Y (bands meet)
axd.plot(xs[keep],M2[keep],c=cB,lw=1.8); axd.plot(xs[keep],M1[keep],c=cA,lw=1.8)
axd.plot([xG],[0],"o",ms=5,mfc="white",mec=cA,mew=1.2,zorder=5)
axd.text(0.47,1.62,r"$|\langle c_2|v_x|v_1\rangle|^2$: $B_{3u}$ pair, allowed",color=cB,fontsize=11,ha="left",va="center")
axd.text(0.03,-0.36,r"$|\langle c_1|v_x|v_1\rangle|^2$: $A_u$ pair, zero along X$-\Gamma$ and at $\Gamma$",color=cA,fontsize=11,ha="left",va="center")
axd.set_ylabel(r"$|\langle c|v_x|v_1\rangle|^2$" "\n" r"(eV$^2$ Å$^2$)",fontsize=11.5); axd.tick_params(axis="y",labelsize=11)
axd.set_ylim(-0.62,np.nanmax(M2[keep])*1.12); axd.set_yticks([0,1,2])
axd.text(0.012,axd.get_ylim()[1]*0.96,"(b)",fontsize=12,va="top")
axd.set_xticks([0,xG,x0]); axd.set_xticklabels(["X",r"$\Gamma$","Y"],fontsize=12)
fig.subplots_adjust(left=0.15,right=0.80,bottom=0.07,top=0.98)
fig.savefig(ROOT/"figures/fig3_bands.png",dpi=300); plt.close(fig)
# ---------------- Fig. 4: symmetry residuals and Wannier-centre relocation
ops=v6.pmmn_ops(L,[0.5,0.5,0]); res={}
for n in v6.OPNAMES[1:]:
    om=v6.orbital_map(ops[n],C,L,orbs); res[n]=v6.op_residual(ops[n],om,Rvec,HR,degen)*1e3
names=["C2x","sy","sz","C2y","C2z","I","sx"]
lbl={"C2x":r"$2_{1x}$","sy":r"$\sigma_y$","sz":r"$n_z$","C2y":r"$2_{1y}$","C2z":r"$C_{2z}$","I":r"$I$","sx":r"$\sigma_x$"}
te=[a[1] for a in atoms if a[0]=='Te']; ta=[a[1] for a in atoms if a[0]=='Ta']
seeds=ta+te+[np.array([0.5,0.5,0.0108789126880766])*L,np.array([0,0,0.9891210873119234])*L]
def wrap(v,c):  # put v into the periodic image nearest to c
    d=v-c; return c+d-np.round(d/L)*L
def yz(r): y=r[1]-L[1] if r[1]>L[1]/2 else r[1]; z=r[2]-L[2] if r[2]>L[2]/2 else r[2]; return y,z
fig=plt.figure(figsize=(7.2,5.2))
gs=fig.add_gridspec(2,1,height_ratios=[1.25,1.1],hspace=0.38)
axb=fig.add_subplot(gs[0]); axc=fig.add_subplot(gs[1],sharex=axb)
for e,r in atoms:
    y,z=yz(r); axb.scatter(y,z,s=75,facecolor="none",edgecolor=COL[e],lw=1.2,zorder=2)
for i in range(nw):
    c=wrap(C[i],seeds[i]); ys,zs=yz(seeds[i]); yc,zc=ys+(c-seeds[i])[1],zs+(c-seeds[i])[2]
    el="Ta" if i<4 else ("Te" if i<14 else "Pd")
    if np.hypot(yc-ys,zc-zs)>0.3:
        f=0.25 if el=="Te" else 1.0          # Te-seeded offsets drawn at 1/4 length: direction only
        ye,ze=ys+f*(yc-ys),zs+f*(zc-zs)
        axb.add_patch(FancyArrowPatch((ys,zs),(ye,ze),arrowstyle="-|>",mutation_scale=9,lw=1.1,color=COL[el],zorder=3,shrinkA=4,shrinkB=0))
    axb.scatter(yc,zc,s=26,marker={"Ta":"o","Te":"s","Pd":"D"}[el],c=COL[el],edgecolor="k",lw=0.4,zorder=4)
axb.axhline(0,c="0.85",lw=0.6,zorder=0)
axb.set_ylabel(r"$z$ (Å)",fontsize=11); axb.set_ylim(-2.7,2.7); axb.tick_params(labelsize=10)
axb.set_title(r"(a) side view ($yz$): Wannier centres relative to their projection seeds",loc="left",fontsize=10.5)
plt.setp(axb.get_xticklabels(),visible=False)
for xm in [0,L[0]/2]: axc.axhline(xm,ls="--",c="0.5",lw=0.8,zorder=0)
for xm,t in [(0,r"$\sigma_x$ ($x=0$)"),(L[0]/2,r"$\sigma_x$ ($x=a/2$)")]:
    axc.text(9.9,xm+0.12,t,fontsize=8.5,color="0.35",va="bottom",ha="right")
for e,r in atoms:
    y,_=yz(r); x=((r[0]+L[0]/4)%L[0])-L[0]/4
    axc.scatter(y,x,s=75,facecolor="none",edgecolor=COL[e],lw=1.2,zorder=2)
for i in range(nw):
    el="Ta" if i<4 else ("Te" if i<14 else "Pd")
    y,_=yz(C[i]); x=((C[i][0]+L[0]/4)%L[0])-L[0]/4
    axc.scatter(y,x,s=26,marker={"Ta":"o","Te":"s","Pd":"D"}[el],c=COL[el],edgecolor="k",lw=0.4,zorder=4)
axc.set_xlabel(r"$y$ (Å)",fontsize=11); axc.set_ylabel(r"$x$ (Å)",fontsize=11); axc.tick_params(labelsize=10)
axc.set_ylim(-1.2,3.1); axc.set_xlim(-10,10)
axc.set_title(r"(b) top view ($yx$): Wannier centres relative to the $\sigma_x$ mirror planes",loc="left",fontsize=10.5)
hs=[]
for e,m,n in [("Ta","o",r"Ta $d_{z^2}$ Wannier centre"),("Te","s",r"Te $p_x$-seeded Wannier centre"),("Pd","D",r"Pd $d_{xz}$ Wannier centre")]:
    hs.append(axc.scatter([],[],marker=m,c=COL[e],edgecolor="k",lw=0.4,s=26,label=n))
for e,n in [("Ta","Ta atom"),("Te","Te atom"),("Pd","Pd atom")]:
    hs.append(axc.scatter([],[],s=70,facecolor="none",edgecolor=COL[e],lw=1.2,label=n))
from matplotlib.lines import Line2D
hs.append(Line2D([0],[0],color=COL["Te"],lw=1.1,marker=">",markersize=5,label=r"Te-seeded offset: direction only, length $\times$0.25"))
hs.append(Line2D([0],[0],color=COL["Pd"],lw=1.1,marker=">",markersize=5,label="Pd offset: true length"))
fig.legend(handles=hs,loc="lower center",ncol=3,fontsize=8.5,frameon=False,bbox_to_anchor=(0.5,-0.10))
fig.savefig(ROOT/"figures/fig4_2d_projections.png",dpi=300,bbox_inches="tight"); plt.close(fig)
print({n:round(v,1) for n,v in res.items()})
