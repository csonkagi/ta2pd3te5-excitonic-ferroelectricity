"""Interactive, rotatable 3D version of Fig. 4 (self-contained HTML, works offline). Requires plotly."""
import numpy as np, plotly.graph_objects as go
from fig4_geometry import ROOT, L, nw, ATOMS, SEEDS, CENTRES, OFFSET, SPREADS, COL, ELW, nearest_sigma_x, bonds

MK = {"Ta": "circle", "Te": "square", "Pd": "diamond"}
LAB = {"Ta": "Ta d<sub>z²</sub>", "Te": "Te p<sub>x</sub>-seeded", "Pd": "Pd d<sub>xz</sub>"}
fig = go.Figure()
for xm in (0.0, L[0]/2):
    y0, y1, z0, z1 = -L[1]/2, L[1]/2, -2.8, 2.8
    fig.add_trace(go.Mesh3d(x=[xm]*4, y=[y0, y1, y1, y0], z=[z0, z0, z1, z1], i=[0, 0], j=[1, 2], k=[2, 3], color="#888888", opacity=0.12,
                            name="σ<sub>x</sub> mirror planes", legendgroup="planes", showlegend=bool(xm == 0), hoverinfo="name"))
bx, by, bz = [], [], []
for r1, r2 in bonds(): bx += [r1[0], r2[0], None]; by += [r1[1], r2[1], None]; bz += [r1[2], r2[2], None]
fig.add_trace(go.Scatter3d(x=bx, y=by, z=bz, mode="lines", line=dict(color="#aaaaaa", width=3), name="bonds < 3.05 Å", visible="legendonly", hoverinfo="skip"))
for e in ("Ta", "Te", "Pd"):
    P = np.array([r for el, r in ATOMS if el == e])
    fig.add_trace(go.Scatter3d(x=P[:, 0], y=P[:, 1], z=P[:, 2], mode="markers", name=f"{e} atom",
                               marker=dict(size=13, color=COL[e], opacity=0.45, line=dict(color=COL[e], width=1)),
                               hovertemplate=f"{e} atom<br>x=%{{x:.2f}} Å<br>y=%{{y:.2f}} Å<br>z=%{{z:.2f}} Å<extra></extra>"))
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
for i in range(nw):
    r = 0.3*np.sqrt(SPREADS[i])
    fig.add_trace(go.Mesh3d(x=(CENTRES[i, 0] + r*np.cos(u)*np.sin(v)).ravel(), y=(CENTRES[i, 1] + r*np.sin(u)*np.sin(v)).ravel(),
                            z=(CENTRES[i, 2] + r*np.cos(v)).ravel(), alphahull=0, color=COL[ELW[i]], opacity=0.15,
                            name="Wannier extent (0.3 × rms radius)", legendgroup="extent", showlegend=bool(i == 0), visible="legendonly", hoverinfo="skip"))
for e in ("Te", "Pd"):
    idx = [i for i in range(nw) if ELW[i] == e]; colr = COL["Te"] if e == "Te" else "#555555"
    lx, ly, lz = [], [], []
    for i in idx: lx += [SEEDS[i, 0], CENTRES[i, 0], None]; ly += [SEEDS[i, 1], CENTRES[i, 1], None]; lz += [SEEDS[i, 2], CENTRES[i, 2], None]
    fig.add_trace(go.Scatter3d(x=lx, y=ly, z=lz, mode="lines", line=dict(color=colr, width=5), name=f"{e}: seed atom → Wannier centre", legendgroup=f"arr{e}", hoverinfo="skip"))
    d = CENTRES[idx] - SEEDS[idx]; un = d/np.linalg.norm(d, axis=1)[:, None]
    fig.add_trace(go.Cone(x=CENTRES[idx, 0], y=CENTRES[idx, 1], z=CENTRES[idx, 2], u=un[:, 0], v=un[:, 1], w=un[:, 2], anchor="tip",
                          sizemode="absolute", sizeref=0.15, colorscale=[[0, colr], [1, colr]], showscale=False, legendgroup=f"arr{e}", showlegend=False, hoverinfo="skip"))
for e in ("Ta", "Te", "Pd"):
    idx = [i for i in range(nw) if ELW[i] == e]
    txt = [f"WF {i+1}: {LAB[e]}<br>spread {SPREADS[i]:.2f} Å² (rms {np.sqrt(SPREADS[i]):.2f} Å)<br>offset from seed {OFFSET[i]:.2f} Å"
           f"<br>distance to σ<sub>x</sub> plane {nearest_sigma_x(CENTRES[i, 0])[1]:.2f} Å" for i in idx]
    fig.add_trace(go.Scatter3d(x=CENTRES[idx, 0], y=CENTRES[idx, 1], z=CENTRES[idx, 2], mode="markers", name=f"{LAB[e]} Wannier centre",
                               marker=dict(size=6, color=COL[e], symbol=MK[e], line=dict(color="black", width=1)), text=txt, hovertemplate="%{text}<extra></extra>"))
rx, ry, rz = [], [], []
for i in range(nw):
    xm, dd = nearest_sigma_x(CENTRES[i, 0])
    if dd > 0.1: rx += [CENTRES[i, 0], xm, None]; ry += [CENTRES[i, 1]]*2 + [None]; rz += [CENTRES[i, 2]]*2 + [None]
fig.add_trace(go.Scatter3d(x=rx, y=ry, z=rz, mode="lines", line=dict(color="#d62728", width=7), name="offset from nearest σ<sub>x</sub> plane", hoverinfo="skip"))
ex, ey, ez = L[0] + 1.0, L[1], 6.0
ASP = {"true proportions": dict(x=ex/ey, y=1.0, z=ez/ey), "x stretched ×3": dict(x=3*ex/ey, y=1.0, z=ez/ey)}
CAMS = {"perspective": dict(eye=dict(x=1.4, y=-1.6, z=0.9)), "side view (yz)": dict(eye=dict(x=2.4, y=0.0, z=0.0)),
        "top view (xy)": dict(eye=dict(x=0.0, y=0.0, z=2.5), up=dict(x=1, y=0, z=0)), "along y (xz)": dict(eye=dict(x=0.0, y=-2.5, z=0.2))}
fig.update_layout(
    title=dict(text="Monolayer Ta<sub>2</sub>Pd<sub>3</sub>Te<sub>5</sub>: maximally localized Wannier centres relative to their projection seeds and the σ<sub>x</sub> mirror planes"
               "<br><sup>Drag to rotate, scroll to zoom, hover for values; click legend entries to show or hide layers (bonds and Wannier extent are hidden initially).</sup>", x=0.01, font=dict(size=15)),
    font=dict(family="Carlito, Calibri, Arial, sans-serif", size=13),
    scene=dict(xaxis=dict(title="x (Å), chain axis", range=[-L[0]/4 - 0.5, 0.75*L[0] + 0.5], tickvals=[0, L[0]/2], ticktext=["0", "a/2"]),
               yaxis=dict(title="y (Å)", range=[-L[1]/2, L[1]/2]), zaxis=dict(title="z (Å)", range=[-3, 3]),
               aspectmode="manual", aspectratio=ASP["true proportions"], camera=CAMS["perspective"]),
    legend=dict(x=0.0, y=0.0, xanchor="left", yanchor="bottom", bgcolor="rgba(255,255,255,0.8)"), margin=dict(l=0, r=0, t=70, b=0), height=780,
    updatemenus=[dict(type="buttons", direction="right", x=0.52, y=1.0, xanchor="left", yanchor="top",
                      buttons=[dict(label=k, method="relayout", args=[{"scene.camera": c}]) for k, c in CAMS.items()]),
                 dict(type="buttons", direction="right", x=0.99, y=1.0, xanchor="right", yanchor="top",
                      buttons=[dict(label=k, method="relayout", args=[{"scene.aspectratio": a}]) for k, a in ASP.items()])])
out = ROOT/"figures/fig4_interactive.html"
fig.write_html(str(out), include_plotlyjs=True, full_html=True, config=dict(displaylogo=False))
print("wrote", out)
