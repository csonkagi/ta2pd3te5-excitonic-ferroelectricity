# Ta₂Pd₃Te₅ excitonic ferroelectricity: data and code

Data, inputs, outputs and code for:

> G. I. Csonka, H. Tang and A. Ruzsinszky, *When can an excitonic insulator be a purely electronic ferroelectric? Symmetry conditions and the chiral, nonpolar order of monolayer Ta₂Pd₃Te₅* (2026).

The repository reproduces every number and figure in the paper from the Wannier90 Hamiltonian of monolayer Ta₂Pd₃Te₅ (r²SCAN+rVV10, 16 Wannier functions).

## Quick start

```bash
pip install -r requirements.txt
bash reproduce.sh          # about 5 minutes on a laptop
```

`reproduce.sh` recomputes the tables in `results/` and redraws every figure in `figures/`.

## Contents

| Path | Contents | Paper |
| --- | --- | --- |
| `data/wannier/wannier90_hr.dat` | 16-orbital tight-binding Hamiltonian (Wannier90 output) | all monolayer results |
| `data/wannier/wannier90.win` | Wannier90 input: projections, windows, k-mesh, atomic positions | Methods II.B |
| `data/wannier/wannier90_wout_final_state.txt` | Final Wannier centres and spreads | Fig. 4, Sec. III.E |
| `data/structures/` | Relaxed monolayer and bulk structures | Figs. 1–2, Methods II.A |
| `data/vasp_inputs/` | VASP input files (INCAR, KPOINTS, POSCAR) | Methods II.A |
| `data/images/` | VESTA renderings used for Figs. 1 and 2 | Figs. 1–2 |
| `data/external/` | Third-party image used for Fig. S1 (CC BY 4.0, see `ATTRIBUTION.md`) | Fig. S1 |
| `code/ta2pd3te5_channel_mf_v6.py` | Symmetry test of the Wannier Hamiltonian and channel-resolved mean field | Methods II.C–D, SI S2 |
| `code/compute_results.py` | Recomputes Table III, Table I (right columns), Γ-point interband dipoles, Wannier-centre analysis, critical couplings | Tables I, III; Figs. 3–5 |
| `code/make_fig3_fig4.py`, `code/make_fig5.py` | Figure 3 (bands, pairings at Γ, interband velocity matrix elements along X–Γ–Y) and Figure 5 | Figs. 3, 5 |
| `code/fig4_geometry.py` | Shared geometry for Fig. 4: atoms, projection seeds, Wannier centres (folded into one cell) | Fig. 4 |
| `code/make_fig4.py` | Static Fig. 4 (3D perspective, matches the default view of the interactive version) | Fig. 4 |
| `code/make_fig4_interactive.py` | Rotatable 3D Fig. 4 as self-contained HTML (`figures/fig4_interactive.html`, requires plotly) | Fig. 4 |
| `code/combine_fig1.py`, `code/combine_fig2.py` | Assemble Figs. 1–2 from the VESTA renderings | Figs. 1–2 |
| `code/build_figS1.py` | Rebuild Fig. S1 on calibrated axes | Fig. S1 |
| `results/` | Computed tables (CSV) and a run log | Tables I, III; Figs. 4–5 |
| `figures/` | Figures as used in the paper, plus `fig4_interactive.html` (open in any browser: rotate, zoom, hover, toggle layers) | Figs. 1–5, S1 |

## Results files

| File | Contents |
| --- | --- |
| `results/symmetry_residuals.csv` | Residual r_g of each nontrivial Pmmn operation (Table III) |
| `results/gamma_eigenvalues.csv` | Γ-point energies and eigenvalues of v₁, c₁, c₂ under the exact operations (Table I) |
| `results/interband_dipoles_gamma.csv` | Squared interband velocity matrix elements at Γ for the c₁–v₁ (Aᵤ) and c₂–v₁ (B₃ᵤ) pairs (Fig. 3b) |
| `results/wannier_centres.csv` | Centre, spread, rms radius, offset from projection seed and distance to the σ_x planes for each Wannier function (Fig. 4) |
| `results/v0c_merged_channels.csv` | Critical couplings V₀ᶜ of the merged channels for r_c = 3.4, 4.0, 4.6 Å (Fig. 5) |
| `results/channel_mf_v6_run_rc3.4.log` | Full run log of `ta2pd3te5_channel_mf_v6.py`, including self-consistent solutions |

## Running the channel analysis directly

```bash
python code/ta2pd3te5_channel_mf_v6.py --hr data/wannier/wannier90_hr.dat \
       --wout data/wannier/wannier90_wout_final_state.txt --rcut 3.4 --V0 0.5 1.0 1.5
```

The program first tests every Pmmn operation on the Hamiltonian, reports which irreducible representations it can separate, then computes critical couplings and self-consistent solutions only for those channels. Run with `--help` for all options.

## Software

DFT: VASP 6 with PAW potentials, r²SCAN + rVV10 (b = 11.95, C = 0.0093), 600 eV cutoff. Wannierization: Wannier90. Analysis: Python 3 with NumPy, SciPy, Matplotlib and Pillow (`requirements.txt`). Figures use the Carlito font if installed and fall back otherwise.

VASP POTCAR files are licensed and are not included; `data/vasp_inputs/README.md` lists the PAW datasets used.

## Licence

Code: MIT (`LICENSE`). Data, results and figures: CC BY 4.0 (`LICENSE-DATA.md`). The third-party image in `data/external/` remains under its original CC BY 4.0 licence with the attribution given there.

## Citation

See `CITATION.cff`. Please cite the paper and this repository.
