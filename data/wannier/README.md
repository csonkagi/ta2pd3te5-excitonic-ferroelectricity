# Wannier90 data (monolayer, r²SCAN)

| File | Description |
| --- | --- |
| `wannier90_hr.dat` | Real-space Hamiltonian, 16 Wannier functions, 145 lattice vectors (all in plane) |
| `wannier90.win` | Input: 4 Ta d_z², 10 Te p_x, 2 Pd d_xz projections; disentanglement window [−5.5, 1.3] eV; frozen window [−3.1, −2.4] eV; 28×5×1 mesh; guiding centres |
| `wannier90_wout_final_state.txt` | Final-state Wannier centres and spreads extracted from `wannier90.wout` |

To add (author): the full `wannier90.wout` and, optionally, `wannier90_band.dat`, `.kpt`, `.gnu`, `.agr`.

Do not use `wannier90_hr_B3u.dat` from earlier work: it was produced by a superseded script with an incorrect energy functional.

Note: the frozen window ends at −2.4 eV, 27 meV below the valence-band maximum (−2.373 eV), so the band-edge states v₁, c₁ and c₂ are fixed by the disentanglement rather than frozen.
