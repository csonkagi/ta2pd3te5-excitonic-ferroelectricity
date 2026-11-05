# The other Wannier models of SI Table S1

Each folder holds `wannier90_hr.dat` (Wannier90 format) and `wannier90_centres.wout` (cell vectors and final centres and spreads, the format `ta2pd3te5_channel_mf_v6.py` reads). Orbital order as in the main model: 4 Ta, 10 Te, 2 Pd functions.

| Folder | Model | Options for v6 and `channel_tr_blocks.py` |
| --- | --- | --- |
| `vasp_projection_only/` | Hong Tang's VASP r²SCAN overlaps, Wannier90 4.0.3 with `num_iter = 0` and the original windows | `--nocc 12 --symtol 0.1` (three operations hold only to 46–77 meV) |
| `vasp_projection_only_frozen_raised/` | as above, frozen window raised to −2.05 eV; band edges exact at Γ | `--nocc 10 --symtol 0.1` (σ_z holds only to 96 meV; the T-even/T-odd coupling in χ is 6%, so the time-reversal split is approximate) |
| `qe_symmetry_adapted/` | Hong Tang's Quantum ESPRESSO r²SCAN run (ONCV PseudoDojo pseudopotentials, 90 Ry, 28×6×1), Wannier90 with site symmetry, no frozen window | `--nocc 12` (all eight operations exact within 9 meV; semimetallic, band edges 30–140 meV off) |

The main model of the paper is in `../wannier/`. Commands that regenerate the tables: `reproduce.sh`.
