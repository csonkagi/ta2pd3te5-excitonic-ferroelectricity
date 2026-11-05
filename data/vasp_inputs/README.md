# VASP inputs

To add (author), one folder per calculation, each with INCAR, KPOINTS and POSCAR:

| Folder | Calculation |
| --- | --- |
| `bulk_scf/` | Bulk r²SCAN+rVV10, Γ-centred 7×28×5 mesh |
| `monolayer_relax/` | Monolayer relaxation, r²SCAN (no dispersion correction) |
| `monolayer_wannier/` | Wannier step on the full 28×5×1 mesh (LWANNIER90 = .TRUE., NBANDS = 144; Wannier90 uses num_bands = 134 after exclude_bands = 1-10); the self-consistent step ran with symmetry (45 irreducible k-points) |
| `bulk_bands_arpes/` | Bulk bands along Y–Γ–Y for Fig. S1 |

Common settings: METAGGA = R2SCAN, LASPH = .TRUE.; monolayer runs ENCUT = 550 eV. Bulk runs add the rVV10 correction: LUSE_VDW = .TRUE., BPARAM = 11.95, CPARAM = 0.0093.

POTCAR files are licensed by VASP Software GmbH and must not be redistributed. List here the PAW datasets used (for example `PAW_PBE Ta_pv`, `PAW_PBE Pd`, `PAW_PBE Te`), with their dates as they appear in the POTCAR TITEL lines.
