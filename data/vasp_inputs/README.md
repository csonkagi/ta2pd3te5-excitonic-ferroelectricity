# VASP inputs

Monolayer calculations of Hong Tang (VASP, PAW, r²SCAN without dispersion correction, ENCUT = 550 eV), one folder per calculation:

| Folder | Calculation | Files in the repository |
| --- | --- | --- |
| `monolayer_relax/` | relaxation of the monolayer | `INCAR`, `KPOINTS`, `POSCAR`, `CONTCAR_1L_r2scan_relaxed` |
| `monolayer_wannier/` | self-consistent step (`INCAR_SCF`, symmetry-reduced: 45 irreducible k-points) and Wannier step (`INCAR_WAN`, full 28×5×1 mesh, LWANNIER90 = .TRUE., NBANDS = 144; Wannier90 uses num_bands = 134 after exclude_bands = 1-10) | `INCAR_SCF`, `INCAR_WAN`, `INCAR` (identical to `INCAR_WAN`), `KPOINTS`, `POSCAR` |
| `monolayer_bandstr/` | DFT band structure along the high-symmetry path | `KPOINTS` (band path); its `INCAR` and `POSCAR` are not available |
| `bulk_scf/` | bulk r²SCAN+rVV10 calculations for Fig. S1 (ENCUT = 550 eV) | not yet available |

Common settings: METAGGA = R2SCAN, LASPH = .TRUE., ENCUT = 550 eV. Bulk runs add the rVV10 correction: LUSE_VDW = .TRUE., BPARAM = 11.95, CPARAM = 0.0093.

Not published: `OUTCAR`, `vasprun.xml`, `WAVECAR` and other output files (see `.gitignore`), and the `POTCAR` files, which are licensed by VASP Software GmbH and must not be redistributed. PAW datasets used (POTCAR TITEL lines): to be added.
