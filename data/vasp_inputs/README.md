# VASP inputs

To add (author), one folder per calculation, each with INCAR, KPOINTS and POSCAR:

| Folder | Calculation |
| --- | --- |
| `bulk_scf/` | Bulk r²SCAN+rVV10, Γ-centred 7×28×5 mesh |
| `monolayer_relax/` | Monolayer relaxation, r²SCAN+rVV10 |
| `monolayer_wannier/` | Self-consistent run on the full 28×5×1 mesh (ISYM = −1, LWANNIER90 = .TRUE., NBANDS = 134) |
| `bulk_bands_arpes/` | Bulk bands along Y–Γ–Y for Fig. S1 |

Common settings: METAGGA = R2SCAN, LUSE_VDW = .TRUE., BPARAM = 11.95, CPARAM = 0.0093, LASPH = .TRUE., ENCUT = 600 eV.

POTCAR files are licensed by VASP Software GmbH and must not be redistributed. List here the PAW datasets used (for example `PAW_PBE Ta_pv`, `PAW_PBE Pd`, `PAW_PBE Te`), with their dates as they appear in the POTCAR TITEL lines.
