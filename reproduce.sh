#!/usr/bin/env bash
# Recompute all tables and redraw all figures of the paper.
set -e
cd "$(dirname "$0")"
python code/compute_results.py
python code/combine_fig1.py
python code/combine_fig2.py
python code/make_fig3_fig4.py
python code/make_fig4.py
python code/make_fig4_interactive.py
python code/make_fig5.py
python code/build_figS1.py
V=code/ta2pd3te5_channel_mf_v6.py; M=data/wannier_models
python code/channel_tr_blocks.py --v6 $V --hr data/wannier/wannier90_hr.dat --wout data/wannier/wannier90_wout_final_state.txt --nocc 12 --rcut 3.4 4.0 4.6 --out results/channel_tr_blocks_mlwf.csv
python code/channel_tr_blocks.py --v6 $V --hr $M/vasp_projection_only/wannier90_hr.dat --wout $M/vasp_projection_only/wannier90_centres.wout --nocc 12 --rcut 3.4 --symtol 0.1 --out results/channel_tr_blocks_projection_only.csv
python code/channel_tr_blocks.py --v6 $V --hr $M/vasp_projection_only_frozen_raised/wannier90_hr.dat --wout $M/vasp_projection_only_frozen_raised/wannier90_centres.wout --nocc 10 --rcut 3.4 --symtol 0.1 --out results/channel_tr_blocks_frozen_raised.csv
python code/channel_tr_blocks.py --v6 $V --hr $M/qe_symmetry_adapted/wannier90_hr.dat --wout $M/qe_symmetry_adapted/wannier90_centres.wout --nocc 12 --rcut 3.4 --out results/channel_tr_blocks_qe_symmetric.csv
echo "Done: results/ and figures/ are up to date."
