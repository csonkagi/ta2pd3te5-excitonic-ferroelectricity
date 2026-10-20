#!/usr/bin/env bash
# Recompute all tables and redraw all figures of the paper.
set -e
cd "$(dirname "$0")"
python code/compute_results.py
python code/combine_fig1.py
python code/combine_fig2.py
python code/make_fig3_fig4.py
python code/make_fig5.py
python code/build_figS1.py
echo "Done: results/ and figures/ are up to date."
