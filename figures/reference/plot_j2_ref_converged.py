#!/usr/bin/env python3
"""Figure H.29 as six separate panels: through-thickness convergence of the C3D8B reference of the J2 interphase (stiff,
alpha = 10, 60 x 60), the change of the layer from n = 5 to 10 and from n = 10 to 20 in e_u, e_[u], e_t, for H = 0.1 and
100 MPa, from appendix_bbar/bbar.json (the data of the ref_converged panel of appendix_reference_figs.py). Same style and
legend as Figure H.28 (refinement_ev60/plot_ev_converged.py). The grey curve is the error of the PSIE against the layer
with n = 20, stiff, from the paper's IJES_submission_Revision/tables/data/j2_scores_H100_sq.json, criterion (ii).
usage: plot_j2_ref_converged.py <out dir>   ->  <out>/H<0.1|100>_{global_relative_error,displacement_jump,traction}.pdf"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.ticker import LogLocator  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
B = json.loads((HERE / "appendix_bbar" / "bbar.json").read_text())
P = json.loads((HERE.parents[1] / "IJES_submission_Revision/tables/data/j2_scores_H100_sq.json").read_text())
PKEYS = ("e_u_global", "e_jump", "e_traction")
plt.rcParams.update({
    "text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"], "font.size": 16,
    "axes.labelsize": 28, "axes.titlesize": 32, "xtick.labelsize": 28, "ytick.labelsize": 28, "legend.fontsize": 18,
    "savefig.format": "pdf", "savefig.bbox": "tight", "figure.dpi": 100,
})
HS = ["0.01", "0.02", "0.04", "0.06", "0.08", "0.1"]
h = np.array([float(x) for x in HS])
PANELS = (("global_relative_error", r"$e_u$"), ("displacement_jump", r"$e_{[\![u]\!]}$"), ("traction", r"$e_t$"))
PAIRS = {"10": ("#4292c6", "o", "--"), "20": ("#08306b", "s", "-")}   # n = 5/10 and n = 10/20, as in Fig. H.28
for H in ("0.1", "100"):
    for k, (stem, title) in enumerate(PANELS):
        fig, ax = plt.subplots(figsize=(7, 6))
        plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.12)
        for key, (c, m, ls) in PAIRS.items():
            ax.loglog(h, [B[f"{x}_{H}"]["succ"][key][k] for x in HS], linestyle=ls, linewidth=2, marker=m,
                      markersize=7, color=c)
        ax.loglog(h, [P[f"psie/stiff/10/{H}/{x}"][PKEYS[k]] for x in HS], linestyle=":", linewidth=2, marker="^",
                  markersize=7, color="0.35")
        ax.set_xlabel(r"$h$ [mm]")
        ax.set_title(title, fontsize=32, pad=15)
        ax.grid(True, alpha=0.3)
        ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
        ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
        fig.savefig(OUT / f"H{H}_{stem}.pdf"); fig.savefig(OUT / f"H{H}_{stem}.png", dpi=150)
        plt.close(fig)
print("wrote", OUT, "| H=100 h=0.1 n=10 vs 20:", B["0.1_100"]["succ"]["20"])
