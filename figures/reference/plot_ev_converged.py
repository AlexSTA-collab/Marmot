#!/usr/bin/env python3
"""Figure H.28 as six separate panels: through-thickness convergence of the elastic and viscoelastic resolved references
(stiff, alpha = 10, 60 x 60), the change of the layer from n = 5 to 10 ('succ/10') and from n = 10 to 20 ('succ/20') in
e_u, e_[u], e_t, from scores_ev60.json (score_ev60.py). One panel per measure and layer, in the style of the convergence
panels of the paper, so that the figure can be set as two captioned rows of subfigures, and a one-row legend. The grey
curve is the error of the interface model (EIE, YIQUAD4) against the layer with n = 20, stiff column of the paper's
tables (IJES_submission_Revision/tables/data/{elastic_table_angle10,visco10_table_angle10}.md), criterion (ii).
usage: plot_ev_converged.py <out dir>   ->  <out>/{elastic,visco}_{global_relative_error,displacement_jump,traction}.pdf,
                                            <out>/legend_ev.pdf"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import LogLocator  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
S = json.loads((HERE / "scores_ev60.json").read_text())
TABLES = HERE.parents[1] / "IJES_submission_Revision/tables/data"


def model_errors(md):
    """{measure index: [stiff YIQUAD4 error over HS]} from a matched-convergence table of the paper."""
    out, k = {}, -1
    for line in (TABLES / md).read_text().splitlines():
        if line.startswith("## "):
            k += 1; out[k] = {}
        elif line.startswith("| 0") and k >= 0:
            c = [s.strip() for s in line.strip("|").split("|")]; out[k][c[0]] = float(c[1])
    return {k: [v[x] for x in HS] for k, v in out.items()}



plt.rcParams.update({
    "text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"], "font.size": 16,
    "axes.labelsize": 28, "axes.titlesize": 32, "xtick.labelsize": 28, "ytick.labelsize": 28, "legend.fontsize": 18,
    "savefig.format": "pdf", "savefig.bbox": "tight", "figure.dpi": 100,
})
HS = ["0.01", "0.02", "0.04", "0.06", "0.08", "0.1"]
h = np.array([float(x) for x in HS])
PANELS = (("global_relative_error", r"$e_u$"), ("displacement_jump", r"$e_{[\![u]\!]}$"), ("traction", r"$e_t$"))
PAIRS = {"succ/10": ("#4292c6", "o", "--", r"$n$: 5/10"),
         "succ/20": ("#08306b", "s", "-", r"$n$: 10/20")}
MODEL = ("0.35", "^", ":", "interface model")
MD = {"elastic": "elastic_table_angle10.md", "visco": "visco10_table_angle10.md"}
for layer in ("elastic", "visco"):
    em = model_errors(MD[layer])
    for k, (stem, title) in enumerate(PANELS):
        fig, ax = plt.subplots(figsize=(7, 6))
        plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.12)
        for key, (c, m, ls, _) in PAIRS.items():
            e = [S[f"{layer}/{x}/{key}"][k] for x in HS]
            ax.loglog(h, e, linestyle=ls, linewidth=2, marker=m, markersize=7, color=c)
        ax.loglog(h, em[k], linestyle=MODEL[2], linewidth=2, marker=MODEL[1], markersize=7, color=MODEL[0])
        ax.set_xlabel(r"$h$ [mm]")
        ax.set_title(title, fontsize=32, pad=15)
        ax.grid(True, alpha=0.3)
        ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
        ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
        fig.savefig(OUT / f"{layer}_{stem}.pdf"); fig.savefig(OUT / f"{layer}_{stem}.png", dpi=150)
        plt.close(fig)
handles = [Line2D([], [], color=c, marker=m, ls=ls, lw=2, ms=9, label=lab) for c, m, ls, lab in (*PAIRS.values(), MODEL)]
fig = plt.figure(figsize=(1.7, 1.1))                    # one column, beside the panels, as the legend of Fig. 19
fig.legend(handles=handles, loc="center", ncol=1, frameon=False, fontsize=13, handlelength=2.0)
fig.savefig(OUT / "legend_ev.pdf"); fig.savefig(OUT / "legend_ev.png", dpi=150)
print("wrote", OUT, "| elastic h=0.1 n=10 vs 20:", S["elastic/0.1/succ/20"])
