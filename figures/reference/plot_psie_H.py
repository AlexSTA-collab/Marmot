#!/usr/bin/env python3
"""Figure I.32 as three separate panels: error of the PSIE against the C3D8B layer with n = 20 over the thickness sweep,
H = 0.1 MPa (dashed, open markers) and H = 100 MPa (solid, filled), stiff family, alpha = 10; data
plastic_bbar/scores_H100_sq.json, as in the ref_psie_H panel of appendix_reference_figs.py. Panel style and side legend
as in Figures H.28 and H.29.
usage: plot_psie_H.py <out dir>  ->  <out>/{global_relative_error,displacement_jump,traction}.pdf, <out>/legend_H.pdf"""
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
SC = json.loads((HERE.parent / "plastic_bbar" / "scores_H100_sq.json").read_text())
plt.rcParams.update({
    "text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"], "font.size": 16,
    "axes.labelsize": 28, "axes.titlesize": 32, "xtick.labelsize": 28, "ytick.labelsize": 28,
    "savefig.format": "pdf", "savefig.bbox": "tight", "figure.dpi": 100,
})
HS = ["0.01", "0.02", "0.04", "0.06", "0.08", "0.1"]
h = np.array([float(x) for x in HS])
C = "#1f77b4"
PANELS = (("global_relative_error", "e_u_global", r"$e_u$"), ("displacement_jump", "e_jump", r"$e_{[\![u]\!]}$"),
          ("traction", "e_traction", r"$e_t$"))
CASES = {"0.1": ("--", "white", r"$H=0.1$~MPa"), "100": ("-", C, r"$H=100$~MPa")}
for stem, key, title in PANELS:
    fig, ax = plt.subplots(figsize=(7, 6))
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.12)
    for H, (ls, fill, _) in CASES.items():
        ax.loglog(h, [SC[f"psie/stiff/10/{H}/{x}"][key] for x in HS], color=C, ls=ls, lw=2, marker="o", markersize=7, mfc=fill)
    ax.set_xlabel(r"$h$ [mm]")
    ax.set_title(title, fontsize=32, pad=15)
    ax.grid(True, alpha=0.3)
    ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
    fig.savefig(OUT / f"{stem}.pdf")
    plt.close(fig)
handles = [Line2D([], [], color=C, ls=ls, lw=2, marker="o", ms=9, mfc=fill, label=lab) for ls, fill, lab in CASES.values()]
fig = plt.figure(figsize=(1.9, 0.9))                    # one column, beside the panels, as in Figs. H.28 and H.29
fig.legend(handles=handles, loc="center", ncol=1, frameon=False, fontsize=13, handlelength=2.0)
fig.savefig(OUT / "legend_H.pdf")
print("wrote", OUT, "| H=100 h=0.01:", SC["psie/stiff/10/100/0.01"])
