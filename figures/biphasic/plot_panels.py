#!/usr/bin/env python3
"""Panels of the biphasic benchmark in the style of the convergence figures of the paper (plot_visco.py).

    python3 plot_panels.py errors.json <figures dir>

Writes <figures dir>/angle_{0,10}/{global_relative_error,displacement_jump,traction}.pdf and legend_biphasic.pdf.
Colour: the model (full, [C] kept / suppressed); line style: the contrast E+/E-.
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator

res = json.loads(Path(sys.argv[1]).read_text())
OUT = Path(sys.argv[2])
H = [0.01, 0.02, 0.04, 0.06, 0.08, 0.1]

plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"],
                     "font.size": 16, "axes.labelsize": 28, "axes.titlesize": 32, "xtick.labelsize": 28,
                     "ytick.labelsize": 28, "legend.fontsize": 18, "savefig.format": "pdf",
                     "savefig.bbox": "tight", "figure.dpi": 100})

COLOR = {"full": "#1f77b4", "supp": "#d62728"}
MARKER = {"full": "o", "supp": "s"}
STYLE = {"c100": "-", "c10": "--"}
LABEL = {"full": r"$[\![C]\!]$ kept", "supp": r"$[\![C]\!]$ removed"}
CASE = {"c100": r"$E^+\!/E^-=100$", "c10": r"$E^+\!/E^-=10$"}
PANELS = {"global": ("global_relative_error", r"$\frac{\|u^\mathrm{h}_i-u^0_i\|_{L^2}}{\|u^0_i\|_{L^2}}$"),
          "e_u": ("displacement_jump",
                  r"$\frac{\|[\![u^\mathrm{h}_i]\!] - [\![u^0_i]\!]\|_{L^2}}{\|u^0_i\|_{L^2}}$"),
          "e_t": ("traction", r"$\frac{\|t^\mathrm{h}_i - t^0_i\|_{L^2}}{\|t^0_i\|_{L^2}}$")}

for angle in (0, 10):
    d = OUT / f"angle_{angle}"
    d.mkdir(parents=True, exist_ok=True)
    for key, (stem, title) in PANELS.items():
        fig, ax = plt.subplots(figsize=(7, 6))
        plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.12)
        for case in ("c100", "c10"):
            for model in ("full", "supp"):
                y = [res[f"{case}/{angle}/{h:g}"][model][key] for h in H]
                plt.loglog(H, y, linestyle=STYLE[case], linewidth=2, marker=MARKER[model], markersize=7,
                           color=COLOR[model])
        plt.xlabel(r"$h$ [mm]")
        plt.title(title, fontsize=32, pad=15)
        plt.grid(True, alpha=0.3)
        ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
        ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
        fig.savefig(d / f"{stem}.pdf")
        plt.close(fig)

fig = plt.figure(figsize=(3.2, 2.6))
handles = [Line2D([], [], color=COLOR[m], marker=MARKER[m], linestyle="-", linewidth=2, markersize=7, label=LABEL[m])
           for m in ("full", "supp")]
handles += [Line2D([], [], color="k", linestyle=STYLE[c], linewidth=2, label=CASE[c]) for c in ("c100", "c10")]
fig.legend(handles=handles, loc="center", frameon=False, fontsize=18)
fig.savefig(OUT / "legend_biphasic.pdf")
print("wrote", OUT)
