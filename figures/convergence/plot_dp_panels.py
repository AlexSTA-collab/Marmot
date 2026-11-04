#!/usr/bin/env python3
"""Convergence panels of the Drucker-Prager interphase with tension cap (paper section on non-associated flow), in the
style of paper_scripts/figures/visco_matched/plot_visco.py: EIE against the matched 20-layer resolved layer, alpha = 0 or 10,
stiff, initial and soft families (colour and marker), phi = 30 with psi = 0 and 10 (line style), and a separate legend.
Scores: ~/psie26/plastic_bbar/dp/fam_<family>_<law>_cap/scores.json, entries eie/10/<h> = [e_u, e_[u], e_t] (score.py).
usage: plot_dp_panels.py <scores dir> <output dir> [angle = 10]"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import LogLocator  # noqa: E402

SRC, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({
    "text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"], "font.size": 16,
    "axes.labelsize": 28, "axes.titlesize": 32, "xtick.labelsize": 28, "ytick.labelsize": 28, "legend.fontsize": 18,
    "savefig.format": "pdf", "savefig.bbox": "tight", "figure.dpi": 100,
})
FAMILY_COLORS = {"stiff": "#1f77b4", "initial": "#ff7f0e", "soft": "#d62728"}
FAMILY_MARKERS = {"stiff": "o", "initial": "s", "soft": "^"}
# the paper reports phi = 30 with psi = 0 and 10; the runs vm (0, 0) and vmna (0, 10) exist as well
LAWS = {"psi0": (r"$\varphi=30^\circ$, $\psi=0^\circ$", "-"), "dp": (r"$\varphi=30^\circ$, $\psi=10^\circ$", "--")}
PANELS = (("global_relative_error", r"$\frac{\|u^\mathrm{h}_i-u^0_i\|_{L^2}}{\|u^0_i\|_{L^2}}$"),
          ("displacement_jump", r"$\frac{\|[\![u^\mathrm{h}_i]\!] - [\![u^0_i]\!]\|_{L^2}}{\|u^0_i\|_{L^2}}$"),
          ("traction", r"$\frac{\|t^\mathrm{h}_i - t^0_i\|_{L^2}}{\|t^0_i\|_{L^2}}$"))
HS = ["0.01", "0.02", "0.04", "0.06", "0.08", "0.1"]
h = np.array([float(x) for x in HS])
ANGLE = sys.argv[3] if len(sys.argv) > 3 else "10"
data = {(f, l): np.array([json.loads((SRC / f"fam_{f}_{l}_cap" / "scores.json").read_text())[f"eie/{ANGLE}/{x}"] for x in HS])
        for f in FAMILY_COLORS for l in LAWS}

for k, (stem, title) in enumerate(PANELS):
    fig, ax = plt.subplots(figsize=(7, 6))
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.12)
    for (f, l), e in data.items():
        plt.loglog(h, e[:, k], linestyle=LAWS[l][1], linewidth=2, marker=FAMILY_MARKERS[f], markersize=7,
                   color=FAMILY_COLORS[f])
    plt.xlabel(r"$h$ [mm]")
    plt.title(title, fontsize=32, pad=15)
    plt.grid(True, alpha=0.3)
    ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
    fig.savefig(OUT / f"{stem}.pdf")
    plt.close(fig)
    print(f"  {stem}: " + "  ".join(f"{f}/{l}={np.polyfit(np.log(h), np.log(e[:, k]), 1)[0]:.2f}"
                                    for (f, l), e in data.items()))

handles = [Line2D([], [], color=FAMILY_COLORS[f], marker=FAMILY_MARKERS[f], ls="", ms=9, label=f.capitalize())
           for f in FAMILY_COLORS]
handles += [Line2D([], [], color="0.25", ls=ls, lw=2, label=lab) for lab, ls in LAWS.values()]
fig = plt.figure(figsize=(14, 0.6))
fig.legend(handles=handles, loc="center", ncol=5, frameon=False, fontsize=18, handlelength=2.4, columnspacing=1.4)
fig.savefig(OUT / "legend_dp.pdf")
print("wrote", OUT)
