#!/usr/bin/env python3
"""Global response of the Drucker-Prager interphase (initial family, alpha = 10, h = 0.01, phi = 30, psi = 10): panel (a) of
the paper's figure on the global response, F / F0 against the applied displacement, and panel (b), the departure of the EIE
from the resolved layer, |F - F0| / max|F0|. Reads global/reaction_dp.json (extract_dp_global.py).
usage: plot_dp_global.py <reaction_dp.json> <output dir>   -> druckerprager_force.pdf, druckerprager_departure.pdf"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.ticker import LogLocator  # noqa: E402

plt.rcParams.update({
    "text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"], "font.size": 16,
    "axes.labelsize": 28, "axes.titlesize": 32, "xtick.labelsize": 28, "ytick.labelsize": 28, "legend.fontsize": 18,
    "savefig.format": "pdf", "savefig.bbox": "tight", "figure.dpi": 100,
})
d = json.loads(Path(sys.argv[1]).read_text())["initial_dp"]
out = Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
uF, FF = np.array(d["full"][0]), np.array(d["full"][1])
uE, FE = np.array(d["EIE"][0]), np.array(d["EIE"][1])
umax, F0 = uF[-1], np.abs(FF).max()
x = uF / umax

fig, ax = plt.subplots(figsize=(6, 6))
plt.subplots_adjust(left=0.2, right=0.95, top=0.88, bottom=0.15)
slope = (FF[1] - FF[0]) / (uF[1] - uF[0]) if len(FF) > 1 else 1.0
ax.plot(x, FF / F0, "-", color="k", lw=2.5, label="resolved layer")
ax.plot(uE / umax, FE / F0, "--", color="#1f77b4", lw=2.5, label="EIE")
ax.plot([0, 1.05], [0, slope * 1.05 * umax / F0], ":", color="0.45", lw=1.5, label="elastic tangent")
ax.set_ylim(-0.05, 1.05)
ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_xlabel(r"$\bar u_x\,/\,\bar u_x^{\max}$")
ax.set_title(r"$\frac{F}{F^0}$", pad=15)
ax.grid(True, alpha=0.3)
ax.legend(loc="upper left", frameon=False, fontsize=14, handlelength=1.8, labelspacing=0.25)
fig.savefig(out / "druckerprager_force.pdf")
plt.close(fig)

fig, ax = plt.subplots(figsize=(6, 6))
plt.subplots_adjust(left=0.2, right=0.95, top=0.88, bottom=0.15)
dep = np.abs(FE - FF) / F0
k = x > 0
ax.semilogy(x[k], dep[k], "-", color="#1f77b4", lw=2.5, label="EIE")
ax.set_xlabel(r"$\bar u_x\,/\,\bar u_x^{\max}$")
ax.set_title(r"$\frac{|F-F^0|}{\max|F^0|}$", pad=15)
ax.grid(True, which="both", alpha=0.3)
ax.legend(loc="lower right", frameon=False, fontsize=14)
fig.savefig(out / "druckerprager_departure.pdf")
print("departure at the end", dep[-1], "F0", F0)
