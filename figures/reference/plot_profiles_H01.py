#!/usr/bin/env python3
"""Figure I.31 as four separate panels: strain of each element row through the thickness of the C3D8B reference, divided
by its mean over the thickness, at H = 0.1 and 100 MPa, for h = 0.1 and 0.01 mm (shear strain and normal strain), with
n = 5, 10, 20; data row_strain/h<h>_H<H>_n<n>.json, as in the ref_profiles_H panel of appendix_reference_figs.py. Set as two
captioned rows with a one-column legend on the right, like Figures H.28 and H.29.
usage: plot_profiles_H01.py <out dir>  ->  <out>/H<H>_h<h>_{gamma,eps}.pdf, <out>/legend_n.pdf"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
ROWS = HERE / "row_strain"
plt.rcParams.update({
    "text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"], "font.size": 16,
    "axes.labelsize": 28, "axes.titlesize": 32, "xtick.labelsize": 28, "ytick.labelsize": 28,
    "savefig.format": "pdf", "savefig.bbox": "tight", "figure.dpi": 100,
})
C = "#1f77b4"
STYLE = {5: ((0, (1, 1.5)), 1.8, 0.45), 10: ((0, (5, 2)), 1.8, 0.7), 20: ("-", 2.0, 1.0)}
PANELS = (("gamma", r"$\gamma_k/\bar\gamma$"), ("eps", r"$\varepsilon_{nn,k}/\bar\varepsilon_{nn}$"))
# one x-range per thickness and strain, shared by the two hardening moduli, as in ref_profiles_H, so that the uniform
# shear strain at H = 100 MPa is not stretched over the whole axis
XLIM = {}
for h in ("0.1", "0.01"):
    for key, _ in PANELS:
        v = [x for H in ("0.1", "100") for n in (5, 10, 20)
             for x in json.loads((ROWS / f"h{h}_H{H}_n{n}.json").read_text())[key]]
        pad = 0.05 * (max(v) - min(v))
        XLIM[h, key] = (min(v) - pad, max(v) + pad)
for H, h in [(H, h) for H in ("0.1", "100") for h in ("0.1", "0.01")]:
    for key, xlabel in PANELS:
        fig, ax = plt.subplots(figsize=(7, 6))
        plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.12)
        for n in (5, 10, 20):
            d = json.loads((ROWS / f"h{h}_H{H}_n{n}.json").read_text())
            ls, lw, al = STYLE[n]
            ax.stairs(d[key], d["zeta"], orientation="horizontal", baseline=None, color=C, ls=ls, lw=lw, alpha=al)
        ax.axvline(1.0, color="0.6", lw=0.8, zorder=0)
        ax.grid(True, alpha=0.3); ax.set_ylim(0, 1); ax.set_xlim(*XLIM[h, key])
        ax.ticklabel_format(axis="x", useOffset=False)
        ax.set_xlabel(xlabel); ax.set_ylabel(r"$\zeta$")
        fig.savefig(OUT / f"H{H}_h{h}_{key}.pdf")
        plt.close(fig)
handles = [Line2D([], [], color=C, ls=STYLE[n][0], lw=STYLE[n][1], alpha=STYLE[n][2], label=rf"$n={n}$") for n in (5, 10, 20)]
fig = plt.figure(figsize=(1.5, 1.1))                    # one column, beside the panels, as in Figs. H.28 and H.29
fig.legend(handles=handles, loc="center", ncol=1, frameon=False, fontsize=13, handlelength=2.0)
fig.savefig(OUT / "legend_n.pdf")
print("wrote", OUT)
