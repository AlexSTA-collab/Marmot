#!/usr/bin/env python3
"""Convergence of the Drucker-Prager references of the paper through the thickness (H = 0.075 G, phi = 30, psi = 0 and 10,
stiff and initial family, alpha = 10, roller from t = 1, C3D8 layer on the 60 x 60 mesh): change from n = 5 to 10 and from
n = 10 to 20 elements over h, and the error of the EIE against the layer with n = 20 (grey, criterion (ii)).
One panel per case and measure: <case>_{global_relative_error,displacement_jump,traction}.{pdf,png}, legend_dp_ref.{pdf,png}.
Data: <sweep json> (~/psie26/tools/plastic/refconv_s75_sweep.py, keys <case>/<h>/<a>-<b>) and the paper's
IJES_submission_Revision/tables/data/dp/<case>.json, keys eie/10/<h>; each entry [e_u, e_[u], e_t].
usage: plot_dp_refconv_s75.py <sweep json> <output dir>"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import LogLocator  # noqa: E402

rc, OUT = json.loads(Path(sys.argv[1]).read_text()), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
DATA = Path(__file__).resolve().parents[2] / "IJES_submission_Revision/tables/data/dp"
plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"],
                     "figure.dpi": 100, "savefig.bbox": "tight", "font.size": 17, "axes.labelsize": 28,
                     "axes.titlesize": 32, "xtick.labelsize": 28, "ytick.labelsize": 28})
HS = ["0.01", "0.02", "0.04", "0.06", "0.08", "0.1"]
x = [float(h) for h in HS]
CASES = ("stiff_psi0", "stiff_dp", "initial_psi0", "initial_dp")
PANELS = (("global_relative_error", r"$e_u$"), ("displacement_jump", r"$e_{[\![u]\!]}$"), ("traction", r"$e_t$"))
CURVES = {"5-10": ("#4292c6", "o", "--", r"$n$: 5/10"), "10-20": ("#08306b", "s", "-", r"$n$: 10/20"),
          "eie": ("0.35", "^", ":", "interface model")}
for c in CASES:
    eie = json.loads((DATA / f"{c}.json").read_text())
    for k, (stem, title) in enumerate(PANELS):
        fig, ax = plt.subplots(figsize=(7, 6))
        plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.12)
        for key, (col, m, ls, _) in CURVES.items():
            y = [eie[f"eie/10/{h}"][k] if key == "eie" else rc[f"{c}/{h}/{key}"][k] for h in HS]
            ax.loglog(x, y, linestyle=ls, linewidth=2, marker=m, markersize=7, color=col)
        ax.set_xlabel(r"$h$ [mm]"); ax.set_title(title, fontsize=32, pad=15); ax.grid(True, alpha=0.3)
        ax.xaxis.set_major_locator(LogLocator(base=10, numticks=6)); ax.yaxis.set_major_locator(LogLocator(base=10, numticks=6))
        for ext in ("pdf", "png"):
            fig.savefig(OUT / f"{c}_{stem}.{ext}", dpi=150)
        plt.close(fig)
handles = [Line2D([], [], color=col, marker=m, ls=ls, lw=2, ms=9, label=lab) for col, m, ls, lab in CURVES.values()]
fig = plt.figure(figsize=(1.7, 1.1))
fig.legend(handles=handles, loc="center", ncol=1, frameon=False, fontsize=13, handlelength=2.0)
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"legend_dp_ref.{ext}", dpi=150)
# the rows of the convergence table: largest 10/20 change over h, its ratio to the 5/10 change, smallest factor below the EIE
for c in CASES:
    eie = json.loads((DATA / f"{c}.json").read_text())
    b = [[rc[f"{c}/{h}/10-20"][k] for h in HS] for k in range(3)]
    a = [[rc[f"{c}/{h}/5-10"][k] for h in HS] for k in range(3)]
    e = [[eie[f"eie/10/{h}"][k] for h in HS] for k in range(3)]
    ratio = [a[k][i] / b[k][i] for k in range(3) for i in range(6)]
    factor = [e[k][i] / b[k][i] for k in range(3) for i in range(6)]
    print(f"{c:13s} max 10-20: " + " ".join(f"{max(b[k]):.1e}" for k in range(3))
          + f"  ratio {min(ratio):.1f}-{max(ratio):.1f}  factor min {min(factor):.0f}"
          + "  per measure " + " ".join(f"{min(e[k][i] / b[k][i] for i in range(6)):.0f}" for k in range(3)))
print("wrote", OUT)
