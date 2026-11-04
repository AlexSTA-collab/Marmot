#!/usr/bin/env python3
"""Figure of the appendix on the in-plane refinement of the PSIE (reviewer comment 2.9): normal traction t_n along the middle
row of interface cells at the end of the loading, stiff family, alpha = 10, h = 0.01, for the in-plane meshes N x N, with the
PSIE (left) and with the squeeze-flow term (middle), and the L2 difference of each profile from the finest mesh (right).
Top row H = 100 MPa, bottom row H = 0.1 MPa.
Data: psie_inplane_data/{psie,psie_squeeze}.json, written by psie_inplane.py (ROOT=psie_inplane / psie_squeeze).
usage: plot_psie_inplane.py <output dir>"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import LogLocator  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"],
                     "figure.dpi": 100, "savefig.bbox": "tight", "font.size": 17})
D = {k: json.loads((HERE / "psie_inplane_data" / f"{k}.json").read_text()) for k in ("psie", "psie_squeeze")}
# meshes are ordered: one hue, light to dark
SHADE = {15: "#9ecae1", 30: "#6baed6", 60: "#2171b5", 120: "#08306b"}
LW = {15: 1.4, 30: 1.6, 60: 1.8, 120: 2.0}
ELEM = {"psie": ("PSIE", "#d62728", "o"), "psie_squeeze": ("PSIE with squeeze-flow term", "#1f77b4", "s")}

fig, axes = plt.subplots(2, 3, figsize=(15, 8.4), gridspec_kw={"width_ratios": [1, 1, 0.8]})
for i, H in enumerate(("100", "0.1")):
    for j, k in enumerate(("psie", "psie_squeeze")):
        ax = axes[i, j]
        for N in sorted(D[k][H], key=int):
            r = D[k][H][N]
            ax.plot(r["s"], r["t_n"], color=SHADE[int(N)], lw=LW[int(N)])
        ax.set_title(rf"{ELEM[k][0]}, $H={H}$~MPa", fontsize=19, pad=8)
        ax.grid(True, alpha=0.3)
        ax.set_ylabel(r"$t_n$ [MPa]")
    lo = min(axes[i, 0].get_ylim()[0], axes[i, 1].get_ylim()[0]); hi = max(axes[i, 0].get_ylim()[1], axes[i, 1].get_ylim()[1])
    for j in (0, 1):
        axes[i, j].set_ylim(lo, hi)
    ax = axes[i, 2]
    for k, (lab, col, mk) in ELEM.items():
        Ns = [int(N) for N in sorted(D[k][H], key=int) if int(N) != 120]
        ax.loglog(Ns, [D[k][H][str(N)]["dtn"] for N in Ns], color=col, marker=mk, ms=8, lw=2, label=lab)
    ax.set_xticks([15, 30, 60]); ax.set_xticklabels(["15", "30", "60"]); ax.minorticks_off()
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=6))
    ax.grid(True, alpha=0.3)
    ax.set_title(r"$\|t_n^N-t_n^{120}\|/\|t_n^{120}\|$" + rf" ($H={H}$~MPa)", fontsize=19, pad=8)
for j in (0, 1):
    axes[1, j].set_xlabel(r"$s$ [mm]")
axes[1, 2].set_xlabel(r"$N$")
handles = [Line2D([], [], color=SHADE[N], lw=LW[N], label=rf"${N}\times{N}$") for N in SHADE]
handles += [Line2D([], [], color=c, marker=m, ms=8, lw=2, label=l) for l, c, m in ELEM.values()]
fig.tight_layout(rect=(0, 0.07, 1, 1))
fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False, fontsize=16)
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"psie_inplane.{ext}", dpi=150)
print("wrote", OUT / "psie_inplane.pdf")
