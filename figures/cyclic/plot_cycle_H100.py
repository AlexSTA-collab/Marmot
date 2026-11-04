#!/usr/bin/env python3
"""H = 100 MPa, 20-layer C3D8B reference: copy of cyclic/plot_cycle_paper.py with the inputs changed.

Panels of the cyclic test of Section 4.4.7 in the typesetting of the other panels of the paper (plot_reaction.py).

    python3 plot_cycle_H100.py <figures dir> [data dir, default cyc100; cyc01 for the H = 0.1 MPa runs]

Writes <figures dir>/{loop_h0.1,loop_h0.01,program_error}.pdf and the PNG of each.
With SIDE_LEGEND=1 (Figure I.34) the panels carry no legend and no thickness label, and <figures dir>/legend_cyclic.pdf
holds one legend in a column for the figure, as in Figures H.28 and I.32.
Reads loop.npz of the four runs (diag/loop2.py) and the error histories of analyze.py.
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = Path(sys.argv[1] if len(sys.argv) > 1 else HERE / "figures")
OUT.mkdir(parents=True, exist_ok=True)
C = HERE / (sys.argv[2] if len(sys.argv) > 2 else "cyc100")   # loop.npz (cyclic/diag/loop2.py), err_<h>.json (cyclic/analyze.py)
REF = {h: C / "ref" / f"h_{h}" for h in ("0.1", "0.01")}
PSIE = {h: C / "psie" / f"h_{h}" for h in ("0.1", "0.01")}
ERR = {h: C / f"err_{h}.json" for h in ("0.1", "0.01")}
U0 = 1.0  # the program is normalised by its own amplitude
SIDE = bool(os.environ.get("SIDE_LEGEND"))

plt.rcParams.update({
    "text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"],
    "font.size": 16, "axes.labelsize": 28, "axes.titlesize": 32, "xtick.labelsize": 28, "ytick.labelsize": 28,
    "legend.fontsize": 18, "savefig.format": "pdf", "savefig.bbox": "tight", "figure.dpi": 100,
})
BLACK, RED = "#111111", "#d62728"


def load(run):
    z = np.load(Path(run) / "loop.npz")
    return z["t"], z["u"], -z["F"]  # reaction positive in the loading direction


def save(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf")
    fig.savefig(OUT / f"{stem}.png", dpi=110)
    plt.close(fig)


for h in ("0.1", "0.01"):
    tr, ur, Fr = load(REF[h])
    tp, up, Fp = load(PSIE[h])
    u0, F0 = np.max(np.abs(ur)), np.max(np.abs(Fr))
    fig, ax = plt.subplots(figsize=(7, 6))
    plt.subplots_adjust(left=0.22, right=0.96, top=0.88, bottom=0.17)
    ax.plot(ur / u0, Fr / F0, "-", color=BLACK, lw=2.5, label="resolved layer")
    ax.plot(up / u0, Fp / F0, ":", color=RED, lw=2.5, label="PSIE")
    print(f"h={h}: F/F0 at the reversals (ref) " + " ".join(f"{np.interp(s, tr, Fr) / F0:+.3f}" for s in (1, 2, 3, 4, 5))
          + f"; max |F_psie - F_ref| / F0 = {np.max(np.abs(Fp - np.interp(tp, tr, Fr))) / F0:.2e}")
    ax.set_xlabel(r"$\bar{u}_x\,/\,\bar{u}^{\max}_x$")
    ax.set_title(r"$\frac{F}{F^0}$", pad=15)
    ax.set_xticks([-1, 0, 1])
    ax.set_yticks([-1, 0, 1])
    # inside the loop when the layer is close to perfect plasticity and the loop is a parallelogram
    if not SIDE:
        if C.name == "cyc100":   # bottom right corner, inside the frame
            ax.text(0.97, 0.04, rf"$h={h}$~mm", transform=ax.transAxes, ha="right", va="bottom", fontsize=24)
        else:
            ax.text(0.40, 0.25, rf"$h={h}$~mm", transform=ax.transAxes, ha="center", va="center", fontsize=24)
    ax.grid(True, alpha=0.3)
    if not SIDE:
        ax.legend(frameon=False, loc="upper left", handlelength=1.2, handletextpad=0.4,
                  borderaxespad=0.2, labelspacing=0.3)
    save(fig, f"loop_h{h}")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 6), sharex=True, gridspec_kw={"height_ratios": [1, 1.25]})
plt.subplots_adjust(left=0.30, right=0.96, top=0.96, bottom=0.17, hspace=0.12)
t = np.linspace(0, 5, 501)
ax1.plot(t, np.interp(t, [0, 1, 2, 3, 4, 5], [0, 1, 0, -1, 0, 1]), "-", color=BLACK, lw=2.5)
ax1.set_ylabel(r"$\bar{u}_x\,/\,\bar{u}^{\max}_x$", fontsize=26)
ax1.set_yticks([-1, 0, 1])
ax1.grid(True, alpha=0.3)
for h, ls, mk in (("0.1", "-", "o"), ("0.01", "--", "s")):
    r = json.loads(ERR[h].read_text())
    tt, e = np.array(r["t"]), np.array(r["e"])
    print(f"h={h}: e_u " + " ".join(f"t={a:.2f}:{b:.2e}" for a, b in zip(tt[1:], e[1:])))
    ax2.semilogy(tt[1:], e[1:], ls, color=RED, lw=2.5, marker=mk, markersize=6, label=rf"$h={h}$")
ax2.set_xlabel(r"$t$")
ax2.set_xticks(range(6))
ax2.set_ylabel(r"$e_u$")
ax2.set_yticks([1e-6, 1e-4])
ax2.grid(True, which="both", alpha=0.3)
if not SIDE:
    ax2.legend(frameon=False, loc="lower right", handlelength=1.8, handletextpad=0.4, borderaxespad=0.2, labelspacing=0.3)
save(fig, "program_error")
if SIDE:
    from matplotlib.lines import Line2D
    hs = [Line2D([], [], color=BLACK, ls="-", lw=2.5, label="resolved layer"),
          Line2D([], [], color=RED, ls=":", lw=2.5, label="PSIE"),
          Line2D([], [], color=RED, ls="-", lw=2.5, marker="o", ms=6, label=r"$e_u$, $h=0.1$~mm"),
          Line2D([], [], color=RED, ls="--", lw=2.5, marker="s", ms=6, label=r"$e_u$, $h=0.01$~mm")]
    fig = plt.figure(figsize=(2.0, 1.4))
    fig.legend(handles=hs, loc="center", ncol=1, frameon=False, fontsize=13, handlelength=2.0)
    fig.savefig(OUT / "legend_cyclic.pdf"); plt.close(fig)
print("wrote", OUT)
