#!/usr/bin/env python3
"""Figures of Section 4.3.3 (biphasic viscoelastic layer), paper style (cf. cyclic/plot_cycle_paper.py).

  histories_strain.pdf   strain of the two halves of the resolved layer
  histories_ratio.pdf    their ratio
  histories_jump.pdf     difference of the shear strains of the bottom (-) and top (+) half: resolved layer, full interface
                         model and the model with [C^ev] = 0 (--supp-hist, the same history of the suppressed model)
  response_creep.pdf     normalised displacement of the loaded face: resolved layer, full model, [C^ev] = 0
  response_error.pdf     displacement error over the hold

The [C^ev] = 0 model is the one with the mean parameters of the two halves; the fitted single material is not drawn.
usage: paper_figures.py <ref hist.json> <model hist.json> <resp h=0.1 json> [<resp h=0.01 json>] <out dir>
                        [--supp-hist=<suppressed model hist.json>]
With SIDE_LEGEND=1 (Figure 14 set as Figure 12) the panels of the figure carry no legend, and <out dir>/legend_bv.pdf
holds one legend in a column for the three panels."""
import json
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SUPP = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--supp-hist=")), None)
args = [a for a in sys.argv[1:] if not a.startswith("--supp-hist=")]
refh, modh, r1 = (json.loads(Path(a).read_text()) for a in args[:3])
r01 = json.loads(Path(args[3]).read_text()) if len(args) == 5 else None
OUT = Path(args[-1]); OUT.mkdir(parents=True, exist_ok=True)
LAB, TIT = float(os.environ.get("LAB", 26)), float(os.environ.get("TIT", 29))
plt.rcParams.update({
    "text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"],
    # placed at \convfigwidth (0.29 textwidth) these sizes match Figure 13a on the page: 8.0 pt labels, 8.9 pt titles
    "font.size": 16, "axes.labelsize": LAB, "axes.titlesize": TIT, "xtick.labelsize": LAB, "ytick.labelsize": LAB,
    "legend.fontsize": 18, "savefig.format": "pdf", "savefig.bbox": "tight", "figure.dpi": 100,
})
BLACK, RED, BLUE, GREEN = "#111111", "#d62728", "#1f77b4", "#2ca02c"
SIDE = bool(os.environ.get("SIDE_LEGEND"))


def save(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf"); fig.savefig(OUT / f"{stem}.png", dpi=110); plt.close(fig)


def panel():
    fig, ax = plt.subplots(figsize=(7, 6))
    plt.subplots_adjust(left=0.24, right=0.96, top=0.88, bottom=0.17)
    ax.grid(True, which="both", alpha=0.3)
    return fig, ax


tb, gb = np.array(refh["layerbottom"]["t"]), np.array(refh["layerbottom"]["gamma_t"])
tt, gt = np.array(refh["layertop"]["t"]), np.array(refh["layertop"]["gamma_t"])
k = tb >= 1.0
# (a) strain of the two halves
fig, ax = panel()
ax.loglog(tb[k], gb[k], "-", color=BLACK, lw=2.5, label="bottom half")
ax.loglog(tt[k], gt[k], "--", color=BLACK, lw=2.5, label="top half")
ax.set_xlabel(r"$t$~[days]"); ax.set_title(r"$\gamma_{tn}$", pad=15)
ax.legend(frameon=False, loc="center right", handlelength=1.8, handletextpad=0.4, borderaxespad=0.2, labelspacing=0.3)
save(fig, "histories_strain")
# (b) ratio
fig, ax = panel()
ax.semilogx(tb[k], gb[k] / gt[k], "-", color=BLACK, lw=2.5)
ax.set_xlabel(r"$t$~[days]"); ax.set_title(r"$\gamma^{-}_{tn}\,/\,\gamma^{+}_{tn}$", pad=15)
save(fig, "histories_ratio")
# (c) difference of the face strains: resolved layer against the interface model (accumulated jump of the normal gradient)
fig, ax = panel()
ax.loglog(tb[k], (gb - gt)[k], "-", color=BLACK, lw=2.5, label="resolved layer")
tm = np.array(modh["t"]); zc = np.abs(np.cumsum(np.array(modh["z_t"])))
km = tm >= 1.0
ax.loglog(tm[km], zc[km], ":", color=RED, lw=2.5, label="full model")
if SUPP:
    sh = json.loads(Path(SUPP).read_text())
    ts = np.array(sh["t"]); zs = np.abs(np.cumsum(np.array(sh["z_t"]))); ks = ts >= 1.0
    ax.loglog(ts[ks], zs[ks], "-.", color=BLUE, lw=2.5, label=r"$[C^{ev}]=0$")
ax.set_xlabel(r"$t$~[days]"); ax.set_title(r"$\big|T_{ik}\,(u^{-}_{k,j}-u^{+}_{k,j})\,N_{jl}\big|$", pad=15)
if not SIDE:
    ax.legend(frameon=False, loc="center right", handlelength=1.8, handletextpad=0.4, borderaxespad=0.2, labelspacing=0.3)
save(fig, "histories_jump")

# (d) creep of the loaded face
fig, ax = panel()
t = np.array(r1["ref"]["t"]); kk = t >= 1.0
u0 = lambda key: np.array(r1[key]["u_face"])[0]
for key, sty, col, lab in (("ref", "-", BLACK, "resolved layer"), ("full", ":", RED, r"full model"),
                           ("supp", "-.", BLUE, r"$[C^{ev}]=0$")):
    u = np.array(r1[key]["u_face"]); ax.semilogx(t[kk], u[kk] / u0("ref"), sty, color=col, lw=2.5, label=lab)
ax.set_xlabel(r"$t$~[days]"); ax.set_title(r"$u^{\mathrm{L}}_x\,/\,u^{\mathrm{L}}_x(0^+)$", pad=15)
if not SIDE:
    ax.legend(frameon=False, loc="upper left", handlelength=1.8, handletextpad=0.4, borderaxespad=0.2, labelspacing=0.3, fontsize=15)
save(fig, "response_creep")
# (e) displacement error over the hold, h = 0.1 (solid, circles) and h = 0.01 (dashed, squares)
from matplotlib.lines import Line2D
fig, ax = panel()
for key, col, lab in (("full", RED, "full model"), ("supp", BLUE, r"$[C^{ev}]=0$")):
    for r, ls, mk in ((r1, "-", "o"), (r01, "--", "s")):
        if r is None:
            continue
        tt_ = np.array(r[key]["t"]); k_ = tt_ >= 1.0
        ax.semilogy(tt_[k_], np.array(r[key]["e_u"])[k_], ls, color=col, lw=2.5, marker=mk, markersize=5)
ax.set_xlabel(r"$t$~[days]"); ax.set_title(r"$e_u$", pad=15)
handles = [Line2D([0], [0], color=RED, lw=2.5), Line2D([0], [0], color=BLUE, lw=2.5),
           Line2D([0], [0], color=BLACK, lw=2.5, marker="o", markersize=5), Line2D([0], [0], color=BLACK, lw=2.5, ls="--", marker="s", markersize=5)]
labels = ["full model", r"$[C^{ev}]=0$", r"$h=0.1$", r"$h=0.01$"]
if not SIDE:
    ax.legend(handles, labels, frameon=False, loc="upper center", ncol=2, handlelength=1.6, handletextpad=0.4, borderaxespad=0.2, labelspacing=0.2, columnspacing=0.8, fontsize=13)
ax.set_ylim(5e-5, 5e1)
save(fig, "response_error")
if SIDE:
    hs = [Line2D([0], [0], color=BLACK, lw=2.5, label="resolved layer"), Line2D([0], [0], color=RED, ls=":", lw=2.5, label="full model"),
          Line2D([0], [0], color=BLUE, ls="-.", lw=2.5, label=r"$[C^{ev}]=0$"),
          Line2D([0], [0], color="0.3", lw=2.5, marker="o", markersize=5, label=r"$h=0.1$~mm"),
          Line2D([0], [0], color="0.3", lw=2.5, ls="--", marker="s", markersize=5, label=r"$h=0.01$~mm")]
    fig = plt.figure(figsize=(1.9, 1.9))
    fig.legend(handles=hs, loc="center", ncol=1, frameon=False, fontsize=13, handlelength=2.0)
    fig.savefig(OUT / "legend_bv.pdf"); plt.close(fig)
print("wrote", OUT)
