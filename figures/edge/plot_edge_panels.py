#!/usr/bin/env python3
"""Figure L.48 as separate panels: e_u, e_[u], e_t against the width of the strip removed at every edge of the
interface, 60 x 60 meshes, alpha = 10, h = 0.1 (solid, filled) and 0.01 mm (dashed, open), up to the margin of the
paper, for the elastic, viscoelastic, biphasic elastic and biphasic viscoelastic interphases. Same data and curves as the
elastic figure of plot_sweep60.py (sweep60/<analysis>.json, entries <series>|<h>|<strip> = [e_u, e_[u], e_t]); one PDF
per analysis and measure in the panel style of the paper, and one legend in a column, so that the figure can be set as
four captioned rows.
usage: python3 plot_edge_panels.py <out dir> [elastic | plastic]  ->  <out>/<analysis>_{global_relative_error,displacement_jump,traction}.pdf,
                                                 <out>/legend_edge.pdf"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
XMAX = 0.15
plt.rcParams.update({
    "text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"], "font.size": 16,
    "axes.labelsize": 28, "axes.titlesize": 32, "xtick.labelsize": 24, "ytick.labelsize": 28,
    "savefig.format": "pdf", "savefig.bbox": "tight", "figure.dpi": 100,
})
FAM = {"stiff": ("#1f77b4", "o", "stiff"), "initial": ("#ff7f0e", "s", "initial"), "soft": ("#d62728", "^", "soft"),
       "full": ("#2ca02c", "o", r"full model"), "supp": ("#9467bd", "D", r"$[C]=0$")}
PANELS = (("global_relative_error", r"$e_u$"), ("displacement_jump", r"$e_{[\![u]\!]}$"), ("traction", r"$e_t$"))
SETS = {"elastic": ("elastic", "visco", "biphasic", "biphasic_visco"),   # Figure L.48
        "plastic": ("eie", "psie", "dp")}                             # Figure L.49: J2 EIE, J2 PSIE, Drucker-Prager EIE
ANALYSES = SETS[sys.argv[2] if len(sys.argv) > 2 else "elastic"]
used = []
for a in ANALYSES:
    d = json.loads((HERE / "sweep60" / f"{a}.json").read_text())
    series = sorted({k.split("|")[0] for k in d}, key=list(FAM).index)
    used += [s for s in series if s not in used]
    for i, (stem, title) in enumerate(PANELS):
        fig, ax = plt.subplots(figsize=(7, 6))
        plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.12)
        for s in series:
            col, mk, _ = FAM[s]
            for h, ls in (("0.1", "-"), ("0.01", "--")):
                pts = sorted((float(k.split("|")[2]), v) for k, v in d.items() if k.startswith(f"{s}|{h}|"))
                pts = [q for q in pts if q[0] <= XMAX + 1e-9]
                if pts:
                    ax.semilogy([q[0] for q in pts], [q[1][i] for q in pts], ls, color=col, marker=mk, ms=7, lw=2,
                                mfc=col if ls == "-" else "none")
        ax.axvline(0.15, color="k", lw=1, ls=":")
        ax.grid(True, alpha=0.3)
        ax.set_xlim(-0.005, XMAX + 0.005)
        ax.set_xlabel(r"strip [mm]")
        ax.set_title(title, fontsize=32, pad=15)
        fig.savefig(OUT / f"{a}_{stem}.pdf")
        plt.close(fig)
hand = [Line2D([], [], color=FAM[s][0], marker=FAM[s][1], lw=2, ms=8, label=FAM[s][2]) for s in used]
hand += [Line2D([], [], color="k", ls="-", marker="o", lw=2, ms=7, label=r"$h=0.1$~mm"),
         Line2D([], [], color="k", ls="--", marker="o", mfc="none", lw=2, ms=7, label=r"$h=0.01$~mm"),
         Line2D([], [], color="k", ls=":", lw=1, label=r"margin")]
fig = plt.figure(figsize=(1.9, 2.6))                    # one column, beside the panels, as in Figs. H.28 and K.47
fig.legend(handles=hand, loc="center", ncol=1, frameon=False, fontsize=13, handlelength=2.0)
fig.savefig(OUT / "legend_edge.pdf")
print("wrote", OUT, "| series:", used)
