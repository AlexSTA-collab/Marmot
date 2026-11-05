#!/usr/bin/env python3
"""The legends of Figures I.39 and I.40 (free edges), one per figure, in the style of plot_edge_panels.py:
  legend_edge_families.pdf  stiff, initial, soft (elastic and viscoelastic layers, Figure I.39)
  legend_edge_biphasic.pdf  mismatch kept / removed (biphasic layers, Figure I.40)
usage: legends_edge.py <out dir>"""
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"],
                     "savefig.format": "pdf", "savefig.bbox": "tight", "figure.dpi": 100})
SERIES = {"families": [("#1f77b4", "o", "stiff"), ("#ff7f0e", "s", "initial"), ("#d62728", "^", "soft")],
          "biphasic": [("#2ca02c", "o", r"$[C^{e(\cdot)}_{ijkl}]$ kept"), ("#9467bd", "D", r"$[C^{e(\cdot)}_{ijkl}]$ removed")]}
for name, ser in SERIES.items():
    hand = [Line2D([], [], color=c, marker=m, lw=2, ms=8, label=lab) for c, m, lab in ser]
    hand += [Line2D([], [], color="k", ls="-", marker="o", lw=2, ms=7, label=r"$h=0.1$~mm"),
             Line2D([], [], color="k", ls="--", marker="o", mfc="none", lw=2, ms=7, label=r"$h=0.01$~mm"),
             Line2D([], [], color="k", ls=":", lw=1, label=r"margin")]
    fig = plt.figure(figsize=(1.9, 2.6))
    fig.legend(handles=hand, loc="center", ncol=1, frameon=False, fontsize=13, handlelength=2.0)
    fig.savefig(OUT / f"legend_edge_{name}.pdf")
    plt.close(fig)
    print("wrote", OUT / f"legend_edge_{name}.pdf")
