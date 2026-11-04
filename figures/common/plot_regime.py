#!/usr/bin/env python3
"""Regime migration figures, in the three-panel layout of the convergence figures.

Two figures, both one panel per material family, so that placement at 0.29\\textwidth gives
them the same height as every other results figure:

  regime_chi/{stiff,initial,soft}.pdf   the jump measured against the membrane deformation
                                        of the layer, ||[u_i]|| / (h ||<u^S_{i,j}>||),
                                        against the step length, one curve per interphase
                                        thickness. Both are L2 over the interface; <u^S_{i,j}>
                                        is the average over the two faces of the symmetric
                                        surface gradient, so the ratio compares the cohesive
                                        contribution with what the layer would carry over its
                                        own thickness if it followed the adherents.

                                        with the two adherent moduli marked. A material
                                        property: no thickness enters it.

plus a shared legend beside each row.

Style follows master_parametric_script.py: the quantity is the panel TITLE, not a y-label;
grid alpha 0.3; no per-panel legend; figsize (7,6) with fixed margins.

    python3 figures/regime/plot_regime.py figures/regime
"""
import pickle
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors                                    # noqa: E402
import numpy as np                                                 # noqa: E402
from matplotlib.lines import Line2D                                # noqa: E402
from matplotlib.ticker import LogLocator, NullFormatter            # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "verification"))
from kelvin_chain import PowerLawChain                             # noqa: E402

# rcParams verbatim from master_parametric_script.py, so these panels typeset
# identically to the convergence figures they sit beside. At 0.29\textwidth the 28 pt
# labels land at 7.3 pt on the page, against an 8 pt caption.
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 16,
    "axes.labelsize": 28,
    "axes.titlesize": 32,
    "xtick.labelsize": 28,
    "ytick.labelsize": 28,
    "legend.fontsize": 18,
    "savefig.format": "pdf",
    "savefig.bbox": "tight",
    "figure.dpi": 100,
})

HERE = Path(__file__).resolve().parent
OUT = Path(sys.argv[1] if len(sys.argv) > 1 else HERE)

FAMILIES = ["stiff", "initial", "soft"]
HEIGHTS = ["0.01", "0.02", "0.04", "0.06", "0.08", "0.1"]
E0 = {"stiff": 4e5, "initial": 4e3, "soft": 4e-5}
E_TOP, E_BOT = 2e5, 4e3
M, N = 5e-3, 0.3

# One panel per family, so the family is what the reader identifies the panel by and
# the panel should carry the family's colour, the same blue, orange and red the
# convergence figures use. Thickness is sequential within a panel, so it runs as a
# gradient of that one colour rather than across a second, unrelated scale: a light
# tint at the thinnest layer to the family colour itself at the coarsest.
FAMCOL = {"stiff": "#1f77b4", "initial": "#ff7f0e", "soft": "#d62728"}


def _ramp(hexcol, n):
    """n shades of one colour, from a light tint to a shade slightly darker than it."""
    base = np.array(colors.to_rgb(hexcol))
    out = []
    for i in range(n):
        t = i / (n - 1)
        if t <= 0.75:                       # white -> the colour itself
            f = 0.78 - 0.78 * (t / 0.75)
            out.append(tuple(base + f * (1.0 - base)))
        else:                               # the colour -> a darker shade of it
            f = (t - 0.75) / 0.25 * 0.30
            out.append(tuple(base * (1.0 - f)))
    return out


HCOL = {fam: dict(zip(HEIGHTS, _ramp(c, len(HEIGHTS)))) for fam, c in FAMCOL.items()}

C = pickle.loads((HERE / "regime_curves.pkl").read_bytes())


def panel(fam, kind, outdir):
    fig, ax = plt.subplots(figsize=(7, 6))
    plt.subplots_adjust(left=0.17, right=0.95, top=0.95, bottom=0.13)

    if kind == "chi":
        ax.set_xscale("log")
        for h in HEIGHTS:
            rec = C.get(("dtfull", fam, h))
            if not rec:
                continue
            dts = sorted(rec)
            ax.plot(dts, [rec[d]["interface"] for d in dts], "-o",
                        color=HCOL[fam][h], lw=2.0, ms=7, mfc="none", mew=2.0)
        ax.set_yscale("log")
        # snap the limits to whole decades: without this the soft family spans less than
        # one decade, its ticks come out as 6x10^8 and the panel is squeezed narrower
        # than its neighbours by the width of those labels
        vals = [v["interface"] for h in HEIGHTS if (r := C.get(("dtfull", fam, h)))
                for v in r.values()]
        lo, hi = np.log10(min(vals)), np.log10(max(vals))
        ax.set_ylim(10.0 ** np.floor(lo), 10.0 ** np.ceil(hi))
        ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=5))
        ax.yaxis.set_minor_formatter(NullFormatter())   # decades only, no 6x10^8 clutter
        ax.set_title(r"$\frac{\|[u_i]\|_{L^2}}{h\,\|\langle u^S_{i,j}\rangle\|_{L^2}}$",
                     pad=15)
    else:
        chain = PowerLawChain(E0[fam], 0.3, M, N)
        grid = np.logspace(-8, 4, 300)
        ax.loglog(grid, [chain.E_ev(d) for d in grid], "-", color="#1f77b4", lw=2.5)
        for E, lab in ((E_TOP, r"$E^{+}$"), (E_BOT, r"$E^{-}$")):
            ax.axhline(E, color="0.55", ls="--", lw=1.5)
            ax.text(2e-8, E * 1.7, lab, fontsize=24, color="0.35")
        ax.set_ylim(1e-6, 1e7)
        ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
        ax.set_title(r"$E^{ev}(\Delta t)$ [MPa]", pad=15)

    ax.set_xlabel(r"$\Delta t$ [days]")
    ax.set_xlim(1e-9, 1e5)
    ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
    ax.grid(True, which="both", alpha=0.3)
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{fam}.pdf")
    plt.close(fig)


def legend(kind, outdir):
    fig = plt.figure(figsize=(2.2, 6))
    if kind == "chi":
        # One legend serves three panels that now each carry their own family colour,
        # so it can only show the ordering: a neutral ramp, light to dark with
        # increasing thickness, which is how each panel is shaded in its own colour.
        grey = _ramp("#3a3a3a", len(HEIGHTS))
        handles = [Line2D([], [], color=grey[i], marker="o", mfc="none", mew=2.0,
                          lw=2.0, ms=7, label=rf"$h={h}$") for i, h in enumerate(HEIGHTS)]
    else:
        handles = [Line2D([], [], color="#1f77b4", lw=2.5, label=r"$E^{ev}$"),
                   Line2D([], [], color="0.55", ls="--", lw=1.5, label="adherents")]
    fig.legend(handles=handles, loc="center", frameon=False)
    fig.savefig(outdir / "common_legend.pdf")
    plt.close(fig)


for kind, sub in (("chi", "regime_chi"),):
    d = OUT / sub
    for fam in FAMILIES:
        panel(fam, kind, d)
    legend(kind, d)
    print(f"wrote {d}/{{{','.join(FAMILIES)}}}.pdf + common_legend.pdf")

print("\nindicator over the sweep (interface model), and worst deviation from the resolved layer:")
for fam in FAMILIES:
    for h in HEIGHTS:
        rec = C.get(("dtfull", fam, h))
        if not rec:
            continue
        dts = sorted(rec)
        dev = [abs(rec[d]["interface"] - rec[d]["full"]) / abs(rec[d]["full"])
               for d in dts if "full" in rec[d] and rec[d]["full"]]
        print(f"  {fam:>8} h={h:<5} {rec[dts[0]]['interface']:.4f} -> "
              f"{rec[dts[-1]]['interface']:.4f}   dev {max(dev)*100 if dev else float('nan'):.2f} %")
