#!/usr/bin/env python3
"""Global response of the J2 interphase at H = 100 MPa against the 20-layer C3D8B resolved layer (stiff, alpha = 10,
h = 0.01): reaction-displacement curve and departure of EIE and PSIE, in the style of paper/figures/reaction/plot_reaction.py
(von Mises block, same sizes and panels). The resolved layer takes more increments than the interface models, so it is
interpolated in the applied displacement at the increments of each interface model. The dotted grey line continues the
elastic tangent of the resolved layer, so that the departure of the hardening branch from it shows at the yield kink. An
interface model that stops before the end of the loading (the EIE at H = 0.1) is drawn up to its last converged increment,
marked by a cross.
usage: plot_reaction_H100.py <reaction_H100.json> <out dir> [legend location of the force panel, default upper left]"""
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

C = {k: (np.array(v[0]), np.array(v[1])) for k, v in json.loads(Path(sys.argv[1]).read_text()).items()}
OUT = Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)
LOC = sys.argv[3] if len(sys.argv) > 3 else "upper left"
plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"],
                     "savefig.format": "pdf", "savefig.bbox": "tight", "figure.dpi": 100})
VM_RC = {"font.size": 16, "axes.labelsize": 28, "axes.titlesize": 32,
         "xtick.labelsize": 28, "ytick.labelsize": 28, "legend.fontsize": 18}
STY = {"full": ("#111111", "-", "resolved layer"), "EIE": ("#1f77b4", "--", "EIE"), "PSIE": ("#d62728", ":", "PSIE")}
sgn = np.sign(C["full"][1][-1])
u0, F0 = np.max(np.abs(C["full"][0])), np.max(np.abs(C["full"][1]))
with plt.rc_context(VM_RC):
    fig, ax = plt.subplots(figsize=(7, 6))
    plt.subplots_adjust(left=0.19, right=0.96, top=0.94, bottom=0.16)
    u, F = C["full"]; s = np.diff(sgn * F) / np.where(np.diff(u) > 0, np.diff(u), np.nan)
    s0 = s[np.isfinite(s)][0]                                   # elastic tangent of the resolved layer
    k = int(np.argmax(np.isfinite(s) & (s < 0.9 * s0)))       # first increment off it: the yield kink
    for tag, (col, ls, lab) in STY.items():
        if tag not in C:
            continue
        u, F = C[tag]
        ax.plot(u / u0, sgn * F / F0, ls, color=col, lw=2.5, label=lab)
        if u[-1] < 0.999 * u0:   # stopped before the end of the loading
            ax.plot(u[-1] / u0, sgn * F[-1] / F0, "x", color=col, ms=7, mew=1.6, zorder=5)
    ylim = ax.get_ylim()
    ue = min(1.0, ylim[1] / (s0 * u0 / F0))
    ax.plot([0, ue], [0, ue * s0 * u0 / F0], ":", color="#7f7f7f", lw=1.8, label="elastic tangent", zorder=1)
    ax.set_ylim(ylim)
    ax.set_xlabel(r"$\bar{u}_x\,/\,\bar{u}^{\max}_x$"); ax.set_title(r"$\frac{F}{F^0}$", pad=15)
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False, loc=LOC, handlelength=1.2, handletextpad=0.4, borderaxespad=0.2, labelspacing=0.3)
    # where the yield kink lands on the saved (tight) page, for the TikZ magnifier of Fig. vm global
    u, F = C["full"]
    fig.canvas.draw(); r = fig.canvas.get_renderer(); bb = fig.get_tightbbox(r)
    px, py = ax.transData.transform((u[k] / u0, sgn * F[k] / F0)) / fig.dpi
    print(f"yield kink at u/umax = {u[k] / u0:.3f}, F/F0 = {sgn * F[k] / F0:.3f}: page fraction "
          f"({(px - bb.x0) / bb.width:.3f}, {(py - bb.y0) / bb.height:.3f})")
    ua = u[k] / u0 + 0.03     # just past the kink, where the branch leaves the tangent
    px, py = ax.transData.transform((ua, np.interp(ua * u0, u, sgn * F) / F0)) / fig.dpi
    print(f"magnifier centre at u/umax = {ua:.3f}: page fraction ({(px - bb.x0) / bb.width:.3f}, {(py - bb.y0) / bb.height:.3f})")
    fig.savefig(OUT / "vonmises_force.pdf"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 6))
    plt.subplots_adjust(left=0.24, right=0.96, top=0.94, bottom=0.16)
    for tag in [t for t in ("EIE", "PSIE") if t in C]:
        u, F = C[tag]
        d = np.abs(F - np.interp(u, *C["full"])) / F0
        ax.semilogy(u[1:] / u0, np.maximum(d[1:], 1e-12), "-", color=STY[tag][0], lw=2.5, label=STY[tag][2])
        if u[-1] < 0.999 * u0:
            ax.semilogy(u[-1] / u0, max(d[-1], 1e-12), "x", color=STY[tag][0], ms=12, mew=2.5)
        print(f"{tag}: max |F - F0| / max|F0| = {d.max():.2e}, at the end {d[-1]:.2e}")
    ax.set_xlabel(r"$\bar{u}_x\,/\,\bar{u}^{\max}_x$"); ax.set_title(r"$\frac{|F-F^{0}|}{\max|F^{0}|}$", pad=15)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, loc="lower right", handlelength=1.2, handletextpad=0.4, borderaxespad=0.2, labelspacing=0.3)
    fig.savefig(OUT / "vonmises_departure.pdf"); plt.close(fig)
