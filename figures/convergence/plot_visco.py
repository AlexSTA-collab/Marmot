#!/usr/bin/env python3
"""Convergence panels for the viscoelastic interphase, EIQUAD4.

Source: the matched-thickness sweep in the PowerLaw study's
proper_convergence/ (alpha = 10) and proper_convergence_angle0/ (alpha = 0),
where the full-Cauchy reference has a physically h-thick layer at every h and
the constitutive thickness equals the kinematic length. EIQUAD4 columns only;
the warping element is not the subject of this paper.

Style is copied verbatim from master_parametric_study's plotting block so the
panels are indistinguishable from the elastic and von Mises ones: the error
expression is the TITLE (fontsize 32, pad 15), not a y-label; grid alpha 0.3;
LogLocator with numticks=6; figsize (7,6) with fixed margins; no per-panel
legend, since common_legend.pdf is placed alongside the row.

Usage:  plot_visco.py <convergence_table.md> <output_dir>
"""
import re, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

# Optional --limits=<json>: fixed ordinate limits per figure directory and panel, keyed by the
# output directory relative to figures/ (figures/axes_limits.json holds the limits the panels
# were printed with before the node-matched recomputation, so the replotted figures keep their
# axes). A limit is widened only where a value would otherwise fall outside it, and that is
# reported.
LIMITS = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--limits=")), None)
sys.argv = [a for a in sys.argv if not a.startswith("--limits=")]
TABLE = Path(sys.argv[1])
OUT = Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)

# Optional third argument: a second table whose range is merged into the ordinate of
# every panel. The single-reference and matched-thickness measurements of the same case
# agree at h = h0 and separate away from it; drawn on axes chosen independently, that
# shared point lands at a different height in each figure and the pair cannot be read
# against each other. Passing the counterpart table puts both on one scale.
PAIR = Path(sys.argv[3]) if len(sys.argv) > 3 else None

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

FAMILY_COLORS  = {"stiff": "#1f77b4", "initial": "#ff7f0e", "soft": "#d62728"}
FAMILY_MARKERS = {"stiff": "o", "initial": "s", "soft": "^"}

# section in convergence_table.md -> (output stem, title expression)
# The jump errors are reported either against the reference field norm and E^ref, or
# against the norm of the reference jump itself. Which one is recorded in the table's
# header, and the axis titles follow it so a figure always states its own normalisation.
_GLOBAL = r"$\frac{\|u^\mathrm{h}_i-u^0_i\|_{L^2}}{\|u^0_i\|_{L^2}}$"
TITLES = {
    "paper": {
        "Displacement jump":
            r"$\frac{\|[\![u^\mathrm{h}_i]\!] - [\![u^0_i]\!]\|_{L^2}}{\|u^0_i\|_{L^2}}$",
        "Traction jump":
            r"$\frac{\|[\![t^\mathrm{h}_i]\!] - [\![t^0_i]\!]\|_{L^2}}{E^\mathrm{ref}}$"},
    "self": {
        "Displacement jump":
            r"$\frac{\|[\![u^\mathrm{h}_i]\!] - [\![u^0_i]\!]\|_{L^2}}"
            r"{\|[\![u^0_i]\!]\|_{L^2}}$",
        "Traction jump":
            r"$\frac{\|[\![t^\mathrm{h}_i]\!] - [\![t^0_i]\!]\|_{L^2}}"
            r"{\|[\![t^0_i]\!]\|_{L^2}}$"},
}
_MODE = "self" if "mode: self" in TABLE.read_text() else "paper"
# The interface traction, relative to the traction being reproduced. Both sections are drawn
# when present, so the jump panels stay reproducible from the same table.
_TRACTION = r"$\frac{\|t^\mathrm{h}_i - t^0_i\|_{L^2}}{\|t^0_i\|_{L^2}}$"
SECTIONS = {
    "Displacement field (bodies)": ("global_relative_error", _GLOBAL),
    "Displacement jump": ("displacement_jump", TITLES[_MODE]["Displacement jump"]),
    "Traction jump": ("traction_jump", TITLES[_MODE]["Traction jump"]),
    "Traction": ("traction", _TRACTION),
}


def parse(text):
    """Pull each family's column out of each markdown section, whatever the element."""
    out = {}
    for sec in SECTIONS:
        m = re.search(rf"##\s*{re.escape(sec)}\s*\n(.*?)(?=\n##|\Z)", text, re.S)
        if not m:
            continue
        lines = [l for l in m.group(1).strip().splitlines() if l.startswith("|")]
        header = [c.strip() for c in lines[0].strip("|").split("|")]
        cols = {name: i for i, name in enumerate(header)}
        data = {}
        for fam in FAMILY_COLORS:
            key = next((k for k in cols if k == fam or k.startswith(f"{fam} ")), None)
            if key is None:
                continue
            hs, es = [], []
            for row in lines[2:]:
                cells = [c.strip() for c in row.strip("|").split("|")]
                if cells[0].lower() == "rate":
                    continue
                hs.append(float(cells[0])); es.append(float(cells[cols[key]]))
            data[fam] = (np.array(hs), np.array(es))
        out[sec] = data
    return out


tables = parse(TABLE.read_text())
pair = parse(PAIR.read_text()) if PAIR else {}

for sec, (stem, title_text) in SECTIONS.items():
    if sec not in tables or not tables[sec]:
        print(f"  missing section: {sec}")
        continue

    fig, ax = plt.subplots(figsize=(7, 6))
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.12)

    for fam in ("stiff", "initial", "soft"):
        if fam not in tables[sec]:
            continue
        h, e = tables[sec][fam]
        order = np.argsort(h); h, e = h[order], e[order]
        plt.loglog(h, e, linestyle="-", linewidth=2,
                   marker=FAMILY_MARKERS[fam], markersize=7,
                   label=fam.capitalize(), color=FAMILY_COLORS[fam])

    if sec in pair and pair[sec]:
        lo = min(float(np.min(v[1])) for v in list(tables[sec].values()) + list(pair[sec].values())
                 if np.min(v[1]) > 0)
        hi = max(float(np.max(v[1])) for v in list(tables[sec].values()) + list(pair[sec].values()))
        ax.set_ylim(lo / 3.0, hi * 3.0)

    if LIMITS:
        import json
        key = OUT.resolve().as_posix().split("/figures/", 1)[1].rstrip("/")
        lo, hi = json.loads(Path(LIMITS).read_text())[key][stem]
        # both tables of a pair, so that the two panels widen together and keep one ordinate
        vals = np.concatenate([v[1] for v in list(tables[sec].values())
                               + list(pair.get(sec, {}).values())])
        vals = vals[vals > 0]
        if vals.min() < lo or vals.max() > hi:
            print(f"  {stem}: limits [{lo:.1e}, {hi:.1e}] widened to hold [{vals.min():.1e}, {vals.max():.1e}]")
            lo, hi = min(lo, vals.min() / 1.5), max(hi, vals.max() * 1.5)
        ax.set_ylim(lo, hi)

    plt.xlabel(r"$h$ [mm]")
    plt.title(title_text, fontsize=32, pad=15)
    plt.grid(True, alpha=0.3)
    ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))

    fig.savefig(OUT / f"{stem}.pdf")
    plt.close(fig)

    rates = {f: np.polyfit(np.log(h), np.log(e), 1)[0]
             for f, (h, e) in tables[sec].items()}
    print(f"  {stem}: " + "  ".join(f"{f}={r:.3f}" for f, r in rates.items()))

print("wrote:", ", ".join(p.name for p in sorted(OUT.glob("*.pdf"))))
