#!/usr/bin/env python3
"""The error over the creep hold (Appendix H) from the 60 x 60 viscoelastic runs of make_jobs60.py, alpha = 10: the
EIE (YIQUAD4) against the resolved layer of the same thickness at every output time.

    e_u(t) = ||u^h(t) - u^0(t)||_L2 / ||u^0(t)||_L2      over both bodies, as e_u of errors_end.measures

The meshes of the two models coincide in the bodies, so u^h is taken at its twin node and the L2 norm is the volume
integral of errors_end.integrate. The two Bochner norms are sup_t e_u and (int_0^T e_u^2 dt)^(1/2), T = 100 days, by
the trapezoidal rule on the output times (the end of the loading step, t = 1e-8, then every day). The undeformed state
at t = 0 is left out.

    python3 bochner60.py <out dir> [h of the history panel, 0.04]
    RERUN60_JOBS=<jobs root> reads the runs from there instead of jobs/ (the shared-node reruns); NO_PLOT=1 stops after
    bochner60.json; BOCHNER_JSON=<file> draws from that file instead of reading the runs.
        ->  bochner_Linf.pdf, bochner_L2t.pdf, error_history.pdf, common_legend.pdf, bochner60.json"""
import json
import os
import re
import sys
from glob import glob
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import LogLocator, NullFormatter  # noqa: E402

HERE = Path(__file__).resolve().parent
JOBS = Path(os.environ.get("RERUN60_JOBS", HERE / "jobs"))
sys.path.insert(0, str(HERE.parent / "refinement")); sys.path.insert(0, str(HERE.parent / "port_2611"))
import errors_end as E  # noqa: E402
from ensight_reader_2611 import read_geo, read_var_node_all_steps  # noqa: E402

OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
HH = sys.argv[2] if len(sys.argv) > 2 else "0.04"
FAMS, HS = ("stiff", "initial", "soft"), ("0.01", "0.02", "0.04", "0.06", "0.08", "0.1")
T = 100.0
COL = {"stiff": "#1f77b4", "initial": "#ff7f0e", "soft": "#d62728"}
MRK = {"stiff": "o", "initial": "s", "soft": "^"}
TRAPZ = getattr(np, "trapezoid", None) or np.trapz
plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"],
                     "font.size": 16, "axes.labelsize": 28, "axes.titlesize": 32, "xtick.labelsize": 28,
                     "ytick.labelsize": 28, "legend.fontsize": 18, "savefig.format": "pdf", "savefig.bbox": "tight",
                     "figure.dpi": 100})
PANEL_BOX = dict(left=0.25, right=0.95, top=0.82, bottom=0.17)   # one box for the three panels of the row


def history_body(run, prefix, side, h):
    d = E.casedir(run, prefix, side, h)
    parts = read_geo(str(d / "geometry.geo_0000"))
    files = sorted(glob(str(d / f"displacement_{side}.var_[0-9][0-9][0-9][0-9]")))
    pn = None
    U = []
    for f in files:
        u = read_var_node_all_steps(f, parts, 3)
        pn = pn or list(u)[0]
        U.append(u[pn][-1])
    case = Path(str(d) + ".case").read_text()
    t = np.array(re.findall(r"[-+]?\d+\.\d+e[-+]\d+", case.split("time values:")[-1]), float)[:len(files)]
    return parts[pn]["xyz"], parts[pn]["blocks"][0][1], np.array(U), t


def history(fam, h):
    ref = JOBS / "visco10" / "ref" / fam / "angle_10" / f"h_{h}"
    ifc = JOBS / "visco10" / "ifc" / fam / "angle_10" / f"h_{h}"
    num = den = 0.0
    for side in ("bottom", "top"):
        xr, cr, ur, tr = history_body(ref, "Cauchy_Full_model", side, h)
        xm, _, um, tm = history_body(ifc, "Cohesive_zone_model", side, h)
        assert len(tr) == len(tm) and np.allclose(tr, tm), (fam, h, len(tr), len(tm))
        w = np.zeros(len(xr))                       # integrate(sq) = w @ sq, errors_end.integrate
        np.add.at(w, cr.ravel(), np.repeat(E.volumes(xr, cr) / cr.shape[1], cr.shape[1]))
        um = um[:, E.matched(xr, xm), :]
        num = num + np.sum((ur - um) ** 2, axis=2) @ w
        den = den + np.sum(ur ** 2, axis=2) @ w
    k = tr > 0.0
    return tr[k], np.sqrt(num[k] / den[k])


def rate(v):
    v = np.asarray(v, float)
    return np.polyfit(np.log([float(h) for h in HS]), np.log(v), 1)[0]


res = json.loads(Path(os.environ["BOCHNER_JSON"]).read_text()) if os.environ.get("BOCHNER_JSON") else {}
for fam in FAMS if not res else ():
    for h in HS:
        t, e = history(fam, h)
        res[f"{fam}/{h}"] = {"t": t.tolist(), "e_u": e.tolist(), "Linf": float(e.max()),
                             "L2t": float(np.sqrt(TRAPZ(e ** 2, t)))}
        print(f"{fam:8s} h={h:5s} steps={len(t)} e(load)={e[0]:.3e} e(T)={e[-1]:.3e} sup={e.max():.3e} "
              f"L2t={res[f'{fam}/{h}']['L2t']:.3e}", flush=True)
for fam in FAMS:
    print(f"{fam:8s} rate sup {rate([res[f'{fam}/{h}']['Linf'] for h in HS]):.3f}  "
          f"L2t {rate([res[f'{fam}/{h}']['L2t'] for h in HS]):.3f}")
(OUT / "bochner60.json").write_text(json.dumps(res, indent=1))
if os.environ.get("NO_PLOT"):
    sys.exit(0)


def save(fig, name):
    with plt.rc_context({"savefig.bbox": "standard"}):   # full 7 x 6 canvas, so the three panels match
        fig.savefig(OUT / name)
    plt.close(fig)


for key, title, name in (("Linf", r"$\sup_t e_u$", "bochner_Linf.pdf"),
                         ("L2t", r"$\left(\int_0^T e_u^2\,\mathrm{d}t\right)^{1/2}$", "bochner_L2t.pdf")):
    fig, ax = plt.subplots(figsize=(7, 6)); plt.subplots_adjust(**PANEL_BOX)
    for fam in FAMS:
        ax.loglog([float(h) for h in HS], [res[f"{fam}/{h}"][key] for h in HS], "-" + MRK[fam], color=COL[fam],
                  lw=2.0, ms=9, mfc="none", mew=2.0)
    ax.set_xlabel(r"$h$ [mm]"); ax.set_title(title, pad=15)
    ax.grid(True, which="both", alpha=0.3)
    ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
    ax.yaxis.set_minor_formatter(NullFormatter())
    save(fig, name)

fig, ax = plt.subplots(figsize=(7, 6)); plt.subplots_adjust(**PANEL_BOX)
for fam in FAMS:
    t, e = (np.array(res[f"{fam}/{HH}"][k]) for k in ("t", "e_u"))
    ax.semilogy(t, e, "-", color=COL[fam], lw=2.4)
    ax.plot(t[0], e[0], MRK[fam], color=COL[fam], ms=11, mfc="white", mew=2.0)
ax.set_xlabel(r"$t$ [days]"); ax.set_title(r"$e_{u}(t)$", pad=15)
ax.set_xlim(-2, 100)
ax.grid(True, which="both", alpha=0.3)
ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
ax.yaxis.set_minor_formatter(NullFormatter())
save(fig, "error_history.pdf")

fig = plt.figure(figsize=(2.6, 6))
fig.legend(handles=[Line2D([], [], color=COL[f], marker=MRK[f], mfc="none", mew=2.0, lw=2.0, ms=9,
                           label=f.capitalize()) for f in FAMS], loc="center", frameon=False)
fig.savefig(OUT / "common_legend.pdf"); plt.close(fig)
print("wrote", OUT)
