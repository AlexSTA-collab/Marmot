#!/usr/bin/env python3
"""In-plane refinement of the EIE, the PSIE without and the PSIE with the squeeze-flow term (Appendix L, reviewer 2.9):
stiff family, alpha = 10, h = 0.01, sigma_y = 5 MPa, H = 100 and 0.1 MPa, interface meshes N x N = 15, 30, 60, 120 with
the bodies refined in proportion (refinement/mesh_if/L-1 ... L2).

Runs (end of the loading):
  EIE           eie_inplane/H_<H>/L<l>     (H_100/L1 is osc100/eie, the same job; H = 0.1 with the reltol Marmot, as the
                                            paper's EIE at H = 0.1)
  PSIE, no SQ   psie_inplane/H_<H>/L<l>
  PSIE, SQ      psie_squeeze_frozen/H_<H>/L<l>   (mu_t fixed within an increment, the element of the paper)

1. Self-convergence: every mesh against the 120 x 120 mesh of the same element, in the measures of the paper
   (refinement/errors_end.py, after port_2611/scripts/node_matched_errors.py), evaluated on the coarser mesh, whose
   nodes and cells are nested in the finer one:
     e_u    ||u^N - u^120|| / ||u^120|| over both bodies, the fine field taken at the coarse nodes;
     e_[u]  W ||[u^N] - [u^120]|| / ||u^120||, at the coarse face nodes, W = sqrt(0.5)/N, edge margin 0.15 mm;
     e_t    ||t^N - t^120|| / ||t^120||, the fine traction averaged over each coarse cell, edge margin 0.15 mm.
2. The interface pressure at the end of the loading on every mesh (the field of Fig. 21 a,b: the PSIE's own nodal
   pressure, for the EIE the pressure -tr(sigma+)/3 of its face stress).

    python3 inplane_study.py <figure dir>   ->  <dir>/H_<H>/{global_relative_error,displacement_jump,traction}.pdf,
                                               <dir>/legend_inplane.pdf, <dir>/contours_H<H>.pdf, inplane_scores.json
"""
import json
import sys
from glob import glob
from pathlib import Path

import matplotlib
import matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import cm, colors  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402
from matplotlib.ticker import LogLocator  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "refinement")); sys.path.insert(0, str(HERE.parent / "port_2611"))
import errors_end as E  # noqa: E402
from ensight_reader_2611 import read_geo  # noqa: E402

OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
LEVELS = {"L-1": 15, "L0": 30, "L1": 60, "L2": 120}
HS = ("100", "0.1")
ELEMENTS = {"eie": "EIE", "psie": r"PSIE, no squeeze flow", "sq": r"PSIE, squeeze flow"}
STYLE = {"eie": ("#1f77b4", "o", "--"), "psie": ("#d62728", "s", ":"), "sq": ("#111111", "D", "-")}  # EIE, PSIE as in Fig. 20
SUF = "0.01"


def run(elem, H, lev):
    if elem == "eie" and H == "100" and lev == "L1":
        return HERE / "osc100" / "eie"
    root = {"eie": "eie_inplane", "psie": "psie_inplane", "sq": "psie_squeeze_frozen"}[elem]
    return HERE / root / f"H_{H}" / lev


# ---------------------------------------------------------------- self-convergence
def scores(coarse, fine, n):
    out = {}
    cb = {s: E.body(coarse, "Cohesive_zone_model", s, SUF) for s in ("bottom", "top")}
    fb = {s: E.body(fine, "Cohesive_zone_model", s, SUF) for s in ("bottom", "top")}
    num = den = 0.0
    for s in ("bottom", "top"):
        (xc, cc, uc), (xf, _, uf) = cb[s], fb[s]
        uf_c = uf[E.matched(xc, xf)]                      # coarse nodes are fine nodes
        num += E.integrate(np.sum((uc - uf_c) ** 2, axis=1), xc, cc)
        den += E.integrate(np.sum(uf_c ** 2, axis=1), xc, cc)
    out["e_u"] = float(np.sqrt(num / den))

    jc, pc = E.jump(cb["bottom"], cb["top"], 0.01, 10)
    jf, pf = E.jump(fb["bottom"], fb["top"], 0.01, 10)
    jf_c = jf[E.matched(pc, pf)]
    k = E.interior(pc, 10, E.MARGIN)
    out["e_jump"] = float(np.sqrt(0.5) / n * np.sqrt(np.sum((jc - jf_c)[k] ** 2)) / np.sqrt(den))

    ccell, tc = E.model_traction(coarse, SUF)
    fcell, tf = E.model_traction(fine, SUF)
    owner = cKDTree(ccell).query(fcell)[1]                # the coarse cell each fine cell lies in
    cnt = np.bincount(owner, minlength=len(ccell))
    assert cnt.min() == cnt.max() == (len(fcell) // len(ccell)), "cells are not nested"
    tavg = np.stack([np.bincount(owner, weights=tf[:, i], minlength=len(ccell)) for i in range(3)], axis=1) / cnt[:, None]
    k = E.interior(ccell, 10, E.MARGIN)
    out["e_t"] = float(np.linalg.norm((tc - tavg)[k]) / np.linalg.norm(tavg[k]))
    return out


res = {}
for H in HS:
    for el in ELEMENTS:
        fine = run(el, H, "L2")
        for lev in ("L-1", "L0", "L1"):
            res[f"{el}/{H}/{LEVELS[lev]}"] = scores(run(el, H, lev), fine, LEVELS[lev])
            print(el, H, LEVELS[lev], "  ".join(f"{k}={v:.3e}" for k, v in res[f"{el}/{H}/{LEVELS[lev]}"].items()), flush=True)
(OUT / "inplane_scores.json").write_text(json.dumps(res, indent=1))

plt.rcParams.update({
    "text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"], "font.size": 16,
    "axes.labelsize": 28, "axes.titlesize": 32, "xtick.labelsize": 28, "ytick.labelsize": 28, "legend.fontsize": 18,
    "savefig.format": "pdf", "savefig.bbox": "tight", "figure.dpi": 100,
})
PANELS = (("global_relative_error", "e_u", r"$\frac{\|u^{N}_i-u^{120}_i\|_{L^2}}{\|u^{120}_i\|_{L^2}}$"),
          ("displacement_jump", "e_jump", r"$\frac{\|[\![u^{N}_i]\!] - [\![u^{120}_i]\!]\|_{L^2}}{\|u^{120}_i\|_{L^2}}$"),
          ("traction", "e_t", r"$\frac{\|t^{N}_i - t^{120}_i\|_{L^2}}{\|t^{120}_i\|_{L^2}}$"))
NS = np.array([15, 30, 60])
ds = 1.0 / NS
for H in HS:
    d = OUT / f"H_{H}"; d.mkdir(exist_ok=True)
    for stem, key, title in PANELS:
        fig, ax = plt.subplots(figsize=(7, 6))
        plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.12)
        for el in ELEMENTS:
            c, m, ls = STYLE[el]
            e = np.array([res[f"{el}/{H}/{n}"][key] for n in NS])
            plt.loglog(ds, e, linestyle=ls, linewidth=2, marker=m, markersize=7, color=c)
            print(f"  H={H} {stem} {el}: slope {np.polyfit(np.log(ds), np.log(e), 1)[0]:.2f}")
        plt.xlabel(r"$\Delta s=1/N$ [mm]")
        plt.title(title, fontsize=32, pad=15)
        plt.grid(True, alpha=0.3)
        ax.set_xticks(ds, [rf"$1/{n}$" for n in NS])
        ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax.set_xlim(ds[-1] / 1.25, ds[0] * 1.25)
        ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
        fig.savefig(d / f"{stem}.pdf")
        plt.close(fig)
handles = [Line2D([], [], color=STYLE[el][0], marker=STYLE[el][1], ls=STYLE[el][2], lw=2, ms=9, label=lab)
           for el, lab in ELEMENTS.items()]
fig = plt.figure(figsize=(2.2, 2.0))
fig.legend(handles=handles, loc="center", ncol=1, frameon=False, fontsize=13, handlelength=1.8)
fig.savefig(OUT / "legend_inplane.pdf")
plt.close(fig)
fig = plt.figure(figsize=(9.0, 0.45))                   # the same, in one row under the panels
fig.legend(handles=handles, loc="center", ncol=3, frameon=False, fontsize=13, handlelength=2.4, columnspacing=2.0)
fig.savefig(OUT / "legend_inplane_row.pdf")
plt.close(fig)

# ---------------------------------------------------------------- pressure contours
import os  # noqa: E402
if os.environ.get("PANELS_ONLY"):
    sys.exit(0)
import pyvista as pv  # noqa: E402

PIX = (1200, 1200)
# fixed colour ranges: at H = 0.1 the squeeze-flow row varies by a few tenths of a percent; it is drawn on the range of
# the PSIE without the term, shifted down so that the uniform field sits on the warm side
CLIM_FIXED = {("0.1", "sq"): (0.19, 0.365)}
plt.rcParams.update({"font.size": 11, "axes.labelsize": 12, "xtick.labelsize": 10, "ytick.labelsize": 10,
                     "figure.dpi": 200, "savefig.dpi": 400})


def pressure_mesh(elem, H, lev):
    d = run(elem, H, lev) / "Cohesive_zone_model_interface_0.01"
    parts = read_geo(str(d / "geometry.geo_0000"))
    pn = next(k for k, v in parts.items() if v["desc"] == "interface-body")
    name = "p_plus" if elem == "eie" else "pressure_interface"
    p = E.read_elem_var(sorted(glob(str(d / f"{name}.var_[0-9][0-9][0-9][0-9]")))[-1], parts, 1)[pn][:, 0]
    m = pv.UnstructuredGrid({pv.CellType.HEXAHEDRON: parts[pn]["blocks"][0][1]}, parts[pn]["xyz"])
    m.cell_data["p"] = p
    return m


def interior_values(mesh, trim=0.18):
    c = mesh.cell_centers().points
    lo, hi = c.min(0), c.max(0)
    span = hi - lo
    keep = np.ones(len(c), bool)
    for k in range(3):
        if span[k] > 1e-12:
            keep &= (c[:, k] > lo[k] + trim * span[k]) & (c[:, k] < hi[k] - trim * span[k])
    return np.asarray(mesh.cell_data["p"]).ravel()[keep]


def render(mesh, clim, n, zoom):
    p = pv.Plotter(off_screen=True, window_size=PIX)
    p.set_background("white")
    p.add_mesh(mesh, scalars="p", cmap="RdBu_r", clim=clim, show_edges=True, edge_color="#5a5a5a",
               line_width=0.35 if n <= 60 else 0.15, show_scalar_bar=False, lighting=False,
               interpolate_before_map=False)
    p.enable_anti_aliasing("ssaa")
    p.view_vector((0.30, -0.62, 0.72), viewup=(0, 0, 1))
    p.camera.zoom(zoom)
    img = p.screenshot(return_img=True)
    p.close()
    return np.asarray(img)


def crop_square(img, frac):
    h, w = img.shape[:2]
    s = int(min(h, w) * frac)
    return img[(h - s) // 2:(h - s) // 2 + s, (w - s) // 2:(w - s) // 2 + s]


for H in HS:
    fig = plt.figure(figsize=(14.0, 10.6))
    gs = fig.add_gridspec(3, 6, width_ratios=[1, 1, 1, 1, 0.035, 0.22], wspace=0.04, hspace=0.12)
    allm = {el: {lev: pressure_mesh(el, H, lev) for lev in LEVELS} for el in ELEMENTS}
    lims = {}
    for group in (("eie",), ("psie",), ("sq",)):  # one colour range per row
        vals = np.concatenate([interior_values(m) for el in group for m in allm[el].values()])
        mid = float(np.median(vals)); amp = float(np.percentile(np.abs(vals - mid), 98))
        for el in group:
            lims[el] = (mid - amp, mid + amp)
    lims.update({e: c for (h, e), c in CLIM_FIXED.items() if h == H})
    for r, el in enumerate(ELEMENTS):
        meshes, clim = allm[el], lims[el]
        for c, (lev, n) in enumerate(LEVELS.items()):
            ax = fig.add_subplot(gs[r, c])
            full = render(meshes[lev], clim, n, 1.0)
            ax.imshow(full); ax.set_axis_off()
            if r == 0:
                ax.set_title(rf"${n}\times{n}$", fontsize=13, pad=4)
            if c == 0:
                ax.text(-0.04, 0.5, ELEMENTS[el], transform=ax.transAxes, rotation=90, ha="right", va="center",
                        fontsize=12)
            zoom_img = crop_square(render(meshes[lev], clim, n, 9.0), 0.55)
            Hh, W = full.shape[:2]
            ax.add_patch(Circle((0.50 * W, 0.55 * Hh), 0.075 * W, fill=False, lw=1.0, ec="#111111"))
            ins = ax.inset_axes([0.62, 0.62, 0.40, 0.40])
            ins.imshow(zoom_img); ins.set_xticks([]); ins.set_yticks([])
            for s in ins.spines.values():
                s.set_linewidth(1.0)
        cax = fig.add_subplot(gs[r, 4])
        fig.colorbar(cm.ScalarMappable(norm=colors.Normalize(*clim), cmap="RdBu_r"), cax=cax)
        cax.set_ylabel(r"$p$ [MPa]", fontsize=12)
    fig.savefig(OUT / f"contours_H{H}.pdf")
    plt.close(fig)
    print("wrote", OUT / f"contours_H{H}.pdf", flush=True)
