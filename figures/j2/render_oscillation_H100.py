#!/usr/bin/env python3
"""H = 100 MPa version of paper/figures/oscillation/render_oscillation.py (Figure 21): the interface pressure field of the
EIE (YIQUAD4) and the PSIE (YIQUAD4_STABP_MINI) at the end of the analysis, stiff family, alpha = 10, h = 0.01, 60 x 60,
from the runs osc100/{eie,psie} (jobs/<element>/stiff/angle_10/H_100/h_0.01 with the extra output p_plus = -tr(sigma+)/3,
mean over the four quadrature points). The PSIE is drawn with its own pressure field, the EIE, which has none, with p_plus. Everything else unchanged, including
panel (c).

    python3 render_oscillation_H100.py <out dir>

Original docstring:
Interface traction/pressure field with and without the MINI bubble.

Renders the two analyses of the controlled ablation -- byte-identical except
for MARMOT_MINI_BUBBLE -- as a two-panel figure with a circular ROI zoom, in
the manner of Fig. 15 of Stathas & Neuner (2026).

PyVista does the 3D rendering off-screen; the composition, the ROI inset and
the colour bar are done in matplotlib so the typography matches the rest of the
paper and nothing looks like a screenshot: no orientation widget, no bounding
box, no default scalar bar, white ground, shared symmetric colour range.
"""
import sys
from pathlib import Path
import pathlib
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, ConnectionPatch
from matplotlib import cm, colors

BASE = Path("/home/alexsta1993/Documents/STATHAS_VIENNA_POSTDOC/Applications_Examples/"
            "EdelweisModels/Cohesive_Zone_models/3D/Three_Body_Cubit/"
            "Comparizon_LinearElastic/New_Functionality/mini_bubble_ablation")
OUT = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
OUT.mkdir(parents=True, exist_ok=True)

FIELD = "pressure_interface"     # the quantity the ablation quantifies
CMAP = "RdBu_r"
PIX = (1500, 1500)
INTERFACE_BLOCK = 1              # the 900-cell interface patch

plt.rcParams.update({
    "text.usetex": True, "font.family": "serif",
    "font.serif": ["Computer Modern Roman"], "font.size": 11,
    "axes.labelsize": 12, "savefig.format": "pdf",
    "savefig.bbox": "tight", "figure.dpi": 200, "savefig.dpi": 600,
})


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "refinement")); sys.path.insert(0, str(HERE.parent / "port_2611"))
import errors_end as E                                  # noqa: E402
from ensight_reader_2611 import read_geo                # noqa: E402
import os
from glob import glob                                   # noqa: E402
OSC = os.environ.get("OSC_DIR", "osc100")   # osc100_sq: the PSIE with the squeeze-flow term
RUN = {"OFF": HERE / OSC / "eie", "ON": HERE / OSC / "psie"}


def load(tag):
    d = RUN[tag] / "Cohesive_zone_model_interface_0.01"
    parts = read_geo(str(d / "geometry.geo_0000"))
    pn = next(k for k, v in parts.items() if v["desc"] == "interface-body")
    # the PSIE has its own pressure field; in the EIE the pressure is that of the face stress
    name = "pressure_interface" if tag == "ON" else "p_plus"
    p = E.read_elem_var(sorted(glob(str(d / f"{name}.var_[0-9][0-9][0-9][0-9]")))[-1], parts, 1)[pn][:, 0]
    m = pv.UnstructuredGrid({pv.CellType.HEXAHEDRON: parts[pn]["blocks"][0][1]}, parts[pn]["xyz"])
    m.cell_data[FIELD] = p
    return m


def render(mesh, clim, zoom=1.0, focus=None):
    """Off-screen render of the interface patch, viewed along the normal."""
    p = pv.Plotter(off_screen=True, window_size=PIX)
    p.set_background("white")
    p.add_mesh(mesh, scalars=FIELD, cmap=CMAP, clim=clim,
               show_edges=True, edge_color="#5a5a5a", line_width=0.35,
               show_scalar_bar=False, lighting=False, interpolate_before_map=False)
    p.enable_anti_aliasing("ssaa")
    p.view_vector((0.30, -0.62, 0.72), viewup=(0, 0, 1))
    p.camera.zoom(zoom)
    if focus is not None:
        p.camera.focal_point = focus
    img = p.screenshot(return_img=True)
    p.close()
    return np.asarray(img)


def crop_square(img, frac=0.34, cx=0.5, cy=0.5):
    h, w = img.shape[:2]
    s = int(min(h, w) * frac)
    x0 = int(np.clip(cx * w - s // 2, 0, w - s))
    y0 = int(np.clip(cy * h - s // 2, 0, h - s))
    return img[y0:y0 + s, x0:x0 + s]


meshes = {t: load(t) for t in ("OFF", "ON")}

def interior(mesh, trim=0.18):
    """Cell values away from the patch boundary.

    The edges of the patch carry a smooth boundary layer an order of magnitude
    larger than the interior field. Scaling the colour map to it would hide the
    alternation that is the subject of the figure, so the range is taken from
    the interior, exactly as the ablation diagnostic trims its sampled row."""
    c = mesh.cell_centers().points
    lo, hi = c.min(0), c.max(0)
    span = hi - lo
    keep = np.ones(len(c), bool)
    for k in range(3):
        if span[k] > 1e-12:
            keep &= (c[:, k] > lo[k] + trim * span[k]) & (c[:, k] < hi[k] - trim * span[k])
    return np.asarray(mesh.cell_data[FIELD]).ravel()[keep]

ref = interior(meshes["OFF"])
mid = float(np.median(ref))
amp = float(np.percentile(np.abs(ref - mid), 98))
clim = (mid - amp, mid + amp)

# One row: the two pressure fields, their colour bar, and the hardening sweep that
# says when the oscillation appears. Keeping them in a single figure is what makes the
# three panels share a height.
fig = plt.figure(figsize=(17.6, 4.5))
# An empty column separates the colour bar from the hardening panel; without it the
# colour-bar label and the panel's tick labels and axis label sit on top of one another.
gs = fig.add_gridspec(1, 5, width_ratios=[1, 1, 0.045, 0.30, 0.95], wspace=0.06)
axes = [fig.add_subplot(gs[0, i]) for i in range(2)]
cax = fig.add_subplot(gs[0, 2])
axh = fig.add_subplot(gs[0, 4])

TITLES = {"OFF": r"(a)  EIE element", "ON": r"(b)  PSIE element"}
for ax, tag in zip(axes, ("OFF", "ON")):
    full = render(meshes[tag], clim, zoom=1.32)
    ax.imshow(full)
    ax.set_axis_off()
    ax.set_title(TITLES[tag], fontsize=13, pad=5, loc="left")

    # circular ROI on the patch, and its magnified copy inset top-right
    zoom_img = crop_square(render(meshes[tag], clim, zoom=9.0), frac=0.55)
    H, W = full.shape[:2]
    roi_xy, roi_r = (0.50 * W, 0.55 * H), 0.075 * W
    ax.add_patch(Circle(roi_xy, roi_r, fill=False, lw=1.2, ec="#111111"))

    ins = ax.inset_axes([0.60, 0.60, 0.42, 0.42])
    ins.imshow(zoom_img)
    ins.set_xticks([]); ins.set_yticks([])
    ins.patch.set_alpha(0)
    for s in ins.spines.values():
        s.set_visible(False)
    clip = Circle((zoom_img.shape[1] / 2, zoom_img.shape[0] / 2),
                  zoom_img.shape[0] / 2, transform=ins.transData)
    ins.images[0].set_clip_path(clip)
    ins.add_patch(Circle((zoom_img.shape[1] / 2, zoom_img.shape[0] / 2),
                         zoom_img.shape[0] / 2, fill=False, lw=1.2, ec="#111111"))
    for ang in (0.65, -0.65):
        ax.add_artist(ConnectionPatch(
            xyA=(roi_xy[0] + roi_r * np.cos(ang), roi_xy[1] - roi_r * np.sin(ang)),
            xyB=(zoom_img.shape[1] / 2 + zoom_img.shape[0] / 2 * np.cos(ang),
                 zoom_img.shape[0] / 2 - zoom_img.shape[0] / 2 * np.sin(ang)),
            coordsA="data", coordsB="data", axesA=ax, axesB=ins,
            lw=0.8, color="#111111"))

cb = fig.colorbar(cm.ScalarMappable(norm=colors.Normalize(*clim), cmap=CMAP),
                  cax=cax)
cb.set_label(r"interface pressure  $p$ [MPa]", labelpad=12)
cb.outline.set_linewidth(0.6)

# --- (c) the hardening sweep -------------------------------------------------------
# 60 x 60 EIE runs, stiff and initial, h = 0.01: hsweep60/eie/<fam>/H_<H> (reltol Marmot), H = 100 from jobs/eie.
# Second gradient of the normal traction along the middle interior row of cells, three trimmed at each end
# (eq. alternation): mean |t(e+1) - 2 t(e) + t(e-1)| / ds^2, final step.
from matplotlib.ticker import LogLocator, NullFormatter
FAM_C = {"stiff": "#1f77b4", "initial": "#ff7f0e"}
FAM_M = {"stiff": "o", "initial": "s"}
HARDENING = ["0.1", "1", "10", "100", "1000"]
NRM = np.array([np.sin(np.radians(10.0)), 0.0, np.cos(np.radians(10.0))])
TRIM = 3


def d2_tn(run):
    d = run / "Cohesive_zone_model_interface_0.01"
    parts = read_geo(str(d / "geometry.geo_0000"))
    pn = next(k for k, v in parts.items() if v["desc"] == "interface-body")
    f = E.read_elem_var(sorted(glob(str(d / "force_interface.var_[0-9][0-9][0-9][0-9]")))[-1], parts, 3)[pn]
    tn = f @ NRM
    c = parts[pn]["xyz"][parts[pn]["blocks"][0][1]].mean(axis=1)
    ux, uy = np.unique(np.round(c[:, 0], 6)), np.unique(np.round(c[:, 1], 6))
    col, row = np.searchsorted(ux, np.round(c[:, 0], 6)), np.searchsorted(uy, np.round(c[:, 1], 6))
    sel = np.where(row == len(uy) // 2)[0]
    sel = sel[np.argsort(col[sel])]
    ds = float(np.mean(np.linalg.norm(np.diff(c[sel], axis=0), axis=1)))
    line = tn[sel][TRIM:len(sel) - TRIM]
    return float(np.mean(np.abs(np.diff(line, n=2)))) / ds ** 2


for fam in ("stiff", "initial"):
    xs, ys = [], []
    for Hh in HARDENING:
        run = (HERE / "jobs" / "eie" / fam / "angle_10" / "H_100" / "h_0.01" if Hh == "100"
               else HERE / "hsweep60" / "eie" / fam / f"H_{Hh}")
        log = run / "run.log"
        if not log.exists() or "Job computation time" not in log.read_text() or "Simulation failed" in log.read_text():
            print("missing", run); continue
        xs.append(float(Hh)); ys.append(d2_tn(run))
        print(f"{fam:8s} H={Hh:6s} |t_n,ss| = {ys[-1]:.4e} MPa/mm^2")
    axh.plot(xs, ys, "-" + FAM_M[fam], color=FAM_C[fam], lw=2.0, ms=7, mfc="none", mew=2.0, label=fam.capitalize())
axh.set_xscale("log"); axh.set_yscale("log")
axh.yaxis.set_major_locator(LogLocator(base=10.0, numticks=5))
axh.yaxis.set_minor_formatter(NullFormatter())
axh.xaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
axh.set_xlabel(r"$H$ [MPa]", fontsize=13, labelpad=6)
# The second gradient of the normal traction along the row in finite
# differences, averaged in absolute value over the row. Section 4.4.6 defines
# the symbol; written out it is too long for a panel title.
axh.set_title(r"(c)  $\overline{|t_{n,ss}|}$ [MPa/mm$^2$]",
              fontsize=13, pad=5, loc="left")
axh.tick_params(labelsize=11)
axh.grid(True, which="both", alpha=0.3)
axh.legend(frameon=False, fontsize=11)

fig.savefig(OUT / "mini_oscillation.pdf")
print("wrote", OUT / "mini_oscillation.pdf")
print(f"clim = {clim[0]:.4f} .. {clim[1]:.4f}")
