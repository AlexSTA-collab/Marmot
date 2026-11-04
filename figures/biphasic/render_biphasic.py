#!/usr/bin/env python3
"""Deformed shape of the biphasic benchmark: resolved layer, [C] kept, [C] removed (Figure 10).

    xvfb-run -a python3 render_biphasic.py <out dir> [case angle h]      (default c100 10 0.1)

Same composition as render_deformed.py (Fig. 12b of the paper): the whole three-body model coloured by
||u||/max||u||, one camera and one magnification for all panels, and a circular magnification of the layer. In the
inset the displacement is magnified ZOOM_FACTOR times, relative to the centre of the interface so that the layer stays in
view, on the same colour scale. The resolved layer is drawn with all its elements. The interface models are drawn as they
are, with the single interface element between the two bodies; nothing is added to it. The inset of the resolved layer
carries the shear strain gamma = 2 eps_tn of each phase (t, n tangent and normal of the interface), read off the
displacement gradient of the layer elements, and the script prints it next to t_t / G of each phase as a check.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
CASE, ANGLE, H = (sys.argv[2], int(sys.argv[3]), sys.argv[4]) if len(sys.argv) > 4 else ("c100", 10, "0.1")
ET, EB = {"c100": (4e5, 4e3), "c10": (4e4, 4e3)}[CASE]
NU = 0.3
G = lambda E: E / (2 * (1 + NU))

# borrow the renderer and the helpers of render_deformed.py, up to its driver loop
src = (HERE.parent / "paper_scripts/figures/deformed/render_deformed.py").read_text()
src = src[:src.index("HI = (4200, 4500)")]
src = src.replace('scalars="umag", cmap=CMAP, clim=clim,', 'scalars=SC[0], cmap=SC[1], clim=clim, nan_color="#d0d0d0",')
sys.argv = [sys.argv[0], str(OUT)]
src = src.replace('edge_color="#2b2b2b", line_width=0.6', 'edge_color=EDGE[0], line_width=EDGE[1]')
CMAP_T = None
SC = ["umag", "coolwarm"]
EDGE = ["#2b2b2b", 0.6]
exec(compile(src, "render_deformed_helpers", "exec"))
import pyvista as pv  # noqa: E402

# The interphase sits at 0.90 ... 0.98 of the largest displacement, the top of the scale, where a linear map renders it as
# one saturated colour. The colour is therefore a power of the displacement, v**GAMMA with v = ||u||/max||u||, which spreads
# the top of the range; the colour bar shows the same map, its ticks at the values of v.
GAMMA = 6.0
CMAP_T = plt.cm.coolwarm

ROOT = HERE / __import__("os").environ.get("ROOT", "runs_as2605") / CASE
HI = (4200, 4500)
ZOOM_FACTOR = float(__import__('os').environ.get('ZOOM_FACTOR', 120))
a = np.radians(ANGLE)
n_vec = np.array([np.sin(a), 0.0, np.cos(a)])
t_vec = np.array([np.cos(a), 0.0, -np.sin(a)])


def view_case(run, stem):
    """EnSight case of the displacement of one output, in the form of the 26.05 outputs (a FILE section), which
    the converted .case of to_2605_layout.py lacks; VTK's reader needs it to find the time steps."""
    d = Path(run) / stem
    var = next(d.glob("displacement_*.var"))
    steps = var.read_bytes().count(b"BEGIN TIME STEP")
    values = "\n".join(f"{float(i):.8e}" for i in range(steps))
    case = Path(run) / f"{stem}.view.case"
    case.write_text(f"FORMAT\ntype: ensight gold\nTIME\ntime set: 2 no description\nnumber of steps: 1\n"
                    f"filename start number: 0\nfilename increment: 1\ntime values: 0.00000000e+00\n"
                    f"time set: 1 no description\nnumber of steps: {steps}\nfilename start number: 0\n"
                    f"filename increment: 1\ntime values: {values}\nFILE\nfile set: 2\nnumber of steps: 1\n"
                    f"file set: 1\nnumber of steps: {steps}\nGEOMETRY\nmodel: 2 2 {stem}/geometry.geo\n"
                    f"VARIABLE\nvector per node: 1 1 {var.stem} {stem}/{var.name}\n")
    return case


def part(model, prefix, name):
    run = ROOT / model / f"angle_{ANGLE}" / f"h_{H}"
    return read_body(view_case(run, f"{prefix}_{name}_{H}"))


def layer_strain(m):
    """gamma = 2 eps_tn per element of a layer half, from the displacement gradient."""
    g = m.compute_derivative(scalars="u", gradient="grad").point_data_to_cell_data()
    grad = np.asarray(g.cell_data["grad"]).reshape(-1, 3, 3)          # grad[k, i, j] = du_i/dx_j
    eps = 0.5 * (grad + grad.transpose(0, 2, 1))
    gam = 2 * np.einsum("i,kij,j->k", t_vec, eps, n_vec)
    m.cell_data["gam"] = gam
    return gam


lam = lambda E: E * NU / ((1 + NU) * (1 - 2 * NU))
Q = lambda E: G(E) * np.eye(3) + (lam(E) + G(E)) * np.outer(n_vec, n_vec)


def traction_vec(model):
    import node_matched_errors as nme
    run = ROOT / model / f"angle_{ANGLE}" / f"h_{H}"
    c, tr = nme.model_traction(run, H)
    return np.asarray(tr[np.argmin(np.linalg.norm(c[:, :2], axis=1))])


def traction_t(model):
    """Tangential traction of the interface element nearest the centre of the interface."""
    import node_matched_errors as nme
    run = ROOT / model / f"angle_{ANGLE}" / f"h_{H}"
    c, tr = nme.model_traction(run, H)
    k = np.argmin(np.linalg.norm(c[:, :2], axis=1))
    return float(tr[k] @ t_vec)


def halves(m, gam_top, gam_bot, g_bot):
    """Split the hexahedra of an interface element body at the mid-surface into two layers carrying their own strain.
    The nodes of the mid-surface carry the mean displacement of the two faces."""
    pts, disp, cells, vals = [], [], [], []
    X, U = np.asarray(m.points), np.asarray(m.point_data["u"])
    for c in m.cells_dict[pv.CellType.HEXAHEDRON]:
        # the element lists the nodes of one face first and those of the other face after, in matching order
        lo, hi = (c[:4], c[4:]) if (X[c[:4]] @ n_vec).mean() < (X[c[4:]] @ n_vec).mean() else (c[4:], c[:4])
        lat = (X[hi] - X[lo]) - ((X[hi] - X[lo]) @ n_vec)[:, None] * n_vec
        assert np.abs(lat).max() < 0.5 * np.abs((X[hi] - X[lo]) @ n_vec).min(), "face nodes not paired"
        hn = ((X[hi] - X[lo]) @ n_vec).mean()
        mid_x, mid_u = 0.5 * (X[lo] + X[hi]), U[lo] + 0.5 * hn * g_bot[None, :]
        for (px, pu), (qx, qu), g in (((X[lo], U[lo]), (mid_x, mid_u), gam_bot),
                                      ((mid_x, mid_u), (X[hi], U[hi]), gam_top)):
            base = len(pts)
            pts += list(px) + list(qx); disp += list(pu) + list(qu)
            cells.append([8] + list(range(base, base + 8))); vals.append(g)
    grid = pv.UnstructuredGrid(np.hstack(cells), np.full(len(cells), pv.CellType.HEXAHEDRON), np.array(pts))
    grid.cell_data["gam"] = np.array(vals)
    grid.point_data["u"] = np.array(disp)
    grid.point_data["umag"] = np.linalg.norm(np.array(disp), axis=1)
    return grid


# --- data -----------------------------------------------------------------------------
res = {k: part("ref", "Cauchy_Full_model", k) for k in ("bottom", "top", "layerbottom", "layertop")}
gam_bot = layer_strain(res["layerbottom"]); gam_top = layer_strain(res["layertop"])
for k in ("bottom", "top"):
    res[k].cell_data["gam"] = np.full(res[k].n_cells, np.nan)
tau = traction_t("full")
print(f"resolved layer: gamma bottom {np.median(gam_bot):.4e}  top {np.median(gam_top):.4e}   "
      f"traction-based: bottom {tau / G(EB):.4e}  top {tau / G(ET):.4e}")

STRAIN = {"ref": (float(np.median(gam_top)), float(np.median(gam_bot)))}
models = {"ref": [res["bottom"], res["layerbottom"], res["layertop"], res["top"]]}
for name in ("full", "supp"):
    p = {k: part(name, "Cohesive_zone_model", k) for k in ("bottom", "top", "interface")}
    tau_m = traction_t(name)
    gt, gb = (tau_m / G(ET), tau_m / G(EB)) if name == "full" else (tau_m / G(0.5 * (ET + EB)),) * 2
    strip = p["interface"]
    STRAIN[name] = (gt, gb)
    # mean shear of the single element, the relative displacement of its faces over the normal separation: the check
    # of what the figure shows, against the mean of the two face strains
    Xs, Us = np.asarray(strip.points), np.asarray(strip.point_data["u"])
    sh = []
    for c in strip.cells_dict[pv.CellType.HEXAHEDRON]:
        lo, hi = (c[:4], c[4:]) if (Xs[c[:4]] @ n_vec).mean() < (Xs[c[4:]] @ n_vec).mean() else (c[4:], c[:4])
        sh.append(((Us[hi] - Us[lo]).mean(axis=0) @ t_vec) / ((Xs[hi] - Xs[lo]).mean(axis=0) @ n_vec))
    print(f"{name}: shear of the interface element (median over elements) {np.median(sh):.4e}, "
          f"mean of the two face strains {0.5 * (gt + gb):.4e}")
    models[name] = [p["bottom"], strip, p["top"]]

# --- scales -----------------------------------------------------------------------------
allm = [m for v in models.values() for m in v]
umax = max(float(np.max(m.point_data["umag"])) for m in allm)
span = float(np.linalg.norm(np.ptp(np.vstack([m.points for m in models["ref"]]), axis=0)))
fac = 0.05 * span / max(umax, 1e-30)
for m in allm:
    m.point_data["umag"] = (m.point_data["umag"] / max(umax, 1e-30)) ** GAMMA
GCLIM = (0.0, 1.05 * float(np.median(gam_bot)))
print(f"umax {umax:.3e} mm, magnification x{fac:.0f}, strain scale 0 .. {GCLIM[1]:.3e}")

order = [("ref", "resolved biphasic layer"), ("full", r"interface, $[C]$ kept"), ("supp", r"interface, $[C]$ removed")]
warped = {k: warp(models[k], fac) for k, _ in order}


def focus_of(k):
    """Centre of the layer on its front face (y minimal), which is the face the camera sees."""
    layer = [m for m in warped[k] if m is not warped[k][0] and m is not warped[k][-1]]
    xyz = np.vstack([m.points for m in layer])
    return (float(xyz[:, 0].mean()), float(xyz[:, 1].min()), float(xyz[:, 2].mean()))


SC[0], SC[1] = "umag", CMAP_T
imgs, cam = {}, None
for k, _ in order:
    img, cam_k, _ = render(warped[k], (0, 1), cam=cam, size=HI)
    cam = cam if cam is not None else cam_k
    imgs[k] = img
HALF = int(0.03 * HI[0])
pts = {k: project(warped[k], (0, 1), cam, focus_of(k), HI) for k, _ in order}
SC[0], SC[1] = "umag", CMAP_T
EDGE[0], EDGE[1] = "#f7f7f7", 2.4        # the mesh of the inset is the point of it: light, heavy lines
zoom = {}
for k, _ in order:
    layer = [m for m in models[k] if m is not models[k][0] and m is not models[k][-1]]
    shift = np.mean(np.vstack([m.point_data["u"] for m in layer]), axis=0)
    c0 = np.mean(np.vstack([m.points + fac * m.point_data["u"] for m in layer]), axis=0)
    meshes = []
    for m in models[k]:
        w = m.copy()
        w.points = w.points + ZOOM_FACTOR * (np.asarray(w.point_data["u"]) - shift)
        if "gam" not in w.cell_data:
            w.cell_data["gam"] = np.full(w.n_cells, np.nan)
        meshes.append(w)
    # the magnified layer is put where the layer of the whole-model view is, so the same crop frames it
    c1 = np.mean(np.vstack([w.points for w in meshes if w in meshes[1:-1] or True][1:-1]), axis=0)
    meshes = [w for w in meshes]
    for w in meshes:
        w.points = w.points + (c0 - c1)
    meshes = [w.extract_surface() for w in meshes]
    zimg, _, _ = render(meshes, (0, 1), cam=cam, size=HI)
    zoom[k] = crop_at(zimg, pts[k], HALF)

fig = plt.figure(figsize=(7.4, 2.9))
gs = fig.add_gridspec(1, 3, wspace=0.03, left=0.01, right=0.89)
for i, (k, ttl) in enumerate(order):
    ax = fig.add_subplot(gs[0, i])
    ax.imshow(imgs[k]); ax.set_axis_off(); ax.set_title(ttl, fontsize=9, pad=3, loc="left")
    ax.add_patch(Circle(pts[k], HALF, fill=False, lw=1.0, ec="#111111"))
    ins = ax.inset_axes([0.50, 0.03, 0.52, 0.52])
    ins.imshow(zoom[k]); ins.set_xticks([]); ins.set_yticks([]); ins.patch.set_alpha(0)
    for s_ in ins.spines.values():
        s_.set_visible(False)
    R = zoom[k].shape[0] / 2
    ins.images[0].set_clip_path(Circle((zoom[k].shape[1] / 2, R), R, transform=ins.transData))
    ins.add_patch(Circle((zoom[k].shape[1] / 2, R), R, fill=False, lw=1.0, ec="#111111"))
    gt, gb = STRAIN[k]
    fmt = lambda v: r"$\gamma=%.1f\cdot10^{%d}$" % (v / 10 ** np.floor(np.log10(v)), np.floor(np.log10(v)))
    for y, v in (((0.80, gt), (0.20, gb)) if k == "ref" else ()):
        ins.text(0.98, y, fmt(v), transform=ins.transAxes, ha="right", va="center", fontsize=7.5,
                 bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
    for ang in (1.0, -0.45):
        ax.add_artist(ConnectionPatch(xyA=(pts[k][0] + HALF * np.cos(ang), pts[k][1] - HALF * np.sin(ang)),
                                      xyB=(zoom[k].shape[1] / 2 + R * np.cos(ang), R - R * np.sin(ang)),
                                      coordsA="data", coordsB="data", axesA=ax, axesB=ins, lw=0.7, color="#111111"))
c1 = fig.add_axes([0.92, 0.12, 0.014, 0.76])
cb = fig.colorbar(cm.ScalarMappable(norm=colors.PowerNorm(GAMMA, 0, 1), cmap=CMAP_T), cax=c1,
                  ticks=[0, 0.8, 0.9, 0.95, 1.0])
cb.set_label(r"$\|u_i\|\,/\,\max\|u_i\|$", labelpad=4); cb.outline.set_linewidth(0.6)
fig.savefig(OUT / "deformed_biphasic.pdf")
fig.savefig(OUT / "deformed_biphasic.png", dpi=220)
print("wrote", OUT / "deformed_biphasic.pdf")
