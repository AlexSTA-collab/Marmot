#!/usr/bin/env python3
"""Deformed shape at the end of the analysis: resolved layer vs interface model (Drucker--Prager, EIE).

Stiff family, alpha = 10 deg, h = h0. Both panels share camera, colour scale and
displacement magnification, so the comparison is of shapes and not of settings.
The magnification is chosen from the model size so that the deformation is
visible; it is reported in the caption.

Composition and typography in matplotlib, as for the other rendered figure, so
nothing reads as a screenshot.
"""
import sys
from pathlib import Path
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.patches import Circle, ConnectionPatch

NF = Path("/home/alexsta1993/Documents/STATHAS_VIENNA_POSTDOC/Applications_Examples/"
          "EdelweisModels/Cohesive_Zone_models/3D/Three_Body_Cubit/"
          "Comparizon_LinearElastic/New_Functionality")
OUT = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
OUT.mkdir(parents=True, exist_ok=True)

CMAP = "coolwarm"    # as Figs. 4-5 of Stathas & Neuner (2026); lighter at the
                     # extremes than RdBu_r, which rendered almost black there

# Drucker-Prager version of the paper's render_deformed.py (Figure 20 c): initial family, alpha = 10, h = 0.01, phi = 30,
# psi = 10, resolved layer (20 C3D8) against the EIE; runs copied to druckerprager/global/initial_dp/{ref,eie}
G = Path(__file__).resolve().parent / "global" / "initial_dp"
CASES = {
    "druckerprager": dict(full=G / "ref", iface=G / "eie", fpre="Cauchy_Full_model", ipre="Cohesive_zone_model",
                          label="Drucker--Prager"),
}

plt.rcParams.update({
    "text.usetex": True, "font.family": "serif",
    "font.serif": ["Computer Modern Roman"], "font.size": 8,
    "axes.labelsize": 9, "savefig.format": "pdf",
    "savefig.bbox": "tight", "figure.dpi": 200,
})


def read_body(case):
    """Final-step mesh of one body, with the displacement attached."""
    if not Path(case).is_file():
        return None
    r = pv.get_reader(str(case))
    r.set_active_time_set(1)
    r.set_active_time_value(r.time_values[-1])
    blocks = [b for b in r.read() if b is not None and b.n_cells > 0
              and any("displacement" in k for k in b.point_data.keys())]
    if not blocks:
        return None
    b = max(blocks, key=lambda x: x.n_cells)
    key = next(k for k in b.point_data.keys() if "displacement" in k)
    b.point_data["u"] = np.asarray(b.point_data[key])
    b.point_data["umag"] = np.linalg.norm(b.point_data["u"], axis=1)
    return b


def assemble(folder, prefix, parts=("bottom", "top", "interface")):
    got = []
    for p in parts:
        m = read_body(Path(folder) / f"{prefix}_{p}_0.01.case")
        if m is not None:
            got.append(m)
    return got


def warp(meshes, factor):
    out = []
    for m in meshes:
        w = m.copy()
        w.points = w.points + factor * np.asarray(w.point_data["u"])
        out.append(w.extract_surface())
    return out


def render(meshes, clim, cam=None, size=(1400, 1500), zoom=None, focus=None):
    """Off-screen render.

    A zoom on the interface must PAN, not rotate: moving only the focal point
    swings the camera about its position and ends up looking through the
    specimen. Position and focal point are therefore translated together, which
    preserves the view direction, and only then is the zoom applied.
    """
    # The default light kit washes the faces out: coolwarm saturates at both ends,
    # so with flat lighting the block reads as two solid patches and its shape is
    # lost. A key light off the camera axis with a weaker fill from the opposite
    # side puts a gradient across each face, and a little specular picks out the
    # edges between them.
    p = pv.Plotter(off_screen=True, window_size=list(size), lighting="none")
    p.set_background("white")
    for m in meshes:
        p.add_mesh(m, scalars="umag", cmap=CMAP, clim=clim,
                   show_edges=True, edge_color="#2b2b2b", line_width=0.6,
                   show_scalar_bar=False, smooth_shading=False,
                   ambient=0.16, diffuse=0.88, specular=0.26, specular_power=18)
    key = pv.Light(light_type="camera light", intensity=1.05)
    key.position = (0.55, 0.35, 1.0)
    fill = pv.Light(light_type="camera light", intensity=0.30)
    fill.position = (-0.75, -0.25, 0.55)
    back = pv.Light(light_type="camera light", intensity=0.20)
    back.position = (0.0, 0.9, -0.6)
    for L in (key, fill, back):
        p.add_light(L)
    p.enable_anti_aliasing("ssaa")
    if cam is None:
        p.view_vector((0.55, -0.95, 0.42), viewup=(0, 0, 1))
        p.camera.zoom(1.28)
    else:
        # camera_position carries only position/focal/up; the zoom lives in the
        # view angle (or the parallel scale), so restoring it alone would let
        # the second render re-fit itself and come out at a different size
        p.camera_position = cam["pos"]
        p.camera.view_angle = cam["angle"]
        p.camera.parallel_projection = cam["par"]
        if cam["par"]:
            p.camera.parallel_scale = cam["scale"]
    if focus is not None:
        pos = np.asarray(p.camera.position)
        fp = np.asarray(p.camera.focal_point)
        tgt = np.asarray(focus, float)
        p.camera.position = tuple(pos + (tgt - fp))     # translate, keep direction
        p.camera.focal_point = tuple(tgt)
    if zoom is not None:
        p.camera.zoom(zoom)
    # after panning/zooming the near plane can cut into the specimen, which
    # renders the interior as a flat colour; recompute it
    p.renderer.reset_camera_clipping_range()

    img = np.asarray(p.screenshot(return_img=True))

    proj = None
    if focus is not None:
        ren = p.renderer
        ren.SetWorldPoint(float(focus[0]), float(focus[1]), float(focus[2]), 1.0)
        ren.WorldToDisplay()
        dx, dy, _ = ren.GetDisplayPoint()
        sy, sx = img.shape[0], img.shape[1]
        proj = (dx * sx / size[0], sy - dy * sy / size[1])
    cam_out = dict(pos=p.camera_position, angle=p.camera.view_angle,
                   par=bool(p.camera.parallel_projection),
                   scale=p.camera.parallel_scale)
    p.close()
    return img, cam_out, proj


def project(meshes, clim, cam, world, size=(1400, 1500)):
    """Screen position of a world point in the rendered view."""
    pl = pv.Plotter(off_screen=True, window_size=list(size))
    pl.add_mesh(meshes[0], show_scalar_bar=False)
    pl.camera_position = cam["pos"]
    pl.camera.view_angle = cam["angle"]
    pl.camera.parallel_projection = cam["par"]
    if cam["par"]:
        pl.camera.parallel_scale = cam["scale"]
    pl.render()
    ren = pl.renderer
    ren.SetWorldPoint(float(world[0]), float(world[1]), float(world[2]), 1.0)
    ren.WorldToDisplay()
    dx, dy, _ = ren.GetDisplayPoint()
    pl.close()
    return (dx, size[1] - dy)


def interface_focus(meshes):
    """A point on the interface.

    When the interphase is present as its own body it is the thinnest of the
    three, so its centre is the natural target; otherwise the interface sits at
    the top of the lowest body.
    """
    if len(meshes) >= 3:
        thin = min(meshes, key=lambda m: m.bounds[5] - m.bounds[4])
        return tuple(thin.center)
    bot = min(meshes, key=lambda m: m.center[2])
    b = bot.bounds
    return ((b[0] + b[1]) / 2, (b[2] + b[3]) / 2, b[5])


def circle_crop(img, frac=0.46):
    h, w = img.shape[:2]
    s = int(min(h, w) * frac)
    y0, x0 = (h - s) // 2, (w - s) // 2
    return img[y0:y0 + s, x0:x0 + s]


def find_transition(img, x_hint, half):
    """Row of the blue-to-red crossing, i.e. the interface, down a column.

    Searching the whole column finds the silhouette against the white ground,
    whose gradient is larger than the interface itself. The search is therefore
    restricted to the middle of the specimen, where the interface must lie, and
    keys on the sign change of red-minus-blue, which for a diverging map is the
    crossing rather than merely a steep edge.
    """
    x = int(np.clip(x_hint, half, img.shape[1] - half - 1))
    col = img[:, max(x - 40, 0):x + 40, :3].astype(float).mean(axis=1)
    solid = np.where(col.mean(axis=1) < 250)[0]
    if len(solid) < 40:
        return img.shape[0] // 2
    lo = int(solid.min() + 0.25 * (solid.max() - solid.min()))
    hi = int(solid.min() + 0.75 * (solid.max() - solid.min()))
    rb = col[lo:hi, 0] - col[lo:hi, 2]
    sign = np.sign(rb - 0.5 * (rb.max() + rb.min()))
    cross = np.where(np.diff(sign) != 0)[0]
    if len(cross):
        return int(lo + cross[len(cross) // 2])
    return int(lo + np.argmax(np.abs(np.gradient(rb))))


def crop_at(img, xy, half):
    """Square crop centred on a projected point, clamped to the image."""
    h, w = img.shape[:2]
    cx, cy = int(round(xy[0])), int(round(xy[1]))
    x0 = int(np.clip(cx - half, 0, max(w - 2 * half, 0)))
    y0 = int(np.clip(cy - half, 0, max(h - 2 * half, 0)))
    return img[y0:y0 + 2 * half, x0:x0 + 2 * half]


HI = (4200, 4500)          # rendered once, large enough that a tight crop
                           # still resolves the five elements through the layer

for tag, cfg in CASES.items():
    full = assemble(cfg["full"], cfg["fpre"])
    iface = assemble(cfg["iface"], cfg["ipre"])
    if not full or not iface:
        print(f"  {tag}: missing outputs, skipped")
        continue

    umax = max(float(np.max(m.point_data["umag"])) for m in full + iface)
    span = float(np.linalg.norm(np.ptp(np.vstack([m.points for m in full]), axis=0)))
    fac = 0.05 * span / max(umax, 1e-30)
    # The companion panels report displacement relative to a reference value, so the
    # deformed shape is coloured by the same kind of quantity: the magnitude divided by
    # its own maximum over the two models. The shared scale is then 0 to 1 by
    # construction and the two renders remain directly comparable.
    for m in full + iface:
        m.point_data["umag"] = m.point_data["umag"] / max(umax, 1e-30)
    clim = (0.0, 1.0)

    wf, wi = warp(full, fac), warp(iface, fac)
    img_i, cam, _ = render(wi, clim, size=HI)
    img_f, _, _ = render(wf, clim, cam=cam, size=HI)

    # The magnified view is a crop of the same render rather than a second
    # camera: no panning, no clipping surprises, and the two panels are
    # guaranteed to show the same region at the same scale.
    pf = project(wf, clim, cam, interface_focus(wf), HI)
    pi = project(wi, clim, cam, interface_focus(wi), HI)
    HALF = int(0.036 * HI[0])
    pf = (pf[0], find_transition(img_f, pf[0], HALF))
    # the interface model: the colour crossing that find_transition keys on lies inside its bottom body, and the projected
    # centre of the interface part lies behind the front face; both renders share camera, size and magnification, and the
    # bodies sit at the same place, so the interface point found on the resolved layer is used for both panels
    pi = pf
    zf, zi = crop_at(img_f, pf, HALF), crop_at(img_i, pi, HALF)

    fig = plt.figure(figsize=(6.2, 2.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.05], wspace=0.02)
    for i, (img, zimg, xy, ttl) in enumerate(
            ((img_f, zf, pf, r"resolved Cauchy layer"),
             (img_i, zi, pi, r"interface model"))):
        ax = fig.add_subplot(gs[0, i])
        ax.imshow(img)
        ax.set_axis_off()
        ax.set_title(ttl, fontsize=9, pad=3, loc="left")
        ax.add_patch(Circle(xy, HALF, fill=False, lw=1.0, ec="#111111"))

        ins = ax.inset_axes([0.58, 0.03, 0.44, 0.44])
        ins.imshow(zimg)
        ins.set_xticks([]); ins.set_yticks([])
        ins.patch.set_alpha(0)
        for s in ins.spines.values():
            s.set_visible(False)
        R = zimg.shape[0] / 2
        ins.images[0].set_clip_path(Circle((zimg.shape[1] / 2, R), R,
                                           transform=ins.transData))
        ins.add_patch(Circle((zimg.shape[1] / 2, R), R, fill=False,
                             lw=1.0, ec="#111111"))
        for ang in (1.0, -0.45):
            ax.add_artist(ConnectionPatch(
                xyA=(xy[0] + HALF * np.cos(ang), xy[1] - HALF * np.sin(ang)),
                xyB=(zimg.shape[1] / 2 + R * np.cos(ang), R - R * np.sin(ang)),
                coordsA="data", coordsB="data", axesA=ax, axesB=ins,
                lw=0.7, color="#111111"))

    cax = fig.add_subplot(gs[0, 2])
    cb = fig.colorbar(cm.ScalarMappable(norm=colors.Normalize(*clim), cmap=CMAP),
                      cax=cax)
    cb.set_label(r"$\|u_i\|\,/\,\max\|u_i\|$", labelpad=5)
    cb.outline.set_linewidth(0.6)
    fig.savefig(OUT / f"deformed_{tag}.pdf")
    plt.close(fig)

    for nm, z in (("(a) zoom", zf), ("(b) zoom", zi)):
        g = z[..., :3].astype(float).mean(axis=2)
        print(f"      {nm}: std={g.std():6.2f}  {'FLAT' if g.std() < 6 else 'ok'}")
    print(f"  {tag}: umax={umax:.4e} mm, magnification x{fac:.0f} -> deformed_{tag}.pdf")
