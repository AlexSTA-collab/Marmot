#!/usr/bin/env python3
"""The three error measures of the convergence figures, read straight off the stored results.

The interface model and the resolved layer of the same thickness carry the same meshes for
the two bodies -- the interface model is the resolved one with the layer replaced. So every
body node of one model has a twin at the same coordinates in the other -- exactly for the
viscoelastic studies, within 7e-6 mm for the elastic and J_2 ones, whose two meshes came from
separate Cubit runs. The interface model's displacement is evaluated at the reference's nodes
with its own trilinear shape functions, which is the nodal value itself where the meshes
coincide, and the two are subtracted node by node. The study's ParaView pipeline instead
resamples one model's field onto the other mesh, locating each node in a cell to a 1e-5
tolerance; where the model error is small (alpha = 0) that interpolation is all it measures.

At the last step of each analysis:

  (a) e_u    ||u^h - u^0||_L2 / ||u^0||_L2 over both bodies. The squared nodal difference is
             averaged onto each hexahedron and multiplied by its volume, as
             Evaluate_relative_error.py integrates it, but with twin nodes subtracted.
  (b) e_[u]  W_AREA ||[u^h] - [u^0]|| / ||u^0||_L2, the jump between the two body faces that
             bound the layer. The upper face is read at the normal correspondent of each
             lower-face node by the study's own trace_mapping.py -- the identity at
             alpha = 0, a bilinear read on the face grid at alpha = 10.
  (c) e_t    ||t^h - t^0|| / ||t^0||: the element's generalizedForce against the resolved
             layer's sigma.n averaged through the thickness, paired by in-plane position,
             as verification/traction_direct/evaluate_traction.py pairs them.

The edge margin (0.15 mm) and W_AREA are those of recompute_jump_errors.py.

--fixed compares every model with the single resolved layer of thickness h0 = 0.01 instead
(Appendix G). The top body of that reference sits h - h0 lower than the model's, so the two
meshes no longer share their nodes: the model's displacement is evaluated at each reference
node by locating the node in the model's elements and inverting their trilinear map, and nodes
that lie outside the model's body (the strip of height h - h0 below its top body) are left out
of both integrals of (a), as the ParaView resampling leaves out the points it cannot locate.
(b) compares each model's own jump with the reference's, paired by the lower-face node, which
the two share; (c) pairs the cells by their sorted in-plane position, as
recompute_traction_errors.py does. At h = h0 the two protocols compare the same pair of
solutions and return the same numbers.

    python3 verification/node_matched_errors.py <case> <angle> <table.md> [--fixed]
    case: elastic | vm_eie | vm_psie | visco
"""
import os
import sys
from pathlib import Path

import numpy as np

NF = Path("/home/alexsta1993/Documents/STATHAS_VIENNA_POSTDOC/Applications_Examples/"
          "EdelweisModels/Cohesive_Zone_models/3D/Three_Body_Cubit/"
          "Comparizon_LinearElastic/New_Functionality")
sys.path.insert(0, str(NF / "Parametric_Study_Marmot-nextv26.05_Von_Mises_StabPMini"
                       / "evaluation_scripts"))
sys.path.insert(0, str(NF / "Parametric_Study_Marmot-nextv26.05_PowerLaw_Warping"
                       / "evaluation_scripts"))
from ensight_reader import read_geo, read_var_node_all_steps, NPE      # noqa: E402
from trace_mapping import interpolate_top_trace_along_normal           # noqa: E402
from scipy.spatial import cKDTree                                      # noqa: E402

DATA = Path(os.environ.get("RESULTS_ROOT", str(NF)))   # port_2611: where the runs to be scored live
EL = DATA / "Parametric_Study_Marmot-nextv26.05_LinearElasticExtended"
VME = DATA / "Parametric_Study_Marmot-nextv26.05_Von_Mises_ExtendedInterfaceElement"
VMP = DATA / "Parametric_Study_Marmot-nextv26.05_Von_Mises_StabPMini"
VE = DATA / "Parametric_Study_Marmot-nextv26.05_PowerLaw_Transition"
VEL = DATA / "Parametric_Study_Marmot-nextv26.05_PowerLaw_Transition_LoadScaled"

# Constants of recompute_jump_errors.py, copied: importing that module runs it, and it
# rewrites verification/convergence_table_angle{10,0}.md.
W_AREA = np.sqrt(1.0 / 30.0 ** 2 * 0.5)
MARGIN = 0.15
HEIGHTS = ["0.01", "0.02", "0.04", "0.06", "0.08", "0.1"]
FAMILIES = ["stiff", "initial", "soft"]
ELEMENT = {"elastic": "YIQUAD4", "vm_eie": "YIQUAD4", "vm_psie": "YIQUAD4_STABP_MINI",
           "visco": "YIQUAD4"}


def model_dir(case, fam, angle, h):
    if case == "elastic":
        return EL / "YIQUAD4_runs" / fam / f"angle_{angle}" / f"h_{h}"
    if case == "vm_eie":
        return VME / "YIQUAD4_runs" / fam / f"angle_{angle}" / f"h_{h}"
    if case == "vm_psie":
        return VMP / fam / f"angle_{angle}" / "Interface_models"
    if angle == 10:
        return VE / "proper_convergence" / fam / "YIQUAD4" / f"h_{h}"
    # alpha = 0 is the load-scaled rerun, with the interface traction written
    return VEL / "proper_convergence_angle0_traction" / fam / "YIQUAD4" / f"h_{h}"


def ref_dir(case, fam, angle, h):
    if case == "visco":
        study = VE if angle == 10 else VEL
        if h == "0.01":
            return study / fam / f"angle_{angle}" / "Full_layer_models"
        return (study / "proper_heff_check" / fam / f"h_{h}" if angle == 10
                else study / "proper_convergence_angle0" / fam / f"h_{h}")
    study = EL if case == "elastic" else VME        # the two J_2 studies share one reference
    if h == "0.01":
        return study / fam / f"angle_{angle}" / "Full_layer_models"
    return study / "matched_reference" / fam / f"angle_{angle}" / f"h_{h}"


def interior(coords, angle, margin):
    """Points further than `margin` from every edge, in the interface's own tangential
    coordinates -- recompute_jump_errors.interior, verbatim."""
    if margin <= 0.0:
        return np.ones(len(coords), dtype=bool)
    a = np.radians(float(angle))
    t1 = np.array([np.cos(a), 0.0, -np.sin(a)])
    t2 = np.array([0.0, 1.0, 0.0])
    keep = np.ones(len(coords), dtype=bool)
    for t in (t1, t2):
        s = coords @ t
        s = s - 0.5 * (s.min() + s.max())
        keep &= np.abs(s) <= np.abs(s).max() - margin
    return keep


# Twins further apart than this are refused. The viscoelastic models share their body
# meshes exactly; the elastic and J_2 references were meshed by separate Cubit runs and
# their nodes sit 5e-6 to 7e-6 mm from the interface model's, against an element size of
# 3.3e-2 mm, so the nearest node is unambiguous.
TWIN_TOL = 1e-4


def matched(xr, xp):
    """Index array taking the interface-model points onto the reference points, one to one."""
    d, idx = cKDTree(xp).query(xr)
    if d.max() > TWIN_TOL or len(np.unique(idx)) != len(idx):
        raise SystemExit(f"no one-to-one twin within {TWIN_TOL:g} mm (worst {d.max():.1e})")
    return idx


# --- reading -------------------------------------------------------------------------
def output(run, side, h):
    d = Path(run) / next(p.name for p in Path(run).iterdir()
                         if p.is_dir() and p.name.endswith(f"_{side}_{h}"))
    return d, read_geo(str(d / "geometry.geo"))


def read_elem_var(path, parts, ncomp):
    """{part: (nelem, ncomp)} of a per-element EnSight Gold variable at its last step."""
    d = open(path, "rb").read()
    i = d.rfind(b"BEGIN TIME STEP") + 160             # marker, description
    out = {}
    while d[i:i + 4] == b"part":
        pn = int(np.frombuffer(d[i + 80:i + 84], "<i4")[0])
        i += 84
        rows = []
        while d[i:i + 80].split(b"\x00")[0].decode().strip() in NPE:
            typ = d[i:i + 80].split(b"\x00")[0].decode().strip()
            ne = len(next(c for t, c in parts[pn]["blocks"] if t == typ))
            rows.append(np.frombuffer(d[i + 80:i + 80 + 4 * ne * ncomp], "<f4")
                        .astype(float).reshape(ncomp, ne).T)
            i += 80 + 4 * ne * ncomp
        out[pn] = np.vstack(rows)
    return out


def body(run, side, h):
    """Coordinates, hexahedra and last-step displacement of one body."""
    d, parts = output(run, side, h)
    u = read_var_node_all_steps(str(d / f"displacement_{side}.var"), parts, 3)
    pn = list(u)[0]
    return parts[pn]["xyz"], parts[pn]["blocks"][0][1], u[pn][-1]


def centroids(parts, pn):
    xyz = parts[pn]["xyz"]
    return np.vstack([xyz[c].mean(axis=1) for _, c in parts[pn]["blocks"]])


def normal(angle):
    a = np.radians(float(angle))
    return np.array([np.sin(a), 0.0, np.cos(a)])


# --- evaluating one model's field at the other model's points ------------------------
def evaluate_at(xyz, conn, u, points, twin):
    """u of a trilinear hexahedral field at `points`, by inverting the element map.

    `twin[i]` is the node of this mesh nearest to points[i]; the candidates are the
    elements sharing that node, and the one whose natural coordinates lie inside [-1, 1]^3
    (or least outside, for a point on the boundary 5e-6 mm off the surface) is used.
    Where the meshes coincide the natural coordinates land on a vertex and this returns
    the nodal value itself.
    """
    node_elems = [[] for _ in range(len(xyz))]
    for e, c in enumerate(conn):
        for n in c:
            node_elems[n].append(e)
    pi, ei = [], []
    for i, n in enumerate(twin):
        for e in node_elems[n]:
            pi.append(i); ei.append(e)
    values, worst = _invert(xyz, conn, u, points, np.array(pi), np.array(ei))
    if worst.max() > 1e-3:
        raise SystemExit(f"a point lies outside every candidate element (by {worst.max():.1e})")
    return values


def locate(xyz, conn, u, points, candidates=12):
    """u at arbitrary `points`, and a mask of the points inside this mesh.

    The candidates for each point are the elements with the nearest centroids."""
    cen = xyz[conn].mean(axis=1)
    _, near = cKDTree(cen).query(points, k=candidates)
    pi = np.repeat(np.arange(len(points)), candidates)
    values, worst = _invert(xyz, conn, u, points, pi, near.ravel())
    return values, worst <= 1e-3


def _invert(xyz, conn, u, points, pi, ei):
    """Evaluate u at points[pi] in elements ei; keep, per point, the least-outside element.
    Returns the values and, per point, how far outside [-1, 1]^3 that element's coordinates lie."""
    X, P = xyz[conn[ei]], points[pi]                           # (k, 8, 3), (k, 3)
    xi = np.zeros((len(pi), 3))
    for _ in range(8):
        N = np.prod(1 + _XI[None] * xi[:, None, :], axis=2) / 8.0             # (k, 8)
        dN = np.stack([_XI[None, :, k] * np.prod([1 + _XI[None, :, m] * xi[:, None, m]
                                                  for m in range(3) if m != k], axis=0) / 8.0
                       for k in range(3)], axis=2)                          # (k, 8, 3)
        r = np.einsum("kn,knl->kl", N, X) - P
        J = np.einsum("knl,knm->klm", X, dN)
        # a far candidate can drive the iterate into a degenerate map; hold it there, since
        # it is never the element a point lies in
        ok = np.abs(np.linalg.det(J)) > 1e-30
        J[~ok] = np.eye(3)
        step = np.linalg.solve(J, r[..., None])[..., 0]
        step[~ok] = 0.0
        xi = np.clip(xi - step, -3.0, 3.0)
    outside = np.abs(xi).max(axis=1) - 1.0
    outside[~np.isfinite(outside)] = np.inf
    order = np.lexsort((outside, pi))                          # per point, least outside first
    first = np.ones(len(order), bool)
    first[1:] = pi[order][1:] != pi[order][:-1]
    k = order[first]                                           # one row per point, in point order
    worst = np.full(len(points), np.inf)
    worst[pi[k]] = outside[k]
    N = np.prod(1 + _XI[None] * xi[k][:, None, :], axis=2) / 8.0
    values = np.zeros((len(points), u.shape[1]))
    values[pi[k]] = np.einsum("kn,knc->kc", N, u[conn[ei[k]]])
    return values, worst


# --- (a) -------------------------------------------------------------------------------
_XI = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)


def volumes(xyz, conn):
    """Hexahedron volumes, 2x2x2 Gauss."""
    X = xyz[conn]
    vol = np.zeros(len(conn))
    for g in _XI / np.sqrt(3.0):
        dN = np.stack([_XI[:, k] * np.prod([1 + _XI[:, m] * g[m] for m in range(3) if m != k],
                                           axis=0) / 8.0 for k in range(3)], axis=1)
        vol += np.abs(np.linalg.det(np.einsum("enk,nl->ekl", X, dN)))
    return vol


def integrate(sq, xyz, conn):
    """Integral of a nodal field: cell average of the nodal values times the cell volume."""
    return float(np.sum(sq[conn].mean(axis=1) * volumes(xyz, conn)))


# --- (b) -------------------------------------------------------------------------------
def face(xyz, point, n, tol=1e-5):
    return np.where(np.abs((xyz - point) @ n) < tol)[0]


def jump(bottom, top, h, angle, on=None):
    """([u] at the lower-face nodes that have a normal correspondent, those nodes).

    With `on` = (reference bottom, reference top), the displacement is this model's but the
    face nodes are the reference's: each body's field is evaluated there first."""
    (xb, cb, ub), (xt, ct, ut) = bottom, top
    n, H = normal(angle), float(h)
    if on is not None:
        (xrb, _, _), (xrt, _, _) = on
        ub = evaluate_at(xb, cb, ub, xrb, matched(xrb, xb))
        ut = evaluate_at(xt, ct, ut, xrt, matched(xrt, xt))
        xb, xt = xrb, xrt
    ib = face(xb, np.zeros(3), n)
    it = face(xt, np.array([0.0, 0.0, H]), n)
    ib = ib[np.lexsort((xb[ib, 1], xb[ib, 0]))]
    valid, _, ut_mapped, _ = interpolate_top_trace_along_normal(xb[ib], xt[it], ut[it], H, angle)
    return ut_mapped - ub[ib][valid], xb[ib][valid]


# --- (c) -------------------------------------------------------------------------------
def by_position(c, v):
    order = np.lexsort((np.round(c[:, 1], 6), np.round(c[:, 0], 6)))
    return c[order], v[order]


def model_traction(run, h):
    d, parts = output(run, "interface", h)
    f = read_elem_var(str(d / "force_interface.var"), parts, 3)
    pn = list(f)[0]
    return by_position(centroids(parts, pn), f[pn])


def layer_traction(run, h, angle):
    """sigma.n of every layer element, averaged over each through-thickness column.
    Marmot writes the stress in Voigt order 11 22 33 12 13 23."""
    d, parts = output(run, "interface", h)
    s = read_elem_var(str(d / "stress_interface.var"), parts, 6)
    pn = list(s)[0]
    s, c, n = s[pn], centroids(parts, pn), normal(angle)
    t = np.stack([s[:, 0] * n[0] + s[:, 3] * n[1] + s[:, 4] * n[2],
                  s[:, 3] * n[0] + s[:, 1] * n[1] + s[:, 5] * n[2],
                  s[:, 4] * n[0] + s[:, 5] * n[1] + s[:, 2] * n[2]], axis=1)
    key = np.stack([np.round(c[:, 0], 5), np.round(c[:, 1], 5)], axis=1)
    uniq, inv = np.unique(key, axis=0, return_inverse=True)
    inv = inv.ravel()
    cc = np.array([c[inv == i].mean(axis=0) for i in range(len(uniq))])
    tt = np.array([t[inv == i].mean(axis=0) for i in range(len(uniq))])
    return by_position(cc, tt)


# --- the three measures ----------------------------------------------------------------
def measures(case, fam, angle, h, fixed=False):
    # at h = h0 the single reference is the matched one: same pair of solutions, same path
    fixed = fixed and h != "0.01"
    h0 = "0.01" if fixed else h
    mdir, rdir = model_dir(case, fam, angle, h), ref_dir(case, fam, angle, h0)
    ref = {s: body(rdir, s, h0) for s in ("bottom", "top")}
    mod = {s: body(mdir, s, h) for s in ("bottom", "top")}

    num = den = den_all = 0.0
    for s in ("bottom", "top"):
        (xr, cr, ur), (xm, cm_, um) = ref[s], mod[s]
        if fixed:
            um_at_ref, inside = locate(xm, cm_, um, xr)
            w = inside.astype(float)
        else:
            um_at_ref, w = evaluate_at(xm, cm_, um, xr, matched(xr, xm)), 1.0
        num += integrate(np.sum((ur - um_at_ref) ** 2, axis=1) * w, xr, cr)
        den += integrate(np.sum(ur ** 2, axis=1) * w, xr, cr)
        den_all += integrate(np.sum(ur ** 2, axis=1), xr, cr)
    out = {"global": np.sqrt(num / den)}

    jr, pr = jump(ref["bottom"], ref["top"], h0, angle)
    if fixed:
        jm, pm = jump(mod["bottom"], mod["top"], h, angle)
        d, idx = cKDTree(pm).query(pr)
        both = d < TWIN_TOL
        jr, pr, jm = jr[both], pr[both], jm[idx[both]]
    else:
        jm, pm = jump(mod["bottom"], mod["top"], h, angle, on=(ref["bottom"], ref["top"]))
        if not np.array_equal(pr, pm):
            raise SystemExit(f"{case} {fam} {angle} {h}: face traces differ")
    k = interior(pr, angle, MARGIN)
    out["e_u"] = W_AREA * np.sqrt(np.sum((jm - jr)[k] ** 2)) / np.sqrt(den_all)

    cm, tm = model_traction(mdir, h)
    cr_, tr = layer_traction(rdir, h0, angle)
    if len(cm) != len(cr_):
        raise SystemExit(f"{case} {fam} {angle} {h}: {len(cm)} interface cells, {len(cr_)} columns")
    if fixed:
        k = interior(cm, angle, MARGIN) & interior(cr_, angle, MARGIN)
    else:
        tm = tm[matched(cr_[:, :2], cm[:, :2])]
        k = interior(cr_, angle, MARGIN)
    out["e_t"] = np.linalg.norm((tm - tr)[k]) / np.linalg.norm(tr[k])
    return out


def rate(v):
    v = np.asarray(v, float)
    hs = np.array([float(h) for h in HEIGHTS])
    ok = np.isfinite(v) & (v > 0)
    return np.polyfit(np.log(hs[ok]), np.log(v[ok]), 1)[0] if ok.sum() > 2 else np.nan


def main():
    fixed = "--fixed" in sys.argv
    args = [a for a in sys.argv[1:] if a != "--fixed"]
    case, angle, out_path = args[0], int(args[1]), Path(args[2])
    res = {}
    for fam in FAMILIES:
        for h in HEIGHTS:
            res[(fam, h)] = m = measures(case, fam, angle, h, fixed)
            print(f"{case} a={angle} {fam:8s} h={h:5s} " +
                  "  ".join(f"{k}={v:.4e}" for k, v in m.items()))

    lines = [f"# {'Single-reference' if fixed else 'Matched-thickness'} convergence, angle {angle}, {case}",
             f"# reference: {'fixed' if fixed else 'matched'}, edge margin: {MARGIN:g}",
             "# node-matched: verification/node_matched_errors.py, no resampling", ""]
    for title, k in (("Displacement field (bodies)", "global"),
                     ("Displacement jump", "e_u"), ("Traction", "e_t")):
        lines += [f"## {title}", "",
                  "| h | " + " | ".join(f"{f} {ELEMENT[case]}" for f in FAMILIES) + " |",
                  "|---|" + "---|" * len(FAMILIES)]
        for h in HEIGHTS:
            lines.append(f"| {h} | " + " | ".join(f"{res[(f, h)][k]:.4e}" for f in FAMILIES) + " |")
        lines.append("| rate | " + " | ".join(f"{rate([res[(f, h)][k] for h in HEIGHTS]):.3f}"
                                             for f in FAMILIES) + " |")
        lines.append("")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines))
    print("wrote", out_path)


if __name__ == "__main__":
    main()
