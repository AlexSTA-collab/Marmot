#!/usr/bin/env python3
"""End-of-analysis errors of the interface model against the resolved layer, for one discretization level.

The three measures of the paper (port_2611/scripts/node_matched_errors.py, whose function bodies are copied below
unchanged), read from 26.11 output (one file per step) and evaluated at the LAST step:
  (a) e_u    ||u^h - u^0||_L2 / ||u^0||_L2 over both bodies (twin nodes subtracted)
  (b) e_[u]  W ||[u^h] - [u^0]|| / ||u^0||_L2, jump between the two body faces bounding the layer, edge margin 0.15 mm.
             W = sqrt(0.5) / N turns the sum over the N x N face grid into the area integral; the paper's
             W_AREA = sqrt(1/30**2 * 0.5) is this at N = 30.
  (c) e_t    ||t^h - t^0|| / ||t^0||, the element force of the interface model against sigma.n of the layer averaged
             through each thickness column, paired by in-plane position, edge margin 0.15 mm.

usage: errors_end.py <model run> <model suffix> <reference run> <reference suffix> <N cells per side> [out.json]
  the model's case directories are Cohesive_zone_model_<side>_<suffix>, the reference's Cauchy_Full_model_<side>_<suffix>
"""
import json
import os
import re
import sys
from glob import glob
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "port_2611"))
from ensight_reader_2611 import NPE, read_geo, read_var_node_all_steps  # noqa: E402
from trace_mapping import interpolate_top_trace_along_normal  # noqa: E402

ANGLE = 10
H0 = float(os.environ.get("GAP", 0.01))  # layer height h; the paper study is h = 0.01
MARGIN = 0.15
TWIN_TOL = 1e-4


def normal(angle):
    a = np.radians(float(angle))
    return np.array([np.sin(a), 0.0, np.cos(a)])


def matched(xr, xp):
    """Index array taking the interface-model points onto the reference points, one to one."""
    d, idx = cKDTree(xp).query(xr)
    if d.max() > TWIN_TOL or len(np.unique(idx)) != len(idx):
        raise SystemExit(f"no one-to-one twin within {TWIN_TOL:g} mm (worst {d.max():.1e})")
    return idx


def centroids(parts, pn):
    xyz = parts[pn]["xyz"]
    return np.vstack([xyz[c].mean(axis=1) for _, c in parts[pn]["blocks"]])


def normal(angle):
    a = np.radians(float(angle))
    return np.array([np.sin(a), 0.0, np.cos(a)])


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


def by_position(c, v):
    order = np.lexsort((np.round(c[:, 1], 6), np.round(c[:, 0], 6)))
    return c[order], v[order]


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



def jump(bottom, top, h, angle, on=None):
    """([u] at the lower-face nodes that have a normal correspondent, those nodes); `on` = the reference's faces."""
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


def casedir(run, prefix, side, suffix):
    return Path(run) / f"{prefix}_{side}_{suffix}"


def last(d, name):
    """The file of the last step, or of the step nearest ERR_T (environment variable) when that is set."""
    files = sorted(glob(str(d / f"{name}.var_[0-9][0-9][0-9][0-9]")))
    tstar = os.environ.get("ERR_T")
    if tstar is None:
        return files[-1]
    case = Path(str(d) + ".case").read_text()
    t = np.array(re.findall(r"[-+0-9.e]+", case.split("time values:")[-1].split("GEOMETRY")[0])[: len(files)], float)
    return files[int(np.argmin(np.abs(t - float(tstar))))]


def body(run, prefix, side, suffix):
    """Coordinates, hexahedra and last-step displacement of one body."""
    d = casedir(run, prefix, side, suffix)
    parts = read_geo(str(d / "geometry.geo_0000"))
    u = read_var_node_all_steps(last(d, f"displacement_{side}"), parts, 3)
    pn = list(u)[0]
    return parts[pn]["xyz"], parts[pn]["blocks"][0][1], u[pn][-1]


def read_elem_var(path, parts, ncomp):
    """{part: (nelem, ncomp)} of a per-element variable of one 26.11 step file."""
    d = open(path, "rb").read()
    i = 80  # description
    out = {}
    while d[i:i + 4] == b"part":
        pn = int(np.frombuffer(d[i + 80:i + 84], "<i4")[0])
        i += 84
        rows = []
        while i < len(d) and d[i:i + 80].split(b"\x00")[0].decode().strip() in NPE:
            typ = d[i:i + 80].split(b"\x00")[0].decode().strip()
            ne = len(next(c for t, c in parts[pn]["blocks"] if t == typ))
            rows.append(np.frombuffer(d[i + 80:i + 80 + 4 * ne * ncomp], "<f4").astype(float).reshape(ncomp, ne).T)
            i += 80 + 4 * ne * ncomp
        out[pn] = np.vstack(rows)
    return out


def interface_part(run, prefix, suffix):
    d = casedir(run, prefix, "interface", suffix)
    return d, read_geo(str(d / "geometry.geo_0000"))


def model_traction(run, suffix):
    d, parts = interface_part(run, "Cohesive_zone_model", suffix)
    f = read_elem_var(last(d, "force_interface"), parts, 3)
    pn = list(f)[0]
    return by_position(centroids(parts, pn), f[pn])


def layer_traction(run, suffix, angle):
    """sigma.n of every layer element, averaged over each through-thickness column (Voigt order 11 22 33 12 13 23)."""
    d, parts = interface_part(run, "Cauchy_Full_model", suffix)
    s = read_elem_var(last(d, "stress_interface"), parts, 6)
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


def measures(mrun, msuffix, rrun, rsuffix, ncells):
    ref = {s: body(rrun, "Cauchy_Full_model", s, rsuffix) for s in ("bottom", "top")}
    mod = {s: body(mrun, "Cohesive_zone_model", s, msuffix) for s in ("bottom", "top")}
    num = den = 0.0
    for s in ("bottom", "top"):
        (xr, cr, ur), (xm, cm_, um) = ref[s], mod[s]
        um_at_ref = evaluate_at(xm, cm_, um, xr, matched(xr, xm))
        num += integrate(np.sum((ur - um_at_ref) ** 2, axis=1), xr, cr)
        den += integrate(np.sum(ur ** 2, axis=1), xr, cr)
    out = {"e_u_global": float(np.sqrt(num / den))}

    jr, pr = jump(ref["bottom"], ref["top"], H0, ANGLE)
    jm, pm = jump(mod["bottom"], mod["top"], H0, ANGLE, on=(ref["bottom"], ref["top"]))
    assert np.array_equal(pr, pm), "face traces differ"
    k = interior(pr, ANGLE, MARGIN)
    w_area = np.sqrt(0.5) / ncells
    out["e_jump"] = float(w_area * np.sqrt(np.sum((jm - jr)[k] ** 2)) / np.sqrt(den))

    cm, tm = model_traction(mrun, msuffix)
    cr_, tr = layer_traction(rrun, rsuffix, ANGLE)
    assert len(cm) == len(cr_), f"{len(cm)} interface cells, {len(cr_)} columns"
    tm = tm[matched(cr_[:, :2], cm[:, :2])]
    k = interior(cr_, ANGLE, MARGIN)
    out["e_traction"] = float(np.linalg.norm((tm - tr)[k]) / np.linalg.norm(tr[k]))
    out["n_face_nodes"] = int(k.sum())
    return out


if __name__ == "__main__":
    mrun, msuffix, rrun, rsuffix, ncells = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
    res = measures(mrun, msuffix, rrun, rsuffix, ncells)
    print(f"model {msuffix} vs reference {rsuffix} (N = {ncells}): " + "  ".join(f"{k}={v:.4e}" if isinstance(v, float) else f"{k}={v}" for k, v in res.items()))
    if len(sys.argv) > 6:
        Path(sys.argv[6]).write_text(json.dumps(res))
