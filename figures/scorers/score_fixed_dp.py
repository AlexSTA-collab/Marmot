#!/usr/bin/env python3
"""Single-reference errors (Appendix J) of the Drucker-Prager interphase: the EIE at every h against the ONE resolved
layer of thickness h0 = 0.01 (20 C3D8, 60 x 60), with the measures of plastic_bbar/score_fixed.py (copied below
unchanged), for the runs of Section 4.5: dp/R075_<family>_<psi0|dp>/{eie,ref}/angle_<a>/h_<h>. At h = h0 the two
protocols compare the same pair of solutions, so those entries are the matched scores of <family dir>/scores.json.
usage: score_fixed_dp.py <family dir> <out.json>      keys eie/<angle>/<h> = [e_u, e_[u], e_t]"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path.home() / "psie26/tools/refinement"))
import errors_end as E  # noqa: E402

root, out = Path(sys.argv[1]), Path(sys.argv[2])
matched = json.loads((root / "scores.json").read_text())
res = json.loads(out.read_text()) if out.exists() else {}
N, H0 = 60, "0.01"


def locate(xyz, conn, u, points, candidates=12):
    cen = xyz[conn].mean(axis=1)
    _, near = cKDTree(cen).query(points, k=candidates)
    pi = np.repeat(np.arange(len(points)), candidates)
    values, worst = E._invert(xyz, conn, u, points, pi, near.ravel())
    return values, worst <= 1e-3


def measures(mrun, h, rrun, angle):
    ref = {s: E.body(rrun, "Cauchy_Full_model", s, H0) for s in ("bottom", "top")}
    mod = {s: E.body(mrun, "Cohesive_zone_model", s, h) for s in ("bottom", "top")}
    num = den = den_all = 0.0
    for s in ("bottom", "top"):
        (xr, cr, ur), (xm, cm_, um) = ref[s], mod[s]
        um_at_ref, inside = locate(xm, cm_, um, xr)
        w = inside.astype(float)
        num += E.integrate(np.sum((ur - um_at_ref) ** 2, axis=1) * w, xr, cr)
        den += E.integrate(np.sum(ur ** 2, axis=1) * w, xr, cr)
        den_all += E.integrate(np.sum(ur ** 2, axis=1), xr, cr)
    o = [float(np.sqrt(num / den))]
    jr, pr = E.jump(ref["bottom"], ref["top"], H0, angle)
    jm, pm = E.jump(mod["bottom"], mod["top"], h, angle)
    d, idx = cKDTree(pm).query(pr)
    both = d < E.TWIN_TOL
    jr, pr, jm = jr[both], pr[both], jm[idx[both]]
    k = E.interior(pr, angle, E.MARGIN)
    o.append(float(np.sqrt(0.5) / N * np.sqrt(np.sum((jm - jr)[k] ** 2)) / np.sqrt(den_all)))
    cm, tm = E.model_traction(mrun, h)
    cr_, tr = E.layer_traction(rrun, H0, angle)
    assert len(cm) == len(cr_), f"{len(cm)} interface cells, {len(cr_)} columns"
    k = E.interior(cm, angle, E.MARGIN) & E.interior(cr_, angle, E.MARGIN)
    o.append(float(np.linalg.norm((tm - tr)[k]) / np.linalg.norm(tr[k])))
    return o


for a in (0, 10):
    rrun = root / "ref" / f"angle_{a}" / f"h_{H0}"
    for h in ("0.01", "0.02", "0.04", "0.06", "0.08", "0.1"):
        key = f"eie/{a}/{h}"
        if key in res:
            continue
        if h == H0:
            res[key] = matched[key]
        else:
            E.ANGLE, E.H0 = a, float(h)
            res[key] = measures(root / "eie" / f"angle_{a}" / f"h_{h}", h, rrun, a)
        print(root.name, key, " ".join(f"{v:.3e}" for v in res[key]), flush=True)
        out.write_text(json.dumps(res, indent=1, sort_keys=True))
