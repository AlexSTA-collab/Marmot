#!/usr/bin/env python3
"""Single-reference errors (Appendix G): EIE and PSIE at every h against the ONE resolved layer of thickness h0 = 0.01
(20 C3D8B through the thickness, 60 x 60), the --fixed path of paper/verification/node_matched_errors.py on the 26.11
output, with the readers and measures of errors_end.py:
  (a) e_u   the model's displacement is evaluated at each reference node by locating the node in the model's elements;
            nodes outside the model's bodies (the strip of height h - h0 below its top body) are left out of both
            integrals.
  (b) e_[u] each model's own jump against the reference's, paired by the lower-face node, which the two share;
            normalised by the reference displacement norm over the whole bodies; W = sqrt(0.5) / N.
  (c) e_t   cells paired by their sorted in-plane position, edge margin on both.
At h = h0 the two protocols compare the same pair of solutions: those entries are taken from the matched scores.
usage: score_fixed.py <jobs dir> <H> <matched scores.json> <out.json>"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path.home() / "psie26/tools/refinement"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "refinement"))
import errors_end as E  # noqa: E402

root, H, matched_scores, out = Path(sys.argv[1]), sys.argv[2], json.loads(Path(sys.argv[3]).read_text()), Path(sys.argv[4])
res = json.loads(out.read_text()) if out.exists() else {}
N, H0 = 60, "0.01"


def locate(xyz, conn, u, points, candidates=12):
    """u at arbitrary `points`, and a mask of the points inside this mesh (node_matched_errors.locate, verbatim)."""
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
    o = {"e_u_global": float(np.sqrt(num / den))}
    jr, pr = E.jump(ref["bottom"], ref["top"], H0, angle)
    jm, pm = E.jump(mod["bottom"], mod["top"], h, angle)
    d, idx = cKDTree(pm).query(pr)
    both = d < E.TWIN_TOL
    jr, pr, jm = jr[both], pr[both], jm[idx[both]]
    k = E.interior(pr, angle, E.MARGIN)
    o["e_jump"] = float(np.sqrt(0.5) / N * np.sqrt(np.sum((jm - jr)[k] ** 2)) / np.sqrt(den_all))
    cm, tm = E.model_traction(mrun, h)
    cr_, tr = E.layer_traction(rrun, H0, angle)
    assert len(cm) == len(cr_), f"{len(cm)} interface cells, {len(cr_)} columns"
    k = E.interior(cm, angle, E.MARGIN) & E.interior(cr_, angle, E.MARGIN)
    o["e_traction"] = float(np.linalg.norm((tm - tr)[k]) / np.linalg.norm(tr[k]))
    return o


for kind in ("psie", "eie"):
    for fam in ("stiff", "initial", "soft"):
        for a in (10, 0):
            rrun = root / "ref" / fam / f"angle_{a}" / f"H_{H}" / f"h_{H0}"
            for h in ("0.01", "0.02", "0.04", "0.06", "0.08", "0.1"):
                key = f"{kind}/{fam}/{a}/{H}/{h}"
                if key in res:
                    continue
                if h == H0:
                    res[key] = matched_scores[key]
                else:
                    E.ANGLE, E.H0 = a, float(h)
                    res[key] = measures(root / kind / fam / f"angle_{a}" / f"H_{H}" / f"h_{h}", h, rrun, a)
                print(key, {x: f"{res[key][x]:.3e}" for x in ("e_u_global", "e_jump", "e_traction")}, flush=True)
                out.write_text(json.dumps(res, indent=1, sort_keys=True))
