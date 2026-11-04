#!/usr/bin/env python3
"""Errors of the full model and of the model with [C^ev] = 0 against the resolved biphasic viscoelastic layer, at chosen times.

usage: score_visco.py <out.json> <case> <angle> <h> [times=1e-8,1,10,100] [margin=0.15]
Measures as in the paper (refinement/errors_end.py): e_u over both bodies, e_[u] the displacement-jump error, e_t the
traction, interior margin as given. The layer traction of the reference is sigma.n of both halves averaged through each
column (six elements, three per phase)."""
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "refinement")); sys.path.insert(0, str(HERE.parent / "port_2611"))
import errors_end as E  # noqa: E402

out, case, angle, h = Path(sys.argv[1]), sys.argv[2], int(sys.argv[3]), sys.argv[4]
times = [float(t) for t in (sys.argv[5] if len(sys.argv) > 5 else "1e-8,1,10,100").split(",")]
E.MARGIN = float(sys.argv[6]) if len(sys.argv) > 6 else 0.15
E.ANGLE, E.H0 = angle, float(h)
ROOT = HERE / "runs" / case
ref = ROOT / "ref" / f"angle_{angle}" / f"h_{h}"


def layer_traction(run):
    cs, ss = [], []
    for half, name in (("layerbottom", "stress_layer_bottom"), ("layertop", "stress_layer_top")):
        d = E.casedir(run, "Cauchy_Full_model", half, h)
        parts = E.read_geo(str(d / "geometry.geo_0000"))
        s = E.read_elem_var(E.last(d, name), parts, 6)
        pn = list(s)[0]
        cs.append(E.centroids(parts, pn)); ss.append(s[pn])
    s, c, n = np.vstack(ss), np.vstack(cs), E.normal(angle)
    t = np.stack([s[:, 0] * n[0] + s[:, 3] * n[1] + s[:, 4] * n[2],
                  s[:, 3] * n[0] + s[:, 1] * n[1] + s[:, 5] * n[2],
                  s[:, 4] * n[0] + s[:, 5] * n[1] + s[:, 2] * n[2]], axis=1)
    key = np.stack([np.round(c[:, 0], 5), np.round(c[:, 1], 5)], axis=1)
    uniq, inv = np.unique(key, axis=0, return_inverse=True)
    inv = inv.ravel()
    cc = np.array([c[inv == i].mean(axis=0) for i in range(len(uniq))])
    tt = np.array([t[inv == i].mean(axis=0) for i in range(len(uniq))])
    return E.by_position(cc, tt)


res = {}
for tstar in times:
    os.environ["ERR_T"] = repr(tstar)
    refb = {s: E.body(ref, "Cauchy_Full_model", s, h) for s in ("bottom", "top")}
    cr_, tr = layer_traction(ref)
    kt = E.interior(cr_, angle, E.MARGIN)
    res[tstar] = {"u_top_ref": float(np.mean(refb["top"][2][:, 0]))}
    for model in ("full", "supp"):
        run = ROOT / model / f"angle_{angle}" / f"h_{h}"
        mod = {s: E.body(run, "Cohesive_zone_model", s, h) for s in ("bottom", "top")}
        num = den = 0.0
        for s in ("bottom", "top"):
            (xr, cr, ur), (xm, cm_, um) = refb[s], mod[s]
            at = E.evaluate_at(xm, cm_, um, xr, E.matched(xr, xm))
            num += E.integrate(np.sum((ur - at) ** 2, axis=1), xr, cr)
            den += E.integrate(np.sum(ur ** 2, axis=1), xr, cr)
        jr, pr = E.jump(refb["bottom"], refb["top"], E.H0, angle)
        jm, pm = E.jump(mod["bottom"], mod["top"], E.H0, angle, on=(refb["bottom"], refb["top"]))
        k = E.interior(pr, angle, E.MARGIN)
        cm, tm = E.model_traction(run, h)
        tm = tm[E.matched(cr_[:, :2], cm[:, :2])]
        res[tstar][model] = {"e_u_global": float(np.sqrt(num / den)),
                             "e_jump": float(np.sqrt(0.5) / 30 * np.linalg.norm((jm - jr)[k]) / np.sqrt(den)),
                             "e_traction": float(np.linalg.norm((tm - tr)[kt]) / np.linalg.norm(tr[kt])),
                             "u_top": float(np.mean(mod["top"][2][:, 0]))}
    print(f"t={tstar:g}  " + "  ".join(f"{m}: e_u={res[tstar][m]['e_u_global']:.3e} e_jump={res[tstar][m]['e_jump']:.3e} "
                                       f"e_t={res[tstar][m]['e_traction']:.3e}" for m in ("full", "supp")), flush=True)
out.write_text(json.dumps(res, indent=1))
