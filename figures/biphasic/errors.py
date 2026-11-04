#!/usr/bin/env python3
"""Errors of the full model and of the model with the mismatch suppressed against the resolved biphasic layer.

    python3 errors.py <out.json>
    ROOT=runs60_as2605 NCELLS=60 python3 errors.py <out.json>     (the 60 x 60 meshes, make_jobs_60.py)

Per case, angle and thickness and for each interface model (full, supp) the two measures of the study, read off
the last step by the study's own node-matched routines (port_2611/scripts/node_matched_errors.py):
  global  ||u^h - u^0|| / ||u^0||        over both bodies
  e_u     the displacement-jump error, normalised as in the study
plus the force on the loaded face, F, of all three models and the ratio F_model / F_ref.
"""
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "port_2611" / "scripts"))
import node_matched_errors as nme  # noqa: E402

ROOT = HERE / os.environ.get("ROOT", "runs_as2605")
# W_AREA = sqrt(0.5) / N turns the sum over the N x N face grid into the area integral (N = 30 in the paper's study)
nme.W_AREA = np.sqrt(0.5) / int(os.environ.get("NCELLS", 30))
HEIGHTS = ["0.01", "0.02", "0.04", "0.06", "0.08", "0.1"]


def force_x(run, side, h):
    d, parts = nme.output(run, side, h)
    r = nme.read_var_node_all_steps(str(d / f"reaction_{side}.var"), parts, 3)
    return float(sum(v[-1][:, 0].sum() for v in r.values()))


def layer_traction(run, h, angle):
    """sigma.n of every layer element of both phases, averaged over each through-thickness column
    (six equal elements, three per phase), paired by in-plane position as in nme.layer_traction."""
    cs, ss = [], []
    for half in ("layerbottom", "layertop"):
        d, parts = nme.output(run, half, h)
        s = nme.read_elem_var(str(d / f"stress_{half[:5]}_{half[5:]}.var"), parts, 6)
        pn = list(s)[0]
        cs.append(nme.centroids(parts, pn)); ss.append(s[pn])
    s, c, n = np.vstack(ss), np.vstack(cs), nme.normal(angle)
    t = np.stack([s[:, 0] * n[0] + s[:, 3] * n[1] + s[:, 4] * n[2],
                  s[:, 3] * n[0] + s[:, 1] * n[1] + s[:, 5] * n[2],
                  s[:, 4] * n[0] + s[:, 5] * n[1] + s[:, 2] * n[2]], axis=1)
    key = np.stack([np.round(c[:, 0], 5), np.round(c[:, 1], 5)], axis=1)
    uniq, inv = np.unique(key, axis=0, return_inverse=True)
    inv = inv.ravel()
    cc = np.array([c[inv == i].mean(axis=0) for i in range(len(uniq))])
    tt = np.array([t[inv == i].mean(axis=0) for i in range(len(uniq))])
    return nme.by_position(cc, tt)


def measures(case, angle, h):
    ref = ROOT / case / "ref" / f"angle_{angle}" / f"h_{h}"
    refb = {s: nme.body(ref, s, h) for s in ("bottom", "top")}
    out = {"F_ref": force_x(ref, "top", h)}
    for model in ("full", "supp"):
        run = ROOT / case / model / f"angle_{angle}" / f"h_{h}"
        mod = {s: nme.body(run, s, h) for s in ("bottom", "top")}
        num = den = den_all = 0.0
        for s in ("bottom", "top"):
            (xr, cr, ur), (xm, cm, um) = refb[s], mod[s]
            at = nme.evaluate_at(xm, cm, um, xr, nme.matched(xr, xm))
            num += nme.integrate(np.sum((ur - at) ** 2, axis=1), xr, cr)
            den += nme.integrate(np.sum(ur ** 2, axis=1), xr, cr)
        jr, pr = nme.jump(refb["bottom"], refb["top"], h, angle)
        jm, pm = nme.jump(mod["bottom"], mod["top"], h, angle, on=(refb["bottom"], refb["top"]))
        k = nme.interior(pr, angle, nme.MARGIN)
        cm, tm = nme.model_traction(run, h)
        cr_, tr = layer_traction(ref, h, angle)
        tm = tm[nme.matched(cr_[:, :2], cm[:, :2])]
        kt = nme.interior(cr_, angle, nme.MARGIN)
        out[model] = {"global": float(np.sqrt(num / den)),
                      "e_t": float(np.linalg.norm((tm - tr)[kt]) / np.linalg.norm(tr[kt])),
                      "e_u": float(nme.W_AREA * np.linalg.norm((jm - jr)[k]) / np.sqrt(den)),
                      "F": force_x(run, "top", h)}
        out[model]["F_ratio"] = out[model]["F"] / out["F_ref"]
    return out


res = {}
for case in ("c10", "c100"):
    for angle in (0, 10):
        for h in HEIGHTS:
            m = measures(case, angle, h)
            res[f"{case}/{angle}/{h}"] = m
            print(f"{case} a={angle:2d} h={h:5s}  F_ref={m['F_ref']:.5e}  "
                  + "  ".join(f"{k}: global={m[k]['global']:.3e} e_u={m[k]['e_u']:.3e} e_t={m[k]['e_t']:.3e} F/Fref={m[k]['F_ratio']:.5f}"
                              for k in ("full", "supp")), flush=True)
Path(sys.argv[1]).write_text(json.dumps(res, indent=1))
