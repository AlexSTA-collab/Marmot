#!/usr/bin/env python3
"""Through-thickness convergence of the elastic and viscoelastic resolved references on the 60 x 60 mesh (jobs/ of make_ev60_jobs.py; n = 20 is the rerun60 reference): the
paper measures e_u, e_[u], e_t (errors_end, interior margin 0.15, N = 60) of the
layer with n elements against the layer with 2n, at the last output (the end of the creep hold in the viscoelastic
case): 'succ/10' = n = 5 against 10, 'succ/20' = n = 10 against 20, 'conv/5' = n = 5 against 20.
usage: score_ev.py <out.json>"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "refinement")); sys.path.insert(0, str(Path.home() / "psie26/tools/refinement"))
import errors_end as E  # noqa: E402

N = 60
N20 = {"elastic": HERE / "jobs" / "elastic" / "n20", "visco": HERE / "jobs" / "visco" / "n20"}


def done(run):
    log = run / "run.log"
    return log.exists() and "Job computation time" in log.read_text() and "Simulation failed" not in log.read_text()


def ref_vs_ref(run, rrun, h):
    A = {s: E.body(str(run), "Cauchy_Full_model", s, h) for s in ("bottom", "top")}
    R = {s: E.body(str(rrun), "Cauchy_Full_model", s, h) for s in ("bottom", "top")}
    num = den = 0.0
    for s in ("bottom", "top"):
        (xr, cr, ur), (xa, ca, ua) = R[s], A[s]
        d = ur - E.evaluate_at(xa, ca, ua, xr, E.matched(xr, xa))
        num += E.integrate(np.sum(d ** 2, axis=1), xr, cr)
        den += E.integrate(np.sum(ur ** 2, axis=1), xr, cr)
    jr, pr = E.jump(R["bottom"], R["top"], E.H0, E.ANGLE)
    ja, _ = E.jump(A["bottom"], A["top"], E.H0, E.ANGLE, on=(R["bottom"], R["top"]))
    k = E.interior(pr, E.ANGLE, E.MARGIN)
    ej = np.sqrt(0.5) / N * np.linalg.norm((ja - jr)[k]) / np.sqrt(den)
    ca_, ta = E.layer_traction(str(run), h, E.ANGLE)
    cr_, tr = E.layer_traction(str(rrun), h, E.ANGLE)
    ta = ta[E.matched(cr_[:, :2], ca_[:, :2])]
    k = E.interior(cr_, E.ANGLE, E.MARGIN)
    return [float(np.sqrt(num / den)), float(ej), float(np.linalg.norm((ta - tr)[k]) / np.linalg.norm(tr[k]))]


out = Path(sys.argv[1])
res = {}
E.ANGLE = 10
for kind in ("elastic", "visco"):
    for h in ("0.01", "0.02", "0.04", "0.06", "0.08", "0.1"):
        E.H0 = float(h)
        run = {n: HERE / "jobs" / kind / f"n{n}" / f"h_{h}" for n in (5, 10)}
        run[20] = N20[kind] / f"h_{h}"
        for key, (a, b) in {"succ/10": (5, 10), "succ/20": (10, 20), "conv/5": (5, 20)}.items():
            if done(run[a]) and done(run[b]):
                res[f"{kind}/{h}/{key}"] = ref_vs_ref(run[a], run[b], h)
                print(kind, h, key, " ".join(f"{v:.3e}" for v in res[f"{kind}/{h}/{key}"]), flush=True)
out.write_text(json.dumps(res, indent=1, sort_keys=True))
