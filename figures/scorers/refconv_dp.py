#!/usr/bin/env python3
"""Through-thickness convergence of the Drucker-Prager resolved layer (paper measures e_u, e_[u], e_t, margin 0.15,
N = 60), alpha = 10: 'succ' = n/2 vs n layers, 'conv' = n vs 20 layers, at h = 0.01 and 0.1.
The 20-layer runs are <case>/ref/, the 5- and 10-layer ones <case>/refconv/n{5,10}/.
usage: refconv_dp.py <case dir> <out.json>"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path.home() / "psie26/tools/refinement"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "refinement"))
import errors_end as E  # noqa: E402

N = 60


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


if __name__ == "__main__":
    root, out = Path(sys.argv[1]), Path(sys.argv[2])
    res = {}
    E.ANGLE = 10
    for h in ("0.01", "0.1"):
        E.H0 = float(h)
        run = {n: root / "refconv" / f"n{n}" / "angle_10" / f"h_{h}" for n in (5, 10)}
        run[20] = root / "ref" / "angle_10" / f"h_{h}"
        for key, (a, b) in {"succ/5-10": (5, 10), "succ/10-20": (10, 20), "conv/5-20": (5, 20)}.items():
            res[f"{h}/{key}"] = ref_vs_ref(run[a], run[b], h)
            print(h, key, " ".join(f"{v:.3e}" for v in res[f"{h}/{key}"]), flush=True)
    out.write_text(json.dumps(res, indent=1, sort_keys=True))
