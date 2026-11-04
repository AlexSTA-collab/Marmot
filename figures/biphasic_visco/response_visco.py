#!/usr/bin/env python3
"""Response of the loaded face and displacement error over the creep time, resolved layer against the interface models.

usage: response_visco.py <out.json> <case> <angle> <h> [models=full,supp,supp_eff]
For every output time: the mean x displacement of the loaded face (top face of the top body) of every model and the displacement
error e_u (paper definition, both bodies, node-matched) of every interface model against the resolved layer."""
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "refinement")); sys.path.insert(0, str(HERE.parent / "port_2611"))
import errors_end as E  # noqa: E402

out, case, angle, h = Path(sys.argv[1]), sys.argv[2], int(sys.argv[3]), sys.argv[4]
models = (sys.argv[5] if len(sys.argv) > 5 else "full,supp,supp_eff").split(",")
TIMES = [1e-8, 1, 2, 3, 5, 7, 10, 15, 20, 30, 40, 60, 80, 100]
ROOT = HERE / os.environ.get("RUNS", "runs") / case
ref = ROOT / "ref" / f"angle_{angle}" / f"h_{h}"
res = {m: {"t": [], "u_face": [], "e_u": []} for m in ["ref"] + models}
for tstar in TIMES:
    os.environ["ERR_T"] = repr(tstar)
    rb = {s: E.body(ref, "Cauchy_Full_model", s, h) for s in ("bottom", "top")}
    xt = rb["top"][0]; ut = rb["top"][2]
    face = xt[:, 2] > xt[:, 2].max() - 1e-6
    res["ref"]["t"].append(tstar); res["ref"]["u_face"].append(float(ut[face, 0].mean())); res["ref"]["e_u"].append(0.0)
    for m in models:
        run = ROOT / m / f"angle_{angle}" / f"h_{h}"
        mb = {s: E.body(run, "Cohesive_zone_model", s, h) for s in ("bottom", "top")}
        num = den = 0.0
        for s in ("bottom", "top"):
            (xr, cr, ur), (xm, cm_, um) = rb[s], mb[s]
            at = E.evaluate_at(xm, cm_, um, xr, E.matched(xr, xm))
            num += E.integrate(np.sum((ur - at) ** 2, axis=1), xr, cr)
            den += E.integrate(np.sum(ur ** 2, axis=1), xr, cr)
        xm_t, um_t = mb["top"][0], mb["top"][2]
        fm = xm_t[:, 2] > xm_t[:, 2].max() - 1e-6
        res[m]["t"].append(tstar); res[m]["u_face"].append(float(um_t[fm, 0].mean())); res[m]["e_u"].append(float(np.sqrt(num / den)))
    print(f"t={tstar:8.3g} " + " ".join(f"{m}: e_u={res[m]['e_u'][-1]:.2e}" for m in models), flush=True)
out.write_text(json.dumps(res, indent=1))
