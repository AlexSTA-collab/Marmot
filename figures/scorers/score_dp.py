#!/usr/bin/env python3
"""Paper error measures (errors_end.measures, margin 0.15, N = 60) of EIE against the 20-layer C3D8B resolved layer of
the same thickness, for the Drucker-Prager jobs of make_dp_jobs.py; adds what is finished to the output json.
Also e_u over the interior: the same integral restricted to the body nodes whose in-plane position, projected along
the normal onto the interface plane, lies 0.15 mm or more inside its edges (the face is 1/cos(alpha) long along t1
and 1 wide in y), the margin the jump and traction measures already use (key eie_int2/...). The keys eie_int/... of
older json files used E.interior on the 3D body nodes, which at alpha = 10 keeps a slab leaning through the body
instead; they are dropped.
usage: score_dp.py <case dir (with eie/ and ref/)> <out.json> [model dir = eie, e.g. psie]"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path.home() / "psie26/tools/refinement"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "refinement"))
import errors_end as E  # noqa: E402

root, out = Path(sys.argv[1]), Path(sys.argv[2])
kind = sys.argv[3] if len(sys.argv) > 3 else "eie"
res = json.loads(out.read_text()) if out.exists() else {}
res = {k: v for k, v in res.items() if "_int/" not in k}


def in_plane_interior(x, a):
    ar = np.radians(a)
    t1 = np.array([np.cos(ar), 0.0, -np.sin(ar)])
    return (np.abs(x @ t1) <= 0.5 / np.cos(ar) - E.MARGIN) & (np.abs(x[:, 1]) <= 0.5 - E.MARGIN)


def eu_interior(mrun, h, rrun, a):
    num = den = 0.0
    for s in ("bottom", "top"):
        (xr, cr, ur), (xm, cm_, um) = E.body(rrun, "Cauchy_Full_model", s, h), E.body(mrun, "Cohesive_zone_model", s, h)
        w = in_plane_interior(xr, a).astype(float)
        um_at_ref = E.evaluate_at(xm, cm_, um, xr, E.matched(xr, xm))
        num += E.integrate(((ur - um_at_ref) ** 2).sum(axis=1) * w, xr, cr)
        den += E.integrate((ur ** 2).sum(axis=1) * w, xr, cr)
    return (num / den) ** 0.5


def done(d):
    log = d / "run.log"
    return log.exists() and "Job computation time" in log.read_text() and "Simulation failed" not in log.read_text()


for a in ("10", "0"):
    for h in ("0.01", "0.02", "0.04", "0.06", "0.08", "0.1"):
        key = f"{kind}/{a}/{h}"
        m, r = root / kind / f"angle_{a}" / f"h_{h}", root / "ref" / f"angle_{a}" / f"h_{h}"
        if key in res or not (done(m) and done(r)):
            continue
        E.ANGLE, E.H0 = int(a), float(h)
        o = E.measures(str(m), h, str(r), h, 60)
        res[key] = [o["e_u_global"], o["e_jump"], o["e_traction"]]
        print(key, " ".join(f"{v:.3e}" for v in res[key]), flush=True)
    for h in ("0.01", "0.02", "0.04", "0.06", "0.08", "0.1"):
        key = f"{kind}_int2/{a}/{h}"
        m, r = root / kind / f"angle_{a}" / f"h_{h}", root / "ref" / f"angle_{a}" / f"h_{h}"
        if key in res or not (done(m) and done(r)):
            continue
        E.ANGLE, E.H0 = int(a), float(h)
        res[key] = eu_interior(str(m), h, str(r), int(a))
        print(key, f"e_u interior {res[key]:.3e}  (whole domain {res[f'{kind}/{a}/{h}'][0]:.3e})", flush=True)
out.write_text(json.dumps(res, indent=1, sort_keys=True))
