#!/usr/bin/env python3
"""Through-thickness convergence of the Drucker-Prager references of the paper (S75: H = 0.075 G, roller from t = 1) over the thickness sweep (paper measures e_u,
e_[u], e_t of refconv_dp.ref_vs_ref, margin 0.15, N = 60), alpha = 10, families stiff and initial, phi 30, psi 0 and 10:
'5-10' = 5 against 10 layers, '10-20' = 10 against 20 layers. The 20-layer runs are S75_<case>/ref/, the 5- and
10-layer ones S75_<case>/refconv/n{5,10}/; keys <case>/<h>/<a>-<b>. Adds what is finished to the output json.
usage: refconv_dp_sweep.py <dp dir> <out.json>"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from refconv_dp import E, ref_vs_ref  # noqa: E402

root, out = Path(sys.argv[1]), Path(sys.argv[2])
res = json.loads(out.read_text()) if out.exists() else {}
E.ANGLE = 10
for f in ("stiff_psi0", "stiff_dp", "initial_psi0", "initial_dp"):
    case = root / f"S75_{f}"
    for h in ("0.01", "0.02", "0.04", "0.06", "0.08", "0.1"):
        E.H0 = float(h)
        run = {n: case / "refconv" / f"n{n}" / "angle_10" / f"h_{h}" for n in (5, 10)}
        run[20] = case / "ref" / "angle_10" / f"h_{h}"
        for a, b in ((5, 10), (10, 20)):
            key = f"{f}/{h}/{a}-{b}"
            done = all("Job computation time" in (run[n] / "run.log").read_text() for n in (a, b)
                       if (run[n] / "run.log").exists()) and all((run[n] / "run.log").exists() for n in (a, b))
            if key in res or not done:
                continue
            res[key] = ref_vs_ref(run[a], run[b], h)
            print(key, " ".join(f"{v:.3e}" for v in res[key]), flush=True)
            out.write_text(json.dumps(res, indent=1, sort_keys=True))
