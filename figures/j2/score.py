#!/usr/bin/env python3
"""Paper error measures (errors_end.measures: e_u, e_[u], e_t, interior margin 0.15, N = 60) of PSIE and EIE against the
20-layer C3D8B resolved layer of the same thickness, for every family, inclination, H and h of make_jobs.py.
usage: score.py <jobs dir> <out.json>"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "refinement"))
import errors_end as E  # noqa: E402

root, out = Path(sys.argv[1]), Path(sys.argv[2])
res = json.loads(out.read_text()) if out.exists() else {}


def done(d):
    log = d / "run.log"   # a failed run also prints its computation time
    return log.exists() and "Job computation time" in log.read_text() and "Simulation failed" not in log.read_text()


for ref in sorted(root.glob("ref/*/angle_*/H_*/h_*")):
    fam, a, H, h = ref.parts[-4], ref.parts[-3][6:], ref.parts[-2][2:], ref.parts[-1][2:]
    if not done(ref):
        continue
    E.ANGLE, E.H0 = int(a), float(h)
    for kind in ("psie", "eie"):
        key = f"{kind}/{fam}/{a}/{H}/{h}"
        m = root / kind / fam / f"angle_{a}" / f"H_{H}" / f"h_{h}"
        if key in res or not done(m):
            continue
        try:
            res[key] = E.measures(str(m), h, str(ref), h, 60)
        except (SystemExit, Exception) as e:   # noqa: BLE001
            res[key] = {"error": str(e)[:200]}
        print(key, res[key], flush=True)
out.write_text(json.dumps(res, indent=1, sort_keys=True))
