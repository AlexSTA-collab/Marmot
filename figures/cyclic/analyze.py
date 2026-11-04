#!/usr/bin/env python3
"""History of the cyclic test: force on the loaded face and displacement error of the PSIE against the resolved layer.

    python3 analyze.py <family> <angle> <h> [out.json]

Reads the 26.11 output of runs/ref and runs/psie directly (one file per step). At the steps the two runs have in common:
  u(t)   the prescribed displacement of the loaded face, along the interface, SCALE * 0.01 * triangle(t) (runs_tri: five steps, reversals at the step ends)
  F(t)   the reaction of the top body along that direction, summed over its nodes (the reaction sign of 26.11 is
         flipped against 26.05; here F is oriented so that it has the sign of u)
  e(t)   ||u_psie - u_ref|| / max_t ||u_ref||, over the nodes of both bodies, with the nodes of the two meshes paired as
         in node_matched_errors.py (one to one, within 1e-4 mm). The norm of the reference is that of the whole history,
         since that of the step vanishes where u does.
"""
import json
import os
import re
import sys
from glob import glob
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "port_2611"))
from ensight_reader_2611 import read_geo, read_var_node_all_steps  # noqa: E402

fam, angle, h = sys.argv[1], int(sys.argv[2]), sys.argv[3]
SCALE = 0.5
U0 = SCALE * 0.01
d_load = np.array([np.cos(np.radians(angle)), 0.0, -np.sin(np.radians(angle))])


def times(case):
    t = Path(case).read_text()
    blocks = re.findall(r"time set: \d+ no description\nnumber of steps: (\d+)\nfilename start number: \d+\nfilename increment: \d+\ntime values: ([^\n]*\n(?:[-+0-9.e]+\n)*)", t)
    n, vals = max(blocks, key=lambda b: int(b[0]))
    return np.array([float(x) for x in vals.split()])


def read(run, prefix, side):
    d = Path(run) / f"{prefix}_{side}_{h}"
    parts = read_geo(str(d / "geometry.geo_0000"))
    out = {}
    for var in ("displacement", "reaction"):
        files = sorted(glob(str(d / f"{var}_{side}.var_*")))
        got = [read_var_node_all_steps(f, parts, 3) for f in files]
        pn = list(got[0])[0]
        out[var] = np.array([g[pn][-1] for g in got])
    return parts[pn]["xyz"], out["displacement"], out["reaction"], times(Path(run) / f"{prefix}_{side}_{h}.case")


R = {}
for model, prefix in (("ref", "Cauchy_Full_model"), ("psie", "Cohesive_zone_model")):
    run = HERE / os.environ.get("RUNS", "runs_tri") / model / fam / f"angle_{angle}" / f"h_{h}"
    R[model] = {s: read(run, prefix, s) for s in ("top", "bottom")}

# increments both runs completed; cutbacks may differ, so match them by time
tr_all = R["ref"]["top"][3][: min(len(R["ref"][s][1]) for s in R["ref"])]
tp_all = R["psie"]["top"][3][: min(len(R["psie"][s][1]) for s in R["psie"])]
ir = [i for i, a in enumerate(tr_all) if np.any(np.abs(tp_all - a) < 1e-9)]
ip = [int(np.argmin(np.abs(tp_all - tr_all[i]))) for i in ir]
for m, ix in (("ref", ir), ("psie", ip)):
    R[m] = {s: (R[m][s][0], R[m][s][1][ix], R[m][s][2][ix], R[m][s][3][ix]) for s in R[m]}
n = len(ir)
tr, tp = R["ref"]["top"][3], R["psie"]["top"][3]
t = tr
u = U0 * np.interp(t, [0, 1, 2, 3, 4, 5], [0, 1, 0, -1, 0, 1])

# force on the loaded face
F = {m: np.array([(R[m]["top"][2][i] @ d_load).sum() for i in range(n)]) for m in R}
sgn = np.sign(np.sum(F["ref"] * u)) or 1.0
F = {m: sgn * F[m] for m in F}

# displacement error
num = np.zeros(n)
for side in ("top", "bottom"):
    xr, ur = R["ref"][side][0], R["ref"][side][1][:n]
    xm, um = R["psie"][side][0], R["psie"][side][1][:n]
    dist, idx = cKDTree(xm).query(xr)
    assert dist.max() < 1e-4 and len(np.unique(idx)) == len(idx), f"{side}: nodes not paired one to one ({dist.max():.1e})"
    num += np.sum((ur - um[:, idx]) ** 2, axis=(1, 2))
den = sum(np.sum(R["ref"][s][1][:n] ** 2, axis=(1, 2)) for s in ("top", "bottom")).max()
e = np.sqrt(num / den)

res = dict(t=t.tolist(), u=u.tolist(), F_ref=F["ref"].tolist(), F_psie=F["psie"].tolist(), e=e.tolist())
print(f"{fam} angle {angle} h={h}: {n} common steps, t up to {t[-1]:.3f}")
print(f"  F peak (ref) {np.abs(F['ref']).max():.4e}   max |F_psie - F_ref| / peak {np.abs(F['psie'] - F['ref']).max() / np.abs(F['ref']).max():.3e}")
print(f"  displacement error: max {e.max():.3e}, at the last step {e[-1]:.3e}")
for k in range(0, n, max(1, n // 12)):
    print(f"   t={t[k]:.3f} u={u[k]:+.3e} F_ref={F['ref'][k]:+.4e} F_psie={F['psie'][k]:+.4e} e={e[k]:.2e}")
if len(sys.argv) > 4:
    Path(sys.argv[4]).write_text(json.dumps(res))
