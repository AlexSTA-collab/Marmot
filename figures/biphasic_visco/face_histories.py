#!/usr/bin/env python3
"""Histories of the two halves of the biphasic viscoelastic layer in the resolved reference: shear and normal strain of the
top and the bottom half and the layer traction, over the creep time.

usage: face_histories.py <out.json> <case> <angle> <h>"""
import json
import os
import re
import sys
from glob import glob
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "refinement")); sys.path.insert(0, str(HERE.parent / "port_2611"))
import errors_end as E  # noqa: E402
from ensight_reader_2611 import read_geo, read_var_node_all_steps  # noqa: E402

out, case, angle, h = Path(sys.argv[1]), sys.argv[2], int(sys.argv[3]), sys.argv[4]
ref = HERE / os.environ.get("RUNS", "runs") / case / "ref" / f"angle_{angle}" / f"h_{h}"
a = np.radians(angle)
n = np.array([np.sin(a), 0.0, np.cos(a)])
t = np.array([np.cos(a), 0.0, -np.sin(a)])


def times(d, nsteps):
    txt = Path(str(d) + ".case").read_text()
    v = re.findall(r"[-+0-9.e]+", txt.split("time values:")[-1].split("GEOMETRY")[0])
    return np.array(v[:nsteps], float)


res = {}
for half, name in (("layerbottom", "displacement_layer_bottom"), ("layertop", "displacement_layer_top")):
    d = E.casedir(ref, "Cauchy_Full_model", half, h)
    parts = read_geo(str(d / "geometry.geo_0000"))
    files = sorted(glob(str(d / f"{name}.var_[0-9][0-9][0-9][0-9]")))
    pn = list(read_var_node_all_steps(files[0], parts, 3))[0]
    xyz = parts[pn]["xyz"]
    s = xyz @ n
    lev = np.unique(np.round(s, 7))
    lo, hi = lev.min(), lev.max()
    cols = {}
    for i, p in enumerate(xyz):
        cols.setdefault((round(p[0] - s[i] * n[0], 5), round(p[1], 5)), []).append(i)
    inner = [i for k, v in cols.items() if abs(k[0]) < 0.3 and abs(k[1]) < 0.3 for i in v]
    inner = np.array(inner)
    glo = inner[np.abs(s[inner] - lo) < 1e-6]
    ghi = inner[np.abs(s[inner] - hi) < 1e-6]
    gam, eps = [], []
    for f in files:
        u = read_var_node_all_steps(f, parts, 3)[pn][-1]
        dut = (u[ghi] @ t).mean() - (u[glo] @ t).mean()
        dun = (u[ghi] @ n).mean() - (u[glo] @ n).mean()
        gam.append(dut / (hi - lo)); eps.append(dun / (hi - lo))
    res[half] = {"t": times(d, len(files)).tolist(), "gamma_t": gam, "eps_n": eps, "thickness": float(hi - lo)}
out.write_text(json.dumps(res, indent=1))
tb, gb = np.array(res["layerbottom"]["t"]), np.array(res["layerbottom"]["gamma_t"])
tt, gt = np.array(res["layertop"]["t"]), np.array(res["layertop"]["gamma_t"])
for tq in (1e-8, 1, 10, 100):
    i = np.argmin(abs(tb - tq))
    print(f"t={tb[i]:9.3g}  shear strain bottom half {gb[i]: .4e}  top half {gt[i]: .4e}  ratio {gb[i]/gt[i] if gt[i] else float('nan'):8.2f}")
