#!/usr/bin/env python3
"""Face histories of the full interface model (run with the extra outputs kink, stress_plus, stress_minus): mean over the
interior columns of the normal-gradient jump z = g+ - g- along the tangent (the difference of the shear strains of the two
faces, gamma+ - gamma-) and of the two face tractions sigma+.n, sigma-.n along the tangent.

usage: model_histories.py <out.json> <run dir> <h> <angle>"""
import json
import re
import sys
from glob import glob
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "refinement")); sys.path.insert(0, str(HERE.parent / "port_2611"))
import errors_end as E  # noqa: E402
from ensight_reader_2611 import read_geo  # noqa: E402

out, run, h, angle = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], int(sys.argv[4])
a = np.radians(angle)
n = np.array([np.sin(a), 0.0, np.cos(a)]); t = np.array([np.cos(a), 0.0, -np.sin(a)])
d = E.casedir(run, "Cohesive_zone_model", "interface", h)
parts = read_geo(str(d / "geometry.geo_0000"))


def series(name, ncomp):
    files = sorted(glob(str(d / f"{name}.var_[0-9][0-9][0-9][0-9]")))
    vals = []
    for f in files:
        v = E.read_elem_var(f, parts, ncomp); pn = list(v)[0]
        vals.append(v[pn])
    pn = list(E.read_elem_var(files[0], parts, ncomp))[0]
    c = E.centroids(parts, pn)
    txt = Path(str(d) + ".case").read_text()
    tt = np.array(re.findall(r"[-+0-9.e]+", txt.split("time values:")[-1].split("GEOMETRY")[0])[: len(files)], float)
    return tt, c, np.array(vals)


tt, c, z = series("kink", 3)
inner = (np.abs(c[:, 0]) < 0.3) & (np.abs(c[:, 1]) < 0.3)
zt = (z[:, inner, :] @ t).mean(axis=1)
zn = (z[:, inner, :] @ n).mean(axis=1)
# the formulation is incremental: the state normalGradientJump is the jump of the last increment, the face histories are its sums
res = {"t": tt.tolist(), "z_t": zt.tolist(), "z_n": zn.tolist(), "z_t_accumulated": np.cumsum(zt).tolist(), "z_n_accumulated": np.cumsum(zn).tolist()}
for nm in ("stress_plus", "stress_minus"):
    _, cc, s = series(nm, 6)
    s = s[:, inner, :]
    tr = np.stack([s[..., 0] * n[0] + s[..., 3] * n[1] + s[..., 4] * n[2],
                   s[..., 3] * n[0] + s[..., 1] * n[1] + s[..., 5] * n[2],
                   s[..., 4] * n[0] + s[..., 5] * n[1] + s[..., 2] * n[2]], axis=-1)
    res[nm + "_t"] = (tr @ t).mean(axis=1).tolist()
    res[nm + "_n"] = (tr @ n).mean(axis=1).tolist()
out.write_text(json.dumps(res, indent=1))
for tq in (1e-8, 1, 10, 100):
    i = int(np.argmin(abs(tt - tq)))
    print(f"t={tt[i]:9.3g}  z_t={zt[i]: .4e}  tau+={res['stress_plus_t'][i]: .4e}  tau-={res['stress_minus_t'][i]: .4e}")
