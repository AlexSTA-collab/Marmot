#!/usr/bin/env python3
"""Through-thickness distribution of the strain in the resolved layer (last step), from the nodal displacements: per interior
column (margin 0.25) and element row k, the tangential shear gamma_k = d u_t / d n and the normal strain eps_k = d u_n / d n,
each divided by its column mean over the thickness (the homogeneous value); eps_over_gamma divides the normal
strain by the column's mean shear instead, which stays defined where the mean normal strain vanishes. Averaged over the columns.
Reports the row profile, the outer-row ratio and where the maximum sits. Mesh-dependent localisation shows as an outer-row
ratio that keeps growing as the rows are refined.
usage: row_strain.py <run dir> <h> [angle=10] [out.json]"""
import sys
from glob import glob
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import errors_end as E  # noqa: E402

run, h = Path(sys.argv[1]), float(sys.argv[2])
E.ANGLE = int(sys.argv[3]) if len(sys.argv) > 3 else 10
E.H0 = h
a = np.radians(E.ANGLE)
t1, n = np.array([np.cos(a), 0, -np.sin(a)]), E.normal(E.ANGLE)
d = Path(glob(str(run / "Cauchy_Full_model_interface_*[!e]"))[0])
parts = E.read_geo(str(d / "geometry.geo_0000"))
ud = E.read_var_node_all_steps(E.last(d, "displacement_interface"), parts, 3)
pk = list(ud)[0]
X, U = parts[pk]["xyz"], ud[pk][-1]
cols, inv = np.unique(np.round(X[:, :2], 5), axis=0, return_inverse=True)   # the layer is swept vertically: a column shares (x, y)
inv = inv.ravel()
keep = E.interior(np.c_[cols, np.zeros(len(cols))], E.ANGLE, 0.25)
G, P, R, Z = [], [], [], None
for c in np.flatnonzero(keep):
    idx = np.flatnonzero(inv == c)
    idx = idx[np.argsort(X[idx, 2])]
    z = X[idx, 2]
    Z = (z - z[0]) / (z[-1] - z[0]) if Z is None else Z
    ut, un = U[idx] @ t1, U[idx] @ n
    g = np.diff(ut) / np.diff(z)
    e = np.diff(un) / np.diff(z)
    G.append(g / (np.sum(g * np.diff(z)) / (z[-1] - z[0])))
    R.append(e / (np.sum(g * np.diff(z)) / (z[-1] - z[0])))   # normal strain over the column's mean shear
    P.append(e / (np.sum(e * np.diff(z)) / (z[-1] - z[0])) if abs(np.sum(e * np.diff(z))) > 1e-14 else e * np.nan)
G, P = np.array(G), np.array(P)
g, p = G.mean(axis=0), np.nanmean(P, axis=0)
nl = len(g)
print(f"{run.name} ({nl} rows, {len(G)} interior columns)")
print(f"  shear  gamma_k / mean, rows bottom->top: " + " ".join(f"{v:.3f}" for v in g))
print(f"         bottom row {g[0]:.3f}  top row {g[-1]:.3f}  centre {g[nl // 2]:.3f}  max at row {int(np.argmax(g))}  "
      f"column-wise max/mean (median over columns) {np.median(G.max(axis=1)):.3f}")
print(f"  normal eps_k / mean, rows bottom->top:   " + " ".join(f"{v:+.2f}" for v in p))
if len(sys.argv) > 4:
    import json
    Path(sys.argv[4]).write_text(json.dumps({"zeta": Z.tolist(), "gamma": g.tolist(), "eps": p.tolist(),
                                             "gamma_q": np.percentile(G, [10, 90], axis=0).tolist(),
                                             "eps_over_gamma": np.mean(R, axis=0).tolist()}))
