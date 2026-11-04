#!/usr/bin/env python3
"""Single-reference tables (the fixed-reference appendix) of the 60 x 60 elastic and viscoelastic analyses (make_jobs60.py):
the EIE at every h against the ONE resolved layer of thickness h0 = 0.01 (20 C3D8 through the layer), with the measures of
plastic_bbar/score_fixed.py (locate / measures, used verbatim), at the last output. At h = h0 both protocols compare the
same pair of solutions, so that entry is errors_end.measures, as in score60.py. Written in the layout read by
paper_scripts/figures/visco_matched/plot_visco.py.
usage: score_fixed60.py <case: elastic|visco10|visco0> <angle> <out.md>"""
import os
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
JOBS = Path(os.environ.get("RERUN60_JOBS", HERE / "jobs"))   # the shared-node reruns: their own root
sys.path.insert(0, str(HERE.parent / "refinement"))
import errors_end as E  # noqa: E402

case, angle, out = sys.argv[1], int(sys.argv[2]), Path(sys.argv[3])
FAMS, HS, N, H0 = ("stiff", "initial", "soft"), ("0.01", "0.02", "0.04", "0.06", "0.08", "0.1"), 60, "0.01"
KEYS = (("Displacement field (bodies)", "e_u_global"), ("Displacement jump", "e_jump"), ("Traction", "e_traction"))
src = (HERE.parent / "plastic_bbar" / "score_fixed.py").read_text()
ns = {"E": E, "np": np, "cKDTree": cKDTree, "N": N, "H0": H0}
exec(src[src.index("def locate"):src.index("for kind in")], ns)   # locate, measures of the J2 appendix
measures = ns["measures"]


def rate(v):
    v = np.asarray(v, float)
    hs = np.array([float(h) for h in HS])
    ok = np.isfinite(v) & (v > 0)
    return np.polyfit(np.log(hs[ok]), np.log(v[ok]), 1)[0] if ok.sum() > 2 else np.nan


E.ANGLE = angle
res = {}
for fam in FAMS:
    r = JOBS / case / "ref" / fam / f"angle_{angle}" / f"h_{H0}"
    for h in HS:
        m = JOBS / case / "ifc" / fam / f"angle_{angle}" / f"h_{h}"
        E.H0 = float(h)
        res[(fam, h)] = E.measures(str(m), h, str(r), h, N) if h == H0 else measures(m, h, r, angle)
        print(f"{case} a={angle} {fam:8s} h={h:5s} " + "  ".join(f"{k}={res[(fam, h)][k]:.4e}" for _, k in KEYS), flush=True)
lines = [f"# Single-reference convergence, angle {angle}, {case}, 60 x 60 mesh, reference h0 = 0.01, 20 C3D8 through the layer",
         f"# reference: fixed (h0 = 0.01), edge margin: {E.MARGIN:g}, plastic_bbar/score_fixed.py measures with N = 60", ""]
for title, k in KEYS:
    lines += [f"## {title}", "", "| h | " + " | ".join(f"{f} YIQUAD4" for f in FAMS) + " |",
              "|---|" + "---|" * len(FAMS)]
    for h in HS:
        lines.append(f"| {h} | " + " | ".join(f"{res[(f, h)][k]:.4e}" for f in FAMS) + " |")
    lines.append("| rate | " + " | ".join(f"{rate([res[(f, h)][k] for h in HS]):.3f}" for f in FAMS) + " |")
    lines.append("")
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text("\n".join(lines))
print("wrote", out)
