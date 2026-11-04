#!/usr/bin/env python3
"""Matched-thickness convergence tables of the 60 x 60 elastic and viscoelastic analyses (make_jobs60.py): the paper
measures (errors_end.measures: e_u, e_[u], e_t, interior margin 0.15, N = 60) of the EIE (YIQUAD4) against the
resolved layer of the same thickness (20 C3D8 through the layer), at the last output (the end of the creep hold for
the viscoelastic cases). Written in the layout of paper/verification/node_matched_errors.py, so that
figures/visco_matched/plot_visco.py draws them unchanged; the rate is the same least-squares slope.
usage: score60.py <case: elastic|visco10|visco0> <angle> <out.md>"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "refinement"))
import errors_end as E  # noqa: E402

case, angle, out = sys.argv[1], int(sys.argv[2]), Path(sys.argv[3])
FAMS, HS = ("stiff", "initial", "soft"), ("0.01", "0.02", "0.04", "0.06", "0.08", "0.1")
KEYS = (("Displacement field (bodies)", "e_u_global"), ("Displacement jump", "e_jump"), ("Traction", "e_traction"))


def done(d):
    log = d / "run.log"
    return log.exists() and "Job computation time" in log.read_text() and "Simulation failed" not in log.read_text()


def rate(v):
    v = np.asarray(v, float)
    hs = np.array([float(h) for h in HS])
    ok = np.isfinite(v) & (v > 0)
    return np.polyfit(np.log(hs[ok]), np.log(v[ok]), 1)[0] if ok.sum() > 2 else np.nan


E.ANGLE = angle
res = {}
for fam in FAMS:
    for h in HS:
        m = HERE / "jobs" / case / "ifc" / fam / f"angle_{angle}" / f"h_{h}"
        r = HERE / "jobs" / case / "ref" / fam / f"angle_{angle}" / f"h_{h}"
        if not (done(m) and done(r)):
            res[(fam, h)] = {k: np.nan for _, k in KEYS}
            print(f"{case} a={angle} {fam} h={h}: missing", flush=True)
            continue
        E.H0 = float(h)
        res[(fam, h)] = E.measures(str(m), h, str(r), h, 60)
        print(f"{case} a={angle} {fam:8s} h={h:5s} " + "  ".join(f"{k}={res[(fam, h)][k]:.4e}" for _, k in KEYS), flush=True)
lines = [f"# Matched-thickness convergence, angle {angle}, {case}, 60 x 60 mesh, reference 20 C3D8 through the layer",
         f"# reference: matched, edge margin: {E.MARGIN:g}, errors_end.measures with N = 60", ""]
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
