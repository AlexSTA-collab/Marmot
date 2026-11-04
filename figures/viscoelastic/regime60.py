#!/usr/bin/env python3
"""Regime indicator ||[u_i]|| / (h ||<u^S_{i,j}>||) of the 60 x 60 step-length sweep (make_regime60.py), Appendix G.

paper/figures/regime/extract_regime.py computes it, unchanged except for the reader: the runs are 26.11 output (one
file per step), read with port_2611/ensight_reader_2611, last step. The soft family is read from the double-precision
displacements of its rerun with the npydump output (regime_dt_npy, npydump/efe_npy.py): its bodies barely deform while
the layer takes the whole displacement (about 50 mm), and in single precision their surface strain (about 1e-6) is
below the resolution of the output (about 4e-4 per node). Then a copy of paper/figures/regime/plot_regime.py
draws the panels from the new curves.

    python3 regime60.py <out dir>   ->  <out dir>/regime_chi/{stiff,initial,soft}.pdf, common_legend.pdf
    REGIME_ROOT=<root> reads the runs from <root>/<family>/h_<h>/dt_<dt> instead (the shared-node reruns, soft with npydump)."""
import importlib.util
import os
import pickle
import subprocess
import sys
from glob import glob
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent.parent / "paper"
OUT = Path(sys.argv[1]).resolve(); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE.parent / "port_2611"))
from ensight_reader_2611 import read_geo, read_var_node_all_steps  # noqa: E402

spec = importlib.util.spec_from_file_location("extract_regime", PAPER / "figures" / "regime" / "extract_regime.py")
R = importlib.util.module_from_spec(spec); spec.loader.exec_module(R)


def body_2611(run, model, which, h):
    d = run / f"{R.PREFIX[model]}_{which}_{h}"
    parts = read_geo(str(d / "geometry.geo_0000"))
    f = sorted(glob(str(d / f"displacement_{which}.var_[0-9][0-9][0-9][0-9]")))[-1]
    u = read_var_node_all_steps(f, parts, 3)
    name = next(p for p, v in parts.items() if v["desc"] == f"{which}-body" and p in u)
    return parts[name]["xyz"], u[name][-1:]


def body_npy(run, model, which, h):
    xyz, u32 = body_2611(run, model, which, h)
    u = np.load(run / f"npy_{which}.npy")
    assert np.array_equal(u.astype(np.float32), u32[-1].astype(np.float32))   # same nodes in the same order
    return xyz, u[None]


R._body = body_2611
C = {}
for fam in R.FAMILIES:
    for h in R.HEIGHTS:
        row = {}
        for dt in R.DT:
            run = (Path(os.environ["REGIME_ROOT"]) if os.environ.get("REGIME_ROOT")
                   else HERE / ("regime_dt_npy" if fam == "soft" else "regime_dt")) / fam / f"h_{h}" / f"dt_{dt:g}"
            R._body = body_npy if fam == "soft" else body_2611
            log = run / "run.log"
            if not (log.exists() and "Job computation time" in log.read_text() and "Simulation failed" not in log.read_text()):
                print("missing", run); continue
            row[dt] = {"interface": float(R.indicator(run, "interface", h)[0][-1])}
        C[("dtfull", fam, h)] = row
        print(f"{fam:8s} h={h:5s} " + " ".join(f"{row[d]['interface']:.4g}" for d in sorted(row)), flush=True)
work = OUT / "_work"; work.mkdir(exist_ok=True)
(work / "regime_curves.pkl").write_bytes(pickle.dumps(C))
src = (PAPER / "figures" / "regime" / "plot_regime.py").read_text()
src = src.replace('sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "verification"))',
                  f'sys.path.insert(0, "{PAPER / "verification"}")')
(work / "plot_regime.py").write_text(src)
subprocess.run([sys.executable, str(work / "plot_regime.py"), str(OUT)], check=True)
