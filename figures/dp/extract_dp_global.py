#!/usr/bin/env python3
"""Reaction-displacement curves of the Drucker-Prager interphase (initial family, alpha = 10, h = 0.01, phi = 30, psi = 10)
for the resolved layer and the EIE, from the runs copied to global/initial_dp/{ref,eie} (workstation
dp/fam_initial_dp_cap/{ref,eie}/angle_10/h_0.01; the stiff family shows no visible yield in the global response). Step 1 (t <= 1) applies the confining pressure, step 2 the motion;
the applied displacement is the mean u_x of the top face, the reaction the sum of its x reactions, both of step 2 only.
usage: extract_dp_global.py <out.json>"""
import json
import sys
from glob import glob
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "port_2611"))
from ensight_reader_2611 import read_geo, read_var_node_all_steps  # noqa: E402


def curve(run, prefix):
    d = run / f"{prefix}_top_0.01"
    case = Path(str(d) + ".case").read_text()
    t = np.array(case.split("time set: 2")[1].split("time values:")[1].split()[:10000], dtype=object)
    parts = read_geo(str(d / "geometry.geo_0000"))
    pn = next(k for k, v in parts.items() if v["desc"] == "top-body")
    x = parts[pn]["xyz"]
    top = x[:, 2] > x[:, 2].max() - 1e-6
    files = sorted(glob(str(d / "displacement_top.var_[0-9][0-9][0-9][0-9]")))
    rfiles = sorted(glob(str(d / "reaction_top.var_[0-9][0-9][0-9][0-9]")))
    times = np.array([float(v) for v in t[:len(files)]])
    u = np.array([read_var_node_all_steps(f, parts, 3)[pn][-1][top, 0].mean() for f in files])
    F = np.array([read_var_node_all_steps(f, parts, 3)[pn][-1][top, 0].sum() for f in rfiles])
    k = times >= 1.0 - 1e-9
    return (u[k] - u[k][0]).tolist(), (F[k] - F[k][0]).tolist(), times[k].tolist()


out = {}
for rule in ("initial_dp",):   # psi = 10; the psi = 0 fields were pruned after scoring
    out[rule] = {"full": curve(HERE / "global" / rule / "ref", "Cauchy_Full_model"),
                 "EIE": curve(HERE / "global" / rule / "eie", "Cohesive_zone_model")}
    for k, (u, F, t) in out[rule].items():
        print(rule, k, len(u), f"u_end {u[-1]:.4e} F_end {F[-1]:.5f}")
Path(sys.argv[1]).write_text(json.dumps(out))
