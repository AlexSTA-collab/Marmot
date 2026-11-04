#!/usr/bin/env python3
"""Tables of the Drucker--Prager interphase (phi = 30, psi = 10, EIE) in the layout read by
paper/figures/visco_matched/plot_visco.py, so that Figures 23, 24, I.40 and I.41 are drawn as the other convergence figures:
  data/dp/dp_eie_table_angle<a>.md         matched reference (data/dp/<family>_dp.json, dp/score_dp.py)
  data/dp/dp_eie_fixed_table_angle<a>.md   single reference h0 = 0.01 (data/dp_fixed/<family>_dp.json, dp/score_fixed_dp.py)
usage: python3 make_dp_figure_tables.py"""
import json
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parent / "data"
FAMS, HS = ("stiff", "initial", "soft"), ("0.01", "0.02", "0.04", "0.06", "0.08", "0.1")
KEYS = ("Displacement field (bodies)", "Displacement jump", "Traction")


def rate(v):
    v, h = np.asarray(v, float), np.array([float(x) for x in HS])
    ok = np.isfinite(v) & (v > 0)
    return np.polyfit(np.log(h[ok]), np.log(v[ok]), 1)[0] if ok.sum() > 2 else np.nan


for src, name, ref in (("dp", "dp_eie_table", "matched"), ("dp_fixed", "dp_eie_fixed_table", "fixed (h0 = 0.01)")):
    for a in ("0", "10"):
        d = {f: json.loads((DATA / src / f"{f}_dp.json").read_text()) for f in FAMS}
        lines = [f"# Drucker--Prager interphase, phi = 30, psi = 10, H = 0.075 G, angle {a}, EIE, 60 x 60 mesh, reference 20 C3D8",
                 f"# reference: {ref}, edge margin 0.15, errors_end.measures with N = 60", ""]
        for k, title in enumerate(KEYS):
            lines += [f"## {title}", "", "| h | " + " | ".join(f"{f} YIQUAD4" for f in FAMS) + " |", "|---|" + "---|" * len(FAMS)]
            for h in HS:
                lines.append(f"| {h} | " + " | ".join(f"{d[f][f'eie/{a}/{h}'][k]:.4e}" for f in FAMS) + " |")
            lines.append("| rate | " + " | ".join(f"{rate([d[f][f'eie/{a}/{h}'][k] for h in HS]):.3f}" for f in FAMS) + " |")
            lines.append("")
        (DATA / "dp" / f"{name}_angle{a}.md").write_text("\n".join(lines))
        print("wrote", DATA / "dp" / f"{name}_angle{a}.md")
