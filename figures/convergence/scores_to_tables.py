#!/usr/bin/env python3
"""Markdown convergence tables (the layout of paper/verification/node_matched/*_table_angle*.md, read by
paper/figures/visco_matched/plot_visco.py) from score.py output: EIE and PSIE against the 20-layer C3D8B resolved layer
of the same thickness, every family, one table per element and inclination, at one hardening modulus.
usage: scores_to_tables.py <scores.json> <H> <out dir> [--fixed]   (--fixed: score_fixed.py output, single reference h0 = 0.01)"""
import json
import sys
from pathlib import Path

import numpy as np

FIXED = "--fixed" in sys.argv
sys.argv = [a for a in sys.argv if a != "--fixed"]
res, H, out = json.loads(Path(sys.argv[1]).read_text()), sys.argv[2], Path(sys.argv[3])
out.mkdir(parents=True, exist_ok=True)
HS = ["0.01", "0.02", "0.04", "0.06", "0.08", "0.1"]
FAMS = ("stiff", "initial", "soft")
ELEMENT = {"eie": "YIQUAD4", "psie": "YIQUAD4_STABP_MINI"}
SECTIONS = (("Displacement field (bodies)", "e_u_global"), ("Displacement jump", "e_jump"), ("Traction", "e_traction"))
for kind, case in (("eie", "vm_eie"), ("psie", "vm_psie")):
    for a in ("0", "10"):
        lines = [f"# {'Single-reference' if FIXED else 'Matched-thickness'} convergence, angle {a}, {case}, H = {H} MPa",
                 f"# reference: {'fixed (h0 = 0.01)' if FIXED else 'matched'}, 20 C3D8B elements through the thickness, "
                 "60 x 60 in plane, edge margin: 0.15",
                 f"# reproducibility/plastic_bbar/{'score_fixed.py' if FIXED else 'score.py'}", ""]
        for sec, key in SECTIONS:
            lines += [f"## {sec}", "", "| h | " + " | ".join(f"{f} {ELEMENT[kind]}" for f in FAMS) + " |",
                      "|---|" + "---|" * len(FAMS)]
            cols = {f: [res[f"{kind}/{f}/{a}/{H}/{h}"][key] for h in HS] for f in FAMS}
            for i, h in enumerate(HS):
                lines.append(f"| {h} | " + " | ".join(f"{cols[f][i]:.4e}" for f in FAMS) + " |")
            x = np.log([float(h) for h in HS])
            lines += ["| rate | " + " | ".join(f"{np.polyfit(x, np.log(cols[f]), 1)[0]:.3f}" for f in FAMS) + " |", ""]
        name = f"{case}_{'fixed_' if FIXED else ''}table_angle{a}.md"
        (out / name).write_text("\n".join(lines))
        print(out / name)
