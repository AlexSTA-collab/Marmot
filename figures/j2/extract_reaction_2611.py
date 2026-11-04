#!/usr/bin/env python3
"""Reaction-displacement curves of the J2 global response (Figure 20, Appendix H) for plot_reaction_H100.py: the mean
displacement and the total reaction of the loaded top face along the shear direction, from the 26.11 output of each run.
The runs are first copied to the 26.05 EnSight layout (port_2611/to_2605_layout.py) and read with curve() of
paper/figures/reaction/extract_reaction.py, the function behind the published curves, used unchanged.
usage: extract_reaction_2611.py <out.json> <angle> <label>=<run dir> [<label>=<run dir> ...]
  labels as plot_reaction_H100.py expects them: full, EIE, PSIE; each run dir holds <prefix>_top_0.01[.case]"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1] / "paper" / "figures"
ns = {"__file__": str(PAPER / "reaction" / "extract_reaction.py")}
src = (PAPER / "reaction" / "extract_reaction.py").read_text()
exec(compile(src[:src.index("\nVM = {")], "extract_reaction.py", "exec"), ns)   # read_geo, curve

out, angle = Path(sys.argv[1]), float(sys.argv[2])
res = {}
with tempfile.TemporaryDirectory() as tmp:
    for arg in sys.argv[3:]:
        label, run = arg.split("=", 1)
        dst = Path(tmp) / label
        subprocess.run([sys.executable, str(HERE.parent / "port_2611" / "to_2605_layout.py"), run, str(dst)],
                       check=True, stdout=subprocess.DEVNULL)
        prefix = "Cauchy_Full_model" if label == "full" else "Cohesive_zone_model"
        u, F = ns["curve"](dst / f"{prefix}_top_0.01", angle)
        # the 26.05 layout repeats the first increment; keep one of each repeated (u, F) pair, then lead with the
        # undeformed state twice, as the curves of the published figures do (plot_reaction_H100.py drops the first)
        keep = [0] + [i for i in range(1, len(u)) if (u[i], F[i]) != (u[i - 1], F[i - 1])]
        u, F = np.concatenate(([0.0], u[keep])), np.concatenate(([0.0], F[keep]))
        res[label] = [u.tolist(), F.tolist()]
        print(f"{label:5s} {len(u)} points, u_max {u.max():.4e}, F_max {F.max():.4e}", flush=True)
out.write_text(json.dumps(res))
print("wrote", out)
