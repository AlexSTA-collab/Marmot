#!/usr/bin/env python3
"""Rows of Table R1 of the reply (convergence of the resolved references through the thickness), from the data of
Figures R1-R3: for every reference the largest change of the layer from n = 10 to n = 20 over the thickness sweep in
e_u, e_[u], e_t, the range of its ratio to the change from n = 5 to n = 10 (criterion (i)), and the smallest factor by
which it lies below the error of the interface model against the layer with n = 20 (criterion (ii), the grey curves).
Pairs where both changes lie below 1e-7 (the floor of the analysis) are left out of the ratio.
usage: refconv_table.py <S75 refconv sweep json>"""
import json
import sys
from pathlib import Path

LP = Path(__file__).resolve().parents[3]
RP, DATA = LP / "reproducibility", LP / "IJES_submission_Revision/tables/data"
HS = ["0.01", "0.02", "0.04", "0.06", "0.08", "0.1"]
FLOOR = 1e-7


def model_md(md):
    out, k = {}, -1
    for line in (DATA / md).read_text().splitlines():
        if line.startswith("## "):
            k += 1; out[k] = {}
        elif line.startswith("| 0") and k >= 0:
            c = [s.strip() for s in line.strip("|").split("|")]; out[k][c[0]] = float(c[1])
    return [[out[k][h] for h in HS] for k in range(3)]


def row(name, case, a, b, e):
    """a, b, e: [measure][h] of the 5/10 change, the 10/20 change and the interface model."""
    big = [max(b[k]) for k in range(3)]
    ratio = [a[k][i] / b[k][i] for k in range(3) for i in range(6) if max(a[k][i], b[k][i]) > FLOOR]
    factor = min(e[k][i] / b[k][i] for k in range(3) for i in range(6))
    f = lambda v: f"${v:.1e}".replace("e-0", r"\cdot10^{-").replace("e-", r"\cdot10^{-") + "}$"
    print(f"    {name:24s} & {case:22s} & " + " & ".join(f(v) for v in big)
          + f" & ${min(ratio):.1f}$--${max(ratio):.1f}$ & ${factor:.0f}$\\\\")


ev = json.loads((RP / "refinement_ev60/scores_ev60.json").read_text())
for layer, md, name in (("elastic", "elastic_table_angle10.md", "elastic, C3D8"), ("visco", "visco10_table_angle10.md", "viscoelastic, C3D8")):
    a = [[ev[f"{layer}/{h}/succ/10"][k] for h in HS] for k in range(3)]
    b = [[ev[f"{layer}/{h}/succ/20"][k] for h in HS] for k in range(3)]
    row(name, "stiff", a, b, model_md(md))
bb = json.loads((RP / "refinement/appendix_bbar/bbar.json").read_text())
j2 = json.loads((DATA / "j2_scores_H100_sq.json").read_text())
K = ("e_u_global", "e_jump", "e_traction")
for H, case in (("0.1", "$H=0.1$~MPa"), ("100", "$H=100$~MPa")):
    a = [[bb[f"{h}_{H}"]["succ"]["10"][k] for h in HS] for k in range(3)]
    b = [[bb[f"{h}_{H}"]["succ"]["20"][k] for h in HS] for k in range(3)]
    e = [[j2[f"psie/stiff/10/{H}/{h}"][K[k]] for h in HS] for k in range(3)]
    row("$J_2$, C3D8B" if H == "0.1" else "", case, a, b, e)
rc = json.loads(Path(sys.argv[1]).read_text())
for c, case in (("stiff_psi0", r"stiff, $\psi=0^\circ$"), ("stiff_dp", r"stiff, $\psi=10^\circ$"),
                ("initial_psi0", r"initial, $\psi=0^\circ$"), ("initial_dp", r"initial, $\psi=10^\circ$")):
    a = [[rc[f"{c}/{h}/5-10"][k] for h in HS] for k in range(3)]
    b = [[rc[f"{c}/{h}/10-20"][k] for h in HS] for k in range(3)]
    eie = json.loads((DATA / "dp" / f"{c}.json").read_text())
    e = [[eie[f"eie/10/{h}"][k] for h in HS] for k in range(3)]
    row("Drucker--Prager, C3D8" if c == "stiff_psi0" else "", case, a, b, e)
