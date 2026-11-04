#!/usr/bin/env python3
"""Convergence tables of the appendix: for every convergence figure of Section 4, the errors (e_u, e_[u], e_t) at
h = 0.1 and 0.01 mm and the slope fitted over the six thicknesses, at alpha = 0 and 10. Data, the files the
convergence figures are drawn from:
  data/elastic_table_angle{0,10}.md, data/visco0_table_angle0.md, data/visco10_table_angle10.md
      (reproducibility/rerun60/tables, written by rerun60/score60.py)
  data/biphasic_errors60.json   (reproducibility/biphasic/errors60.json, drawn by biphasic/plot_panels.py)
  data/j2_scores_H100_sq.json   (reproducibility/plastic_bbar/scores_H100_sq.json)
  data/dp/<family>_<rule>.json  (workstation ~/psie26/plastic_bbar/dp/fam_<family>_<rule>_cap/scores.json)
A rate is left out (--) where the six errors lie below 1e-7, i.e. at the solver tolerance or at round-off, and for the
soft Drucker--Prager layer at alpha = 0, which stays elastic.
usage: python3 make_rate_tables.py   (writes rates_{elastic,biphasic,visco,j2,dp}.tex)"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
HS = ("0.01", "0.02", "0.04", "0.06", "0.08", "0.1")
FAMS = ("stiff", "initial", "soft")
FLOOR = 1e-7


def sci(v):
    m, e = f"{v:.1e}".split("e")
    return f"${m}\\cdot10^{{{int(e)}}}$"


def cells(v, floor=False):
    r = np.polyfit(np.log([float(h) for h in HS]), np.log(v), 1)[0]
    return [sci(v[-1]), sci(v[0]), "--" if floor or max(v) < FLOOR else f"${r:.2f}$"]


HEAD = (r"\multicolumn{3}{c}{$e_u$} & \multicolumn{3}{c}{$e_{[u]}$} & \multicolumn{3}{c}{$e_t$}\\"
        "\n    %s\n"
        r"    %s & $h=0.1$ & $h=0.01$ & rate & $h=0.1$ & $h=0.01$ & rate & $h=0.1$ & $h=0.01$ & rate\\")


def table(labels, rows):
    n = len(labels)
    rules = "".join(f"\\cmidrule(lr){{{n + 1 + 3 * i}-{n + 3 + 3 * i}}}" for i in range(3))
    out = ["  \\begin{tabular}{" + "l" * n + "rrr" * 3 + "}", "    \\toprule",
           "    " + "& " * n + HEAD % (rules, " & ".join(labels)), "    \\midrule"]
    for r in rows:
        out.append("    \\midrule" if r is None else "    " + " & ".join(r) + r"\\")
    out += ["    \\bottomrule", "  \\end{tabular}"]
    return "\n".join(out) + "\n"


def md(path):
    """the three measures of a score60.py table, {measure index: {family: [six values]}}"""
    out, sec = {}, -1
    for line in path.read_text().splitlines():
        if line.startswith("## "):
            sec += 1
            out[sec] = {f: [] for f in FAMS}
        elif sec >= 0 and line.startswith("| 0"):
            c = [x.strip() for x in line.strip("|").split("|")]
            for f, v in zip(FAMS, c[1:]):
                out[sec][f].append(float(v))
    return out


def deg(a):
    return f"${a}^\\circ$"


# linear elastic and viscoelastic interphase, EIE
for name, files in (("elastic", ("elastic_table_angle0.md", "elastic_table_angle10.md")),
                    ("visco", ("visco0_table_angle0.md", "visco10_table_angle10.md"))):
    rows = []
    for a, f in zip((0, 10), files):
        if rows:
            rows.append(None)
        d = md(DATA / f)
        for fam in FAMS:
            r = [deg(a) if fam == "stiff" else "", fam]
            for k in range(3):
                # the soft elastic layer: displacements at the floor, traction falling to round-off -- no rate
                r += cells(d[k][fam], floor=(name == "elastic" and fam == "soft"))
            rows.append(r)
    (HERE / f"rates_{name}.tex").write_text(table([r"$\alpha$", "family"], rows))

# biphasic elastic interphase, mismatch kept (full) and removed
d = json.loads((DATA / "biphasic_errors60.json").read_text())
rows = []
for a in (0, 10):
    if rows:
        rows.append(None)
    for case, ratio in (("c100", "100"), ("c10", "10")):
        for model, label in (("full", "kept"), ("supp", "removed")):
            r = [deg(a) if (case, model) == ("c100", "full") else "", f"${ratio}$" if model == "full" else "", label]
            for k in ("global", "e_u", "e_t"):    # e_u of errors60.json is the displacement jump
                r += cells([d[f"{case}/{a}/{h}"][model][k] for h in HS])
            rows.append(r)
(HERE / "rates_biphasic.tex").write_text(table([r"$\alpha$", r"$E^{+}/E^{-}$", r"$[C^{e}_{ijkl}]$"], rows))

# J_2 interphase, H = 100 MPa, EIE and PSIE with the squeeze-flow term
d = json.loads((DATA / "j2_scores_H100_sq.json").read_text())
rows = []
for a in (0, 10):
    for el, name in (("eie", "EIE"), ("psie", "PSIE")):
        if rows:
            rows.append(None)
        for fam in FAMS:
            r = [deg(a) if (el, fam) == ("eie", "stiff") else "", name if fam == "stiff" else "", fam]
            for k in ("e_u_global", "e_jump", "e_traction"):
                r += cells([d[f"{el}/{fam}/{a}/100/{h}"][k] for h in HS], floor=(fam == "soft" and (k != "e_traction" or el == "psie")))
            rows.append(r)
(HERE / "rates_j2.tex").write_text(table([r"$\alpha$", "element", "family"], rows))

# Drucker--Prager interphase, EIE, the two flow rules the paper reports
RULES = (("psi0", r"$30^\circ,0^\circ$"), ("dp", r"$30^\circ,10^\circ$"))
rows = []
for a in (0, 10):
    for fam in FAMS:
        if rows:
            rows.append(None)
        for rule, name in RULES:
            s = json.loads((DATA / "dp" / f"{fam}_{rule}.json").read_text())
            r = [deg(a) if (fam, rule) == ("stiff", "psi0") else "", name, fam]
            for i in range(3):
                r += cells([s[f"eie/{a}/{h}"][i] for h in HS], floor=(fam == "soft"))   # the soft layer stays elastic: no rate
            rows.append(r)
(HERE / "rates_dp.tex").write_text(table([r"$\alpha$", r"$\varphi,\psi$", "family"], rows))
print("wrote rates_elastic.tex, rates_biphasic.tex, rates_visco.tex, rates_j2.tex, rates_dp.tex")
