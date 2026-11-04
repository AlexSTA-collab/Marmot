#!/usr/bin/env python3
"""Creep curves and deformed shapes of the viscoelastic analyses (Figure 12) from the 60 x 60 runs of make_jobs60.py.

The paper's own scripts draw them, unchanged except for where they read: paper/figures/reaction/{extract,plot}_reaction.py
for the creep of the loaded face at h = 0.01 and both inclinations, paper/figures/deformed/render_deformed.py for the
deformed configuration, stiff family, alpha = 10. They read the 26.05 EnSight layout, so the h = 0.01 runs are first
copied to that layout (port_2611/to_2605_layout.py) under as2605/.

    xvfb-run -a python3 visco_global60.py <out dir>     ->  viscoelastic_creep.pdf, deformed_viscoelastic.pdf
    RERUN60_JOBS=<jobs root> reads the runs from there (the shared-node reruns), converted to <root>_as2605."""
import os
import pickle
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
JOBS = Path(os.environ.get("RERUN60_JOBS", HERE / "jobs"))
AS2605 = HERE / "as2605" if JOBS == HERE / "jobs" else Path(str(JOBS) + "_as2605")
ROOT = HERE.parent
PAPER = ROOT.parent / "paper" / "figures"
OUT = Path(sys.argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
FAMS = ("stiff", "initial", "soft")
CASE = {0: "visco0", 10: "visco10"}

# 1. 26.05 layout of the h = 0.01 runs
for a, case in CASE.items():
    for kind in ("ifc", "ref"):
        for fam in FAMS:
            src = JOBS / case / kind / fam / f"angle_{a}" / "h_0.01"
            dst = AS2605 / case / kind / fam / f"angle_{a}" / "h_0.01"
            if not (dst / "run.log").exists():
                subprocess.run([sys.executable, str(ROOT / "port_2611" / "to_2605_layout.py"), str(src), str(dst)],
                               check=True, stdout=subprocess.DEVNULL)

# 2. creep curves: extract_reaction.curve() on the converted runs; the von Mises curves of the pickle are kept
ns = {"__file__": str(PAPER / "reaction" / "extract_reaction.py")}
src = (PAPER / "reaction" / "extract_reaction.py").read_text()
exec(compile(src[:src.index("\nVM = {")], "extract_reaction.py", "exec"), ns)
C = pickle.load(open(PAPER / "reaction" / "reaction_curves.pkl", "rb"))
C = {k: v for k, v in C.items() if k[0] != "ve"}
for a, case in CASE.items():
    for fam in FAMS:
        for tag, kind, pre in (("EIE", "ifc", "Cohesive_zone_model"), ("full", "ref", "Cauchy_Full_model")):
            c = ns["curve"](AS2605 / case / kind / fam / f"angle_{a}" / "h_0.01" / f"{pre}_top_0.01", a)
            assert c is not None, (a, fam, tag)
            C[("ve", a, fam, tag)] = c
            print(f"a={a:<2} {fam:<8} {tag:<4} steps={len(c[0])}  u: {c[0][1]:.4e} -> {c[0][-1]:.4e}  "
                  f"ratio {c[0][-1] / c[0][1]:.4f}", flush=True)
work = OUT / "_work"
work.mkdir(exist_ok=True)
pickle.dump(C, open(work / "reaction_curves.pkl", "wb"))
plot = (PAPER / "reaction" / "plot_reaction.py").read_text()
(work / "plot_reaction.py").write_text(plot)
subprocess.run([sys.executable, str(work / "plot_reaction.py"), str(OUT)], check=True)

# 3. deformed shapes, viscoelastic case only. render_deformed reads time set 1, which the converted .case files give
# the single geometry step; a .view.case per output in the 26.05 form (as biphasic/render_biphasic.view_case) fixes it
def view_case(run, stem):
    d = Path(run) / stem
    var = next(d.glob("displacement_*.var"), None)
    if var is None:      # an output without displacement: read_body skips it, as for the 26.05 runs
        return
    steps = var.read_bytes().count(b"BEGIN TIME STEP")
    values = "\n".join(f"{float(i):.8e}" for i in range(steps))
    (Path(run) / f"{stem}.view.case").write_text(
        f"FORMAT\ntype: ensight gold\nTIME\ntime set: 2 no description\nnumber of steps: 1\n"
        f"filename start number: 0\nfilename increment: 1\ntime values: 0.00000000e+00\n"
        f"time set: 1 no description\nnumber of steps: {steps}\nfilename start number: 0\n"
        f"filename increment: 1\ntime values: {values}\nFILE\nfile set: 2\nnumber of steps: 1\n"
        f"file set: 1\nnumber of steps: {steps}\nGEOMETRY\nmodel: 2 2 {stem}/geometry.geo\n"
        f"VARIABLE\nvector per node: 1 1 {var.stem} {stem}/{var.name}\n")


ve = AS2605 / "visco10"
for kind, pre in (("ref", "Cauchy_Full_model"), ("ifc", "Cohesive_zone_model")):
    for p in ("bottom", "top", "interface"):
        view_case(ve / kind / "stiff" / "angle_10" / "h_0.01", f"{pre}_{p}_0.01")
rd = (PAPER / "deformed" / "render_deformed.py").read_text()
k = rd.count('f"{prefix}_{p}_0.01.case"')
assert k == 1, k
rd = rd.replace('f"{prefix}_{p}_0.01.case"', 'f"{prefix}_{p}_0.01.view.case"')
rd = rd.replace('full=NF / "Parametric_Study_Marmot-nextv26.05_PowerLaw_Transition"\n                  / "stiff/angle_10/Full_layer_models"',
                f'full=Path("{ve}/ref/stiff/angle_10/h_0.01")')
rd = rd.replace('iface=NF / "Parametric_Study_Marmot-nextv26.05_PowerLaw_Transition"\n                   / "proper_convergence/stiff/YIQUAD4/h_0.01"',
                f'iface=Path("{ve}/ifc/stiff/angle_10/h_0.01")')
assert rd.count(str(ve)) == 2
rd = rd.replace("for tag, cfg in CASES.items():", "for tag, cfg in [(k, v) for k, v in CASES.items() if k == 'viscoelastic']:")
# The interface output of these runs carries no nodal displacement, so read_body skips the interface and the magnified
# jump opens a gap between the bodies through which the top face of the lower body shows. The interface elements share
# their nodes with the two body faces: draw them from their geometry with the displacement of the coinciding body node.
k = rd.count("\ndef warp(")
assert k == 1, k
rd = rd.replace("\ndef warp(", '''
_assemble = assemble
PORT = "''' + str(ROOT / "port_2611") + '''"


def assemble(folder, prefix, parts=("bottom", "top", "interface")):
    got = _assemble(folder, prefix, parts)
    case = Path(folder) / f"{prefix}_interface_0.01.case"
    if len(got) == 3 or not case.is_file():
        return got
    from scipy.spatial import cKDTree
    sys.path.insert(0, PORT)
    from ensight_reader_2611 import read_geo
    part = next(v for v in read_geo(str(Path(folder) / f"{prefix}_interface_0.01" / "geometry.geo")).values()
                if v["desc"] == "interface-body")
    m = pv.UnstructuredGrid({pv.CellType.HEXAHEDRON: part["blocks"][0][1]}, part["xyz"])
    pts = np.vstack([b.points for b in got])
    u = np.vstack([np.asarray(b.point_data["u"]) for b in got])
    d, i = cKDTree(pts).query(m.points)
    assert d.max() < 1e-6, d.max()
    m.point_data["u"] = u[i]
    m.point_data["umag"] = np.linalg.norm(u[i], axis=1)
    return got + [m]


def warp(''')
(work / "render_deformed.py").write_text(rd)
subprocess.run([sys.executable, str(work / "render_deformed.py"), str(OUT)], check=True)
for f in OUT.glob("vonmises_*.pdf"):
    f.unlink()
print("wrote", OUT / "viscoelastic_creep.pdf", OUT / "deformed_viscoelastic.pdf")
