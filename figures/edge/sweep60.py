#!/usr/bin/env python3
"""Paper error measures against the width of the strip removed at every edge of the interface, 60 x 60 meshes, alpha = 10,
h = 0.1 and 0.01. Strips of 0..8, 10..12 and 14..24 cells (margin (k - 1/2)/60, which removes exactly k rows of nodes,
up to 0.4 mm, leaving a 0.2 x 0.2 mm interior) and the paper margin 0.15 itself in place of 9 cells, so that this point is the value of the paper.
usage: sweep60.py <analysis> <out.json>
  elastic, visco       rerun60/jobs/{elastic,visco10}/{ifc,ref}            errors_end.measures (EIE)
  eie, psie            plastic_bbar/jobs/{eie,ref}, jobs_sq/psie, H = 100   errors_end.measures
  dp                   ~/psie26/plastic_bbar/dp/S75_*_dp (workstation, roller, H = 0.075 G, EIE on the reference's own body nodes), reference R075_*_dp/ref  errors_end.measures
  biphasic             biphasic/runs60_as2605/c100/{full,supp,ref}          biphasic/errors.py measures (node_matched_errors)
  biphasic_visco       biphasic_visco/runs60/creep/{full,supp,ref}, last output, score_visco.py measures
The roots of the shared-node reruns are given by EDGE_EV_JOBS (elastic/visco: <root>/{elastic,visco10}/{ifc,ref}),
EDGE_PSIE and EDGE_EIE (<root>/<fam>/angle_10/H_<EDGE_H, default 100>/h_*), EDGE_J2_REF (the J2 references, same layout), EDGE_BIPHASIC (<root>/{full,supp,ref}) and EDGE_BIPHASIC_VISCO."""
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPRO = Path.home() / "psie26" if (Path.home() / "psie26").exists() else HERE.parents[2] / "reproducibility"
REPRO = Path(os.environ.get("EDGE_REPRO", REPRO))   # the tree holding refinement/, port_2611/, biphasic/, biphasic_visco/
sys.path.insert(0, str(REPRO / "refinement")); sys.path.insert(0, str(Path.home() / "psie26/tools/refinement"))
import errors_end as E  # noqa: E402

CELL = 1.0 / 60.0
STRIPS = ([(k * CELL, 0.0 if k == 0 else (k - 0.5) * CELL) for k in range(9)] + [(0.15, 0.15)]
          + [(k * CELL, (k - 0.5) * CELL) for k in (10, 11, 12, 14, 16, 18, 20, 22, 24)])
HS = ("0.1", "0.01")
what, out = sys.argv[1], Path(sys.argv[2])
res = json.loads(out.read_text()) if out.exists() else {}


def save(key, val):
    res[key] = val
    print(key, val, flush=True)
    out.write_text(json.dumps(res, indent=1))


def sweep_measures(series, mod_ref):
    """series: labels; mod_ref(label, h) -> (model dir, reference dir)"""
    for s in series:
        for h in HS:
            m, r = mod_ref(s, h)
            E.ANGLE, E.H0 = 10, float(h)
            for x, mg in STRIPS:
                key = f"{s}|{h}|{x:.4f}"
                if key in res:
                    continue
                E.MARGIN = mg
                o = E.measures(str(m), h, str(r), h, 60)
                save(key, [o["e_u_global"], o["e_jump"], o["e_traction"]])


if what in ("elastic", "visco"):
    J = Path(os.environ.get("EDGE_EV_JOBS", REPRO / "rerun60" / "jobs")) / {"elastic": "elastic", "visco": "visco10"}[what]
    sweep_measures(("stiff", "initial", "soft"), lambda f, h: (J / "ifc" / f / "angle_10" / f"h_{h}",
                                                              J / "ref" / f / "angle_10" / f"h_{h}"))
elif what in ("eie", "psie"):
    B = REPRO / "plastic_bbar"
    M = {"eie": Path(os.environ.get("EDGE_EIE", B / "jobs" / "eie")), "psie": Path(os.environ.get("EDGE_PSIE", B / "jobs_sq" / "psie"))}[what]
    R = Path(os.environ.get("EDGE_J2_REF", B / "jobs" / "ref"))
    HH = "H_" + os.environ.get("EDGE_H", "100")   # EDGE_H=0.1: the J2 analyses of the revised paper
    sweep_measures(("stiff", "initial", "soft"), lambda f, h: (M / f / "angle_10" / HH / f"h_{h}",
                                                              R / f / "angle_10" / HH / f"h_{h}"))
elif what == "dp":
    D = Path.home() / "psie26/plastic_bbar/dp"
    sweep_measures(("stiff", "initial"), lambda f, h: (D / f"S75_{f}_dp/eie/angle_10/h_{h}",
                                                      D / f"R075_{f}_dp/ref/angle_10/h_{h}"))
elif what == "biphasic":
    sys.path.insert(0, str(REPRO / "port_2611" / "scripts"))
    import node_matched_errors as nme  # noqa: E402
    nme.W_AREA = np.sqrt(0.5) / 60
    ROOT = Path(os.environ.get("EDGE_BIPHASIC", REPRO / "biphasic" / "runs60_as2605" / "c100"))
    src = (REPRO / "biphasic" / "errors.py").read_text()
    ns = {"nme": nme, "np": np, "Path": Path}
    exec(src[src.index("def force_x"):src.index("def measures")], ns)   # force_x, layer_traction of the study
    layer_traction = ns["layer_traction"]
    for h in HS:
        ref = ROOT / "ref" / "angle_10" / f"h_{h}"
        refb = {s: nme.body(ref, s, h) for s in ("bottom", "top")}
        cr_, tr = layer_traction(ref, h, 10)
        jr, pr = nme.jump(refb["bottom"], refb["top"], h, 10)
        for model in ("full", "supp"):
            run = ROOT / model / "angle_10" / f"h_{h}"
            mod = {s: nme.body(run, s, h) for s in ("bottom", "top")}
            num = den = 0.0
            for s in ("bottom", "top"):
                (xr, cr, ur), (xm, cm, um) = refb[s], mod[s]
                at = nme.evaluate_at(xm, cm, um, xr, nme.matched(xr, xm))
                num += nme.integrate(np.sum((ur - at) ** 2, axis=1), xr, cr)
                den += nme.integrate(np.sum(ur ** 2, axis=1), xr, cr)
            jm, _ = nme.jump(mod["bottom"], mod["top"], h, 10, on=(refb["bottom"], refb["top"]))
            cm, tm = nme.model_traction(run, h)
            tm = tm[nme.matched(cr_[:, :2], cm[:, :2])]
            for x, mg in STRIPS:
                k, kt = nme.interior(pr, 10, mg), nme.interior(cr_, 10, mg)
                save(f"{model}|{h}|{x:.4f}", [float(np.sqrt(num / den)),
                                              float(nme.W_AREA * np.linalg.norm((jm - jr)[k]) / np.sqrt(den)),
                                              float(np.linalg.norm((tm - tr)[kt]) / np.linalg.norm(tr[kt]))])
elif what == "biphasic_visco":
    ROOT = Path(os.environ.get("EDGE_BIPHASIC_VISCO", REPRO / "biphasic_visco" / "runs60" / "creep"))
    src = (REPRO / "biphasic_visco" / "score_visco.py").read_text()
    for h in HS:
        E.ANGLE, E.H0 = 10, float(h)
        ns = {"E": E, "np": np, "h": h, "angle": 10}
        exec(src[src.index("def layer_traction"):src.index("res = {}")], ns)   # layer_traction of the study
        ref = ROOT / "ref" / "angle_10" / f"h_{h}"
        refb = {s: E.body(ref, "Cauchy_Full_model", s, h) for s in ("bottom", "top")}
        cr_, tr = ns["layer_traction"](ref)
        jr, pr = E.jump(refb["bottom"], refb["top"], E.H0, 10)
        for model in ("full", "supp"):
            run = ROOT / model / "angle_10" / f"h_{h}"
            mod = {s: E.body(run, "Cohesive_zone_model", s, h) for s in ("bottom", "top")}
            num = den = 0.0
            for s in ("bottom", "top"):
                (xr, cr, ur), (xm, cm_, um) = refb[s], mod[s]
                at = E.evaluate_at(xm, cm_, um, xr, E.matched(xr, xm))
                num += E.integrate(np.sum((ur - at) ** 2, axis=1), xr, cr)
                den += E.integrate(np.sum(ur ** 2, axis=1), xr, cr)
            jm, _ = E.jump(mod["bottom"], mod["top"], E.H0, 10, on=(refb["bottom"], refb["top"]))
            cm, tm = E.model_traction(run, h)
            tm = tm[E.matched(cr_[:, :2], cm[:, :2])]
            for x, mg in STRIPS:
                k, kt = E.interior(pr, 10, mg), E.interior(cr_, 10, mg)
                save(f"{model}|{h}|{x:.4f}", [float(np.sqrt(num / den)),
                                              float(np.sqrt(0.5) / 60 * np.linalg.norm((jm - jr)[k]) / np.sqrt(den)),
                                              float(np.linalg.norm((tm - tr)[kt]) / np.linalg.norm(tr[kt]))])
