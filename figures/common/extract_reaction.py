#!/usr/bin/env python3
"""Total reaction on the loaded face, for comparing global stiffness.

The loaded face is the top of the upper body; its reaction is summed over all
its nodes and projected on the loading direction, which is inclined with the
model by the interface angle. Follows the extraction already used in
analyse_bubble_ablation.py so the numbers are directly comparable.
"""
import math, sys, pickle
from pathlib import Path
import numpy as np

NF = Path("/home/alexsta1993/Documents/STATHAS_VIENNA_POSTDOC/Applications_Examples/"
          "EdelweisModels/Cohesive_Zone_models/3D/Three_Body_Cubit/"
          "Comparizon_LinearElastic/New_Functionality")
sys.path.insert(0, str(NF / "Parametric_Study_Marmot-nextv26.05_Von_Mises_GaussLobatto"
                       / "evaluation_scripts"))
from ensight_reader import read_geo, read_var_node_all_steps   # noqa: E402


def curve(outdir, angle, which="top"):
    """(time index, mean face displacement, total face reaction) along the load."""
    d = Path(outdir)
    if not (d / "geometry.geo").is_file():
        return None
    n = np.array([math.cos(math.radians(angle)), 0.0, -math.sin(math.radians(angle))])
    parts = read_geo(str(d / "geometry.geo"))
    disp = read_var_node_all_steps(str(d / f"displacement_{which}.var"), parts, 3)
    reac = read_var_node_all_steps(str(d / f"reaction_{which}.var"), parts, 3)
    tgt = next(p for p, v in parts.items() if v["desc"] == f"{which}-body" and p in disp)
    xyz = parts[tgt]["xyz"]
    face = np.where(xyz[:, 2] > xyz[:, 2].max() - 1e-9)[0]
    u = (disp[tgt][:, face, :] @ n).mean(axis=1)
    F = -(reac[tgt][:, face, :] @ n).sum(axis=1)
    return u, F


VM = {
    "full":     "Parametric_Study_Marmot-nextv26.05_Von_Mises_ExtendedInterfaceElement/"
                "{fam}/angle_{a}/Full_layer_models/Cauchy_Full_model_top_0.01",
    "EIE":      "Parametric_Study_Marmot-nextv26.05_Von_Mises_ExtendedInterfaceElement/"
                "YIQUAD4_runs/{fam}/angle_{a}/h_0.01/Cohesive_zone_model_top_0.01",
    "PSIE":     "Parametric_Study_Marmot-nextv26.05_Von_Mises_StabPMini/"
                "{fam}/angle_{a}/Interface_models/Cohesive_zone_model_top_0.01",
}
# The viscoelastic results come from the transition study (m = 5e-3). The Warping study
# (m = 2.5e-1) is superseded: its creep compliance is large enough that the interphase
# tangent starts below both adherents, so no regime migration is observable in it.
VE_STUDY = "Parametric_Study_Marmot-nextv26.05_PowerLaw_Transition"
VE = (VE_STUDY + "/proper_convergence{suf}/{fam}/YIQUAD4/h_0.01/Cohesive_zone_model_top_0.01")
VE_FULL = (VE_STUDY + "/{fam}/angle_{a}/Full_layer_models/Cauchy_Full_model_top_0.01")

out = {}
print("=== von Mises, h = 0.01 ===")
for a in (0, 10):
    for fam in ("stiff", "initial", "soft"):
        for tag, tpl in VM.items():
            c = curve(NF / tpl.format(fam=fam, a=a), a)
            if c is None:
                continue
            out[("vm", a, fam, tag)] = c
            print(f"  a={a:<3} {fam:<8} {tag:<8} steps={len(c[0]):<4} "
                  f"u_end={c[0][-1]:.4e}  F_end={c[1][-1]:.6e}")

print("\n=== viscoelastic (power law), h = 0.01 ===")
for a in (0, 10):
    suf = "" if a == 10 else "_angle0"
    for fam in ("stiff", "initial", "soft"):
        c = curve(NF / VE.format(suf=suf, fam=fam), a)
        if c is not None:
            out[("ve", a, fam, "EIE")] = c
            print(f"  a={a:<3} {fam:<8} iface  steps={len(c[0]):<4} "
                  f"u: {c[0][0]:.3e}->{c[0][-1]:.3e}   F: {c[1][0]:.4e}->{c[1][-1]:.4e}")
        cf = curve(NF / VE_FULL.format(fam=fam, a=a), a)
        if cf is not None:
            out[("ve", a, fam, "full")] = cf
            print(f"  a={a:<3} {fam:<8} full   steps={len(cf[0]):<4} "
                  f"u: {cf[0][0]:.3e}->{cf[0][-1]:.3e}   F: {cf[1][0]:.4e}->{cf[1][-1]:.4e}")

pickle.dump(out, open(Path(__file__).with_name("reaction_curves.pkl"), "wb"))
print("\nsaved reaction_curves.pkl")
