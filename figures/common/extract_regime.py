#!/usr/bin/env python3
"""Regime indicator for the viscoelastic interphase.

The claim under test is the one made after equation (viscoelastic limits): because the
adherents are elastic and the interphase is not, the stiffness contrast that selects the
soft or the stiff regime is a function of the observation time rather than a fixed
property of the problem.

The indicator is the share of the assembly's relative motion that the interface takes up,

    chi = < [u]_x >  /  < u_x >_loaded face ,

with the far face fully fixed, so chi -> 0 is a layer stiff enough that the deformation
lives in the bodies (stiff regime, membrane-dominated) and chi -> 1 is a layer that takes
all of it (soft regime, cohesive-dominated). It is measured identically on the interface
model and on the resolved Cauchy layer; for the latter the "jump" is the relative
displacement across the resolved layer.

Two experiments are read:

  dt_sweep/  single-increment loading at step lengths spanning 1e-8 to 1e4 days. This
             sweeps the algorithmic tangent C^ev(Delta t) across both limits of
             (viscoelastic limits) with no creep accumulated, so it isolates the tangent.
  creep/     the standard protocol, instantaneous load then a 100-day hold at fixed
             Delta t = 1 day. Here the tangent is constant and the migration is driven by
             the accumulated compliance instead.

Writes regime_curves.pkl.
"""
import pickle
import sys
from pathlib import Path

import numpy as np

NF = Path("/home/alexsta1993/Documents/STATHAS_VIENNA_POSTDOC/Applications_Examples/"
          "EdelweisModels/Cohesive_Zone_models/3D/Three_Body_Cubit/"
          "Comparizon_LinearElastic/New_Functionality")
STUDY = NF / "Parametric_Study_Marmot-nextv26.05_PowerLaw_Transition"
RT = STUDY / "dt_full"

FAMILIES = ["stiff", "initial", "soft"]
HEIGHTS = ["0.01", "0.02", "0.04", "0.06", "0.08", "0.1"]
DT = [1e-8, 1e-6, 1e-4, 1e-2, 1e0, 1e2, 1e4]

sys.path.insert(0, str(NF / "Parametric_Study_Marmot-nextv26.05_Von_Mises_GaussLobatto"
                       / "evaluation_scripts"))
from ensight_reader import read_geo, read_var_node_all_steps      # noqa: E402

PREFIX = {"interface": "Cohesive_zone_model", "full": "Cauchy_Full_model"}
TOL = 1e-9


def _surface(xyz, side):
    """Indices of the interface-facing surface: the extreme z of each (x, y) column.

    The interface is inclined, so this is not a plane of constant z; taking the extreme
    per column follows it. The mesh is structured, so the columns are exact.
    """
    key = np.round(xyz[:, :2] / TOL).astype(np.int64)
    order = {}
    for i, k in enumerate(map(tuple, key)):
        z = xyz[i, 2]
        if k not in order or (z > xyz[order[k], 2] if side == "max" else z < xyz[order[k], 2]):
            order[k] = i
    return order


def _body(run, model, which, h):
    """(xyz, displacement over all steps) for one body of one run."""
    d = run / f"{PREFIX[model]}_{which}_{h}"
    parts = read_geo(str(d / "geometry.geo"))
    disp = read_var_node_all_steps(str(d / f"displacement_{which}.var"), parts, 3)
    name = next(p for p, v in parts.items() if v["desc"] == f"{which}-body" and p in disp)
    return parts[name]["xyz"], disp[name]


def times(run, model, h):
    """Step times, read from the EnSight case file rather than assumed."""
    txt = (run / f"{PREFIX[model]}_top_{h}.case").read_text()
    block = txt.split("time set: 1")[1].split("time values:")[1]
    vals = []
    for tok in block.split():
        try:
            vals.append(float(tok))
        except ValueError:
            break
    return np.array(vals)


def _grid(xyz, idx):
    """Order a structured interface surface into an (ni, nj) index grid."""
    key = np.round(xyz[idx][:, :2] / TOL).astype(np.int64)
    xs = np.unique(key[:, 0])
    ys = np.unique(key[:, 1])
    g = np.full((len(xs), len(ys)), -1, dtype=int)
    xi = {v: i for i, v in enumerate(xs)}
    yi = {v: j for j, v in enumerate(ys)}
    for n, (kx, ky) in enumerate(map(tuple, key)):
        g[xi[kx], yi[ky]] = n
    return g


def tangential_sym_grad(xyz, u, idx):
    """Symmetric part of the surface gradient of u on a planar interface surface.

    The two faces of the layer are planar and meshed on a structured grid, so the
    position is affine in the grid indices. With a_1, a_2 the (constant) covariant basis
    vectors dr/di and dr/dj, and b_a their dual (b_a . a_c = delta_ac, both in the tangent
    plane), the surface gradient is

        G_ij = sum_a (du_i/dxi_a) (b_a)_j ,

    which is the part of u_{i,j} that survives differentiation along the interface --- the
    membrane deformation. The complementary through-thickness part is what the jump
    carries. Returns u^S_{ij} = (G_ij + G_ji)/2 at every node of the surface.
    """
    g = _grid(xyz, idx)
    P = xyz[idx][g]                              # (ni, nj, 3) positions
    U = u[:, idx, :][:, g, :]                    # (nsteps, ni, nj, 3) displacements
    a1 = (P[1:, :, :] - P[:-1, :, :]).reshape(-1, 3).mean(axis=0)
    a2 = (P[:, 1:, :] - P[:, :-1, :]).reshape(-1, 3).mean(axis=0)
    A = np.column_stack([a1, a2])                                  # 3 x 2
    B = A @ np.linalg.inv(A.T @ A)                                 # dual basis, 3 x 2
    d1 = np.gradient(U, axis=1)                  # du/dxi_1
    d2 = np.gradient(U, axis=2)                  # du/dxi_2
    G = (d1[..., :, None] * B[:, 0][None, None, None, :]
         + d2[..., :, None] * B[:, 1][None, None, None, :])        # (nsteps, ni, nj, 3, 3)
    return 0.5 * (G + np.swapaxes(G, -1, -2))


def indicator(run, model, h):
    """Jump measured against the membrane deformation of the layer.

    Returns  ||[u_i]|| / (h ||<u^S_{i,j}>||),  L2 over the interface at every step: the
    cohesive contribution divided by the deformation the layer would carry over its own
    thickness if it followed the average surface strain. Below one the membrane term
    dominates, above one the cohesive term does.
    """
    xyz_b, u_b = _body(run, model, "bottom", h)
    xyz_t, u_t = _body(run, model, "top", h)
    top_of_bottom = _surface(xyz_b, "max")
    bot_of_top = _surface(xyz_t, "min")
    shared = sorted(set(top_of_bottom) & set(bot_of_top))
    if not shared:
        raise RuntimeError(f"no paired interface nodes in {run}")
    ib = [top_of_bottom[k] for k in shared]
    it = [bot_of_top[k] for k in shared]

    jump = u_t[:, it, :] - u_b[:, ib, :]                       # (nsteps, npair, 3)
    Sb = tangential_sym_grad(xyz_b, u_b, ib)
    St = tangential_sym_grad(xyz_t, u_t, it)
    Savg = 0.5 * (Sb + St)                                     # <u^S_{i,j}>

    num = np.sqrt((jump ** 2).sum(axis=(1, 2)))
    den = float(h) * np.sqrt((Savg ** 2).sum(axis=(1, 2, 3, 4)))
    with np.errstate(divide="ignore", invalid="ignore"):
        return num / den, num, den, len(shared)


def main():
    out = {}
    print("=== tangent sweep: one increment, over families and thicknesses ===")
    print(f"{'family':>8} {'h':>6} " + "".join(f"{d:>10.0e}" for d in DT) + "   max dev vs resolved")
    for fam in FAMILIES:
        for h in HEIGHTS:
            row, devs = {}, []
            for dt in DT:
                pair = {}
                for model in ("interface", "full"):
                    run = RT / model / fam / f"h_{h}" / f"dt_{dt:g}"
                    if not (run / f"{PREFIX[model]}_top_{h}").is_dir():
                        continue
                    c, j, tot, n = indicator(run, model, h)
                    pair[model] = float(c[-1])
                if pair:
                    row[dt] = pair
                    if len(pair) == 2 and pair["full"]:
                        devs.append(abs(pair["interface"] - pair["full"]) / abs(pair["full"]))
            if row:
                out[("dtfull", fam, h)] = row
                vals = "".join(f"{row[d]['interface']:>10.4f}" if d in row else f"{'--':>10}"
                               for d in DT)
                print(f"{fam:>8} {h:>6} {vals}   {max(devs)*100 if devs else float('nan'):8.2f} %")

    print("\n=== creep hold, same grid ===")
    for fam in FAMILIES:
        for h in HEIGHTS:
            iface = STUDY / "proper_convergence" / fam / "YIQUAD4" / f"h_{h}"
            ref = (STUDY / fam / "angle_10" / "Full_layer_models" if h == "0.01"
                   else STUDY / "proper_heff_check" / fam / f"h_{h}")
            rec = {}
            for model, run in (("interface", iface), ("full", ref)):
                if not (run / f"{PREFIX[model]}_top_{h}").is_dir():
                    continue
                c, j, tot, n = indicator(run, model, h)
                rec[model] = dict(chi=c, t=times(run, model, h))
            if rec:
                out[("creep", fam, h)] = rec
                c = rec["interface"]["chi"]
                print(f"  {fam:>8} h={h:<5} chi: {c[1]:.4f} -> {c[-1]:.4f}")

    (Path(__file__).parent / "regime_curves.pkl").write_bytes(pickle.dumps(out))
    print(f"\nwrote {Path(__file__).parent / 'regime_curves.pkl'}")


if __name__ == "__main__":
    main()
