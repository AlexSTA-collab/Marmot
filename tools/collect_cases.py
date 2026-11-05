#!/usr/bin/env python3
"""Collect the job file of every analysis behind the figures and tables of the IJES revision from the workstation tree
(~/psie26) into cases/, with the mesh includes rewritten to @IJES_ROOT@/meshes/<mesh set>/angle_<a>/height_<h>/<file>,
and write cases/index.tsv. Run on the workstation; the paths below are where the analyses ran.
usage: collect_cases.py <psie26 root> <out dir>"""
import csv
import json
import re
import sys
from pathlib import Path

W, OUT = Path(sys.argv[1]), Path(sys.argv[2])
MESH = [  # (path fragment as included, mesh set in the package)
    ("refinement/mesh/L1nt20", "L1nt20"), ("refinement/mesh/L1nt10", "L1nt10"), ("refinement/mesh/L1nt5", "L1nt5"),
    ("refinement/mesh/L1", "L1"), ("refinement/mesh_if/L-1", "if_Lm1"), ("refinement/mesh_if/L0", "if_L0"),
    ("refinement/mesh_if_flat/L-1", "flat_Lm1"), ("refinement/mesh_if_flat/L0", "flat_L0"),
    ("refinement/mesh_if_flat/L1", "flat_L1"), ("refinement/mesh_if_flat/L2", "flat_L2"),
    ("refinement/mesh_if/L1", "if_L1"), ("refinement/mesh_if/L2", "if_L2"), ("plastic_bbar/dp/S_mesh", "S_mesh"),
    ("mesh/full_layer_60", "full_layer_60")]
INC = re.compile(r"input=(\S*?)(" + "|".join(re.escape(f) for f, _ in MESH) + r")/(angle_\d+/height_[0-9.]+/\S+)")
SETS = dict(MESH)
FAMS, HS = ("stiff", "initial", "soft"), ("0.01", "0.02", "0.04", "0.06", "0.08", "0.1")
used = set()


def rewrite(text, src):
    def sub(m):
        used.add(SETS[m.group(2)]); return f"input=@IJES_ROOT@/meshes/{SETS[m.group(2)]}/{m.group(3)}"
    t, k = INC.subn(sub, text)
    rest = [l for l in t.splitlines() if "input=" in l and "@IJES_ROOT@" not in l]
    assert k and not rest, (src, rest[:2])
    return t


rows = []


def add(case, study, law, fam, angle, h, model, src, reference="", feeds="", env=""):
    src = W / src
    if not (src / "job.inp").exists():
        print("missing", src); return
    dst = OUT / case; dst.mkdir(parents=True, exist_ok=True)
    (dst / "job.inp").write_text(rewrite((src / "job.inp").read_text(), src))
    rows.append(dict(case=case, study=study, law=law, family=fam, angle=angle, h=h, model=model,
                     job=f"{case}/job.inp", reference=reference, env=env, feeds=feeds, source=str(src.relative_to(W))))


# elastic and viscoelastic, shared-node interface mesh, references with 20 C3D8 through the layer
for study, angles, feeds in (("elastic", (0, 10), "Figs 6-7, Table J.2, Figs H.29-H.30, I.39"),
                             ("visco0", (0,), "Fig 11, Table J.4, Fig 13, Fig H.31"),
                             ("visco10", (10,), "Fig 12, Table J.4, Fig 13, Fig G.28, Fig H.32, I.39")):
    for f in FAMS:
        for a in angles:
            for h in HS:
                r = f"{study}/ref/{f}/a{a}/h{h}"
                add(r, study, study.rstrip("010"), f, a, h, "ref", f"smesh_refs/{study}/ref/{f}/angle_{a}/h_{h}", feeds=feeds)
                add(f"{study}/eie/{f}/a{a}/h{h}", study, study.rstrip("010"), f, a, h, "eie",
                    f"smesh/{study}/ifc/{f}/angle_{a}/h_{h}", r, feeds)
# J2 at H = 0.1 MPa: references (C3D8B), EIE (relative tolerance of the return mapping) and PSIE on the shared-node mesh
J2F = "Figs 15-16, 18-20, Table J.5, Figs H.33-H.36, I.41"
for f in FAMS:
    for a in (0, 10):
        for h in HS:
            r = f"j2/ref/{f}/a{a}/H0.1/h{h}"
            add(r, "j2", "j2", f, a, h, "ref", (f"smesh_refs/plastic_bbar/jobs/ref/{f}/angle_{a}/H_0.1/h_{h}" if (f, a) == ("stiff", 10)
                                                 else f"h01/ref/{f}/angle_{a}/h_{h}"), feeds=J2F)
            add(f"j2/eie/{f}/a{a}/H0.1/h{h}", "j2", "j2", f, a, h, "eie", f"h01/eie/{f}/angle_{a}/h_{h}", r, J2F,
                env="MARMOT_VONMISES_RELTOL=1")
            add(f"j2/psie/{f}/a{a}/H0.1/h{h}", "j2", "j2", f, a, h, "psie", (f"smesh_plain/psie/{f}/angle_{a}/H_0.1/h_{h}"
                                                                          if (f, a) == ("stiff", 10) else f"h01/psie/{f}/angle_{a}/h_{h}"), r, J2F)
# Figure 21: the EIE run of j2/eie/stiff/a10/H0.1/h0.01 with the pressure of its face stress as an extra output
add("j2/pressure/eie", "j2", "j2", "stiff", 10, "0.01", "eie", "fig21_h01/eie", "j2/ref/stiff/a10/H0.1/h0.01", "Fig 21",
    env="MARMOT_VONMISES_RELTOL=1")
for h in ("0.1", "0.01"):
    r = f"cyclic/ref/H0.1/h{h}"
    add(r, "cyclic", "j2", "stiff", 10, h, "ref", f"optA/cyclic/ref/H_0.1/h_{h}", feeds="Fig 23")
    add(f"cyclic/psie/H0.1/h{h}", "cyclic", "j2", "stiff", 10, h, "psie", f"smesh_plain/cyclic/psie/H_0.1/h_{h}", r, "Fig 23")
for H in ("0.1",):   # through-thickness convergence of the J2 reference (reply Table R1)
    for n in (5, 10, 20):
        for d, h in (("bb_ref", "0.1"), ("bb01_ref", "0.01"), ("bbhm02_ref", "0.02"), ("bbhm04_ref", "0.04"),
                     ("bbhm06_ref", "0.06"), ("bbhm08_ref", "0.08")):
            add(f"j2_refconv/n{n}/H{H}/h{h}", "j2_refconv", "j2", "stiff", 10, h, f"ref_n{n}", f"optA/bb_refs/{d}{n}_H{H}",
                feeds="reply Table R1")
# Drucker-Prager: EIE on the shared-node mesh, references with 20 C3D8, and the n = 5, 10 references of the reply
for f in FAMS:
    for rule, psi in (("dp", 10),):
        for a in (0, 10):
            for h in HS:
                r = f"dp/ref/{f}_{rule}/a{a}/h{h}"
                add(r, "dp", "dp", f"{f}_psi{psi}", a, h, "ref", f"plastic_bbar/dp/R075_{f}_{rule}/ref/angle_{a}/h_{h}",
                    feeds="Figs 24-26, Table J.7, Figs H.37-H.38, I.41")
                add(f"dp/eie/{f}_{rule}/a{a}/h{h}", "dp", "dp", f"{f}_psi{psi}", a, h, "eie",
                    f"plastic_bbar/dp/S75_{f}_{rule}/eie/angle_{a}/h_{h}", r, "Figs 24-26, Table J.7, Figs H.37-H.38, I.41")
        if f != "soft":
            for n in (5, 10):
                for h in HS:
                    add(f"dp_refconv/n{n}/{f}_{rule}/h{h}", "dp_refconv", "dp", f"{f}_psi{psi}", 10, h, f"ref_n{n}",
                        f"plastic_bbar/dp/S75_{f}_{rule}/refconv/n{n}/angle_10/h_{h}", feeds="reply Table R1")
# biphasic elastic and viscoelastic layers
for c in ("c10", "c100"):
    for a in (0, 10):
        for h in HS:
            r = f"biphasic/ref/{c}/a{a}/h{h}"
            add(r, "biphasic", "biphasic", c, a, h, "ref", f"optA/biphasic/{c}/ref/angle_{a}/h_{h}", feeds="Figs 8-10, Table J.3, Fig I.40")
            for m in ("full", "supp"):
                add(f"biphasic/{m}/{c}/a{a}/h{h}", "biphasic", "biphasic", c, a, h, m,
                    f"optA/biphasic/{c}/{m}/angle_{a}/h_{h}", r, "Figs 8-10, Table J.3, Fig I.40")
for h in ("0.1", "0.01"):
    r = f"biphasic_visco/ref/a10/h{h}"
    add(r, "biphasic_visco", "biphasic_visco", "c100", 10, h, "ref", f"optA/biphasic_visco/ref/angle_10/h_{h}", feeds="Fig 14, Fig I.40")
    for m in ("full", "supp"):
        add(f"biphasic_visco/{m}/a10/h{h}", "biphasic_visco", "biphasic_visco", "c100", 10, h, m,
            f"optA/biphasic_visco/{m}/angle_10/h_{h}", r, "Fig 14, Fig I.40")
# regime of the viscoelastic interface model against the step length
for f in FAMS:
    for d in sorted((W / "optA/regime" / f).glob("h_*/dt_*")):
        h, dt = d.parent.name[2:], d.name[3:]
        add(f"regime/{f}/h{h}/dt{dt}", "regime", "visco", f, 10, h, "eie", str(d.relative_to(W)), feeds="Fig F.27")
# in-plane refinement at h = 0.01, H = 0.1 (Figure 22, Table J.6): each element on the flat interface meshes N x N, against its
# own solution on the 120 x 120 mesh
LEVEL = {"L-1": "15", "L0": "30", "L1": "60", "L2": "120"}
for lev, n in LEVEL.items():
    add(f"inplane/eie/H0.1/{lev}", "inplane", "j2", "stiff", 10, "0.01", "eie", f"optC/eie/H_0.1/{lev}", "", "Fig 22, Table J.6",
        env="MARMOT_VONMISES_RELTOL=1")
    add(f"inplane/psie/H0.1/{lev}", "inplane", "j2", "stiff", 10, "0.01", "psie", f"optC/psie/H_0.1/{lev}", "", "Fig 22, Table J.6")
# computational cost (Table 1), stiff family, alpha = 10, h = 0.01, each run alone with four threads: the elastic, viscoelastic
# and Drucker-Prager rows are the runs of 27 and 31 October 2026 (bench/; the elastic and viscoelastic interface models ran on
# the interface mesh mesh_if/L1 of that date, kept for this table only), the J2 rows the runs at H = 0.1 of 4-5 November (bench01/)
for row, models, root in (("elastic", ("resolved", "EIE"), "bench"), ("visco", ("resolved", "EIE"), "bench"),
                          ("dp", ("resolved", "EIE"), "bench"), ("plastic", ("resolved", "EIE", "PSIE"), "bench01"),
                          ("cyclic", ("resolved", "PSIE"), "bench01")):
    for m in models:
        add(f"bench/{row}/{m}", "bench", row, "stiff", 10, "0.01", m, f"{root}/{row}/{m}", feeds="Table 1",
            env="MARMOT_VONMISES_RELTOL=1" if (row, m) == ("plastic", "EIE") else "")
# through-thickness convergence of the elastic and viscoelastic references (reply)
for law in ("elastic", "visco"):
    for d in sorted((W / "refinement_ev60/jobs" / law).glob("*/*")):
        if (d / "job.inp").exists():
            add(f"ev_refconv/{law}/{d.parent.name}/{d.name}", "ev_refconv", law, "stiff", 10, d.name.split("_")[-1],
                d.parent.name, str(d.relative_to(W)), feeds="reply Table R1")

OUT.mkdir(parents=True, exist_ok=True)
with open(OUT / "index.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
print(len(rows), "cases; mesh sets", sorted(used))
(OUT / "mesh_sets.json").write_text(json.dumps(sorted(used)))
