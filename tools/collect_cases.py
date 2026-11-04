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
for study, angles, feeds in (("elastic", (0, 10), "Figs 6-7, Table L.3, Figs I.32-I.33, Fig K.45"),
                             ("visco0", (0,), "Fig 11, Table L.5, Fig 13, Fig I.34"),
                             ("visco10", (10,), "Fig 12, Table L.5, Fig 13, Fig G.27, Fig I.35, Fig K.45")):
    for f in FAMS:
        for a in angles:
            for h in HS:
                r = f"{study}/ref/{f}/a{a}/h{h}"
                add(r, study, study.rstrip("010"), f, a, h, "ref", f"smesh_refs/{study}/ref/{f}/angle_{a}/h_{h}", feeds=feeds)
                add(f"{study}/eie/{f}/a{a}/h{h}", study, study.rstrip("010"), f, a, h, "eie",
                    f"smesh/{study}/ifc/{f}/angle_{a}/h_{h}", r, feeds)
# J2: references (C3D8B), PSIE with the squeeze-flow term on the shared-node mesh, EIE on the interface mesh of the paper
for f in FAMS:
    for a in (0, 10):
        for H in ("100", "0.1") if (f, a) == ("stiff", 10) else ("100",):   # H = 0.1 only for the stiff family at 10 degrees
            for h in HS:
                r = f"j2/ref/{f}/a{a}/H{H}/h{h}"
                add(r, "j2", "j2", f, a, h, "ref", f"smesh_refs/plastic_bbar/jobs/ref/{f}/angle_{a}/H_{H}/h_{h}",
                    feeds="Figs 15-16, 18-21, Table L.6, Figs H.29-H.30, I.36-I.39, K.46")
                add(f"j2/psie/{f}/a{a}/H{H}/h{h}", "j2", "j2", f, a, h, "psie", f"smesh/j2/psie/{f}/angle_{a}/H_{H}/h_{h}", r,
                    "Figs 18-21, Table L.6, Figs H.29-H.30, I.38-I.39, K.46")
                if H == "100":   # at H = 0.1 the EIE of the paper is the run with the relative tolerance below
                    add(f"j2/eie/{f}/a{a}/H{H}/h{h}", "j2", "j2", f, a, h, "eie", f"optA/j2_eie_oldmesh/{f}/angle_{a}/H_{H}/h_{h}", r,
                        "Figs 15-16, 20, Table L.6, Figs I.36-I.37, K.46")
add("j2/eie_reltol/stiff/a10/H0.1/h0.01", "j2", "j2", "stiff", 10, "0.01", "eie", "optA/j2_eie_reltol/stiff/angle_10/H_0.1/h_0.01",
    "j2/ref/stiff/a10/H0.1/h0.01", "Fig H.30", env="MARMOT_VONMISES_RELTOL=1")
add("j2/oscillation/eie", "j2", "j2", "stiff", 10, "0.01", "eie", "optA/osc100_oldmesh/eie", "j2/ref/stiff/a10/H100/h0.01", "Fig 21")
for H in ("100", "0.1"):
    for h in ("0.1", "0.01"):
        r = f"cyclic/ref/H{H}/h{h}"
        feeds = "Fig 22, Table M.8" if H == "100" else "Fig H.31"
        add(r, "cyclic", "j2", "stiff", 10, h, "ref", f"optA/cyclic/ref/H_{H}/h_{h}", feeds=feeds)
        add(f"cyclic/psie/H{H}/h{h}", "cyclic", "j2", "stiff", 10, h, "psie", f"optA/cyclic/psie/H_{H}/h_{h}", r, feeds)
for H in ("100", "0.1"):
    for n in (5, 10, 20):
        for d, h in (("bb_ref", "0.1"), ("bb01_ref", "0.01"), ("bbhm02_ref", "0.02"), ("bbhm04_ref", "0.04"),
                     ("bbhm06_ref", "0.06"), ("bbhm08_ref", "0.08")):
            add(f"j2_refconv/n{n}/H{H}/h{h}", "j2_refconv", "j2", "stiff", 10, h, f"ref_n{n}", f"optA/bb_refs/{d}{n}_H{H}",
                feeds="Fig H.28; reply Table R1, Figure R2")
# Drucker-Prager: EIE on the shared-node mesh, references with 20 C3D8, and the n = 5, 10 references of the reply
for f in FAMS:
    for rule, psi in (("psi0", 0), ("dp", 10)):
        for a in (0, 10):
            for h in HS:
                r = f"dp/ref/{f}_{rule}/a{a}/h{h}"
                add(r, "dp", "dp", f"{f}_psi{psi}", a, h, "ref", f"plastic_bbar/dp/R075_{f}_{rule}/ref/angle_{a}/h_{h}",
                    feeds="Figs 23-25, Table L.7, Figs I.40-I.41, K.46")
                add(f"dp/eie/{f}_{rule}/a{a}/h{h}", "dp", "dp", f"{f}_psi{psi}", a, h, "eie",
                    f"plastic_bbar/dp/S75_{f}_{rule}/eie/angle_{a}/h_{h}", r, "Figs 23-25, Table L.7, Figs I.40-I.41, K.46")
        if f != "soft":
            for n in (5, 10):
                for h in HS:
                    add(f"dp_refconv/n{n}/{f}_{rule}/h{h}", "dp_refconv", "dp", f"{f}_psi{psi}", 10, h, f"ref_n{n}",
                        f"plastic_bbar/dp/S75_{f}_{rule}/refconv/n{n}/angle_10/h_{h}", feeds="reply Table R1, Figure R3")
# biphasic elastic and viscoelastic layers
for c in ("c10", "c100"):
    for a in (0, 10):
        for h in HS:
            r = f"biphasic/ref/{c}/a{a}/h{h}"
            add(r, "biphasic", "biphasic", c, a, h, "ref", f"optA/biphasic/{c}/ref/angle_{a}/h_{h}", feeds="Figs 8-10, Table L.4, Fig K.45")
            for m in ("full", "supp"):
                add(f"biphasic/{m}/{c}/a{a}/h{h}", "biphasic", "biphasic", c, a, h, m,
                    f"optA/biphasic/{c}/{m}/angle_{a}/h_{h}", r, "Figs 8-10, Table L.4, Fig K.45")
for h in ("0.1", "0.01"):
    r = f"biphasic_visco/ref/a10/h{h}"
    add(r, "biphasic_visco", "biphasic_visco", "c100", 10, h, "ref", f"optA/biphasic_visco/ref/angle_10/h_{h}", feeds="Fig 14, Fig K.45")
    for m in ("full", "supp"):
        add(f"biphasic_visco/{m}/a10/h{h}", "biphasic_visco", "biphasic_visco", "c100", 10, h, m,
            f"optA/biphasic_visco/{m}/angle_10/h_{h}", r, "Fig 14, Fig K.45")
# regime of the viscoelastic interface model against the step length
for f in FAMS:
    for d in sorted((W / "optA/regime" / f).glob("h_*/dt_*")):
        h, dt = d.parent.name[2:], d.name[3:]
        add(f"regime/{f}/h{h}/dt{dt}", "regime", "visco", f, 10, h, "eie", str(d.relative_to(W)), feeds="Fig F.26")
# in-plane refinement of the interface model at h = 0.01 against the 60 x 60 reference (Appendix J): the interface meshes are
# the subject of this study, so these jobs keep the interface meshes mesh_if/L<level>; L3 was not used
LEVEL = {"L-1": "15", "L0": "30", "L1": "60", "L2": "120"}
for H in ("100", "0.1"):
    r = f"j2/ref/stiff/a10/H{H}/h0.01"
    for lev, n in LEVEL.items():
        eie = f"plastic_bbar/eie_inplane/H_{H}/{lev}" if not (H == "100" and lev == "L1") else "optA/osc100_oldmesh/eie"
        add(f"inplane/eie/H{H}/{lev}", "inplane", "j2", "stiff", 10, "0.01", "eie", eie, r, "Figs J.42-J.44, Table J.2")
        add(f"inplane/psie/H{H}/{lev}", "inplane", "j2", "stiff", 10, "0.01", "psie", f"plastic_bbar/psie_inplane/H_{H}/{lev}", r,
            "Figs J.42-J.44, Table J.2")
        add(f"inplane/psie_sq/H{H}/{lev}", "inplane", "j2", "stiff", 10, "0.01", "psie",
            f"repro_test/laptop_jobs/psie_squeeze_frozen/H_{H}/{lev}", r, "Figs J.42-J.44, Table J.2")
# computational cost (Table 1): the timing runs, stiff family, alpha = 10, h = 0.01
for d in sorted(list((W / "bench").glob("*/*")) + list((W / "bench_sq").glob("*/*"))):
    if (d / "job.inp").exists():
        add(f"bench/{d.parent.parent.name}/{d.parent.name}/{d.name}", "bench", d.parent.name, "stiff", 10, "0.01", d.name,
            str(d.relative_to(W)), feeds="Table 1")
# through-thickness convergence of the elastic and viscoelastic references (reply)
for law in ("elastic", "visco"):
    for d in sorted((W / "refinement_ev60/jobs" / law).glob("*/*")):
        if (d / "job.inp").exists():
            add(f"ev_refconv/{law}/{d.parent.name}/{d.name}", "ev_refconv", law, "stiff", 10, d.name.split("_")[-1],
                d.parent.name, str(d.relative_to(W)), feeds="reply Table R1, Figure R1")

OUT.mkdir(parents=True, exist_ok=True)
with open(OUT / "index.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
print(len(rows), "cases; mesh sets", sorted(used))
(OUT / "mesh_sets.json").write_text(json.dumps(sorted(used)))
