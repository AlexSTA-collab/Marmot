#!/usr/bin/env python3
"""results/paper_scores.json: the error measures the paper reports for every interface-model case that `ijes score` can
recompute (errors_end.measures: elastic, viscoelastic, J2, Drucker-Prager), keyed by case, [e_u, e_[u], e_t]. Read from
the data behind the convergence figures and the tables of Appendix L (IJES_submission_Revision/tables/data).
usage: paper_scores.py <tables/data dir> <cases/index.tsv> <out.json>"""
import csv
import json
import sys
from pathlib import Path

DATA, INDEX, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
HS = ("0.01", "0.02", "0.04", "0.06", "0.08", "0.1")


def md(name):
    """{(family, h): [e_u, e_[u], e_t]} from a matched-convergence table (one section per measure)."""
    out, k, fams = {}, -1, []
    for line in (DATA / name).read_text().splitlines():
        if line.startswith("## "):
            k += 1
        elif line.startswith("| h |"):
            fams = [c.split()[0] for c in line.strip("|").split("|")[1:]]
        elif line.startswith("| 0") and k >= 0:
            c = [s.strip() for s in line.strip("|").split("|")]
            for f, v in zip(fams, c[1:]):
                out.setdefault((f, c[0]), [None] * 3)[k] = float(v)
    return out


paper = {}
tab = {("elastic", 0): md("elastic_table_angle0.md"), ("elastic", 10): md("elastic_table_angle10.md"),
       ("visco0", 0): md("visco0_table_angle0.md"), ("visco10", 10): md("visco10_table_angle10.md")}
j2 = json.loads((DATA / "j2_scores_H100_sq.json").read_text())
K = ("e_u_global", "e_jump", "e_traction")
with open(INDEX) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        c, s, f, a, h, m = r["case"], r["study"], r["family"], int(r["angle"]), r["h"], r["model"]
        if not r["reference"]:
            continue
        if s in ("elastic", "visco0", "visco10"):
            paper[c] = tab[(s, a)].get((f, h))
        elif s == "j2" and c.startswith(("j2/eie/", "j2/psie/")):
            H = c.split("/")[4][1:]
            key = f"{m}/{f}/{a}/{H}/{h}"
            if key in j2:
                paper[c] = [j2[key][k] for k in K]
        elif s == "dp":
            fam, rule = c.split("/")[2].split("_")
            d = json.loads((DATA / "dp" / f"{fam}_{rule}.json").read_text())
            paper[c] = d.get(f"eie/{a}/{h}")
paper = {k: v for k, v in paper.items() if v}
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(paper, indent=1, sort_keys=True))
print(len(paper), "cases with a paper value")
