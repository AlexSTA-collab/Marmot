"""Command line of the reproducibility package.

  ijes list [--study S]                 the cases, one per job file of the paper, with the figure or table they feed
  ijes run CASE [CASE ...] [-j N] [-t T] run the job of each case (and nothing else) under work/<case>
  ijes score CASE                       error measures of an interface-model case against its reference case,
                                        next to the value the paper reports
  ijes sweep STUDY [-j N] [-t T]        run every case of a study, references first, then score them
  ijes fetch meshes|outputs [--study S] download the meshes, or the archived outputs of the paper into the work
                                        directory (Zenodo)

A case is a row of cases/index.tsv. Its job file includes the meshes as @IJES_ROOT@/meshes/..., filled in when the
job is copied to $IJES_WORK/<case> (default ./work), so it runs where the package is unpacked once the meshes are
fetched.
"""
import argparse
import csv
import json
import os
import subprocess
import sys
import tarfile
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(os.environ.get("IJES_ROOT", Path(__file__).resolve().parents[1]))
WORK = Path(os.environ.get("IJES_WORK", "work")).resolve()
ZENODO = json.loads((ROOT / "zenodo.json").read_text()) if (ROOT / "zenodo.json").exists() else {}
PAPER = json.loads((ROOT / "results" / "paper_scores.json").read_text())
SCORED = ("elastic", "visco0", "visco10", "j2", "dp")   # the studies measured with errors_end.measures


def cases():
    with open(ROOT / "cases" / "index.tsv") as f:
        return {r["case"]: r for r in csv.DictReader(f, delimiter="\t")}


def workdir(c):
    return WORK / c


def run_one(c, threads):
    row = cases()[c]
    src, dst = ROOT / "cases" / row["job"], workdir(c)
    dst.mkdir(parents=True, exist_ok=True)
    # the includes of the packaged job files start with @IJES_ROOT@/meshes/
    (dst / "job.inp").write_text(src.read_text().replace("@IJES_ROOT@", str(ROOT)))
    log = dst / "run.log"
    if log.exists() and "Job computation time" in log.read_text(errors="ignore"):
        return c, "already done"
    # the env column holds the settings a case was run with, e.g. MARMOT_VONMISES_RELTOL=1
    env = dict(os.environ, **dict(kv.split("=", 1) for kv in row.get("env", "").split(";") if kv))
    rc = subprocess.run(["ijes-run", str(dst), str(threads)], env=env).returncode
    log = (dst / "run.log").read_text(errors="ignore")
    ok = rc == 0 and "Job computation time" in log and "Simulation failed" not in log
    return c, "done" if ok else f"FAILED (see {dst / 'run.log'})"


def score_one(c):
    sys.path.insert(0, str(ROOT / "ijes" / "score"))
    import errors_end as E  # noqa: E402

    row = cases()[c]
    ref = row["reference"]
    if not ref:
        raise SystemExit(f"{c} is a reference; score an interface-model case")
    if row["study"] not in SCORED:
        raise SystemExit(f"{c}: the {row['study']} study has its own measures; see figures/README.md")
    E.ANGLE, E.H0 = int(row["angle"]), float(row["h"])
    o = E.measures(str(workdir(c)), row["h"], str(workdir(ref)), row["h"], 60)
    paper = PAPER.get(c)
    out = {"case": c, "reference": ref, **{k: o[k] for k in ("e_u_global", "e_jump", "e_traction")},
           "paper [e_u, e_[u], e_t]": paper}
    print(json.dumps(out))
    return out


def fetch(what, study):
    if not ZENODO:
        raise SystemExit("zenodo.json is missing: the record is not published yet")
    names = [n for n in ZENODO["files"][what] if study is None or n.startswith(study)]
    # meshes.tar.zst holds meshes/<set>/..., unpacked into the package; <study>.tar.zst holds <case>/... (every file of the
    # run as it ran), unpacked into the work directory, where `ijes score` reads it
    for n in names:
        url, dst = f"{ZENODO['base']}/{n}?download=1", (ROOT if what == "meshes" else WORK) / n
        dst.parent.mkdir(parents=True, exist_ok=True)
        print("fetching", n, flush=True)
        urllib.request.urlretrieve(url, dst)
        with tarfile.open(dst) as t:   # zstd: tarfile of Python 3.14
            t.extractall(dst.parent, filter="data")
        dst.unlink()


def main(argv=None):
    p = argparse.ArgumentParser(prog="ijes", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    s = p.add_subparsers(dest="cmd", required=True)
    a = s.add_parser("list"); a.add_argument("--study")
    a = s.add_parser("run"); a.add_argument("case", nargs="+"); a.add_argument("-j", type=int, default=1); a.add_argument("-t", type=int, default=4)
    a = s.add_parser("score"); a.add_argument("case", nargs="+")
    a = s.add_parser("sweep"); a.add_argument("study"); a.add_argument("-j", type=int, default=1); a.add_argument("-t", type=int, default=4)
    a = s.add_parser("fetch"); a.add_argument("what", choices=("meshes", "outputs")); a.add_argument("--study")
    a = p.parse_args(argv)
    if a.cmd == "list":
        for c, r in cases().items():
            if a.study is None or r["study"] == a.study:
                print(f"{c:48s} {r['model']:6s} -> {r['feeds']}")
    elif a.cmd == "run":
        with ThreadPoolExecutor(a.j) as ex:
            for c, status in ex.map(lambda c: run_one(c, a.t), a.case):
                print(c, status, flush=True)
    elif a.cmd == "score":
        for c in a.case:
            score_one(c)
    elif a.cmd == "sweep":
        rows = [c for c, r in cases().items() if r["study"] == a.study]
        refs = [c for c in rows if not cases()[c]["reference"]]
        mods = [c for c in rows if cases()[c]["reference"]]
        for batch in (refs, mods):
            with ThreadPoolExecutor(a.j) as ex:
                for c, status in ex.map(lambda c: run_one(c, a.t), batch):
                    print(c, status, flush=True)
        res = [score_one(c) for c in mods]
        (WORK / f"{a.study}_scores.json").write_text(json.dumps(res, indent=1))
    elif a.cmd == "fetch":
        fetch(a.what, a.study)


if __name__ == "__main__":
    main()
