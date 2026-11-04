#!/usr/bin/env python3
"""Add the 26.05 EnSight layout next to the 26.11 per-step files, so the paper's scripts read 26.11 runs unchanged.

    python3 to_2605_layout.py <26.11 run dir> <dst run dir>

The 26.11 output stays untouched; the 26.05-layout copy of a run goes to <dst run dir> (job.inp and run.log are
linked there). port_2611 keeps these copies in runs_2611_as2605/, mirroring runs_2611/.

26.05 wrote one file per variable, "C Binary", then per step  BEGIN TIME STEP | step data | END TIME STEP
(80-byte records), and one geometry.geo built the same way. 26.11 writes one file per step
(<name>.var_0000, ...; geometry.geo_0000) without the markers. The derived files are written to the destination under the 26.05 names; the .case is rewritten
the same way (result time set with the first increment repeated, one .var per variable).

Two differences of the 26.11 output are undone so the unmodified scripts give comparable numbers:
  * 26.05 wrote the first increment twice (time 1e-8 appears twice in its .case); the first step is
    repeated here too, so step indices and step counts agree.
  * the reaction of 26.11 is the negative of 26.05's (checked on every body, see FINDINGS.md); the nodal
    reaction vectors are multiplied by -1.
Per-element variables (stress, generalizedForce) are copied unchanged.
"""
import re
import struct
import sys
from pathlib import Path

import numpy as np

REC = 80


def rec(s):
    return s.ljust(REC, "\0").encode() if len(s) < REC else s.encode()[:REC]


def flip_node_vector(blob, nnodes_by_part):
    """Negate every float of a per-node vector variable step (description, then parts)."""
    out = bytearray(blob[:REC])
    i = REC
    while i < len(blob) and blob[i:i + 4] == b"part":
        out += blob[i:i + REC + 4]
        pn = struct.unpack("<i", blob[i + REC:i + REC + 4])[0]
        i += REC + 4
        out += blob[i:i + REC]            # "coordinates"
        i += REC
        n = nnodes_by_part[pn]
        a = np.frombuffer(blob[i:i + 12 * n], dtype="<f4")
        out += (-a).astype("<f4").tobytes()
        i += 12 * n
    return bytes(out)


def parts_nodes(geo):
    d = geo.read_bytes()
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from ensight_reader_2611 import read_geo
    return {pn: len(p["xyz"]) for pn, p in read_geo(str(geo)).items()}


def rewrite_case(text):
    """26.05 form of a 26.11 .case: the result time set repeats the first increment."""
    blocks = re.split(r"(?=time set:)", text)
    out = []
    for b in blocks:
        m = re.match(r"(time set: \d+ no description\nnumber of steps: )(\d+)(\nfilename start number: 0\nfilename increment: 1\ntime values: )([^\n]*\n)((?:[-+0-9.e]+\n)*)", b)
        if m and int(m.group(2)) > 1:
            vals = (m.group(4) + m.group(5)).split()
            vals = vals[:2] + vals[1:]
            b = m.group(1) + str(len(vals)) + m.group(3) + vals[0] + "\n" + "\n".join(vals[1:]) + "\n" + b[m.end():]
        out.append(b)
    return re.sub(r"\.(var|geo)_\*\*\*\*", r".\1", "".join(out))


def convert(run, dst):
    run, dst = Path(run), Path(dst)
    dst.mkdir(parents=True, exist_ok=True)
    for name in ("job.inp", "run.log"):
        if (run / name).exists() and not (dst / name).exists():
            (dst / name).symlink_to((run / name).resolve())
    for c in run.glob("*.case"):
        (dst / c.name).write_text(rewrite_case(c.read_text()))
    for src in sorted(p for p in run.iterdir() if p.is_dir() and (p / "geometry.geo_0000").is_file()):
        d = dst / src.name
        d.mkdir(exist_ok=True)
        nn = parts_nodes(src / "geometry.geo_0000")
        geo = (src / "geometry.geo_0000").read_bytes()
        (d / "geometry.geo").write_bytes(geo[:REC] + rec("BEGIN TIME STEP") + geo[REC:] + rec("END TIME STEP"))
        names = sorted({re.sub(r"_\d{4}$", "", f.name) for f in src.glob("*.var_[0-9][0-9][0-9][0-9]")})
        for name in names:
            steps = sorted(src.glob(f"{name}_[0-9][0-9][0-9][0-9]"))
            blobs = [s.read_bytes() for s in steps]
            if name.startswith("reaction_"):
                blobs = [flip_node_vector(b, nn) for b in blobs]
            blobs = blobs[:2] + blobs[1:]                  # 26.05 wrote the first increment twice
            body = b"".join(rec("BEGIN TIME STEP") + b + rec("END TIME STEP") for b in blobs)
            (d / name).write_bytes(rec("C Binary") + body)


if __name__ == "__main__":
    convert(sys.argv[1], sys.argv[2])
