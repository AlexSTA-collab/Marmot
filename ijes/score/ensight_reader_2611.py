"""EnSight Gold binary reader: full part table with node ids, and nodal
variables over every time step.

The existing ens.py in this scratchpad reads element variables at a single step;
the force-displacement curves need nodal variables across the whole history and
the NSET parts (to pick out the loaded face), so this adds both.
"""
import numpy as np
import struct

NPE = {"hexa8": 8, "quad4": 4, "penta6": 6, "tetra4": 4, "bar2": 2, "point": 1}


def s80(d, i):
    return d[i:i + 80].split(b"\x00")[0].decode("latin1").strip(), i + 80


def read_geo(path):
    """-> {part_number: {desc, ids, xyz, blocks}}"""
    d = open(path, "rb").read()
    k = d.find(b"BEGIN TIME STEP")
    i = k + 80 if k != -1 else 80   # 26.11 writes one geo per step, without the marker
    _, i = s80(d, i)
    _, i = s80(d, i)
    nidm, i = s80(d, i)
    eidm, i = s80(d, i)
    node_ids_given = nidm.startswith("node id given")
    elem_ids_given = eidm.startswith("element id given")

    parts = {}
    while i < len(d):
        tag, j = s80(d, i)
        if not tag.startswith("part"):
            break
        i = j
        pn = struct.unpack("<i", d[i:i + 4])[0]
        i += 4
        desc, i = s80(d, i)
        _coord_tag, i = s80(d, i)
        nn = struct.unpack("<i", d[i:i + 4])[0]
        i += 4
        if node_ids_given:
            ids = np.frombuffer(d[i:i + 4 * nn], dtype="<i4").astype(int)
            i += 4 * nn
        else:
            ids = np.arange(1, nn + 1)
        xyz = np.frombuffer(d[i:i + 4 * 3 * nn], dtype="<f4").astype(float).reshape(3, nn).T
        i += 4 * 3 * nn

        blocks = []
        while i < len(d):
            t, j = s80(d, i)
            key = t.split()[0] if t else ""
            if key not in NPE:
                break
            i = j
            ne = struct.unpack("<i", d[i:i + 4])[0]
            i += 4
            if elem_ids_given:
                i += 4 * ne
            k = NPE[key]
            conn = np.frombuffer(d[i:i + 4 * ne * k], dtype="<i4").astype(int).reshape(ne, k) - 1
            i += 4 * ne * k
            blocks.append((key, conn))

        parts[pn] = dict(desc=desc, ids=ids, xyz=xyz, blocks=blocks)
    return parts


def read_var_node_all_steps(path, parts, ncomp=3):
    """-> {part_number: array (nsteps, nnodes, ncomp)} for a per-node variable."""
    d = open(path, "rb").read()
    marker = b"BEGIN TIME STEP"
    starts = []
    k = d.find(marker)
    if k == -1:    # 26.11: one file per step, no marker
        starts = [-80]
    while k != -1:
        starts.append(k)
        k = d.find(marker, k + 1)

    out = {}
    for step, s in enumerate(starts):
        i = s + 80
        _desc, i = s80(d, i)
        while i < len(d):
            tag, j = s80(d, i)
            if not tag.startswith("part"):
                break
            i = j
            pn = struct.unpack("<i", d[i:i + 4])[0]
            i += 4
            _coord_tag, i = s80(d, i)
            nn = len(parts[pn]["xyz"])
            a = np.frombuffer(d[i:i + 4 * nn * ncomp], dtype="<f4").astype(float)
            a = a.reshape(ncomp, nn).T
            i += 4 * nn * ncomp
            out.setdefault(pn, []).append(a)
    return {pn: np.array(v) for pn, v in out.items()}
