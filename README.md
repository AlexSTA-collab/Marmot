# Reproducibility package: An inelastic imperfect interface model derived from the underlying interphase

International Journal of Engineering Science, manuscript JENS-D-26-01640 (revision, November 2026).

Branch `ijes-reproducibility` of https://github.com/AlexSTA-collab/Marmot. It shares no history with the Marmot sources: the
solver is built from the tags below, and this branch holds the job files, the scorer, the scores and the figure scripts.

Every analysis of the revised manuscript can be rerun from this package, and every convergence figure and table rebuilt
from its scores. It holds

- the **job file of every analysis** of the paper, interface models and resolved references alike (`cases/`), indexed
  by the figure or table it feeds (`cases/index.tsv`);
- the **scorer** that computes the error measures of the paper, e_u, e_[u] and e_t (`ijes/score/errors_end.py`,
  margin 0.15 mm, N = 60 cells per side);
- the **scores** behind the figures and tables (`results/`) and the scripts that draw them (`figures/`);
- a **Dockerfile** that builds the solver exactly as it was run: the conda environment of the workstation
  (`env/linux-64.explicit.txt`, 309 packages), Marmot (tag `ijes-revision-2026-11c`) and EdelweissFE (tag `ijes-revision-2026-11`), and the
  plug-ins of Marmot `plugins/ijes-revision` (C3D8B, the Drucker--Prager material, the PSIE).

The meshes and the archived outputs of all analyses are on Zenodo (one zstd-compressed tar per study, see below) and are
fetched by `ijes fetch`.

## Software

| component | source |
|---|---|
| Marmot | https://github.com/AlexSTA-collab/Marmot, tag `ijes-revision-2026-11c` |
| plug-ins | the same tag, `plugins/ijes-revision/` |
| EdelweissFE | https://github.com/AlexSTA-collab/EdelweissFE, tag `ijes-revision-2026-11` |
| environment | `env/linux-64.explicit.txt` (Python 3.14 free-threading, conda-forge) |

## Quick start

```bash
# the solver (Docker or Apptainer)
docker pull ghcr.io/alexsta-collab/marmot:ijes-reproducibility-2026-11
docker run --rm -it -v $PWD/work:/work ghcr.io/alexsta-collab/marmot:ijes-reproducibility-2026-11 bash
#   apptainer run docker://ghcr.io/alexsta-collab/marmot:ijes-reproducibility-2026-11

# inside the container
ijes fetch meshes                       # once
ijes list --study elastic               # the cases of a study and the figure each one feeds
ijes run elastic/ref/stiff/a10/h0.04 elastic/eie/stiff/a10/h0.04 -j 2
ijes score elastic/eie/stiff/a10/h0.04  # e_u, e_[u], e_t, and the value the paper reports
ijes fetch outputs --study j2            # or score the archived outputs of the paper without rerunning
ijes sweep elastic -j 4                 # a whole thickness sweep: references, models, scores
```

`ijes-run <job dir> [threads]` runs any job directory the way the paper's analyses were run.

## What each study costs

See `cases/index.tsv` for the cases. The resolved references with twenty elements through the layer on the 60 x 60
mesh are the expensive part; the interface models run in seconds to minutes. Memory and run times measured on the
workstation of the paper (Intel Xeon Platinum 8558P, 4 threads per job, one job at a time) are in `results/costs/bench.log`
(Table 1, study `bench`).

## Layout

| path | what |
|---|---|
| `cases/` | job files, `index.tsv` (case, study, law, family, angle, h, model, job, reference, feeds, paper value) |
| `ijes/` | command line (`ijes`), the scorer and its EnSight reader, `efe_run.py` (roller tolerance of the paper, npydump output) |
| `results/` | scores behind every figure and table |
| `figures/` | scripts that draw the figures from `results/` |
| `env/`, `Dockerfile`, `docker/` | the solver |

## The analyses

| study | cases | what |
|---|---|---|
| elastic, visco0, visco10 | 72, 36, 36 | EIE against the resolved layer (20 C3D8) of the same thickness, three families, six thicknesses |
| biphasic, biphasic_visco | 72, 6 | biphasic layers, the mismatch kept (`full`) and removed (`supp`) |
| regime | 126 | viscoelastic interface model against the step length |
| j2 | 109 | J2 layer at H = 0.1 MPa: references (20 C3D8B), EIE and PSIE; `j2/pressure/eie` is the EIE run of Figure 21 with the pressure of its face stress |
| cyclic | 4 | reversing path at H = 0.1 MPa, PSIE and reference, h = 0.1 and 0.01 |
| inplane | 8 | in-plane refinement, EIE and PSIE on the 15, 30, 60, 120 meshes of the interface (flat interface meshes) |
| dp | 72 | Drucker--Prager layer, phi = 30, psi = 10, H = 0.075 G: EIE and references (20 C3D8) |
| ev_refconv, j2_refconv, dp_refconv | 36, 18, 24 | references with 5 and 10 elements through the layer (reply, Table R1) |
| bench | 11 | the timing runs of Table 1: elastic, viscoelastic and Drucker--Prager of 27 and 31 October 2026, J2 at H = 0.1 of 4-5 November 2026 |

The interface models use the shared-node interface mesh `S_mesh` (the body nodes of the reference mesh), except the
in-plane study, which uses its own flat interface meshes, and the October timing runs of the elastic and viscoelastic
interface models in Table 1, which ran on the interface mesh `if_L1` of that date (included for that table only). The `env` column of `cases/index.tsv` holds the settings a case
ran with: the EIE of the J2 layer at H = 0.1 MPa uses `MARMOT_VONMISES_RELTOL=1` (relative tolerance of the return
mapping, tag `ijes-revision-2026-11c`).

## Zenodo deposit

`meshes.tar.zst` holds `meshes/<set>/angle_<a>/height_<h>/...`, the meshes the job files include. `<study>.tar.zst` holds
`<case>/...`, every file of each run of the study as it ran (EnSight output, `run.log`, the job file with the paths of the
workstation). `SHA256SUMS` lists the checksums. Unpack by hand with `tar --zstd -xf <file>` or let `ijes fetch` do it.
