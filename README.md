# Reproducibility package: a generalised interface model for inelastic interphases (IJES, revision 2026-11)

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
  plug-ins of Marmot `plugins/ijes-revision` (C3D8B, the Drucker--Prager material, the PSIE variants).

The meshes (~300 MB) and the archived outputs of all analyses are on Zenodo and are fetched by `ijes fetch`.

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
ijes run elastic/stiff/a10/h0.04/ref elastic/stiff/a10/h0.04/eie -j 2
ijes score elastic/stiff/a10/h0.04/eie  # e_u, e_[u], e_t, and the value the paper reports
ijes sweep elastic -j 4                 # a whole thickness sweep: references, models, scores
```

`ijes-run <job dir> [threads]` runs any job directory the way the paper's analyses were run.

## What each study costs

See `cases/index.tsv` for the cases. The resolved references with twenty elements through the layer on the 60 x 60
mesh are the expensive part; the interface models run in seconds to minutes. Memory and run times measured on the
workstation of the paper (Intel Xeon, 4 threads per job) are listed in `results/costs.tsv`.

## Layout

| path | what |
|---|---|
| `cases/` | job files, `index.tsv` (case, study, law, family, angle, h, model, job, reference, feeds, paper value) |
| `ijes/` | command line (`ijes`), the scorer and its EnSight reader, `efe_run.py` (roller tolerance of the paper, npydump output) |
| `results/` | scores behind every figure and table |
| `figures/` | scripts that draw the figures from `results/` |
| `env/`, `Dockerfile`, `docker/` | the solver |
