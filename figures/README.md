# Figures and tables: script, data, cases

Every figure and table with numbers, the script that drew it, the scores it reads (under `results/`) and the cases it comes
from (studies of `cases/index.tsv`; `ijes list --study <study>`). The scripts are the ones that drew the published
figures, copied unchanged; some still name the directories they ran in, while the data they read are all in `results/`.
Numbers as in the revised manuscript; "reply" marks the figures and the table of the reply to Reviewer 3.

| item | script | data (`results/`) | study |
|---|---|---|---|
| Figs 6-7, Table L.3 | `convergence/plot_visco.py`, `convergence/make_rate_tables.py` | `convergence/elastic_table_angle*.md` | elastic |
| Figs 8-9, Table L.4 | `biphasic/plot_panels.py`, `convergence/make_rate_tables.py` | `biphasic/errors60.json` (scorer `biphasic/errors.py`) | biphasic |
| Fig 10 | `biphasic/render_biphasic.py` | EnSight output of `biphasic/*/c100/a10/h0.1` | biphasic |
| Figs 11-12, Table L.5 | `convergence/plot_visco.py`, `convergence/make_rate_tables.py` | `convergence/visco*_table_angle*.md` | visco0, visco10 |
| Fig 13 | `viscoelastic/visco_global60.py` (with `common/extract_reaction.py`, `common/render_deformed.py`) | EnSight output, h = 0.01 | visco0, visco10 |
| Fig 14 | `biphasic_visco/paper_figures.py` (SIDE_LEGEND=1) | `biphasic_visco/*.json` (from `face_histories.py`, `model_histories.py`, `response_visco.py`) | biphasic_visco |
| Figs 15-16, 18-19, Table L.6 | `convergence/plot_visco.py` on the tables of `convergence/scores_to_tables.py` | `convergence/j2_scores_H100_sq.json` (scorer `j2/score.py`) | j2 |
| Fig 20, Fig H.30 | `j2/plot_reaction_H100.py`, `j2/render_deformed_H100.py` | `global/reaction_H100_smesh.json`, `global/reaction_H0p1_smesh.json` (from `j2/extract_reaction_2611.py`) | j2 |
| Fig 21 | `j2/render_oscillation_H100.py` | EnSight output of `j2/oscillation/eie`, `j2/psie/stiff/a10/H100/h0.01` | j2 |
| Fig 22, Fig H.31 | `cyclic/plot_cycle_H100.py` | `cyclic/H100`, `cyclic/H0.1` (`loop.npz` from `cyclic/loop2.py`, `err_*.json` from `cyclic/analyze.py`) | cyclic |
| Figs 23-24, Table L.7 | `convergence/plot_dp_panels.py`, `convergence/make_rate_tables.py` | `convergence/dp/*.json` (scorer `scorers/score_dp.py`) | dp |
| Fig 25 | `dp/plot_dp_global.py`, `dp/extract_dp_global.py`, `dp/render_deformed_dp.py` | EnSight output of `dp/*/initial_dp/a10/h0.01` | dp |
| Fig F.26 | `viscoelastic/regime60.py` (with `common/extract_regime.py`, `common/plot_regime.py`) | EnSight and npydump output | regime |
| Fig G.27 | `viscoelastic/bochner60.py` | `regime/bochner60.json` | visco10 |
| Fig H.28 | `reference/plot_profiles_H01.py` | `profiles/*.json` (from `scorers/row_strain.py`) | j2_refconv |
| Fig H.29 | `reference/plot_psie_H.py` | `convergence/j2_scores_H100_sq.json` | j2 |
| Figs I.32-I.35 | `convergence/plot_visco.py` | `single_reference/{elastic,visco*}_fixed_table_angle*.md` (from `viscoelastic/score_fixed60.py`) | elastic, visco0, visco10 |
| Figs I.36-I.39 | `convergence/plot_visco.py` on `convergence/scores_to_tables.py --fixed` | `single_reference/j2_scores_fixed_H100.json` (from `j2/score_fixed.py`) | j2 |
| Figs I.40-I.41 | `convergence/plot_dp_panels.py` | `single_reference/dp/*.json` (from `scorers/score_fixed_dp.py`) | dp |
| Figs J.42-J.44, Table J.2 | `inplane/inplane_study.py`, `inplane/plot_psie_inplane.py` | `inplane/inplane_scores.json` | inplane |
| Figs K.45-K.46 | `edge/plot_edge_panels.py` | `edge/*.json` (from `edge/sweep60.py`) | elastic, visco10, biphasic, biphasic_visco, j2, dp |
| Table 1 | typed from the logs | `costs/bench*.log` (wall time, peak memory; increments and iterations from each run.log) | bench |
| Table M.8 | typed from the logs | increments and cutbacks of the run.log of each case | j2, cyclic, dp |
| reply Table R1, Figs R1-R3 | `reference/refconv_table.py`, `reference/plot_ev_converged.py`, `reference/plot_j2_ref_converged.py`, `reference/plot_dp_refconv_s75.py` | `reference_convergence/*.json` (scorers `scorers/score_ev60.py`, `scorers/refconv_s75_sweep.py`) | ev_refconv, j2_refconv, dp_refconv |

Shared readers and scorers are in `scorers/` and in `ijes/score/` (`errors_end.py`, the measures of the paper).
