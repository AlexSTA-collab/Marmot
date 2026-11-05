# Figures and tables: script, data, cases

Every figure and table with numbers, the script that drew it, the scores it reads (under `results/`) and the cases it comes
from (studies of `cases/index.tsv`; `ijes list --study <study>`). The scripts are the ones that drew the published
figures; some still name the directories they ran in, while the data they read are all in `results/`. Numbers as in the
revised manuscript; "reply" marks the table of the reply to Reviewer 3.

All convergence panels are drawn by `convergence/plot_visco.py <table> <out>/figures/<key> [<pair table>]
--limits=axes_limits.json`: the ordinates are fixed per figure, as in the original submission, and the EIE and the PSIE panel
of the J2 interphase share one ordinate (the pair table).

| item | script | data (`results/`) | study |
|---|---|---|---|
| Figs 6-7, Table J.2 | `convergence/plot_visco.py`, `convergence/make_rate_tables.py` | `convergence/elastic_table_angle*.md` | elastic |
| Figs 8-9, Table J.3 | `biphasic/plot_panels.py`, `convergence/make_rate_tables.py` | `biphasic/errors60.json` (scorer `biphasic/errors.py`) | biphasic |
| Fig 10 | `biphasic/render_biphasic.py` | EnSight output of `biphasic/*/c100/a10/h0.1` | biphasic |
| Figs 11-12, Table J.4 | `convergence/plot_visco.py`, `convergence/make_rate_tables.py` | `convergence/visco*_table_angle*.md` | visco0, visco10 |
| Fig 13 | `viscoelastic/visco_global60.py` (with `common/extract_reaction.py`, `common/render_deformed.py`) | EnSight output, h = 0.01 | visco0, visco10 |
| Fig 14 | `biphasic_visco/paper_figures.py` (SIDE_LEGEND=1) | `biphasic_visco/*.json` (from `face_histories.py`, `model_histories.py`, `response_visco.py`) | biphasic_visco |
| Figs 15-16, 18-19, Table J.5 | `convergence/plot_visco.py` on the tables of `convergence/scores_to_tables.py <scores> 0.1`, `convergence/make_rate_tables.py` | `convergence/j2_scores_H01.json` (scorer `j2/score.py`) | j2 |
| Fig 20 | `j2/extract_reaction_2611.py`, `j2/plot_reaction_H100.py` (legend lower right), `j2/render_deformed_H100.py` (GLOBAL_DIR) | `global/reaction_H0p1_smesh.json`; EnSight output of `j2/{ref,eie,psie}/stiff/a10/H0.1/h0.01` | j2 |
| Fig 21 | `j2/render_oscillation_H01.py` | EnSight output of `j2/pressure/eie` and `j2/psie/stiff/a10/H0.1/h0.01` | j2 |
| Fig 22, Table J.6 | `inplane/inplane_study.py` (INPLANE_HS=0.1, INPLANE_ELEMS=eie,psie) | `inplane/inplane_scores.json` | inplane |
| Fig 23 | `cyclic/loop2.py`, `cyclic/analyze.py`, `cyclic/plot_cycle_H100.py <out> <data dir>` | `cyclic/H0.1` (`loop.npz`, `err_*.json`) | cyclic |
| Figs 24-25, Table J.7 | `convergence/plot_visco.py` on the tables of `convergence/make_dp_figure_tables.py`, `convergence/make_rate_tables.py` | `convergence/dp/*_dp.json` (scorer `scorers/score_dp.py`) | dp |
| Fig 26 | `dp/plot_dp_global.py`, `dp/extract_dp_global.py`, `dp/render_deformed_dp.py` | EnSight output of `dp/*/initial_dp/a10/h0.01` | dp |
| Fig F.27 | `viscoelastic/regime60.py` (with `common/extract_regime.py`, `common/plot_regime.py`) | EnSight and npydump output | regime |
| Fig G.28 | `viscoelastic/bochner60.py` | `regime/bochner60.json` | visco10 |
| Figs H.29-H.32 | `convergence/plot_visco.py` | `single_reference/{elastic,visco*}_fixed_table_angle*.md` (from `viscoelastic/score_fixed60.py`) | elastic, visco0, visco10 |
| Figs H.33-H.36 | `convergence/plot_visco.py` on `convergence/scores_to_tables.py <scores> 0.1 --fixed` | `single_reference/j2_scores_fixed_H01.json` (from `j2/score_fixed.py`) | j2 |
| Figs H.37-H.38 | `convergence/plot_visco.py` on the tables of `convergence/make_dp_figure_tables.py` | `single_reference/dp/*_dp.json` (from `scorers/score_fixed_dp.py`) | dp |
| Figs I.39-I.41 | `edge/sweep60.py` (J2: EDGE_H=0.1), `edge/plot_edge_panels.py`, `edge/legends_edge.py` | `edge/*.json` | elastic, visco10, biphasic, biphasic_visco, j2, dp |
| Table 1 | typed from the logs | `costs/bench.log` (wall time, peak memory, increments, iterations, cutbacks; one job at a time, four threads) | bench |
| reply Table R1 | `reference/refconv_table.py` | `reference_convergence/*.json` (scorers `scorers/score_ev60.py`, `scorers/refconv_s75_sweep.py`) | ev_refconv, j2_refconv, dp_refconv |

Shared readers and scorers are in `scorers/` and in `ijes/score/` (`errors_end.py`, the measures of the paper).
