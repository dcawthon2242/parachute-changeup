# parachute-changeup

Pitching research in R and Python on Statcast (MLB 2015–2026), D1 TrackMan,
and broadcast video. The flagship program asks whether a **changeup that spins
on the same axis as the pitcher's four-seamer** ("parachute" changeup) beats
industry stuff models on whiff because deception is a property of the
fastball–changeup *pair*, not of the pitch alone. Around it sit the programs
that fed it or grew out of it: velocity separation, tunneling, swing timing,
arm-slot reconstruction, and a computer-vision pipeline that measures where a
pitcher stands on the rubber.

> **Parachute status (locked spec, arm-clustered): +1.78 ± 0.97 whiff points
> above model, p = 0.067. Not established. Negative on run value. Do not ship
> as a pitch-design recommendation.** The confirmatory test is pre-registered
> in [`parachute_precommit.md`](data/statcast_model/article_assets/parachute_precommit.md)
> and waits on D1 TrackMan 2022 and 2026.

Start here:

| Document | Role |
|---|---|
| [`docs/parachute-changeup-README.md`](docs/parachute-changeup-README.md) | **The case for pair deception.** Advocacy piece with video, exemplar pitchers, and the honest ledger. |
| [`HANDOFF.md`](HANDOFF.md) | Full technical handoff: spec, models, mechanism battery, pitfalls, script catalog, reproduction order. |
| [`CLAUDE.md`](CLAUDE.md) | Operating rules (MUST / MUST NOT) for anyone or any agent working in this repo. |
| [`parachute_decision_memo.md`](data/statcast_model/article_assets/parachute_decision_memo.md) | Scientific status of record. Wins any disagreement with the article draft. |
| [`parachute_useful.md`](data/statcast_model/article_assets/parachute_useful.md) | Prevalence and run-value decomposition. |
| [`parachute_article.md`](data/statcast_model/article_assets/parachute_article.md) | Narrative draft (Aug 25, 2026). Optimistic. **Not confirmatory.** |
| [`figure_captions.md`](data/statcast_model/article_assets/figure_captions.md) | Captions for every figure in `article_assets/`. |
| [`CCAM-Tunneling-Project/README.md`](CCAM-Tunneling-Project/README.md) | The 2022–2024 tunneling metric and swing-decision model that preceded all of this. |

---

## Quick start

```bash
git clone <this repo> parachute-changeup && cd parachute-changeup

# R (data.table, lightgbm, ggplot2, bit64, mgcv, ...)
Rscript install.R

# Python (video / CV pipeline, NCAA scrapers, clip GIFs)
python3 -m venv .venv-cv
.venv-cv/bin/pip install -r baseball/requirements-cv.txt
```

Everything runs from the repo root and addresses data as `data/...`:

```bash
Rscript baseball/spec_lock.R                       # the confirmatory number
Rscript baseball/bin_roster_detail.R searched      # Cease 2026 / Skubal 2021 smoke test
.venv-cv/bin/python baseball/parachute_clips.py --pitcher 656302 --season 2026 --n 3 --pair
```

Raw Statcast, `.rds` caches, and video frames are **not in git** (several GB).
`data/README.md` lists what is vendored, what is linked, and how to rebuild the
rest. The reproduction order below assumes a fresh clone.

---

## Reproduction (parachute program)

Caches are skipped if the RDS exists; delete it or set `REFIT=1` where the
script honours it.

1. `Rscript baseball/scrape_statcast_multi.R 2020 2021 2022 2023 2024 2025 2026`
   Direct Savant CSV, 4-day windows, resume-safe chunks. Hours. **Never `baseballr`.**
2. Place `data/active_spin/active_spin_YYYY.csv` and `data/savant/arsenal_YYYY.csv` (already vendored).
3. `Rscript baseball/parachute_extended_build.R FF` → `data/statcast_model/parachute_ff.rds`
   Four-seam anchor only. `primary` exists for historical comparison and is not allowed for confirmatory work.
4. `Rscript baseball/active_spin_merge.R` (loop must be `2020:2026`) → `active_spin_long.rds`
5. `Rscript baseball/whiff_tjstuff.R` → `whiff_tjstuff.rds` (residuals r1–r4; use **r4** for anything gated on arm)
6. `Rscript baseball/whiff_parallel.R` → `mlb_whiff_locaware.rds`
7. NCAA: raw TrackMan → `ncaa_03_build.R` → `ncaa_04` … `ncaa_10` until `ncaa_whiff_resid.rds`, `ncaa_spineff_pairs.rds`, `ncaa_armangle.rds`, `ncaa_usage.rds` exist. `library(bit64)` **before** `as.character` on `PitcherId`.
8. **`Rscript baseball/spec_lock.R`** → `locked_spec.rds`. Match the table in `HANDOFF.md` §3 before doing anything else.
9. `mech_01` … `mech_11` regenerate the decision memo; `mech_12`–`15` reproduce the slot autopsy.
10. Do not run a new gate search. `optimal_gates.R` is historical; if you demo it, quote the permutation p (≈ .22).

Expected scale of `parachute_ff.rds`: ~300k–465k regular-season changeups with a valid FF anchor.

---

## Repo layout

```
parachute-changeup/
  README.md                 this file
  CLAUDE.md                 operating rules (HANDOFF §0)
  HANDOFF.md                full technical handoff
  install.R                 R package install
  baseball/                 all analysis scripts, flat (293 files; catalog below)
    requirements-cv.txt     Python deps for the video / CV pipeline
    tjstuff_run_values.csv  run-value lookup used by the tjStuff+ port
  CCAM-Tunneling-Project/   2022–2024 tunneling metric + swing-decision model (own README)
  docs/
    parachute-changeup-README.md   the case for pair deception (video, exemplars, ledger)
  data/
    README.md               tracked vs untracked, rebuild instructions
    active_spin/            vendored
    savant/                 vendored arsenal usage
    mlb_stats/, ncaa_rosters/, pitcher_heights.csv   vendored
    statcast_model/
      article_assets/       memos, captions, figures, ext_*.csv, clips/   (tracked)
      *.rds, tunnel_location/                                            (ignored)
    statcast_YYYY/, rubber/, obm/, swing_timing/, open_command/, ncaa/   (ignored)
```

`baseball/` is deliberately flat. `HANDOFF.md` §9 sketches numbered subfolders,
but every script addresses siblings as `baseball/<name>` and data as `data/...`
relative to the repo root, so moving files would break them. The grouping
lives in the catalog instead.

---

## Script catalog

One line per family. Names are prefixes in `baseball/` unless noted.

### Parachute changeup (flagship)

| Scripts | Role |
|---|---|
| `spec_lock.R` | **The frozen confirmatory spec.** Only allowed bin definition for confirmatory claims. |
| `parachute_extended_build.R`, `active_spin_merge.R`, `whiff_tjstuff.R`, `whiff_parallel.R`, `parachute_ff_anchor.R` | Build the FF-anchored changeup table and out-of-fold whiff residuals (r1–r4, location-aware). |
| `bin_roster_detail.R`, `bin_gate_attrition.R`, `bin_pitcher_counts.R`, `core_plus_eff_floor.R` | Rosters and attrition for Core / Wide / mean / searched / eff90 bins. |
| `mech_01_usage` … `mech_15_slot_ladder.R` | Mechanism and falsification battery behind the decision memo (dose, sequencing, tunnel, bands, spin-rate control, era holdout, within-pitcher, splitter, multiplicity, RV decomposition, D1 transfer, miss distance, slot). |
| `optimal_gates.R`, `arm_gate_vs_descriptive.R`, `velo_corr_searched_bin.R`, `arm_slot_sweep.R`, `axis_sweep_arm40.R` | The Aug 11 gate search and its autopsy. **Historical; do not publish the winner.** |
| `parachute_*.R` (roster, slot, spin_efficiency, highspin, within_pitcher, run_value, gb_check, filtered, above_avg, extended_test) | Pre-lock discovery variants. Label historical. |
| `deception_changeups_2026.R`, `find_changeups_2026.R`, `pitch_pair_*.R`, `matched_spin_overperformance.R`, `whiff_matched_spin.R` | Jul 31–Aug 6 discovery: turnover vs parachute, matched-spin residuals. |
| `axis_collapse.R`, `axis_after_tunnel.R`, `axis_sd_changeup.R`, `axis_signflip_validate.R`, `axis_vs_break.R`, `axis_gap_highspin.R`, `zone_axis_decomposition.R` | Spin-axis diagnostics, including the Figure 13 retraction of the spin-similarity finding. |
| `spin_axis_atlas.R`, `spin_axis_audit.R`, `spin_axis_provenance.R`, `spin_examples_and_offspeed.R`, `same_axis_sliders.R`, `movement_only_spin_test.R`, `breaking_spin_recheck.R`, `active_spin_offspeed_recheck.R` | Arsenal-wide axis geometry and measurement provenance. |
| `kick_change_*.R`, `spin_lookalike_kickchange_2026.R` | Kick-change foil: low-spin offspeed underperforms stuff. |
| `cue_comparison.R`, `cue_decomposition.R`, `cue_outlier_pitchers.R`, `fig11_*.R`, `fig3_miss_drivers.R`, `article_visuals.R`, `mlb_chase_parachute.R`, `mlb_topthird.R`, `d1_resid_diagnose.R` | Figures and article support. |
| `parachute_clips.py` | Savant broadcast clips → GIFs (and FF/CH side-by-sides) for the exemplar pitches. |

### Velocity separation (absorbed parachute as a special case)

| Scripts | Role |
|---|---|
| `velo_sep_01_blind.R` … `velo_sep_46_axis_matched_drop.R` | Do stuff models underrate big-separation changeups? Blind test, anchors, within-pitcher, D1 replication, seam-shifted wake, portability, archetypes, feature absorption, RV conversion, usage, teachability. `45`/`46` are the post-lock seam-shift hypothesis in the article draft. |
| `velo_sep_by_type.R`, `velo_bin_interaction.R`, `velo_feature_attribution.R`, `velo_gate_grid.R`, `velo_resid_robustness.R`, `velo_gap_failures.R`, `whiff_velo.R`, `whiff_resid_vs_velo.R`, `chase_velo_slope.R` | Separation vs residual across types, robustness, attribution. |
| `velo_gap_tunnel_tradeoff.R`, `pair_velo_gap_path_match.R`, `poor_tunnel_high_gap.R`, `high_gap_*.R`, `high_velo_gap_timing_misses.R`, `gap_*.R` | The gap-versus-tunnel tradeoff and where big gaps leak value. |
| `hard_slider_vs_gap*.R`, `hard_slider_same_axaz.R`, `slider_gap_by_velo.R`, `why_hard_sliders_stuff.R` | 98/90 vs 98/86 slider question. |

### Tunneling

| Scripts | Role |
|---|---|
| `CCAM-Tunneling-Project/` | Original lattice tunnel metric (200 timestamps release → 150 ms before plate), swing-decision LightGBM, chase-above-expected. Own README. |
| `swing_decision_rv_2026.R` | Out-of-sample rebuild of the swing-decision run value (train 2023–25, apply 2026). |
| `angular_tunnel_*.R` | Tunneling as the hitter's eye receives it (angular, not field-coordinate); grading, validation, leaderboards, two-strike nulls. |
| `tunnel_location_build.R`, `tunnel_location_*.R`, `tunnel_pairing_maps.R`, `tunnel_case_study.R`, `tunnel_command.R`, `tunnel_delta_test.R`, `breaking_location_tunnel.R`, `tight_tunnel_overperform.R`, `path_ratio_correlates.R` | Sequential FB → breaking-ball pairs; location dominates, tunnel adds 2–3%. |
| `location_approach_test.R`, `whiff_chase_surfaces.R`, `combo_performance_channels.R` | Location / approach-angle surfaces. |

### Miss-distance grade model

`miss_grade_features.R`, `miss_grade_train.R`, `miss_grade_augmented.R`, `miss_grade_grid.R`, `miss_grade_putaway.R`, `miss_grade_gainers*.R`, `miss_grade_apply_insights.R`, `miss_decomp_2026.R`, `miss_distance_vs_xwobacon.R`, `risers_non2k_offspeed.R`, `driver_gain_all_types.R`, `remaining_residual.R`, `count_scope_test.R`.
Target is bat-tracking `miss_distance` (2023 All-Star break onward), 2026 holdout.

### Stuff models

`tjstuff_plus_v3.R` (Nestico's tjStuff+ v3.0 port, `tjstuff_run_values.csv`), `stuff_platoon_features.R`, `stuff_platoon_train.R` (LHH/RHH split models; what release and rubber position add on top of shape), `kick_stuff_resid_2026.R`.

### Swing timing, tscore / bscore, contact quality

`swing_timing_fetch.R`, `swing_timing_probe.R`, `swing_timing_*.R`, `timing_*.R`, `tscore_plus.R`, `tscore_vs_zone*.R`, `tscore_whiffs.R`, `tscore_inzone_contact.R`, `whiff_plus.R`, `barrel_score.R`, `barrel_score_plus.R`, `barrel_vs_hardhit_dimensions.R`, `contact_targets*.R`, `tdev_vs_xwobacon_by_pitch.R`, `flail_target_2026.R`, `power_projection.R`.
Timing disruption from batter-relative intercept depth; converts to runs per 100 swings.

### Command, ERA, odds and ends

`command_join.R` (OpenCommand inferred targets → command metrics), `physics_fip.R`, `physics_fip_era.R`, `kershaw_era_compare.R`, `changeup_platoon_2026.R`, `changeup_run_value_2026.R`, `sweeper_platoon_2026.R`, `sweeper_platoon_quality_2026.R`, `mlb_chase_decomp.R`, `scrape_scope.R`.

### Arm slot

`arm_angle_reconstruct.R`, `arm_angle_reconstruct2.R` (arm angle from release point + listed height, to extend seasons), `arm_angle_tunnel.R`, `arm_distribution.R`, `proxy_slot_feasibility.R`, `height_value_test.R`, `build_d1_heights.R`, `batter_silhouette.R`.

### Rubber position (computer vision)

`rubber_01_playids.py` → `rubber_01b_target_games.R` → `rubber_02_select_clips.R` → `rubber_02b_fetch_targets.R` → `rubber_03_fetch_frames.py` → `rubber_04*.py` (measure, label pack, label server, frame score, contact sheet, park triage, CV measure, keypoint, train pack, pose) → `rubber_05_calibrate.R`, `rubber_05b_label_qc.R` → `rubber_06_natural_experiment.R`, `rubber_07_counterfactual.R`, `rubber_07b_channel_decomp.R`, `rubber_08_handedness_pitchtype.R` → `rubber_09_annotate_examples.py`, `rubber_10_paper_figures.py`, `rubber_11_leaderboard.py`.
Where a pitcher sets up on the rubber, measured from broadcast frames, and what it is worth. Deps in `requirements-cv.txt` (opencv, torch, imageio-ffmpeg).

### Hitter eye point (OpenBiomechanics)

`obm_probe.py`, `obm_eye_point.py`, `obm_hitter_eye_points.py`. Average L/R hitter eye position from C3D trials in Statcast feet, for the angular tunnel work. Needs `ezc3d`.

### NCAA D1 port

`ncaa_00_probe.R` … `ncaa_13_sequencing.R` (TrackMan axis comparability, inferred efficiency, parachute port, arm proxy, chase, sequencing), `scrape_ncaa_rosters.py`, `sidearm_discover.py`, `sidearm_probe.py`, `sidearm_scrape.py`, `sidearm_seeds.py` (roster scrapes via school athletics sites).
TrackMan `SpinAxis` is a movement axis, not Hawk-Eye's imaged axis; the matched cut is ≈ 7.1°, not 10°.

### Ingest

`scrape_statcast_multi.R`, `scrape_statcast_2026.R` (direct Savant CSV), `active_spin_merge.R`, `swing_timing_fetch.R`.

---

## Rules that cannot be relaxed

Short form of `CLAUDE.md`:

- Four-seam anchor only. No sinker/cutter fallback.
- Arm-aware residual (`r4`) for anything gated on arm.
- `library(bit64)` before `as.character` on NCAA `PitcherId`.
- `spec_lock.R` is the only bin definition for confirmatory claims.
- Arm-clustered SEs when pooling repeated pitcher-seasons.
- `rv = -delta_run_exp`; positive favours the pitcher.
- Round active-spin gaps to 4 decimals before threshold tests.
- No searching for a spec after a null. No publishing the gate-search winner. No cherry-picked rosters as evidence. No `baseballr` scrapes. No imputing arsenal usage.

---

## Provenance

Code and memos were moved here on 2026-09-09 from a mixed work/baseball
workspace where `baseball/` had never been committed, so this repo's first
commit is the first version control these scripts have had. The source
conversations are cited in `HANDOFF.md`. `CCAM-Tunneling-Project/` is copied
in flat; its own 16-commit history remains at
https://github.com/dcawthon2242/CCAM-Tunneling-Project.
