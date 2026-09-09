# Parachute changeup research — Claude Code handoff

Self-contained brief for standing up this work in a **new repo**. Copy this
file in as `HANDOFF.md` (and paste §0 into `CLAUDE.md`). Do not treat the
current workspace as a clean source tree: it is a mixed MapQuest + baseball
workspace. The baseball work lives under `baseball/` and
`data/statcast_model/` in `/Users/damon.cawthon/System1MapquestWorkspace`.

**Source conversation:** [Parachute changeup deception](4625b8bb-28de-4d17-acec-0337795b4859)
(Fri Jul 31 → Tue Aug 25, 2026). Older truncated snapshot:
[Changeup finder scrape](ddcd4b7b-d037-4351-b368-5e09268c70c4).

---

## 0. Operating rules for the new-repo agent

Read this section first. It overrides the article draft.

### Status of record (do not argue with this)

| Document | Role |
|---|---|
| `parachute_decision_memo.md` | **Scientific status.** Locked cell is unexplained and not established. Do not build product on the bin. |
| `parachute_precommit.md` | **Frozen confirmatory rule** for D1 2022 and 2026. Cannot be rewritten after those seasons are seen. |
| `parachute_useful.md` | Prevalence (publishable as description) + RV decomposition (does not support "useful"). |
| `parachute_article.md` | Narrative draft, Aug 25. Optimistic. **Not confirmatory.** |

The locked arm-clustered effect is **+1.78 ± 0.97 whiff points, p = 0.067**.
Mechanisms that would attribute this to "hitter mistakes the pitch for a
fastball" are null. Held-out MLB 2024–2026 is slightly negative. The next
honest step is the pre-committed D1 2022/2026 run, not another gate search.

### MUST

- Anchor spin-axis gap on a **four-seam fastball only**. No sinker/cutter fallback.
- Score an arm-slot gate against an **arm-aware** residual (`r4` or equivalent). Scoring arm as a "deception" feature against an arm-blind model rediscovers a missing-feature artifact.
- Load `library(bit64)` **before** `as.character` on NCAA `PitcherId`.
- Treat `spec_lock.R` as the only allowed bin definition for confirmatory claims.
- Report **arm-clustered** SEs as primary when pooling repeated pitcher-seasons.
- Sign convention: `rv = -delta_run_exp`. Positive run value favours the **pitcher**.
- Round active-spin gaps to 4 decimals before threshold tests.

### MUST NOT

- Search for an approach that proves a predetermined conclusion after a null. That is how Figure 13's spin-similarity finding was born.
- Publish tuned thresholds from `optimal_gates.R` (or any new 4D grid) as "the" parachute definition. The Aug 11 search failed its permutation null (p ≈ .22) and split-half (86% shrinkage).
- Rebuild Figure 9's cherry-picked-per-bin pitcher list as evidence.
- Use `baseballr::statcast_search()` against current Savant. Direct CSV only.
- Impute missing arsenal usage. Missing row = drop; blank cell = 0.
- Equate early branding (Ragans, Hellickson) with locked-bin membership. Ragans often fails axis ≤ 10°.
- Confuse **whiff residual** with **usefulness**. Locked bin is negative on run value.
- Treat `parachute_article.md` as the result of the confirmatory program.
- Mix **symmetric efficiency-gap** gates (`|as_ch − as_fb| ≤ 0.10`) with **floor** gates (`both ≥ .85`). They select different rosters.
- Use TrackMan `SpinAxis` as if it were Hawk-Eye imaged axis. At D1 it is approximately the movement axis.

### What this repo is for

A standalone R research repo that:

1. Scrapes / stores Statcast 2020–2026 + Savant active-spin + arsenal usage.
2. Builds a four-seam-anchored changeup table and out-of-fold whiff residuals.
3. Reproduces the **frozen** matched-axis bin and the mechanism / falsification battery.
4. Optionally ports the same spec to D1 TrackMan with the documented substitutions.
5. Keeps discovery scripts, but labels them historical.

It is **not** a pitch-design product, a live leaderboard, or a claim that parachute changeups are a proven pitch.

---

## 1. What a parachute changeup is

Two changeup archetypes, named by the user on Aug 3, 2026:

| Archetype | Spin vs four-seam | Typical names | Claim |
|---|---|---|---|
| **Turnover** | Large axis gap (pronated, different clock face) | Matt Boyd, Jason Alexander, Alex Vesia (~32°) | Looks like a changeup out of the hand. Stuff models see the shape. |
| **Parachute** | Small axis gap (spins like the heater), slower, often high slot, then drops | Jeremy Hellickson (historical), Cole Ragans (early branding), **2026 Dylan Cease**, **2021 Tarik Skubal** | Looks like a fastball for ~15 feet. Stuff models grade it as a slow fastball with no shape of its own. Hypothesis: it beats those models on whiff because deception is a *pair* property, not an object property. |

Early 2026 movement filter that started the work (not the locked definition):

- velocity gap vs FB > 10 mph
- |HB gap| < 3 in
- |IVB gap| < 8 in

That filter found Cantillo, Ragans, Ribalta, Sullivan and already lost on run
value (−0.83 RV/100 vs league ~0). Shape ≠ value. The research then moved to
spin-axis match as the discriminator.

**Physical claim, in one sentence:** a changeup whose imaged spin axis matches
the four-seamer is indistinguishable from the heater until velocity (and,
if present, seam-shifted wake) shows up after the swing decision.

**What the data actually supports:** a rare shape (~1.8% of four-seam-primary
changeup seasons) with a fragile positive whiff residual that does not
survive arm-clustering at p < .05, does not convert to run value, and has
no demonstrated perception mechanism.

---

## 2. Intellectual history (so you do not re-walk the traps)

### 2.1 Jul 31 – Aug 3: discovery

- Scrape 2026 Statcast. Movement-shape CH finder.
- Aug 3: turnover vs parachute named. Vesia classified as turnover (axis 31.8°).
- `deception_changeups_2026.R` / `pitch_pair_deception_2026.R`: stuff-only whiff model, residual vs spin-axis relationship.
- Headline that survived: **matched spin alone is worthless**. Matched axis + a big velo/spin kill is the engine. Same axis, small kill → *underperforms*.

### 2.2 Aug 3–6: miss distance, tunneling, specification shopping

- Miss-distance grade model (`miss_grade_*.R`), TJStuff+-style features, 2023H2–2025 train / 2026 holdout.
- Adding `axis_diff` as a model feature did **not** recover parachute overperformance in a gated all-types model — the model learned "bigger spin *difference* → more miss" (splitters, kick-changes) and buried Cease/Vesia.
- Correct proof method: residual over a shape model that **never sees** axis gap, then ask whether matched-axis pitches have positive residuals.
- **Figure 9:** user asked for 10 pitchers per spin-sim bin chosen to prove the point. Do not treat as evidence.
- **Figure 11b spin-similarity:** later **retracted by Figure 13** (`axis_collapse.R`). Putting location in the model collapses breaking-ball axis correlations. Changeup axis never significant (best p ≈ .093).
- Kick changes (Muñoz-type: low spin, near-zero IVB) underperform stuff. Foil, not cousin.

### 2.3 Aug 11: four-seam anchor, Core bin, unsupervised gates

This is the day the user originally asked to recover.

- **Primary-FB fallback retracted.** A sinker's seam-shifted wake decouples spin axis from movement. SI-anchored seasons were 9.4% of the pool but 15% of Core; mean axis gap 18.4° vs 22.3° FF. The old bin was partly a sinkerballer detector.
- Rebuild 2020–2026 with `parachute_extended_build.R FF` → `parachute_ff.rds`.
- Core (`parachute_ff_anchor.R`): FF only, axis ≤ 10°, |active-spin gap| ≤ 0.10, arm ≥ 44°, ≥ 60 CH swings. Mean residual **+0.48**. Drop Cease 2026 and Skubal 2021 → **−0.95**. Shape bin ≠ quality bin.
- User asked: *"What are the optimal gates for whiff% for the stats I gave you above"* (arm, axis, FF active spin, CH active spin). See **§5**.

### 2.4 Aug 12–18: D1 port, lock, mechanisms

- TrackMan axis is not Hawk-Eye axis. Literal 10° does not transfer.
- `spec_lock.R` freezes the cell. `mech_01`–`mech_11` test mechanisms and falsifiers.
- Decision memo: unexplained, not established. Pre-commit written **before** D1 2022/2026.

### 2.5 Aug 19–25: article vs memo split

- Further gate shopping (Rodón/Roark exclusion, "optimized" 17-season set) is reverse-engineered. Not locked.
- Article draft (Aug 25) claims a teachable, run-positive pitch once axis match is paired with seam-shifted wake. Memo still says wait. **Ship the memo as ground truth; keep the article as a draft that overreaches.**

---

## 3. Frozen confirmatory specification

**Source of truth:** `baseball/spec_lock.R` writes `data/statcast_model/locked_spec.rds`.
Everything confirmatory reads the bin from that file and nothing else.

```
changeup active spin   >= 0.85
four-seam active spin  >= 0.85
spin-axis gap vs FF    <= 10 degrees          # unsigned circular
arm slot               >= league's own 67th percentile
                          (recomputed AFTER the population filter)
qualifying             >= 40 changeup swings in the pitcher-season
population filter      ff_use >= si_use
                          missing arsenal ROW dropped, not imputed
                          blank cells treated as 0 usage
                          filter applies to pool AND comparison group
```

**Why a percentile arm gate, not 44°:** D1 arm is a regression estimate,
SD ≈ 68% of MLB's measured spread. Fixed 44° selects the top 17% of D1 vs
the top 27% of MLB. A percentile asks the same question of both leagues.

**Why FF-primary:** axis is measured against the four-seamer. For a
sinker-primary pitcher the hitter's expectation is set by the sinker, so
the measurement is the wrong reference. Honesty note: this exclusion was
informed by the first MLB look (6 of 27 bin seasons were sinker-primary;
dropping them moved 0–10 from +2.05 to +3.00). Locked before new D1
seasons. **Not** a clean out-of-sample rule for the MLB rerun.

**Residual:** season-level velocity-adjusted, both leagues:
`y = resid(lm(y_raw ~ vs))` within league. D1 needed it (aggregation
artifact: season-level residual correlated with mean velo sep at r = +.062
even though pitch-level was ~0). MLB does not need it; applied anyway so
both sides are treated identically.

`y_raw` = 100 × mean pitch-level out-of-fold whiff residual from a
location-aware model (`mlb_whiff_locaware.rds` / `ncaa_whiff_resid.rds`).

**Estimator that decides:** inverse-variance pool of the two league means,
**arm-clustered** standard error. College careers overlap; treating
repeated seasons as independent overstates precision. Season-level p is
secondary.

### Current locked numbers

| | seasons | arms | effect (whiff pp above model) | p |
|---|---|---|---|---|
| MLB 2020–2026 | 20 | 14 | +3.07 ± 1.38 | 0.038 |
| D1 2023–2025 | 14 | 13 | +1.15 ± 1.95 | 0.564 |
| Pooled, season-level | 34 | 27 | +2.43 ± 1.12 | 0.031 |
| **Pooled, one row per arm (decides)** | **27** | **27** | **+1.78 ± 0.97** | **0.067** |

Heterogeneity Q = 0.65 (p = 0.421).

MLB locked arms (`parachute_useful.md`): Kikuchi (3), Cease, Ray, Cabrera,
Skubal (2 each); Sulser, Underwood Jr., May, Smith, Jax, Detmers,
Naughton, Jones, Shuster (1). Prevalence: **1.8% of 1,132** four-seam-primary
MLB changeup seasons; **1.7% of 820** D1.

### Pre-commit for D1 2022 and 2026 (already locked)

When those seasons are added, **same specification, no further tuning**:

- **confirm** if pooled p < 0.05 **and** point estimate stays above **+1.2**
- **kill** if point estimate falls below **+0.8**
- **unresolved** otherwise; wait for MLB 2027

Use the arm-clustered estimator as primary. Extra seasons from arms already
in 2023–2025 are not independent.

### Historical bins (do not mix)

| Name | Rule | Role |
|---|---|---|
| **Core** | axis ≤ 10, \|as_gap\| ≤ 0.10, arm ≥ 44, ≥ 60 swings, FF only | Aug 11 descriptive. Symmetric efficiency *gap*, not a floor. |
| **Wide** | axis ≤ 15, \|as_gap\| ≤ 0.15, arm ≥ 42 | Looser shell. |
| **mean** | FF act > league mean, CH act > league mean, axis ≤ 10, arm ≥ 44 | A-priori round cuts. |
| **searched** | arm ≥ 35, axis ≤ 8, FF ≥ .85, CH ≥ .89 | `optimal_gates.R` winner. Inspection only. |
| **eff90 / eff90axis** | both ≥ .90, arm ≥ 44, ± axis ≤ 10 | MLB roster of the D1 .90 screen. |
| **locked** | §3 | Confirmatory. Efficiency *floors* + percentile arm + FF-primary + 40 swings. |

`bin_roster_detail.R` implements Core/Wide/mean/searched/eff90 via argv:
`Rscript baseball/bin_roster_detail.R searched`

---

## 4. The unsupervised gate search (Aug 11, 2:21 PM)

This is the piece the user asked to recover from transcripts.

**Ask:** optimal gates on arm slot, spin-axis gap, four-seam active spin,
changeup active spin, maximizing **whiff residual**.

**Script:** `baseball/optimal_gates.R`
**Figure:** `data/statcast_model/article_assets/fig37_optimal_gates.png`
**Roster:** `Rscript baseball/bin_roster_detail.R searched`
**Velo check:** `baseball/velo_corr_searched_bin.R`

### Setup

- Population: pitcher-seasons from `whiff_tjstuff.rds` merged to
  `active_spin_long.rds`, **≥ 60 changeup swings**.
- Outcome `w = 100 * mean(r4)` — residual from **tjStuff+ v3.0 + arm + location**.
  Arm is in the model because the search includes an arm gate.
- `axis_diff` is **excluded** from the residual model (otherwise the bin
  definition is in the predictor).
- Grid: arm `seq(25, 55, 2.5)`, axis `seq(6, 30, 2)`, FF/CH active spin
  `seq(.85, .96, .02)` → **6,084** combinations. Minimum bin size **20**.
- Score: mean(bin) − mean(everyone else), unweighted across pitcher-seasons.
- `set.seed(7)`.

### Winner (nominal)

```
arm slot           >= 35°
axis gap           <= 8°
FF active spin     >= 0.85
CH active spin     >= 0.89
```

20 pitcher-seasons, **+2.89** whiff points above model, 95% CI [+0.85, +4.92],
nominal p = 0.0077.

Top of roster (arm-aware residual, raw rates in `bin_roster_detail.R searched`):

| Pitcher | Season | CH | Swings | Arm | Axis | Act FF/CH | Velo | Sep | Spin | Whiff | Expected | Above | RV/100 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Dylan Cease | 2026 | 226 | 84 | 59.7 | 4.8 | .91 / .92 | 82.8 | 15.0 | 1768 | 53.6 | 42.5 | **+11.0** | +0.82 |
| Tarik Skubal | 2021 | 305 | 112 | 60.0 | 3.7 | .94 / .98 | 82.1 | 12.0 | 1441 | 49.1 | 38.8 | **+10.3** | +0.13 |

Then Osvaldo Bido 2025 (+6.5), Robbie Ray 2025 (+6.1), Angel Zerpa 2023 (+5.9), …

### Why that number is not real

**Permutation null (300 shuffles of the outcome, identical search):**
median winner **+2.32**, 90th **+3.36**, max **+4.99**.
Search-corrected p ≈ **0.22**. The bar to clear is +2.32, not zero.

**Split-half (300 splits, minn=10 on the train half):**
in-sample +3.60 vs held-out **+0.50**. **86% shrinkage.** Held-out positive
in 67% of splits; 5th–95th **−1.94 to +2.78**.

**Threshold stability:** arm wanders 25–52.5°, axis 8–22°. Active-spin
thresholds collapse to the grid floor — the optimizer is discarding both
efficiency gates. It only ever reaches for arm and axis, and even those
are not sharply located.

**Reference gates that day (no searching):**

| Gate | Seasons | Whiff above model | p |
|---|---|---|---|
| All four above league average | 155 | −0.12 | .82 |
| Spin floors + axis ≤ 10 | 35 | +1.08 | .24 |
| Spin floors + axis ≤ 10 + arm ≥ 44 | 8 | +4.46 | .046 |
| Searched optimum | 20 | +2.89 | .008 |

The 8-season cell is tiny and not search-corrected. Do not promote it.

**Velo-sep vs residual inside the searched bin:** appears strongly
positive until Cease and Skubal are dropped; then it nearly vanishes.
Spearman << Pearson. Leverage, not a law.

**What to tell a reader:** the *direction* of the gradient (tighter axis,
higher slot → more residual) has been consistent. The *location of a cliff*
is not identifiable at this sample size. Publish round a-priori cuts, not
the search winner.

---

## 5. Models

### 5.1 Feature construction (`parachute_extended_build.R`)

```r
cmean <- function(a) { r <- a*pi/180; (atan2(mean(sin(r)), mean(cos(r)))*180/pi) %% 360 }
circd <- function(a,b) { z <- abs(a-b) %% 360; pmin(z, 360-z) }
```

- Fastball reference: pitcher-season, **≥ 50** pitches of type FF (or FF/SI/FC if `ANCHOR=primary` — do not use primary for confirmatory work).
- `fb_axis` = circular mean of `spin_axis`.
- `axis_diff = circd(spin_axis, fb_axis)` — unsigned, 0–180.
- `speed_diff = release_speed - fb_velo` (negative = slower CH). **`velo_sep = -speed_diff`**.
- `az_diff = az - fb_az`. **`kill = -az_diff`** (IVB-kill proxy).
- `rv = -delta_run_exp`.
- Swings: `swinging_strike`, `swinging_strike_blocked`, `foul`, `foul_tip`, `hit_into_play`, plus bunts.
- Whiff: swinging strikes + foul tip + missed bunt. **Not** fouls.
- VAA/HAA from Statcast kinematics at plate front (`yf = 17/12`).
- `sax/cax = sin/cos(spin_axis)` stored; `atan2(sax,cax)` recovers degrees.
- LHP mirroring on x-features when fitting tjStuff-style models: `x0, ax, ax_diff, spin_axis → 360-axis`.

Filter: `game_type == "R"`, finite kinematics, unique on `game_pk, at_bat_number, pitch_number`.

### 5.2 Whiff residual used by the gate search (`whiff_tjstuff.R`)

Fit on changeup **swings** from `parachute_ff.rds`. LightGBM binary, 4-fold OOF,
12% validation split for early stopping.

| Residual | Features | Use |
|---|---|---|
| r1 | Nestico tjStuff+ v3.0 eleven: speed, spin, extension, ax, az, x0, z0, spin_axis, speed_diff, ax_diff, az_diff (LHP-mirrored) | "above stuff" in the public sense |
| r2 | r1 + `arm_angle` | |
| r3 | r1 + location (`plate_x_arm, plate_z, z_rel_bot, z_rel_top, VAA, HAA_in, stand_R, balls, strikes`) | |
| **r4** | r1 + arm + location | **Default for any analysis that gates on arm** |

**Excluded from all four:** `axis_diff`.

LightGBM params (shared across this project):

```
objective = binary, metric = binary_logloss
learning_rate = 0.06, num_leaves = 31, min_data_in_leaf = 300
feature_fraction = 0.8, bagging_fraction = 0.8, bagging_freq = 1
nrounds = 1500, early_stopping = 50
```

RobustScaler is omitted: trees split on order statistics.

### 5.3 Locked confirmatory residual (`whiff_parallel.R` → `mlb_whiff_locaware.rds`)

Location-aware, includes spin-axis sine/cosine among features (so the
*bin's* axis *gap* is still a between-pitcher contrast, not a pitch-level
feature the model can use to eat the effect). 4-fold OOF. `r_all = whiff - p`.

D1 analogue: `ncaa_whiff_resid.rds`.

### 5.4 Miss-distance grade (historical, adjacent)

`miss_grade_features.R` / `miss_grade_train.R` / `miss_grade_augmented.R`.

- Target: `miss_distance` on competitive swings; contact imputed to 0; untracked-whiff dropped.
- Train 2023-07-14+ through 2025; **2026 holdout** (not k-fold for architecture selection).
- Bat tracking starts 2023 All-Star break. Locked bin is 2020–2026, so miss contrasts are a **subset** (`mech_12`).
- Base features are TJStuff+ shape, **no location**. Location dominates RMSE if added (~0.89" drop on breaking balls) and erases tunneling's marginal value.
- Offspeed augmentation (`axis_diff`, spin-efficiency gap) helps Approach A (all-types, pitch type unknown) mostly as a **pitch-type proxy**, not as deception.

Do not use miss-distance as the confirmatory target unless you are explicitly
reproducing the 2023H2–2026 subset analysis.

### 5.5 Shape+location OOF for RV / GB (`parachute_ff_anchor.R`)

Same 4-fold LightGBM. Targets: `rv` (regression), `whiff` (binary, swings),
`ground_ball` (binary, BIP). **Excludes `axis_diff` and usually `arm_angle`.**
Produces `parachute_ff_resid.rds`.

---

## 6. Mechanism and falsification battery

Scripts `baseball/mech_01_usage.R` … `mech_15_slot_ladder.R`.
Results: `parachute_decision_memo.md`.

| ID | Question | Result |
|---|---|---|
| M1 | Does more four-seam usage raise the bin residual? (dose) | Pooled interaction −0.075 per usage point, p = .47. D1 significantly *negative* (−0.48, p = .007). Leagues disagree. |
| M2 | Do bin changeups overperform more when the previous pitch was a four-seamer? | +2.34 ± 2.62, p = .37. Underpowered (80% power needs ~7.3 points). Changeups after a four-seamer underperform for *everyone* (−1.72, p < .001). |
| M3 | Does tunneling concentrate in the bin? | Bin × tunnel +0.01 ± 1.14, p = .99. General tunnel slope is real: **+0.47 whiff points per SD**, p = .009, both leagues, location already in the model. It is not a parachute finding. |
| Band | Is 30–45° the real cell? | No. 0–10 is the only band positive in both leagues. |
| F1 | Spin-*rate* matching as quality confound | **−3.85 ± 0.73**, p < .001. Opposite of a quality story. Keep this negative control. |
| F2 | Era holdout | 2020–2023 in-sample +4.83 (p = .003); **2024–2026 −0.31** (p = .92). Looks like a threshold fitted to noise. |
| F3 | Within-pitcher | 12 MLB arms on both sides of the bin: +2.50 ± 2.26, p = .29. Cease +18.5 to Jones −9.3. |
| Alt | Handedness, count, release scatter | Handedness flat. Scatter 0.206 vs 0.212; controlling it does not move the bin. Scatter itself predicts (−0.42 residual per SD, p = .003) — another finding that does **not** need the bin. |
| Splitter | Same gates, FS | 4 MLB seasons, 0 D1. Underpowered, not a null. MLB FS only 2023–2026; D1 FS pool tiny. |
| RV | `mech_10_rv_decomp.R` | Locked bin **−0.41 RV/100** vs rest of FF-primary pool. Whiff share 17.0 vs 16.3 (+0.04 contribution). BIP is the hole (−4.93 vs −3.30). xwOBA-implied still −0.23. Expected contact is the wrong sign for a luck story. |
| Multiplicity | `mech_09` | Locked cell survives its own 6-cell arm×axis family (p = .011) and does **not** survive being the best of five axis bands (p = .125). D1's best cell in the same family is a *different* one. |
| Slot autopsy | `mech_14`, `mech_15` | Axis × slot "survived" then died after residual model / arm source / FF filter / velo adjustment / clustering changes. Do not revive it without a pre-registered spec. |

**Findings that do not depend on the bin being real** (act on these if anything):

1. Better tunneling → more whiffs (+0.47 / SD, p = .009), both leagues, pitch-level.
2. Tighter release scatter → more whiffs (−0.42 / SD, p = .003).
3. Changeups after a four-seamer underperform (−1.72, p < .001) for everyone.
4. Tight spin-*rate* matching underperforms (−3.85, p < .001).

### Article draft's later mechanism (not locked)

`parachute_article.md` and `velo_sep_45` / `velo_sep_46` argue:

- Tight axis match without seam-shifted wake: misses bats, loses runs on fly balls (backspin + slow + high slot).
- Axis match **+** seam deviation (≥ ~7.8°) **+** velo sep (≥ ~8.7 mph): the version that converts.
- Drop-beyond-spin (the "parachute down" story) is **null**. Seam shift in the working cell is largely horizontal; direction-agnostic magnitude, not depth.
- Conjunction table (axis-matched, 187 seasons): only the high-sep × high-seam-dev corner is alive (+4.45, p = .0001). Three other corners sit at zero.

Treat this as a **hypothesis generated after the lock**, not as the confirmatory result. If you implement it, put it in `exploratory/` and do not silently replace `spec_lock.R`.

---

## 7. D1 / NCAA transfer

Raw files (outside the current workspace, typically `~/Downloads/`):

- `pbp23tm.csv`, `D1TM24.csv`, `D1TM25.csv`
- Filter `Level == "D1"`
- Build: `ncaa_03_build.R` → `ncaa_d1_pitches.rds`

**The measurement problem:** TrackMan `SpinAxis` ≈ direction of movement
(~5° of break for changeups). Hawk-Eye `spin_axis` is an imaged axis. A
literal 10° cut does not mean the same thing.

**When they converge:** inferred efficiency ≥ .95 on both pitches. Movement
vs imaged gap correlation 0.46 → 0.79; median |diff| 20.5° → 8.9°.

**Correct transfer (often skipped in early D1 reruns):** percentile-match
D1's axis cut to MLB's ≤10° share of the .95 pool → about **7.1°**, not 10°.
`mech_11_d1_p95.R`: D1 .95 + matched axis = **+5.58 ± 3.77 on 5 arms, p = .14**.
Underpowered. MLB at the same .95 floor and literal 10° is +0.55 on 12 seasons.

**D1 cannot test seam-shifted wake.** Movement axis ≡ SpinAxis ⇒ SSW ≈ 0
by construction.

**Arm:** reconstructed (`ncaa_10_armangle.R`, `arm_hat`). Compressed
distribution → percentile gate, not 44°.

**IDs:** `PitcherId` is `integer64`. Without `bit64`, `as.character` yields
strings like `"3.37562447470698e-318"` that still join *within one session*
and silently fail against a script that loaded `bit64`. Mandatory.

**D1 2022 and 2026 are the confirmatory data.** They are not in the current
run. The decision rule is already written. Do not peek and retune.

---

## 8. Data sources

### 8.1 Statcast pitch-level (MLB)

**Do not use `baseballr`.** `baseballr` 1.6.0 hardcodes ~92 column names;
Savant now returns ~119 and name assignment fails.

Scripts: `baseball/scrape_statcast_multi.R`, `baseball/scrape_statcast_2026.R`.

```
Rscript baseball/scrape_statcast_multi.R 2020 2021 2022 2023 2024 2025 2026
```

- Direct `https://baseballsavant.mlb.com/statcast_search/csv?...&type=details`
- 4-day windows, 5s sleep, 4 retries, User-Agent required
- Resume-safe chunk files under `data/statcast_<YEAR>/chunks/`
- Output: `data/statcast_<YEAR>/statcast_<YEAR>_all.csv`
- `hfGT=R|PO|S` then filter `game_type == "R"` in analysis
- 2020 special-cased (60-game season starting ~Jul 23) so empty April–July windows are not treated as throttles
- `miss_distance` and bat tracking arrive automatically in `type=details` from 2023-07-14

Columns the parachute build requires: `game_pk, at_bat_number, pitch_number,
game_date, game_type, pitcher, player_name, pitch_type, p_throws, stand,
balls, strikes, release_speed, release_spin_rate, release_extension,
release_pos_x/y/z, spin_axis, arm_angle, plate_x/z, sz_bot/top, vx0/vy0/vz0,
ax/ay/az, description, events, bb_type, launch_speed/angle,
estimated_woba_using_speedangle, delta_run_exp`. Plus `miss_distance` if you
rebuild the grade model.

Hawk-Eye `arm_angle` exists from **2020**, not 2023. Active-spin leaderboards
also from 2020. That is why the confirmatory window is 2020–2026.

### 8.2 Measured active spin

`data/active_spin/active_spin_2020.csv` … `2026.csv`.

Savant **measured** active-spin % leaderboards, wide (one row per pitcher,
columns per pitch type: `active_spin_fourseam`, `_changeup`, …). Values are
percent; `active_spin_merge.R` divides by 100.

Long form: `data/statcast_model/active_spin_long.rds`
(`pitcher, season, pitch_type, active_spin`).

**Bug to fix on rebuild:** committed `active_spin_merge.R` loops `2023:2026`
but the RDS on disk is 2020–2026. Change the loop to `2020:2026`. Copy CU
onto KC/CS as the script already does.

There is no in-repo downloader. Hand-pull from Baseball Savant or write one.
Do not go back to inferred Magnus `spin_eff` as the bin gate; it was
unreliable. Inferred efficiency is only for the D1 port.

### 8.3 Arsenal usage (FF-primary filter)

`data/savant/arsenal_2020.csv` … `2026.csv`. Columns include `n_ff`, `n_si`.
`spec_lock.R` binds these. Blank → 0. Missing row → drop.

### 8.4 NCAA TrackMan

See §7. Not in git. Document the expected paths in a `.env.example` or
`data/ncaa/README.md`.

---

## 9. Suggested repo layout

```
parachute-changeup/
  CLAUDE.md                         # paste §0
  HANDOFF.md                        # this file
  README.md                         # short: what / status / how to run
  renv.lock                         # pin R pkgs
  baseball/
    00_scrape/
      scrape_statcast_multi.R
      scrape_statcast_2026.R
    01_build/
      parachute_extended_build.R    # ANCHOR=FF
      active_spin_merge.R           # fix year loop
      whiff_tjstuff.R
      whiff_parallel.R
      parachute_ff_anchor.R         # optional Core/Wide residuals
    02_lock/
      spec_lock.R                   # THE spec
      bin_roster_detail.R
    03_mechanisms/
      mech_01_usage.R … mech_15_slot_ladder.R
    04_ncaa/                        # ncaa_00 … ncaa_13
    05_historical/                  # do not use for claims
      deception_changeups_2026.R
      pitch_pair_*.R
      optimal_gates.R
      arm_gate_vs_descriptive.R
      velo_corr_searched_bin.R
      parachute_*.R (pre-lock variants)
      miss_grade_*.R
      axis_collapse.R
    06_exploratory/                 # post-lock hypotheses
      velo_sep_12_parachute.R
      velo_sep_45_sep_vs_seam.R
      velo_sep_46_axis_matched_drop.R
  data/
    README.md                       # what is gitignored vs vendored
    statcast_YYYY/                  # gitignore CSVs (huge)
    active_spin/*.csv               # small; vendor
    savant/arsenal_*.csv            # small; vendor
    ncaa/                           # gitignore raw TrackMan
    statcast_model/                 # rds caches; gitignore or Git LFS
      article_assets/               # memos + captions CAN be committed
  docs/
    parachute_decision_memo.md
    parachute_precommit.md
    parachute_useful.md
    parachute_article.md            # label DRAFT / NOT CONFIRMATORY
    figure_captions.md
```

Copy scripts from the current workspace paths listed in §12. Do not copy
MapQuest, clustering, or rubber/timing work unless asked.

**R packages:** `data.table`, `lightgbm`, `ggplot2`, `curl`, `bit64`.
Optional: `renv`. Python is not required except `scrape_ncaa_rosters.py` if
you rebuild D1 names.

---

## 10. Reproduction order

Caches are skipped if the RDS exists. Delete or set `REFIT=1` /
`Sys.getenv("REFIT")` where scripts honor it.

1. Scrape Statcast 2020–2026. Hours, resume-safe.
2. Place `data/active_spin/active_spin_YYYY.csv` and `data/savant/arsenal_YYYY.csv`.
3. `Rscript baseball/parachute_extended_build.R FF`
   → `data/statcast_model/parachute_ff.rds` + `parachute_ff_fbtype.csv`
4. Fix and run `active_spin_merge.R` (years 2020:2026)
   → `active_spin_long.rds`
5. `Rscript baseball/whiff_tjstuff.R` → `whiff_tjstuff.rds` (r1–r4)
6. `Rscript baseball/whiff_parallel.R` → `mlb_whiff_locaware.rds`
7. Optional historical: `optimal_gates.R` (permutation 300 × 6084 is slow;
   tens of minutes to hours). Then `bin_roster_detail.R searched`.
8. NCAA: raw CSVs → `ncaa_03_build.R` → `ncaa_04` … `ncaa_10` as needed
   until `ncaa_whiff_resid.rds`, `ncaa_spineff_pairs.rds`, `ncaa_armangle.rds`,
   `ncaa_usage.rds` exist.
9. **`Rscript baseball/spec_lock.R`** → `locked_spec.rds`. This is the
   confirmatory number. Match the memo table before proceeding.
10. `mech_01` through `mech_11` (then 12–15 if reproducing the autopsy).
11. Do not run a new gate search. If you must demo `optimal_gates.R`, keep
    its output in `05_historical/` and quote the permutation p.

Expected `parachute_ff.rds` scale: on the order of **~300k–465k** regular-season
changeups 2020–2026 with a valid FF anchor.

---

## 11. Key numbers cheat sheet

**Searched-optimum (inspection only):** Cease 2026 +11.0, Skubal 2021 +10.3
whiff points over r4; 20 seasons +2.89 nominal / permutation p ≈ .22.

**Core (FF, axis≤10, |as_gap|≤.10, arm≥44, 60 swings):** ~21 seasons, mean
+0.48; without Cease+Skubal −0.95; RV/100 negative.

**Locked:** §3 table. Deciding number +1.78, p = .067.

**RV locked bin:** −0.41 / 100 vs rest; expected −0.23. Raw whiff share only
+0.7 points. The model residual and the raw rate can both be true: the model
already prices location and shape, so +3 MLB residual means "more miss than
this stuff should," not "an elite whiff pitch in raw rate."

**Holdout era:** +4.83 (2020–23) → −0.31 (2024–26).

**Article inversion table (different archetypes, 75-swing floor, not locked):**
tight axis match predicted 30.3% / actual 41.3%. Power changeup (big axis)
underperforms. Do not cite as the confirmatory result.

---

## 12. Script catalog (current workspace paths)

Prefix: `/Users/damon.cawthon/System1MapquestWorkspace/`

### Scrape / ingest

| Path | Purpose |
|---|---|
| `baseball/scrape_statcast_multi.R` | Multi-year Savant CSV scrape |
| `baseball/scrape_statcast_2026.R` | 2026-only ancestor |
| `baseball/active_spin_merge.R` | Wide leaderboards → long RDS (fix year loop) |

### Build / residual

| Path | Purpose |
|---|---|
| `baseball/parachute_extended_build.R` | 2020–2026 CH table; `FF` or `primary` |
| `baseball/whiff_tjstuff.R` | r1–r4; gate-search residual |
| `baseball/whiff_parallel.R` | MLB location-aware residual for lock |
| `baseball/parachute_ff_anchor.R` | Core/Wide; FF-only residuals |

### Lock / roster / search

| Path | Purpose |
|---|---|
| `baseball/spec_lock.R` | Frozen cell |
| `baseball/bin_roster_detail.R` | mean / searched / core / wide / eff90 |
| `baseball/optimal_gates.R` | 4D search + permutation + split-half |
| `baseball/arm_gate_vs_descriptive.R` | Arm as gate vs description |
| `baseball/velo_corr_searched_bin.R` | Velo sep vs r4 in searched bin |

### Mechanisms

`baseball/mech_01_usage.R` M1
`baseball/mech_02_sequencing.R` M2
`baseball/mech_03_tunnel.R` M3
`baseball/mech_04_band.R` 30–45 band / sinker-primary
`baseball/mech_05_falsify.R` F1–F3
`baseball/mech_06_alt.R` platoon / count / scatter
`baseball/mech_07_figure.R` decision figure
`baseball/mech_08_splitter.R` FS analog
`baseball/mech_09_multiplicity.R` best-of-k
`baseball/mech_10_rv_decomp.R` RV channels
`baseball/mech_11_d1_p95.R` D1 axis transfer
`baseball/mech_12_miss_gate.R` miss distance on locked recipe
`baseball/mech_13_miss_sig.R` miss significance
`baseball/mech_14_slot_interaction.R` slot re-test
`baseball/mech_15_slot_ladder.R` why slot died

### NCAA

`baseball/ncaa_00_probe.R` … `ncaa_13_sequencing.R`
(axis comparability, inferred efficiency, parachute port, arm proxy, chase, sequencing)

### Historical discovery (label, do not claim)

| Path | Purpose |
|---|---|
| `baseball/deception_changeups_2026.R` | First parachute vs turnover stats |
| `baseball/pitch_pair_deception_2026.R` | Matched/mirror spin vs stuff residual |
| `baseball/pitch_pair_deception_v2_2026.R` | FF-anchor, kill, VAA |
| `baseball/pitch_pair_all_types_2026.R` | Rule across type pairs |
| `baseball/parachute_run_value.R` | Early RV / GB hypotheses |
| `baseball/parachute_within_pitcher.R` | Within-arm Δaxis |
| `baseball/parachute_roster_*.R` | Pre-FF-primary rosters |
| `baseball/parachute_highspin*.R` | Floor-not-gap; missing-feature demo |
| `baseball/parachute_extended_test.R` | 2020–22 OOS of an *earlier* frozen target (GB +6.9) — that target **did not replicate** |
| `baseball/axis_collapse.R` | Figure 13 retraction |
| `baseball/miss_grade_*.R` | Miss-distance grade model |
| `baseball/find_changeups_2026.R` | Original movement filter |

### Post-lock exploratory (velo-sep lineage)

`baseball/velo_sep_12_parachute.R` — parachute traits inside high-sep strata (Bonferroni).
`baseball/velo_sep_45_sep_vs_seam.R` — sep vs seam-shift horse race.
`baseball/velo_sep_46_axis_matched_drop.R` — conjunction / drop-beyond-spin null.

Earlier `velo_sep_01`–`44` are a broader separation research program that
**absorbed** parachute as a special case. Copy only if the new repo is
meant to continue that program, not just the parachute lock.

### Status documents (commit these)

| Path | Role |
|---|---|
| `data/statcast_model/article_assets/parachute_decision_memo.md` | Status of record |
| `data/statcast_model/article_assets/parachute_precommit.md` | Decision rule |
| `data/statcast_model/article_assets/parachute_useful.md` | Prevalence + RV |
| `data/statcast_model/article_assets/parachute_article.md` | Draft, overreaches |
| `data/statcast_model/article_assets/figure_captions.md` | Captions for the *deception/tunneling* article, mixed with parachute |

### Canvases (optional; Cursor canvas format)

Under `~/.cursor/projects/Users-damon-cawthon-System1MapquestWorkspace/canvases/`:

- `parachute-changeups-2020-2026.canvas.tsx` — FF-anchored roster, Elite/Core/Wide
- `parachute-teachability.canvas.tsx` — within-arm teachability (`velo_sep_41`–`44`)
- `matched-axis-changeup-final.canvas.tsx` — discrepancy-gated battery

---

## 13. Pitfalls (full list)

1. **integer64 PitcherId** — §7. Failures are silent.
2. **Sinker fallback** — matched axis vs SI ≠ matched look. FF only.
3. **Arm-blind scoring of an arm gate** — missing-feature artifact. Use r4.
4. **TrackMan axis ≠ Hawk-Eye axis.**
5. **Inferred vs measured active spin** — measured for MLB gates.
6. **Float boundary on as_gap** — round to 4 decimals.
7. **Gate-search overfitting** — permutation + split-half required; still do not publish the winner as the spec.
8. **Miss distance is 2023+ only.**
9. **Slot interaction is not robust** across residual definitions (`mech_15`).
10. **Fail-closed usage** — never impute sinker-primary back in.
11. **D1 arm proxy compressed** — percentile, not 44°.
12. **Sequencing file limits** — `parachute_rv.rds` is CH+FB types only, so prev-pitch match fails for breaking-ball predecessors (`mech_02`).
13. **No baseballr scrape.**
14. **`active_spin_merge.R` year loop** is wrong in the committed file.
15. **RV sign** — Statcast `delta_run_exp` is batter-positive. We flip. Nestico's write-up uses the other convention. Getting this backwards inverts every usefulness claim.
16. **Whiff definition** — swinging strike including blocked + foul tip. Fouls are swings, not whiffs.
17. **Circular statistics** — never take `mean(spin_axis)` arithmetically.
18. **Primary vs FF `fb_type` join** — downstream efficiency filter must use the same anchor type the axis was measured against (`parachute_ff_fbtype.csv`).
19. **2023 July-14 cutoff** — early scripts dropped 2023 H1 for arm-angle reasons that no longer apply. Full 2023 belongs in the lock.
20. **Article vs memo** — if they disagree, memo wins.
21. **Cherry-picked exemplars** — Medina / Heuer tiny samples; SI-anchored Cease in early 2026 pair tables is not the FF-anchored Cease in the lock.
22. **Ragans / Hellickson branding ≠ locked membership.**
23. **Excluding named pitchers then searching a gate** (Rodón, Roark, Skubal '22) is reverse engineering.
24. **Aggregation artifact** — season-level residual can correlate with velo sep even when pitch-level does not. That is why the lock velocity-adjusts `y_raw`.

---

## 14. README blurb (paste into the new repo)

> Matched-axis ("parachute") changeups: changeups whose imaged spin axis
> matches the pitcher's four-seamer, thrown from a high slot with high
> active spin on both pitches. Industry stuff models grade this shape as a
> slow fastball. This repo measures whether it beats those models on whiff
> and whether that residual is useful in run value.
>
> **Current status (locked spec, arm-clustered): +1.78 whiff points, p = 0.067.
> Not established. Do not ship as a pitch-design recommendation.** The next
> test is pre-registered: add D1 TrackMan 2022 and 2026 under
> `docs/parachute_precommit.md` and confirm / kill / wait.
>
> An unsupervised grid search over arm, axis, and active-spin floors (Aug 11
> 2026) put 2026 Dylan Cease and 2021 Tarik Skubal on top of a 20-season
> bin. That search does not beat a permutation null. The locked definition
> is round a-priori cuts, not the search winner.

---

## 15. First tasks for Claude Code in the new repo

Do these in order. Stop after (4) and show the locked table before writing
any article prose.

1. Copy `spec_lock.R`, the three status memos, `parachute_extended_build.R`,
   `whiff_tjstuff.R`, `whiff_parallel.R`, `active_spin_merge.R`, scrape
   scripts, and `bin_roster_detail.R`. Write `CLAUDE.md` from §0.
2. Vendor active-spin and arsenal CSVs. Gitignore raw Statcast.
3. Reproduce `parachute_ff.rds` → residuals → `locked_spec.rds`. Diff the
   locked table against §3. If it does not match, stop and diagnose
   (anchor, year loop, bit64, usage filter, residual cache).
4. Reproduce `bin_roster_detail.R searched` and confirm Cease 2026 / Skubal
   2021 sit at +11.0 / +10.3 on r4. This is a smoke test of the residual
   pipeline, not a license to adopt the searched gates.
5. Copy `mech_01`–`mech_11` and regenerate the decision memo numbers.
6. Only then: NCAA path, canvases, or the exploratory seam-shift hypothesis.

If the user asks to "prove parachute changeups work," refuse the framing.
Point at the pre-commit. Offer to implement D1 2022/2026 under the frozen
spec, or to write a descriptive prevalence piece that does not claim
performance.
