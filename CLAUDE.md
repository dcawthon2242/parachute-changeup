# CLAUDE.md

Project instructions for Claude Code in this repo. Section 0 below is copied
verbatim from `HANDOFF.md`; the handoff is the long-form reference and wins on
any detail not covered here.

## Layout note

The handoff (§9) proposes numbered subfolders under `baseball/`. This repo keeps
`baseball/` **flat** on purpose: every script addresses data as `data/...` and
sibling scripts as `baseball/...` relative to the repo root, and reorganising
would break those paths. The grouping in §9 is documented as a catalog in
`docs/REPO_GUIDE.md` instead. Run everything from the repo root:
`Rscript baseball/<script>.R`, `.venv-cv/bin/python baseball/<script>.py`.

Large data (`data/statcast_YYYY`, `data/rubber`, `data/obm`, `data/swing_timing`,
`data/open_command`, `*.rds` caches) is gitignored. See `data/README.md`.

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

