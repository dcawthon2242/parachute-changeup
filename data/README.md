# data/

What lives here, what is tracked, and how to rebuild what is not.

## Tracked (vendored, small)

| Path | What | Source |
|---|---|---|
| `active_spin/active_spin_2020..2026.csv` | Savant **measured** active-spin leaderboards, wide, percent | Hand-pulled from Baseball Savant; no downloader in repo |
| `savant/arsenal_2020..2026.csv` | Arsenal usage (`n_ff`, `n_si`, …) for the four-seam-primary filter | Baseball Savant |
| `pitcher_heights.csv` | Listed heights for arm-angle reconstruction | MLB Stats API |
| `mlb_stats/pitching_2024..2026.json` | Season pitching lines for FIP/ERA comparisons | MLB Stats API |
| `ncaa_rosters/` | D1 roster scrapes for NCAA name resolution | `baseball/scrape_ncaa_rosters.py`, `sidearm_*.py` |
| `statcast_model/article_assets/` | Memos, figure captions, PNG figures, `ext_*.csv` tables, two embedded GIFs, clip manifests | Written by the analysis scripts |
| `statcast_model/*.csv` (small) | Roster tables, fbtype joins, misc outputs | Written by the analysis scripts |

## Not tracked (large; symlinked or hard-linked to the original workspace)

| Path | Size | Rebuild |
|---|---|---|
| `statcast_2015 … statcast_2026/` | ~0.3–0.9 GB each | `Rscript baseball/scrape_statcast_multi.R 2020 2021 … 2026` (direct Savant CSV, resume-safe chunks) |
| `statcast_model/*.rds` (104 caches) | 2.4 GB | Rerun the producing script; most skip if the RDS exists, set `REFIT=1` to force |
| `statcast_model/tunnel_location/` | | `baseball/tunnel_location_build.R` |
| `rubber/` | 3.2 GB | `baseball/rubber_01_playids.py` → `rubber_03_fetch_frames.py` → … |
| `obm/` | 1.2 GB | OpenBiomechanics hitting C3D trials; `baseball/obm_*.py` |
| `swing_timing/` | 377 MB | `baseball/swing_timing_fetch.R` and downstream |
| `open_command/` | 444 MB | tomdoyo/open-command releases, per-year `targets.csv.gz` + `pbp_info.csv.gz` |
| `ncaa/` | | Raw D1 TrackMan CSVs (`pbp23tm.csv`, `D1TM24.csv`, `D1TM25.csv`); scripts currently read them from `~/Downloads/` |

On the machine this repo was created on, the untracked paths are **symlinks**
into `/Users/damon.cawthon/System1MapquestWorkspace/data/` (directories) or
**hard links** (`.rds` files), so the scripts run unchanged and share caches
with the old workspace. On a fresh clone they are simply absent: rebuild them
in the order given in `README.md` § Reproduction.

## Rules that live here

- Missing arsenal **row** → drop the pitcher-season. Blank **cell** → 0 usage. Never impute.
- Measured active spin for MLB gates. Inferred efficiency only for the D1 port.
- `active_spin_merge.R` must loop `2020:2026` (the committed file once said `2023:2026`).
