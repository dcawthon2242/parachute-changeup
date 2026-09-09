#!/usr/bin/env Rscript

# PHASE 1 - Sequential primary-fastball -> breaking-ball tunnel pairs.
#
# Builds the pair-level table behind the "tunnel location sensitivity" work:
#   pitch N   = the pitcher's PRIMARY fastball (most-thrown of FF/SI/FC that season)
#   pitch N+1 = a breaking ball (SL, ST, or CU with KC folded in), 2-strike count
#   consecutive pitch_number inside the same plate appearance
#
# For every pair we keep both plate locations, the pairwise tunnel metric, and the
# miss-distance target on the breaking ball.
#
# The tunnel metric reuses the EXACT definition already used by the article's model
# (baseball/miss_grade_features.R): a weighted 3D integration of the distance between
# the two reconstructed trajectories over the hitter's pre-commitment window, divided
# by the pitches' separation at the plate. Lower ratio = tighter tunnel.

suppressPackageStartupMessages({ library(data.table) })

MDIR      <- file.path("data", "statcast_model")
OUT_RDS   <- file.path(MDIR, "tunnel_pairs.rds")
SEASONS   <- 2023:2026

FASTBALLS <- c("FF","SI","FC")
BRK_MAP   <- c(SL="SL", ST="ST", CU="CU", KC="CU")   # KC folded into CU
WHIFF     <- c("swinging_strike","swinging_strike_blocked","foul_tip")
SWING     <- c(WHIFF,"foul","hit_into_play")

# tunnel-metric constants, identical to baseball/miss_grade_features.R
WX <- 1.2; WY <- 0.4; WZ <- 1.4; REACT <- 0.150; NSTEP <- 40
yf <- 17/12

NEED <- c("pitcher","player_name","pitch_type","game_type","description","game_date",
          "game_pk","at_bat_number","pitch_number","balls","strikes","miss_distance",
          "p_throws","stand","plate_x","plate_z",
          "release_pos_x","release_pos_y","release_pos_z","vx0","vy0","vz0","ax","ay","az")

pos <- function(r, v, a, t) r + v*t + 0.5*a*t^2

load_season <- function(yr) {
  f <- file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr))
  if (!file.exists(f)) { message(sprintf("MISSING season file: %s", f)); return(NULL) }
  d <- fread(f, select = NEED, showProgress = FALSE)
  d[, season := yr]
  d <- d[game_type == "R" & !is.na(vx0) & pitch_type != "" & !is.na(release_pos_y)]
  # bat tracking (miss_distance) began at the 2023 All-Star break
  if (yr == 2023L) { d[, game_date := as.Date(game_date)]; d <- d[game_date >= as.Date("2023-07-14")] }

  # ---- primary fastball for the season: most-thrown of FF/SI/FC ----
  fbc <- d[pitch_type %in% FASTBALLS, .N, by = .(pitcher, pitch_type)][N >= 50]
  prim <- fbc[order(pitcher, -N)][, .SD[1], by = pitcher][, .(pitcher, prim_fb = pitch_type)]
  d <- merge(d, prim, by = "pitcher", all.x = TRUE)

  # ---- lag the previous pitch inside the plate appearance ----
  setorder(d, game_pk, at_bat_number, pitch_number)
  d[, t_plate := (-vy0 - sqrt(vy0^2 - 2*ay*(release_pos_y - yf)))/ay]
  d[, t_react := pmax(t_plate - REACT, 0.05)]
  lagcols <- c("release_pos_x","release_pos_y","release_pos_z","vx0","vy0","vz0",
               "ax","ay","az","t_react","pitch_type","pitch_number","plate_x","plate_z")
  for (cc in lagcols) d[, (paste0("p_", cc)) := shift(get(cc)), by = .(game_pk, at_bat_number)]

  # ---- qualifying pairs: any fastball immediately followed by a breaking ball, 2 strikes ----
  # setup_type is retained so specific pairings (SI->ST, FF->CU, ...) can be pulled out;
  # is_primary_setup flags the subset where the setup was the pitcher's primary fastball.
  d[, brk_type := BRK_MAP[pitch_type]]
  q <- d[!is.na(brk_type) &
         !is.na(p_pitch_number) & (pitch_number - p_pitch_number == 1) &
         p_pitch_type %in% FASTBALLS &
         strikes == 2L]
  if (!nrow(q)) return(NULL)
  q[, setup_type := p_pitch_type]
  q[, is_primary_setup := !is.na(prim_fb) & p_pitch_type == prim_fb]

  # ---- pairwise tunnel: integrate weighted 3D separation over the pre-commit window ----
  Tmax <- pmax(q$t_react, q$p_t_react)
  acc  <- numeric(nrow(q))
  for (k in 1:NSTEP) {
    tk <- (k - 0.5)/NSTEP * Tmax
    dx <- pos(q$release_pos_x, q$vx0, q$ax, tk) - pos(q$p_release_pos_x, q$p_vx0, q$p_ax, tk)
    dy <- pos(q$release_pos_y, q$vy0, q$ay, tk) - pos(q$p_release_pos_y, q$p_vy0, q$p_ay, tk)
    dz <- pos(q$release_pos_z, q$vz0, q$az, tk) - pos(q$p_release_pos_z, q$p_vz0, q$p_az, tk)
    acc <- acc + sqrt(WX*dx^2 + WY*dy^2 + WZ*dz^2) * (Tmax/NSTEP)
  }
  q[, tunnel    := acc]
  q[, plate_sep := pmax(sqrt((plate_x - p_plate_x)^2 + (plate_z - p_plate_z)^2), 0.1)]
  q[, path_ratio := tunnel / plate_sep]

  # ---- miss target on the breaking ball ----
  q[, is_swing := description %in% SWING]
  q[, is_whiff := description %in% WHIFF]
  q <- q[!(is_swing & is_whiff & is.na(miss_distance))]   # whiff w/o tracked miss = unknown
  q[is_swing == TRUE & is_whiff == FALSE, miss_distance := 0]   # contact -> 0 miss
  q[is_swing == FALSE, miss_distance := NA_real_]               # takes have no miss

  out <- q[, .(season, pitcher, player_name, p_throws, stand, brk_type,
               setup_type, is_primary_setup,
               game_pk, at_bat_number, pitch_number,
               fb_plate_x = p_plate_x, fb_plate_z = p_plate_z,
               bb_plate_x = plate_x,   bb_plate_z = plate_z,
               tunnel, plate_sep, path_ratio, is_swing, is_whiff, miss_distance)]
  message(sprintf("  season %d: %d qualifying FB->BB pairs (%.1f%% swings, %.1f%% whiff|swing)",
                  yr, nrow(out), 100*mean(out$is_swing),
                  100*mean(out$is_whiff[out$is_swing])))
  out
}

message("Building primary-FB -> breaking-ball 2-strike pairs...")
pairs <- rbindlist(lapply(SEASONS, load_season), use.names = TRUE, fill = TRUE)
pairs <- pairs[is.finite(path_ratio)]

# ---- well-tunneled gate: tightest tercile of path_ratio within pairing x pitcher hand ----
# Gated per setup->breaking pairing because tunnel tightness is not comparable across
# pairings (an SI->ST pair separates very differently from an FF->CU pair).
pairs[, tun_cut := quantile(path_ratio, 1/3, na.rm = TRUE),
      by = .(setup_type, brk_type, p_throws)]
pairs[, well_tunneled := path_ratio <= tun_cut]

saveRDS(pairs, OUT_RDS)

cat("\n================ PAIR SUMMARY ================\n")
cat(sprintf("total pairs: %d   (%s)\n", nrow(pairs), OUT_RDS))
cat("\nby setup -> breaking pairing:\n")
print(pairs[, .(pairs = .N, swings = sum(is_swing),
                swing_pct = round(100*mean(is_swing),1),
                mean_miss_on_swing = round(mean(miss_distance[is_swing]), 3),
                med_path_ratio = round(median(path_ratio), 3)),
            by = .(setup_type, brk_type)][order(setup_type, brk_type)])
cat("\nprimary-setup subset (the Phase 2 population):\n")
print(pairs[is_primary_setup == TRUE, .(pairs = .N, swings = sum(is_swing)),
            by = .(brk_type, p_throws)][order(brk_type, p_throws)])
cat("\nby season:\n")
print(pairs[, .(pairs = .N, swings = sum(is_swing)), by = season][order(season)])
cat("\n2026 pitcher x type coverage (swings, for the per-pitcher dashboard):\n")
cov <- pairs[season == 2026 & is_swing == TRUE, .N, by = .(pitcher, brk_type)]
print(cov[, .(pitcher_type_units = .N,
              ge20 = sum(N >= 20), ge40 = sum(N >= 40)), by = brk_type][order(brk_type)])
