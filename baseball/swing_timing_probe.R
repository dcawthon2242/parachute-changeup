#!/usr/bin/env Rscript

# Probe 1: what swing-tracking / timing fields are actually populated, and on which
# pitch outcomes. The published timing axis is the batter-relative intercept DEPTH
# (y, inches toward the pitcher); miss_distance is the scalar bat-to-ball miss.
# Question being set up: is there a mistiming band that yields weak contact rather
# than a whiff, and is that band worth more to the pitcher than the whiff?

suppressPackageStartupMessages({ library(data.table) })

cols <- c("game_type","game_year","description","events","bb_type","balls","strikes",
          "launch_speed","launch_angle","delta_run_exp","bat_speed","swing_length",
          "miss_distance","attack_angle","attack_direction","swing_path_tilt",
          "intercept_ball_minus_batter_pos_x_inches",
          "intercept_ball_minus_batter_pos_y_inches")

probe_season <- function(yr) {
  f <- file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr))
  if (!file.exists(f)) { cat(sprintf("\n[%d] missing %s\n", yr, f)); return(NULL) }
  hdr <- names(fread(f, nrows = 0))
  have <- intersect(cols, hdr)
  dt <- fread(f, showProgress = FALSE, select = have)
  dt <- dt[game_type == "R"]

  cat(sprintf("\n=== %d : %s rows (regular season) ===\n", yr, format(nrow(dt), big.mark=",")))
  missing_cols <- setdiff(cols, hdr)
  if (length(missing_cols)) cat("  columns ABSENT from feed:", paste(missing_cols, collapse=", "), "\n")

  swing <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip","hit_into_play",
             "foul_bunt","missed_bunt","bunt_foul_tip")
  sw <- dt[description %in% swing]
  cat(sprintf("  swings: %s\n", format(nrow(sw), big.mark=",")))
  for (v in intersect(c("bat_speed","swing_length","miss_distance","attack_angle",
                        "attack_direction","swing_path_tilt",
                        "intercept_ball_minus_batter_pos_x_inches",
                        "intercept_ball_minus_batter_pos_y_inches"), have)) {
    cat(sprintf("    %-46s %5.1f%% of swings populated\n", v,
                100 * mean(is.finite(sw[[v]]))))
  }
  invisible(sw)
}

sw26 <- probe_season(2026)
invisible(probe_season(2025))
invisible(probe_season(2024))

if (!is.null(sw26)) {
  cat("\n=== 2026: is miss_distance / intercept published on CONTACT, or only on whiffs? ===\n")
  print(sw26[, .(n = .N,
                 pct_miss_distance = round(100*mean(is.finite(miss_distance)),1),
                 pct_intercept_y   = round(100*mean(is.finite(
                   intercept_ball_minus_batter_pos_y_inches)),1),
                 pct_bat_speed     = round(100*mean(is.finite(bat_speed)),1)),
             by = description][order(-n)])

  cat("\n=== 2026: intercept depth (y, in) and miss_distance by batted-ball type ===\n")
  bip <- sw26[description == "hit_into_play" & bb_type != ""]
  print(bip[, .(n = .N,
                depth_y   = round(mean(intercept_ball_minus_batter_pos_y_inches, na.rm=TRUE),1),
                side_x    = round(mean(intercept_ball_minus_batter_pos_x_inches, na.rm=TRUE),1),
                miss_dist = round(mean(miss_distance, na.rm=TRUE),2),
                bat_speed = round(mean(bat_speed, na.rm=TRUE),1),
                ev        = round(mean(launch_speed, na.rm=TRUE),1),
                la        = round(mean(launch_angle, na.rm=TRUE),1)),
            by = bb_type][order(-n)])

  cat("\n=== 2026: same, for whiffs vs fouls (for contrast) ===\n")
  print(sw26[, .(n = .N,
                 depth_y   = round(mean(intercept_ball_minus_batter_pos_y_inches, na.rm=TRUE),1),
                 miss_dist = round(mean(miss_distance, na.rm=TRUE),2),
                 bat_speed = round(mean(bat_speed, na.rm=TRUE),1)),
             by = .(grp = fifelse(description %in% c("swinging_strike","swinging_strike_blocked"),
                                  "whiff", fifelse(description=="foul","foul",
                                  fifelse(description=="hit_into_play","in_play","other"))))])
}
