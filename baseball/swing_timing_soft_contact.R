#!/usr/bin/env Rscript

# Is there a level of MISTIMING that a pitcher should prefer to a whiff?
#
# The published pitch-level timing axis is the batter-relative intercept DEPTH
# (intercept_ball_minus_batter_pos_y_inches): how far out front, in inches, the bat
# crossed the ball. Larger = met the ball earlier in its flight = EARLY; smaller =
# deeper = LATE. This is the only timing field populated on CONTACT (the scalar
# miss_distance is whiffs-only), so it is the axis used here.
#
# Raw depth confounds pitch location, velocity and the hitter's own stance, so a
# residualized "timing deviation" (tdev) is built alongside the raw version:
#   tdev = depth - E[depth | location, velo, pitch group] - batter's own mean residual
# Positive tdev = early relative to that hitter's norm; negative = late.

suppressPackageStartupMessages({ library(data.table); library(splines) })

SEASONS <- c(2025, 2026)
cols <- c("game_year","game_type","batter","pitcher","stand","pitch_type","description",
          "events","bb_type","balls","strikes","outs_when_up","on_1b","plate_x","plate_z",
          "sz_top","sz_bot","release_speed","launch_speed","launch_angle",
          "estimated_woba_using_speedangle","delta_run_exp","bat_speed","swing_length",
          "attack_angle","attack_direction","miss_distance",
          "intercept_ball_minus_batter_pos_x_inches",
          "intercept_ball_minus_batter_pos_y_inches")

dt <- rbindlist(lapply(SEASONS, function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress = FALSE, select = cols)))
setnames(dt, c("intercept_ball_minus_batter_pos_y_inches",
               "intercept_ball_minus_batter_pos_x_inches"), c("depth","side"))

SWING <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip","hit_into_play")
sw <- dt[game_type == "R" & description %in% SWING &
         is.finite(depth) & is.finite(plate_x) & is.finite(plate_z) &
         is.finite(release_speed) & is.finite(delta_run_exp)]

FB <- c("FF","SI","FC"); BR <- c("SL","ST","CU","KC","SV","CS"); OS <- c("CH","FS","FO")
sw[, pgrp := fifelse(pitch_type %in% FB, "FB",
             fifelse(pitch_type %in% BR, "BR", fifelse(pitch_type %in% OS, "OS", NA_character_)))]
sw <- sw[!is.na(pgrp)]

# batter-handed horizontal location (positive = inside) and height within the zone
sw[, px_bat := fifelse(stand == "R", -plate_x, plate_x)]
sw[, pz_rel := (plate_z - sz_bot) / pmax(sz_top - sz_bot, 0.1)]

# ---- outcome flags ---------------------------------------------------------
sw[, whiff    := description %in% c("swinging_strike","swinging_strike_blocked","foul_tip")]
sw[, is_foul  := description == "foul"]
sw[, inplay   := description == "hit_into_play"]
sw[, gb       := inplay & bb_type == "ground_ball"]
sw[, weak_gb  := gb & is.finite(launch_speed) & launch_speed < 85]
sw[, barrel   := inplay & is.finite(launch_speed) & is.finite(launch_angle) &
                 launch_speed >= 98 & launch_angle >= 26 & launch_angle <= 30]
sw[, pit_rv   := -delta_run_exp]           # run value from the PITCHER's side

one <- c("field_out","force_out","sac_fly","sac_bunt","fielders_choice_out","other_out","strikeout")
sw[, outs_made := fifelse(grepl("triple_play", events), 3L,
                  fifelse(grepl("double_play", events), 2L,
                  fifelse(events %in% one, 1L, 0L)))]
sw[, dp_opp := (on_1b > 0 & !is.na(on_1b)) & outs_when_up < 2]

cat(sprintf("Competitive swings %d-%d with a published intercept depth: %s\n",
            min(SEASONS), max(SEASONS), format(nrow(sw), big.mark = ",")))
cat(sprintf("Sanity: mean pitcher RV on whiffs = %+.4f, on balls in play = %+.4f\n",
            sw[whiff == TRUE, mean(pit_rv)], sw[inplay == TRUE, mean(pit_rv)]))

# ---- residualized timing deviation ----------------------------------------
fit <- lm(depth ~ ns(px_bat, 5) * pgrp + ns(pz_rel, 5) + ns(release_speed, 4) + stand,
          data = sw)
cat(sprintf("Depth model (location + velo + pitch group): R2 = %.3f, residual SD = %.2f in\n",
            summary(fit)$r.squared, sd(residuals(fit))))
sw[, r1 := residuals(fit)]
sw[, nbat := .N, by = batter]
sw[, tdev := r1 - mean(r1), by = batter]
sw <- sw[nbat >= 200]

# ---- the curve -------------------------------------------------------------
BR_EDGES <- c(-Inf, -12, -9, -6, -3, 0, 3, 6, 9, 12, Inf)
LABS <- c("<-12 (very late)","-12..-9","-9..-6","-6..-3","-3..0 (~on time)",
          "0..3 (~on time)","3..6","6..9","9..12","12+ (very early)")
sw[, bin := cut(tdev, BR_EDGES, labels = LABS)]

summarise <- function(d) d[, .(
  n          = .N,
  swing_pct  = round(100 * .N / nrow(d), 1),
  whiff      = round(100 * mean(whiff), 1),
  foul       = round(100 * mean(is_foul), 1),
  inplay     = round(100 * mean(inplay), 1),
  gb         = round(100 * mean(gb), 1),
  weak_gb    = round(100 * mean(weak_gb), 1),
  barrel     = round(100 * mean(barrel), 1),
  ev_con     = round(mean(launch_speed[inplay], na.rm = TRUE), 1),
  xwobacon   = round(mean(estimated_woba_using_speedangle[inplay], na.rm = TRUE), 3),
  outs_sw    = round(mean(outs_made), 3),
  pit_rv     = round(mean(pit_rv), 4)
), by = bin][order(bin)]

cat("\n############ ALL COUNTS: outcome mix and pitcher run value by timing deviation ############\n")
cat("(tdev in inches; negative = LATE / deeper contact, positive = EARLY / out front)\n")
print(summarise(sw))

for (lab in c("0-1 strikes", "2 strikes")) {
  d <- if (lab == "2 strikes") sw[strikes == 2] else sw[strikes < 2]
  cat(sprintf("\n############ %s (n=%s) ############\n", lab, format(nrow(d), big.mark=",")))
  print(summarise(d))
}

cat("\n############ Double plays: does mistimed contact buy the extra out? ############\n")
cat("(runner on 1B, fewer than 2 outs)\n")
print(sw[dp_opp == TRUE, .(n = .N,
      whiff = round(100*mean(whiff),1),
      gb = round(100*mean(gb),1),
      dp = round(100*mean(grepl("double_play", events)),2),
      outs_sw = round(mean(outs_made),3),
      pit_rv = round(mean(pit_rv),4)), by = bin][order(bin)])

cat("\n############ Where is each objective maximised? ############\n")
obj <- function(d, lab) {
  s <- summarise(d)[n >= 2000]
  cat(sprintf("  %-12s  best RV bin: %-18s (RV %+.4f)   |  best whiff bin: %-18s (whiff %.1f%%)   |  best weak-GB bin: %-18s (%.1f%%)\n",
    lab,
    as.character(s$bin[which.max(s$pit_rv)]),  max(s$pit_rv),
    as.character(s$bin[which.max(s$whiff)]),   max(s$whiff),
    as.character(s$bin[which.max(s$weak_gb)]), max(s$weak_gb)))
}
obj(sw, "all counts"); obj(sw[strikes < 2], "0-1 strikes"); obj(sw[strikes == 2], "2 strikes")

cat("\n############ Cross-check on RAW contact depth (no adjustment) ############\n")
sw[, dbin := cut(depth, c(-Inf, 12, 18, 24, 27, 30, 33, 36, 42, Inf))]
print(sw[, .(n = .N, whiff = round(100*mean(whiff),1), gb = round(100*mean(gb),1),
             weak_gb = round(100*mean(weak_gb),1), barrel = round(100*mean(barrel),1),
             ev_con = round(mean(launch_speed[inplay], na.rm=TRUE),1),
             xwobacon = round(mean(estimated_woba_using_speedangle[inplay], na.rm=TRUE),3),
             pit_rv = round(mean(pit_rv),4)), by = dbin][order(dbin)])

out <- file.path("data","statcast_model","swing_timing_soft_contact.csv")
fwrite(summarise(sw), out)
cat(sprintf("\nWrote %s\n", out))
