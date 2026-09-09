#!/usr/bin/env Rscript

# REBUILD THE PARACHUTE DATASET OVER 2020-2026.
#
# Statcast publishes per-pitch arm_angle from 2020, not 2023 as assumed earlier, so the three
# Hawk-Eye seasons before the current window can use the real measurement rather than the
# release-point reconstruction, which was only good enough to recover two thirds of high-slot
# pitchers. Active-spin leaderboards go back to 2020 as well. That makes 2020-2022 a clean
# extension on all three of the bin's criteria.
#
# The point of the extension is not a smaller p-value on the same data. The bin was found in
# 2023-2026, so 2020-2022 is untouched territory and can be scored as a genuine out-of-sample
# replication against a target that is already written down: +6.9 percentage points of ground
# balls above a shape-and-location model, and no reliable whiff or run-value edge.
#
# This script only builds and caches. The test is in parachute_extended_test.R.

#
# ANCHOR. Pass "FF" to require a genuine four-seamer as the reference pitch, or "primary" for
# the older behaviour of falling back to a sinker or cutter when the pitcher has no four-seam.
# The distinction matters because a sinker's seam-shifted wake gives it movement its spin axis
# does not predict, so "matched spin axis" against a sinker is not the same measurement as it
# is against a four-seamer - and sinker-anchored seasons show a mean axis gap of 18.4 degrees
# against 22.3 for four-seam-anchored ones, which means they are over-selected into any bin
# defined on that gap.

suppressPackageStartupMessages({ library(data.table) })
set.seed(1); options(width = 205)
MDIR <- "data/statcast_model"
ANCHOR <- { a <- commandArgs(TRUE); if (length(a)) a[1] else "FF" }
stopifnot(ANCHOR %in% c("FF","primary"))
FASTBALLS <- if (ANCHOR == "FF") "FF" else c("FF","SI","FC"); yf <- 17/12
OUT <- file.path(MDIR, if (ANCHOR == "FF") "parachute_ff.rds" else "parachute_extended.rds")
YEARS <- 2020:2026
message(sprintf("anchor: %s (reference pitch types: %s)", ANCHOR, paste(FASTBALLS, collapse=",")))

WANT <- c("game_pk","at_bat_number","pitch_number","game_date","game_type","pitcher",
          "player_name","pitch_type","p_throws","stand","balls","strikes","release_speed",
          "release_spin_rate","release_extension","release_pos_x","release_pos_y",
          "release_pos_z","spin_axis","arm_angle","plate_x","plate_z","sz_bot","sz_top",
          "vx0","vy0","vz0","ax","ay","az","description","events","bb_type","launch_speed",
          "launch_angle","estimated_woba_using_speedangle","delta_run_exp")

read_year <- function(y) {
  f <- sprintf("data/statcast_%d/statcast_%d_all.csv", y, y)
  if (!file.exists(f)) { message(sprintf("  %d: MISSING, skipped", y)); return(NULL) }
  have <- intersect(WANT, names(fread(f, nrows = 1)))
  k <- fread(f, select = have, showProgress = FALSE)
  miss <- setdiff(WANT, have); for (m in miss) k[, (m) := NA]
  k[, season := y]
  message(sprintf("  %d: %s pitches%s", y, format(nrow(k), big.mark = ","),
                  if (length(miss)) sprintf(" (missing: %s)", paste(miss, collapse=",")) else ""))
  k
}
message("loading seasons:")
d <- rbindlist(lapply(YEARS, read_year), use.names = TRUE, fill = TRUE)
d <- d[game_type == "R" & !is.na(vx0) & pitch_type != "" & !is.na(release_pos_y)]
d <- unique(d, by = c("game_pk","at_bat_number","pitch_number"))
# 2023 is only usable from mid-July in this project because the arm-angle work started there;
# with 2020-2022 now present there is no reason to keep that restriction, so the full 2023 is in.
message(sprintf("total %s regular-season pitches, %d seasons",
                format(nrow(d), big.mark = ","), uniqueN(d$season)))

cmean <- function(a) { r <- a*pi/180; (atan2(mean(sin(r)), mean(cos(r)))*180/pi) %% 360 }
circd <- function(a,b) { z <- abs(a-b) %% 360; pmin(z, 360-z) }
fb <- d[pitch_type %in% FASTBALLS, .(nfb = .N, fb_velo = mean(release_speed, na.rm=TRUE),
        fb_ax = mean(ax, na.rm=TRUE), fb_az = mean(az, na.rm=TRUE),
        fb_arm = mean(arm_angle, na.rm=TRUE),
        fb_axis = cmean(spin_axis[!is.na(spin_axis)])),
        by = .(pitcher, season, pitch_type)][nfb >= 50]
fb[, rk := match(pitch_type, FASTBALLS)]
fb <- fb[order(pitcher, season, rk)][, .SD[1], by = .(pitcher, season)]
setnames(fb, "pitch_type", "fb_type")
d <- merge(d, fb[, .(pitcher, season, fb_type, fb_velo, fb_ax, fb_az, fb_axis, fb_arm)],
           by = c("pitcher","season"), all.x = TRUE)
d[, `:=`(speed_diff = release_speed - fb_velo, ax_diff = ax - fb_ax, az_diff = az - fb_az,
         axis_diff = circd(spin_axis, fb_axis), arm_diff = abs(arm_angle - fb_arm))]

d[, vyf := -sqrt(pmax(vy0^2 - 2*ay*(50-yf), 0))][, tf := (vyf - vy0)/ay]
d[, `:=`(VAA = atan2(vz0 + az*tf, vyf)*180/pi, HAA = atan2(vx0 + ax*tf, vyf)*180/pi,
         below_zone = as.integer(plate_z < sz_bot), z_rel_bot = plate_z - sz_bot,
         z_rel_top = plate_z - sz_top, sax = sin(spin_axis*pi/180), cax = cos(spin_axis*pi/180))]
d[, `:=`(bs = fifelse(stand == "R", -1, 1), ps = fifelse(p_throws == "R", -1, 1))]
d[, `:=`(plate_x_in = bs*plate_x, HAA_in = bs*HAA, plate_x_arm = ps*plate_x,
         stand_R = as.integer(stand == "R"), throws_R = as.integer(p_throws == "R"),
         same_hand = as.integer(stand == p_throws), rv = -delta_run_exp)]
d <- d[is.finite(rv)]

sw   <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip","hit_into_play",
          "foul_bunt","missed_bunt","bunt_foul_tip")
miss <- c("swinging_strike","swinging_strike_blocked","foul_tip","missed_bunt")
d[, `:=`(is_swing = description %in% sw, whiff = as.integer(description %in% miss),
         is_bip = !is.na(launch_speed) & bb_type != "")]

KEEP <- c("pitcher","player_name","season","pitch_type","p_throws","stand","balls","strikes",
          "release_speed","release_spin_rate","release_extension","release_pos_x",
          "release_pos_z","sax","cax","ax","az","arm_angle","arm_diff","speed_diff","ax_diff",
          "az_diff","axis_diff","plate_x","plate_z","below_zone","VAA","HAA","z_rel_bot",
          "z_rel_top","plate_x_in","plate_x_arm","HAA_in","stand_R","throws_R","same_hand",
          "rv","is_swing","whiff","is_bip","bb_type","launch_speed",
          "estimated_woba_using_speedangle")
saveRDS(d[pitch_type == "CH", ..KEEP], OUT)
# The efficiency filter compares the changeup's active spin to the anchor fastball's, so the
# downstream join needs to know which fastball that was.
fwrite(unique(fb[, .(pitcher, season, fb_type)]),
       file.path(MDIR, sprintf("parachute_%s_fbtype.csv", if (ANCHOR == "FF") "ff" else "extended")))
message(sprintf("\nwrote %s: %s changeups", OUT, format(sum(d$pitch_type == "CH"), big.mark = ",")))
print(d[pitch_type == "CH", .(changeups = .N, with_arm = sum(is.finite(arm_angle)),
        with_axis = sum(is.finite(axis_diff))), by = season][order(season)])
