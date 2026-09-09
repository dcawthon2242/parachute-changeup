#!/usr/bin/env Rscript

# Build the modeling table for the miss-distance pitch grade model.
#
# Target : miss_distance on competitive swings (bunts excluded)  = E[miss | swing]
# Train  : 2023 (H2, via non-NA target) + 2024 + 2025
# Holdout: 2026
#
# Features:
#   TJStuff+ v3.0 shape set: release_speed, release_spin_rate, release_extension,
#     release_pos_x, release_pos_z, spin_axis(sin/cos), ax, az,
#     speed_diff, ax_diff, az_diff  (diffs vs the pitcher's primary fastball, per season)
#   Novel (for the augmentation test, carried but toggled later):
#     axis_diff, spin_eff, spin_eff_diff  (per-pitch, vs primary FB)
#     path_ratio                          (arsenal-level: mean path_to_location_ratio
#                                          on breaking/offspeed-after-FB sequences)
#
# Output: data/statcast_model/miss_grade_data.rds

suppressPackageStartupMessages({ library(data.table) })

OUT_DIR <- file.path("data", "statcast_model")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SEASONS_TRAIN <- c(2023L, 2024L, 2025L)
SEASON_TEST   <- 2026L
FASTBALLS <- c("FF","SI","FC")
BREAKING  <- c("SL","ST","CU","KC","SV","CS")
OFFSPEED  <- c("CH","FS","FO")

# Statcast only records miss_distance on WHIFFS (NA on all contact). To realize
# the E[miss | swing] target the model treats contact as a ~0 miss: competitive
# swings = whiffs + fouls + balls-in-play (bunts excluded), whiff misses kept,
# contact misses imputed to 0. Whiffs whose bat-tracking miss wasn't captured are
# dropped (unknown miss, can't assume 0).
WHIFF <- c("swinging_strike","swinging_strike_blocked","foul_tip")
SWING <- c(WHIFF,"foul","hit_into_play")

NEED <- c("pitcher","player_name","pitch_type","game_type","description","game_date",
  "game_pk","at_bat_number","pitch_number","balls","strikes","miss_distance",
  "release_speed","release_spin_rate","release_extension",
  "release_pos_x","release_pos_y","release_pos_z","spin_axis",
  "vx0","vy0","vz0","ax","ay","az","plate_x","plate_z")

g <- 32.174; yf <- 17/12
cmean <- function(a){a<-a[!is.na(a)]; if(!length(a)) return(NA_real_); r<-a*pi/180; ((atan2(mean(sin(r)),mean(cos(r)))*180/pi)+360)%%360}
circd <- function(a,b){d<-abs(a-b)%%360; pmin(d,360-d)}

load_season <- function(yr) {
  f <- file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr))
  if (!file.exists(f)) { message(sprintf("MISSING season file: %s", f)); return(NULL) }
  d <- fread(f, showProgress = FALSE)
  keep <- intersect(NEED, names(d))
  d <- d[, ..keep]
  d[, season := yr]
  d <- d[game_type == "R" & !is.na(vx0) & pitch_type != "" & !is.na(release_pos_y)]
  # Bat tracking (miss_distance) only began at the 2023 All-Star break; restrict
  # 2023 to its tracked half so contact and whiff swings come from the same period.
  if (yr == 2023L) { d[, game_date := as.Date(game_date)]; d <- d[game_date >= as.Date("2023-07-14")] }

  # spin efficiency (Nathan Magnus decomposition), per pitch
  d[, tt := (-vy0 - sqrt(vy0^2 - 2*ay*(50 - yf)))/ay]
  d[, `:=`(vxm = vx0 + ax*tt/2, vym = vy0 + ay*tt/2, vzm = vz0 + az*tt/2)]
  d[, vmag := sqrt(vxm^2 + vym^2 + vzm^2)]
  d[, dotp := (ax*vxm + ay*vym + (az+g)*vzm)/vmag^2]
  d[, amag := sqrt((ax-dotp*vxm)^2 + (ay-dotp*vym)^2 + ((az+g)-dotp*vzm)^2)]
  d[, spin_eff := fifelse(!is.na(release_spin_rate) & release_spin_rate > 0, amag/(vmag*release_spin_rate), NA_real_)]

  # ---- primary fastball anchor per pitcher x season ----
  fb <- d[pitch_type %in% FASTBALLS,
    .(nfb = .N, fb_velo = mean(release_speed, na.rm=TRUE),
      fb_ax = mean(ax, na.rm=TRUE), fb_az = mean(az, na.rm=TRUE),
      fb_axis = cmean(spin_axis), fb_eff = mean(spin_eff, na.rm=TRUE)),
    by = .(pitcher, pitch_type)]
  fb[, rank := fifelse(pitch_type=="FF",1L, fifelse(pitch_type=="SI",2L,3L))]
  fb <- fb[nfb >= 50][order(pitcher, rank)][, .SD[1], by = pitcher][
    , .(pitcher, fb_velo, fb_ax, fb_az, fb_axis, fb_eff)]
  d <- merge(d, fb, by = "pitcher", all.x = TRUE)
  d[, `:=`(speed_diff = release_speed - fb_velo,
           ax_diff = ax - fb_ax, az_diff = az - fb_az,
           axis_diff = circd(spin_axis, fb_axis),
           spin_eff_diff = abs(spin_eff - fb_eff))]

  # ---- path_to_location_ratio (tunnel vs previous pitch), arsenal-level ----
  WX <- 1.2; WY <- 0.4; WZ <- 1.4; REACT <- 0.150; NSTEP <- 40
  setorder(d, game_pk, at_bat_number, pitch_number)
  d[, t_plate := (-vy0 - sqrt(vy0^2 - 2*ay*(release_pos_y - yf)))/ay]
  d[, t_react := pmax(t_plate - REACT, 0.05)]
  lagcols <- c("release_pos_x","release_pos_z","release_pos_y","vx0","vy0","vz0",
               "ax","ay","az","t_react","pitch_type","pitch_number","plate_x","plate_z")
  for (c in lagcols) d[, (paste0("p_",c)) := shift(get(c)), by = .(game_pk, at_bat_number)]
  d[, has_prev := !is.na(p_pitch_number) & (pitch_number - p_pitch_number == 1)]
  pos <- function(r,v,a,t) r + v*t + 0.5*a*t^2
  d[, tunnel := NA_real_]
  idx <- which(d$has_prev)
  if (length(idx)) {
    Tmax <- pmax(d$t_react[idx], d$p_t_react[idx]); acc <- numeric(length(idx))
    for (k in 1:NSTEP) {
      tk <- (k-0.5)/NSTEP * Tmax
      dx <- pos(d$release_pos_x[idx], d$vx0[idx], d$ax[idx], tk) - pos(d$p_release_pos_x[idx], d$p_vx0[idx], d$p_ax[idx], tk)
      dy <- pos(d$release_pos_y[idx], d$vy0[idx], d$ay[idx], tk) - pos(d$p_release_pos_y[idx], d$p_vy0[idx], d$p_ay[idx], tk)
      dz <- pos(d$release_pos_z[idx], d$vz0[idx], d$az[idx], tk) - pos(d$p_release_pos_z[idx], d$p_vz0[idx], d$p_az[idx], tk)
      acc <- acc + sqrt(WX*dx^2 + WY*dy^2 + WZ*dz^2) * (Tmax/NSTEP)
    }
    d$tunnel[idx] <- acc
  }
  d[, location := pmax(sqrt((plate_x - p_plate_x)^2 + (plate_z - p_plate_z)^2), 0.1)]
  d[, path_to_location_ratio := tunnel / location]
  d[, after_fb := p_pitch_type %in% FASTBALLS]

  # arsenal-level mean ratio per pitcher x pitch_type on secondary-after-FB seqs
  pr <- d[after_fb == TRUE & is.finite(path_to_location_ratio) &
          pitch_type %in% c(BREAKING, OFFSPEED),
          .(path_ratio = mean(path_to_location_ratio), npr = .N),
          by = .(pitcher, pitch_type)][npr >= 10]
  d <- merge(d, pr[, .(pitcher, pitch_type, path_ratio)], by = c("pitcher","pitch_type"), all.x = TRUE)

  d[, grp := fifelse(pitch_type %in% FASTBALLS, "fastball",
             fifelse(pitch_type %in% BREAKING, "breaking",
             fifelse(pitch_type %in% OFFSPEED, "offspeed", "other")))]
  d[, is_swing := description %in% SWING]
  d[, is_whiff := description %in% WHIFF]
  d <- d[is_swing == TRUE & !(is_whiff & is.na(miss_distance))]  # drop whiffs w/o tracked miss
  d[is.na(miss_distance), miss_distance := 0]                     # contact -> 0 miss

  out <- d[, .(season, pitcher, player_name, pitch_type, grp, game_pk, at_bat_number, pitch_number,
      balls, strikes, miss_distance, is_whiff,
      release_speed, release_spin_rate, release_extension, release_pos_x, release_pos_z,
      spin_axis, ax, az, speed_diff, ax_diff, az_diff,
      axis_diff, spin_eff, spin_eff_diff, path_ratio)]
  message(sprintf("  season %d: %d competitive swings (%.1f%% whiffs), mean miss=%.2f",
                  yr, nrow(out), 100*mean(out$is_whiff), mean(out$miss_distance)))
  out
}

all_rows <- rbindlist(lapply(c(SEASONS_TRAIN, SEASON_TEST), load_season), use.names = TRUE, fill = TRUE)
all_rows[, sax := sin(spin_axis*pi/180)]
all_rows[, cax := cos(spin_axis*pi/180)]
all_rows[, set := fifelse(season == SEASON_TEST, "holdout", "train")]

saveRDS(all_rows, file.path(OUT_DIR, "miss_grade_data.rds"))
message(sprintf("\nWrote %s : %d rows (train=%d, holdout=%d)",
  file.path(OUT_DIR,"miss_grade_data.rds"), nrow(all_rows),
  sum(all_rows$set=="train"), sum(all_rows$set=="holdout")))
message("Rows missing primary-FB diffs (no qualifying FB): ",
  all_rows[is.na(speed_diff), .N])
