#!/usr/bin/env Rscript

# ANGULAR TUNNELING: the fastball/breaking-ball tunnel measured as the hitter's
# eye actually receives it, rather than as a weighted 3D distance in field
# coordinates.
#
# The existing metric (driver_gain_all_types.R) integrates the separation between
# two trajectories under hand-tuned weights WX=1.2, WY=0.4, WZ=1.4, which stand in
# for perceptual salience. With a hitter eye point in the same frame -- recovered
# from OpenBiomechanics in baseball/obm_eye_point.py and height-adjusted per batter
# in baseball/obm_hitter_eye_points.py -- that weighting can be replaced by the
# real quantity: the angle the two pitches subtend at the eye.
#
# THE BALL IS A DISC, NOT A POINT
# -------------------------------
# A baseball is 2.9 in across, so at the commit point it subtends about 0.6 deg --
# far from negligible. Each pitch is therefore modelled as an angular disc:
#
#   angular radius   rho(t) = asin(r_ball / d(t))        d = eye-to-ball distance
#   angular gap      theta(t) = angle between the two eye->ball unit vectors
#   EDGE-to-edge gap g_pos(t) = theta(t) - rho_A(t) - rho_B(t)
#
# g_pos <= 0 means the two discs overlap on the retina: the hitter cannot separate
# them at all. This is the "two sides of the baseball" formulation.
#
# DEPTH CUE
# ---------
# Angular diameter also encodes distance, so two pitches at the same flight time
# that are at different depths look different sizes even if their centres coincide:
#
#   g_diam(t) = |2*rho_A(t) - 2*rho_B(t)|
#
# ACUITY THRESHOLD (the tuned parameter)
# --------------------------------------
# Neither cue is usable below the hitter's discrimination threshold tau. A pitch
# pair is DISTINGUISHABLE at time t when g_pos(t) > tau OR g_diam(t) > tau. The
# tunnel then ends at the first such t, and later is better. tau is swept rather
# than assumed: foveal acuity is ~0.017 deg but dynamic acuity against a target
# moving this fast is far coarser, so the useful value is an empirical question.
# The sweep is scored by how well the resulting metric predicts whiffs.

suppressPackageStartupMessages(library(data.table))
options(width = 210)

FASTBALLS <- c("FF","SI","FC"); yf <- 17/12
REACT <- 0.150; NSTEP <- 80
R_BALL <- (2.9/2)/12                    # baseball radius, feet
TAUS <- c(0.01, 0.02, 0.05, 0.10, 0.20, 0.40, 0.80)   # degrees
DEG <- 180/pi

cols <- c("game_pk","at_bat_number","pitch_number","pitch_type","game_type","pitcher",
          "player_name","p_throws","batter","stand","balls","strikes","plate_x","plate_z",
          "release_pos_x","release_pos_y","release_pos_z","vx0","vy0","vz0","ax","ay","az",
          "description","zone","miss_distance","delta_pitcher_run_exp")

cf <- list.files(file.path("data","statcast_2026","chunks"), pattern = "^chunk_.*csv$", full.names = TRUE)
k <- rbindlist(lapply(cf, function(f) fread(f, select = cols, showProgress = FALSE)))
k <- k[game_type == "R" & !is.na(vx0) & !is.na(release_pos_y) & pitch_type != "" &
         stand %in% c("L","R") & p_throws %in% c("L","R")]
k <- unique(k, by = c("game_pk","at_bat_number","pitch_number"))
setorder(k, game_pk, at_bat_number, pitch_number)

k[, t_plate := (-vy0 - sqrt(vy0^2 - 2*ay*(release_pos_y - yf)))/ay]
k[, t_react := pmax(t_plate - REACT, 0.05)]
lagc <- c("release_pos_x","release_pos_y","release_pos_z","vx0","vy0","vz0","ax","ay","az",
          "t_react","pitch_type","pitch_number")
for (cc in lagc) k[, (paste0("p_", cc)) := shift(get(cc)), by = .(game_pk, at_bat_number)]
k[, consecutive := !is.na(p_pitch_number) & (pitch_number - p_pitch_number == 1L)]
k[, after_fb := consecutive & p_pitch_type %in% FASTBALLS]

# ---- hitter eye points, height-adjusted per batter ------------------------
eye <- fread("data/obm/hitter_eye_points_2026_by_batter.csv")
eye <- eye[, .(batter, stand, eye_x = eye_x_obm, eye_y, eye_z)]
before <- nrow(k)
k <- merge(k, eye, by = c("batter","stand"), all.x = TRUE)
# batters without a height fall back to the side average
side_avg <- fread("data/obm/hitter_eye_points_summary.csv")
for (s in c("L","R")) {
  a <- side_avg[side == s]
  k[stand == s & is.na(eye_x), `:=`(eye_x = a$eye_x_obm, eye_y = a$eye_y, eye_z = a$eye_z)]
}
setorder(k, game_pk, at_bat_number, pitch_number)
cat(sprintf("pitches: %s (eye point attached to %.1f%%)\n",
            format(nrow(k), big.mark = ","), 100*mean(!is.na(k$eye_x))))

# ---- restrict to the pairs we care about: breaking ball after a fastball ---
d <- k[after_fb == TRUE & pitch_type %in% c("ST","SL","CU","KC","CH","FS")]
cat(sprintf("secondary-after-fastball pairs: %s\n", format(nrow(d), big.mark = ",")))

# ---- angular geometry over the flight -------------------------------------
pos <- function(r, v, a, t) r + v*t + 0.5*a*t^2
n <- nrow(d)
Tmax <- pmax(d$t_react, d$p_t_react)

# first step (per pair) at which the pair becomes distinguishable, per tau
first_pos  <- matrix(NA_integer_, n, length(TAUS))
first_any  <- matrix(NA_integer_, n, length(TAUS))
overlap_steps <- integer(n)      # steps where the discs physically overlap

for (s in 1:NSTEP) {
  tk <- (s - 0.5)/NSTEP * Tmax
  ax_ <- pos(d$release_pos_x, d$vx0, d$ax, tk) - d$eye_x
  ay_ <- pos(d$release_pos_y, d$vy0, d$ay, tk) - d$eye_y
  az_ <- pos(d$release_pos_z, d$vz0, d$az, tk) - d$eye_z
  bx_ <- pos(d$p_release_pos_x, d$p_vx0, d$p_ax, tk) - d$eye_x
  by_ <- pos(d$p_release_pos_y, d$p_vy0, d$p_ay, tk) - d$eye_y
  bz_ <- pos(d$p_release_pos_z, d$p_vz0, d$p_az, tk) - d$eye_z

  da <- sqrt(ax_^2 + ay_^2 + az_^2); db <- sqrt(bx_^2 + by_^2 + bz_^2)
  cosang <- (ax_*bx_ + ay_*by_ + az_*bz_)/(da*db)
  cosang <- pmin(1, pmax(-1, cosang))
  theta <- acos(cosang) * DEG
  rho_a <- asin(pmin(1, R_BALL/da)) * DEG
  rho_b <- asin(pmin(1, R_BALL/db)) * DEG

  g_pos  <- theta - rho_a - rho_b          # edge-to-edge angular gap
  g_diam <- abs(2*rho_a - 2*rho_b)         # angular-size (depth) difference

  overlap_steps <- overlap_steps + (g_pos <= 0)
  for (j in seq_along(TAUS)) {
    tau <- TAUS[j]
    hit <- is.na(first_pos[, j]) & (g_pos > tau)
    first_pos[hit, j] <- s
    hit2 <- is.na(first_any[, j]) & ((g_pos > tau) | (g_diam > tau))
    first_any[hit2, j] <- s
  }
}

# Convert "first distinguishable step" into a break fraction of the decision
# window. 1 = stayed together all the way to the commit point (perfect tunnel).
for (j in seq_along(TAUS)) {
  set(d, j = sprintf("brk_pos_%03d", round(TAUS[j]*100)),
      value = ifelse(is.na(first_pos[, j]), 1, (first_pos[, j] - 0.5)/NSTEP))
  set(d, j = sprintf("brk_any_%03d", round(TAUS[j]*100)),
      value = ifelse(is.na(first_any[, j]), 1, (first_any[, j] - 0.5)/NSTEP))
}
set(d, j = "overlap_frac", value = overlap_steps/NSTEP)

wd <- c("swinging_strike","swinging_strike_blocked","foul_tip","missed_bunt","bunt_foul_tip")
sd <- c(wd, "foul","hit_into_play","hit_into_play_score","hit_into_play_no_out","foul_bunt")
d[, `:=`(swing = description %in% sd, whiff = description %in% wd)]

cat("\n=== 1. HOW LONG THE PAIR STAYS INDISTINGUISHABLE, by pitch type ===\n")
cat("    (break fraction: share of the decision window survived before the hitter\n")
cat("     can separate the two pitches; 1.0 = tunneled all the way to commit)\n\n")
bcols <- grep("^brk_any_", names(d), value = TRUE)
print(d[, c(.(n = .N, overlap = round(mean(overlap_frac), 3)),
            lapply(.SD, function(z) round(mean(z), 3))),
        by = pitch_type, .SDcols = bcols][order(-n)])

cat("\n  positional cue alone (no depth cue):\n")
pcols <- grep("^brk_pos_", names(d), value = TRUE)
print(d[, c(.(n = .N), lapply(.SD, function(z) round(mean(z), 3))),
        by = pitch_type, .SDcols = pcols][order(-n)])

## ---- 2. which tau carries the most signal? ------------------------------
# Score each threshold by how well the break fraction separates whiffs from
# non-whiffs on competitive swings, at the pitcher x type level (the unit
# fig11_rebuild.R established as the honest one).
cat("\n\n=== 2. CHOOSING tau: correlation of break fraction with whiff rate ===\n")
cat("    (pitcher x pitch-type level, min 40 swings; + = better tunnel goes with more whiffs)\n\n")
sw <- d[swing == TRUE]
res <- rbindlist(lapply(c(bcols, pcols), function(v) {
  a <- sw[, .(n = .N, wh = mean(whiff), br = mean(get(v))), by = .(pitcher, pitch_type)][n >= 40]
  a <- a[is.finite(br) & is.finite(wh)]
  rbindlist(lapply(c("ST","SL","CU","CH","FS"), function(ty) {
    b <- a[pitch_type == ty]
    if (nrow(b) < 25) return(NULL)
    ct <- suppressWarnings(cor.test(b$br, b$wh, method = "spearman"))
    data.table(metric = v, ptype = ty, n_arms = nrow(b),
               rho = unname(ct$estimate), p = ct$p.value)
  }))
}))
print(dcast(res, metric ~ ptype, value.var = "rho")[
  , lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)])
cat("\n  p-values:\n")
print(dcast(res, metric ~ ptype, value.var = "p")[
  , lapply(.SD, function(z) if (is.numeric(z)) signif(z, 2) else z)])

# For reference: the OLD metric on the identical rows.
cat("\n  for reference, the existing weighted-distance path ratio on the same rows:\n")
WX <- 1.2; WY <- 0.4; WZ <- 1.4
acc <- numeric(n)
for (s in 1:NSTEP) {
  tk <- (s - 0.5)/NSTEP * Tmax
  dx <- pos(d$release_pos_x,d$vx0,d$ax,tk) - pos(d$p_release_pos_x,d$p_vx0,d$p_ax,tk)
  dy <- pos(d$release_pos_y,d$vy0,d$ay,tk) - pos(d$p_release_pos_y,d$p_vy0,d$p_ay,tk)
  dz <- pos(d$release_pos_z,d$vz0,d$az,tk) - pos(d$p_release_pos_z,d$p_vz0,d$p_az,tk)
  acc <- acc + sqrt(WX*dx^2 + WY*dy^2 + WZ*dz^2) * (Tmax/NSTEP)
}
d[, tunnel_old := acc]
sw <- d[swing == TRUE]
old <- rbindlist(lapply(c("ST","SL","CU","CH","FS"), function(ty) {
  b <- sw[pitch_type == ty, .(n = .N, wh = mean(whiff), tv = mean(tunnel_old)),
          by = .(pitcher, pitch_type)][n >= 40]
  if (nrow(b) < 25) return(NULL)
  ct <- suppressWarnings(cor.test(-b$tv, b$wh, method = "spearman"))
  data.table(ptype = ty, n_arms = nrow(b), rho_neg_tunnel = round(unname(ct$estimate), 3),
             p = signif(ct$p.value, 2))
}))
print(old)

saveRDS(d, file.path("data","statcast_model","angular_tunnel_2026.rds"))
cat("\nwrote data/statcast_model/angular_tunnel_2026.rds\n")
