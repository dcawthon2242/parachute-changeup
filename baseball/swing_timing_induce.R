#!/usr/bin/env Rscript

# Pitcher-facing swing timing: which pitches induce tied-up / flail / over / under
# swings, and what is each category actually worth in contact quality?
#
# Unit of analysis is the cell (pitcher, pitch type, batter side, season). That is
# the grain Savant publishes, and it is forced on us: the nine timing categories
# are defined on the BARREL-to-ball distance at closest approach
#   tied-up / centered / flail   centered = barrel within +/-4 in horizontally
#   late / on-time / early       on time  = barrel within +/-7 milliseconds
#   over / lined-up / under      lined up = barrel within +/-2 in vertically
# and none of those components appears in the per-pitch feed, which carries only
# the scalar miss_distance (whiffs only) plus the ball's position relative to the
# batter's body. So the categories come from the leaderboard and everything else
# -- outcomes, pitch characteristics -- is aggregated from per-pitch data to match.
#
# Two things are deliberately kept separate:
#   1. What each category is WORTH  (category rates -> xwOBACON, run value)
#   2. What INDUCES each category   (pitch characteristics -> category rates)
# Both are run within pitch type and batter side, because raw cross-sectional
# correlations are dominated by the fact that changeups differ from fastballs.
#
# Inputs : data/swing_timing/swing_timing_pitcher_<season>.csv  (swing_timing_fetch.R)
#          data/statcast_<season>/statcast_<season>_all.csv
# Output : data/swing_timing/swing_timing_cells.csv

suppressPackageStartupMessages({ library(data.table) })

SEASONS <- 2023:2026
MIN_SWINGS <- 50L    # binomial noise on a 9-way split gets ugly below this
MIN_BIP    <- 15L    # xwOBACON floor; cells below this are kept but not scored
OUT_DIR <- file.path("data", "swing_timing")

RATES <- c("tied_up_percent", "centered_percent", "flailed_percent",
           "early_percent", "on_time_percent", "late_percent",
           "over_percent", "lined_up_percent", "under_percent")

# ---------------------------------------------------------------- leaderboard --
lb <- rbindlist(lapply(SEASONS, function(y)
  fread(file.path(OUT_DIR, sprintf("swing_timing_pitcher_%d.csv", y)),
        showProgress = FALSE)), fill = TRUE)
# Trap in the pitcher view of this leaderboard. Splitting on bat_side returns BOTH
# a `bat_side` and a `bat_side_formatted`, and despite the name the formatted one
# is not the batter's side at all -- it is the PITCHER's throwing hand (verified
# at 100% agreement against p_throws, and its 74/26 R/L marginal is the league
# split of right-handed pitchers, not of batters). The raw `bat_side` is the
# actual split key and sits at the expected ~50/50. Joining on the formatted
# column silently pairs each leaderboard row with the wrong platoon side, which
# is easy to miss because the join still succeeds for every row.
setnames(lb, c("bat_side_formatted", "api_pitch_type"), c("p_throws_lb", "ptype"))
lb <- lb[n_swings >= MIN_SWINGS]

# ------------------------------------------------------------------ per pitch --
cols <- c("game_date", "game_year", "game_type", "pitcher", "player_name", "stand",
          "p_throws", "pitch_type", "description", "events",
          "estimated_woba_using_speedangle", "delta_run_exp", "launch_speed",
          "release_speed", "release_spin_rate", "release_extension",
          "release_pos_x", "release_pos_z", "arm_angle", "pfx_x", "pfx_z",
          "plate_x", "plate_z", "sz_top", "sz_bot",
          "vx0", "vy0", "vz0", "ax", "ay", "az", "bat_speed")

pp <- rbindlist(lapply(SEASONS, function(y) {
  x <- fread(file.path("data", sprintf("statcast_%d", y),
                       sprintf("statcast_%d_all.csv", y)),
             select = cols, showProgress = FALSE)
  # Swing timing begins in the second half of 2023, so the 2023 leaderboard
  # covers only part of the season. Joining a full-season outcome onto a
  # half-season category rate would be a silent mismatch, so 2023 per-pitch data
  # is cut to the window where bat tracking actually exists.
  if (y == 2023L) x <- x[game_date >= x[is.finite(bat_speed), min(game_date)]]
  x
}))

SWING <- c("swinging_strike", "swinging_strike_blocked", "foul", "foul_tip",
           "hit_into_play")
WHIFF <- c("swinging_strike", "swinging_strike_blocked", "foul_tip")
# The leaderboard folds curve variants into CU.
pp[, ptype := fifelse(pitch_type %in% c("KC", "CS"), "CU", pitch_type)]
pp <- pp[game_type == "R" & ptype != ""]

# geometry: approach angles at the front of the plate
pp[, vyf := -sqrt(pmax(vy0^2 - 2 * ay * (50 - 17/12), 1e-6))]
pp[, tt  := (vyf - vy0) / ay]
pp[, vaa := atan2(vz0 + az * tt, abs(vyf)) * 180 / pi]
pp[, haa := atan2(vx0 + ax * tt, abs(vyf)) * 180 / pi]
# Sign horizontal quantities from the batter's point of view so that platoon
# splits are comparable: positive = toward the batter's inside.
pp[, haa_bat := fifelse(stand == "R", -haa, haa)]
pp[, pfx_x_bat := fifelse(stand == "R", -pfx_x, pfx_x)]
pp[, plate_x_bat := fifelse(stand == "R", -plate_x, plate_x)]
pp[, pz_rel := (plate_z - sz_bot) / pmax(sz_top - sz_bot, 0.1)]
pp[, pit_rv := -delta_run_exp]          # sign so that positive favours the pitcher

sw <- pp[description %in% SWING]
agg <- sw[, .(
  n_pp     = .N,
  p_throws = p_throws[1],
  whiff_pp = mean(description %in% WHIFF),
  n_bip    = sum(description == "hit_into_play"),
  xwobacon = mean(estimated_woba_using_speedangle[description == "hit_into_play"],
                  na.rm = TRUE),
  ev       = mean(launch_speed[description == "hit_into_play"], na.rm = TRUE),
  rv_swing = mean(pit_rv, na.rm = TRUE),
  velo     = mean(release_speed, na.rm = TRUE),
  spin     = mean(release_spin_rate, na.rm = TRUE),
  ext      = mean(release_extension, na.rm = TRUE),
  rel_x    = mean(release_pos_x, na.rm = TRUE),
  rel_z    = mean(release_pos_z, na.rm = TRUE),
  arm_ang  = mean(arm_angle, na.rm = TRUE),
  hb       = mean(pfx_x_bat, na.rm = TRUE),
  ivb      = mean(pfx_z, na.rm = TRUE),
  vaa      = mean(vaa, na.rm = TRUE),
  haa      = mean(haa_bat, na.rm = TRUE),
  loc_x    = mean(plate_x_bat, na.rm = TRUE),
  loc_z    = mean(pz_rel, na.rm = TRUE)
), by = .(id = pitcher, year = game_year, ptype, bat_side = stand)]

# pitch-level velo relative to the pitcher's own fastball, the usual "velo kill"
fb <- pp[pitch_type %in% c("FF", "SI"), .(fb_velo = mean(release_speed, na.rm = TRUE)),
         by = .(id = pitcher, year = game_year)]
agg <- merge(agg, fb, by = c("id", "year"), all.x = TRUE)
agg[, velo_kill := fb_velo - velo]

d <- merge(lb, agg, by = c("id", "year", "ptype", "bat_side"))
d[, platoon := fifelse(bat_side == p_throws, "same", "opp")]
# The leaderboard's throwing hand must agree with the per-pitch feed, or the join
# has drifted onto the wrong pitcher.
stopifnot(all(d$p_throws_lb == d$p_throws))

# QC gate. Savant's api_pitch_type and the raw feed's pitch_type disagree on a
# small tail of cells (a pitcher whose sweepers are logged as sliders, say). Those
# cells join a leaderboard row to the wrong pitch's outcomes, so they are dropped
# on count agreement rather than trusted. The band is calibrated on whiff-rate
# agreement, which is published on both sides and so acts as an independent check:
# inside [0.90, 1.20] agreement is cor 0.995 / mean abs diff 0.008, and it decays
# sharply outside. Per-pitch counts run about 3% above Savant's because Savant
# excludes bunts, so the band is deliberately asymmetric about 1.
d[, cnt_ratio := n_pp / n_swings]
n_pre <- nrow(d)
d <- d[cnt_ratio %between% c(0.90, 1.20)]
cat(sprintf("cells: %d leaderboard >=%d swings -> %d joined -> %d after QC (dropped %d on pitch-type disagreement)\n",
            nrow(lb), MIN_SWINGS, n_pre, nrow(d), n_pre - nrow(d)))
cat(sprintf("  whiff rate agreement: cor = %.4f, mean abs diff = %.4f\n",
            cor(d$whiff_rate, d$whiff_pp), mean(abs(d$whiff_rate - d$whiff_pp))))
cat(sprintf("  seasons %s | pitchers %d | swings %s\n\n",
            paste(range(d$year), collapse = "-"), uniqueN(d$id),
            format(sum(d$n_swings), big.mark = ",")))

# Within-cell-group deviations. Everything below is asked within pitch type,
# batter side and season, so that "changeups flail more than fastballs" does not
# masquerade as a finding about individual pitches.
dem <- function(v) v - mean(v, na.rm = TRUE)
# Approach angles are held out of the joint model on purpose. VAA is very nearly an
# exact function of release height, induced break, velocity and plate height, so
# fitting it alongside them is double-counting: with VAA and HAA in, VIF reaches
# 116 for VAA and 87 for release height and the coefficients become unstable
# linear combinations rather than effects (the first cut of this produced betas
# above 1.0 with cancelling signs). Dropping the two derived angles takes max VIF
# to 3.0. They are still reported below as marginal associations, where they are
# interpretable and are what a pitching coach actually talks about.
FEATS <- c("velo", "velo_kill", "spin", "ext", "rel_x", "rel_z", "arm_ang",
           "hb", "ivb", "loc_x", "loc_z")
DERIVED <- c("vaa", "haa")
for (v in c(RATES, FEATS, DERIVED, "xwobacon", "whiff_pp", "rv_swing", "ev"))
  d[, (paste0("d_", v)) := dem(get(v)), by = .(ptype, bat_side, year)]

# ------------------------------------------------- 1. reliability of the rates --
cat("=========== 1. Is each category a repeatable property of the pitch? ===========\n")
cat("Year-over-year correlation for the same pitcher/pitch/batter-side cell,\n")
cat("pitch-type, batter-side and season effects removed. whiff_rate is the yardstick.\n\n")
a <- copy(d); a[, cell := paste(id, ptype, bat_side, sep = "_")]
b <- copy(a)[, year := year - 1]
kp <- c("cell", "year", paste0("d_", c(RATES, "whiff_pp")))
pr <- merge(a[, ..kp], b[, ..kp], by = c("cell", "year"), suffixes = c("", ".n"))
rel <- rbindlist(lapply(paste0("d_", c(RATES, "whiff_pp")), function(v) {
  x <- pr[[v]]; y <- pr[[paste0(v, ".n")]]; ok <- is.finite(x) & is.finite(y)
  data.table(metric = sub("^d_", "", v), r_yoy = round(cor(x[ok], y[ok]), 3),
             n_pairs = sum(ok))
}))
print(rel[order(-r_yoy)])

# --------------------------------------------- 2. what each category is worth --
cat("\n=========== 2. What is each category worth? ===========\n")
cat(sprintf("Cell-level, within pitch type/side/season, weighted by sample (n_bip>=%d\n", MIN_BIP))
cat("for xwOBACON). Slope is per 10 percentage points of the category rate.\n\n")
sc <- d[n_bip >= MIN_BIP & is.finite(xwobacon)]
worth <- rbindlist(lapply(RATES, function(v) {
  f1 <- lm(d_xwobacon ~ 0 + get(paste0("d_", v)), data = sc, weights = sc$n_bip)
  f2 <- lm(d_whiff_pp ~ 0 + get(paste0("d_", v)), data = d, weights = d$n_swings)
  f3 <- lm(d_rv_swing ~ 0 + get(paste0("d_", v)), data = d, weights = d$n_swings)
  data.table(category = sub("_percent", "", v),
             xwobacon_per_10pp = round(coef(f1)[1] * 0.10, 4),
             t_xw = round(summary(f1)$coefficients[1, 3], 1),
             whiff_per_10pp = round(coef(f2)[1] * 0.10, 4),
             rv_swing_per_10pp = round(coef(f3)[1] * 0.10, 4))
}))
print(worth[order(xwobacon_per_10pp)])
cat("\n(negative xwobacon_per_10pp = the pitcher gets WEAKER contact as the rate rises)\n")

# ------------------------------------------------- 3. what induces each category --
cat("\n=========== 3. What pitch characteristics induce each category? ===========\n")
cat("Standardized betas, joint within-pitch-type regression weighted by swings.\n")
cat("Only |beta| >= 0.08 shown; blank means no material association. Max VIF 3.0.\n\n")
TARGETS <- c("tied_up_percent", "flailed_percent", "over_percent",
             "under_percent", "late_percent", "early_percent")
rhs <- paste(paste0("d_", FEATS), collapse = " + ")
keep <- paste0("d_", c(FEATS, DERIVED))
dj <- d[complete.cases(d[, ..keep])]
ind <- rbindlist(lapply(TARGETS, function(v) {
  fit <- lm(as.formula(sprintf("scale(d_%s) ~ 0 + %s", v, rhs)),
            data = dj, weights = dj$n_swings)
  cf <- coef(fit) * sapply(FEATS, function(f) sd(dj[[paste0("d_", f)]], na.rm = TRUE))
  out <- as.data.table(as.list(round(cf, 2)))
  setnames(out, FEATS); cbind(category = sub("_percent", "", v), out)
}))
show <- copy(ind)
for (f in FEATS) show[abs(get(f)) < 0.08, (f) := NA]
print(show)

cat("\nMarginal correlations for the derived approach angles, which cannot be\n")
cat("identified jointly alongside the inputs they are computed from:\n")
mg <- rbindlist(lapply(TARGETS, function(v) {
  r <- sapply(DERIVED, function(f)
    round(cor(dj[[paste0("d_", v)]], dj[[paste0("d_", f)]], use = "complete.obs"), 2))
  cbind(category = sub("_percent", "", v), as.data.table(as.list(r)))
}))
print(mg)

fwrite(d, file.path(OUT_DIR, "swing_timing_cells.csv"))
cat(sprintf("\nwrote %s (%d cells, %d cols)\n",
            file.path(OUT_DIR, "swing_timing_cells.csv"), nrow(d), ncol(d)))
