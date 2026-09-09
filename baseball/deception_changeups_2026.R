#!/usr/bin/env Rscript

# Identify "same-spin parachute" changeups (spin axis mirrors the four-seam,
# but the pitch is much slower and drops) and quantify the traits that make
# them deceptive in ways spin/velocity-based "stuff" models tend to miss.
#
# Two changeup archetypes:
#   * Turnover  (e.g. Boyd, Alexander): pronated, spin axis differs from the FB.
#   * Parachute (e.g. Hellickson, Ragans): spins on ~the same axis as the FB,
#     so it looks identical out of hand, then parachutes under the barrel.
# The discriminator is the CH-vs-FF spin-axis difference: small = parachute.
#
# "Proof of deception beyond stuff" is built from three families of stats:
#   1. Looks-like-the-FB (perception inputs stuff models underweight):
#        spin-axis diff, spin-rate diff, release-point diff, extension diff.
#   2. Whiff/chase outcomes driven by misreads: Whiff%, Chase%, CSW%.
#   3. Contact suppression: xwOBAcon, GB%, avg exit velo.
# The final section fits a simple "stuff" model (results ~ physical traits) and
# reports each pitch's residual: outperforming its own stuff = deception value.

suppressPackageStartupMessages({
  library(data.table)
})

in_csv  <- file.path("data", "statcast_2026", "statcast_2026_all.csv")
out_csv <- file.path("data", "statcast_2026", "deception_changeups_2026.csv")

REGULAR_ONLY  <- TRUE
MIN_CH        <- 25     # min changeups thrown
MIN_FF        <- 25     # min four-seamers thrown
VELO_DIFF_MIN <- 10     # (ff_velo - ch_velo) must exceed this
AXIS_DIFF_MAX <- 20     # CH-vs-FF spin-axis diff (deg) below this = "same spin"

# ---- Load -----------------------------------------------------------------

dt <- data.table::fread(in_csv, showProgress = FALSE, select = c(
  "pitch_type", "game_type", "pitcher", "player_name", "p_throws",
  "release_speed", "release_spin_rate", "spin_axis", "pfx_x", "pfx_z",
  "release_pos_x", "release_pos_z", "release_extension", "arm_angle",
  "description", "zone", "bb_type", "launch_speed",
  "estimated_woba_using_speedangle", "delta_run_exp"))

if (REGULAR_ONLY) dt <- dt[game_type == "R"]
dt <- dt[pitch_type %in% c("FF", "CH") & !is.na(pfx_x) & !is.na(pfx_z)]

# ---- Circular helpers for spin axis (degrees) -----------------------------

cmean <- function(a) {
  a <- a[!is.na(a)]
  r <- a * pi / 180
  ang <- atan2(mean(sin(r)), mean(cos(r))) * 180 / pi
  (ang + 360) %% 360
}
circdiff <- function(a, b) {
  d <- abs(a - b) %% 360
  pmin(d, 360 - d)
}

# ---- Outcome flags (changeups) --------------------------------------------

swing_desc <- c("swinging_strike", "swinging_strike_blocked", "foul", "foul_tip",
                "hit_into_play", "foul_bunt", "missed_bunt", "bunt_foul_tip")
whiff_desc <- c("swinging_strike", "swinging_strike_blocked", "foul_tip", "missed_bunt")
batted     <- c("fly_ball", "ground_ball", "line_drive", "popup")

dt[, is_swing := description %in% swing_desc]
dt[, is_whiff := description %in% whiff_desc]
dt[, is_called := description == "called_strike"]
dt[, out_zone := !is.na(zone) & zone >= 11]
dt[, in_zone  := !is.na(zone) & zone <= 9]
dt[, is_batted := bb_type %in% batted]

# ---- Per-pitcher summaries -------------------------------------------------

ch <- dt[pitch_type == "CH", .(
  ch_n        = .N,
  ch_velo     = mean(release_speed, na.rm = TRUE),
  ch_spin     = mean(release_spin_rate, na.rm = TRUE),
  ch_axis     = cmean(spin_axis),
  ch_ivb      = mean(pfx_z) * 12,
  ch_hb       = mean(pfx_x) * 12,
  ch_relx     = mean(release_pos_x, na.rm = TRUE),
  ch_relz     = mean(release_pos_z, na.rm = TRUE),
  ch_ext      = mean(release_extension, na.rm = TRUE),
  arm_angle   = mean(arm_angle, na.rm = TRUE),
  # outcomes
  swings      = sum(is_swing),
  whiffs      = sum(is_whiff),
  called      = sum(is_called),
  oz_pitches  = sum(out_zone),
  oz_swings   = sum(out_zone & is_swing),
  batted_n    = sum(is_batted),
  gb          = sum(bb_type == "ground_ball", na.rm = TRUE),
  xwobacon    = mean(estimated_woba_using_speedangle[!is.na(estimated_woba_using_speedangle)]),
  avg_ev      = mean(launch_speed[is_batted & !is.na(launch_speed)]),
  rv100       = -sum(delta_run_exp, na.rm = TRUE) / .N * 100
), by = .(pitcher, player_name, p_throws)]

ff <- dt[pitch_type == "FF", .(
  ff_n    = .N,
  ff_velo = mean(release_speed, na.rm = TRUE),
  ff_spin = mean(release_spin_rate, na.rm = TRUE),
  ff_axis = cmean(spin_axis),
  ff_ivb  = mean(pfx_z) * 12,
  ff_relx = mean(release_pos_x, na.rm = TRUE),
  ff_relz = mean(release_pos_z, na.rm = TRUE),
  ff_ext  = mean(release_extension, na.rm = TRUE)
), by = pitcher]

m <- merge(ch, ff, by = "pitcher")
m <- m[ch_n >= MIN_CH & ff_n >= MIN_FF]

# ---- Derived deception metrics --------------------------------------------

m[, velo_diff := ff_velo - ch_velo]
m[, axis_diff := circdiff(ch_axis, ff_axis)]                 # small = same spin as FB
m[, spin_diff := ff_spin - ch_spin]
m[, ivb_drop  := ff_ivb - ch_ivb]                            # + = CH drops vs FB
m[, rel_diff  := sqrt((ch_relx - ff_relx)^2 + (ch_relz - ff_relz)^2) * 12]  # inches
m[, ext_diff  := ff_ext - ch_ext]                            # ft; ~0 = same extension

m[, whiff_pct := whiffs / swings * 100]
m[, chase_pct := oz_swings / oz_pitches * 100]
m[, csw_pct   := (called + whiffs) / ch_n * 100]
m[, gb_pct    := gb / batted_n * 100]

# ---- Simple "stuff" model: do these beat what their raw stuff predicts? ----
# Fit results as a function of physical traits only (velocity, movement, spin,
# release), weighted by CH count, across ALL qualified changeups. The residual
# (actual - predicted) is the value NOT explained by stuff = deception proxy.

stuff_feat <- ~ ch_velo + velo_diff + ch_ivb + ch_hb + ch_spin + rel_diff + ext_diff
fit_resid <- function(target) {
  d <- as.data.frame(m[is.finite(get(target))])
  f <- update(stuff_feat, paste(target, "~ ."))
  mod <- lm(f, data = d, weights = ch_n)
  # Predict via coefficients to avoid predict() re-evaluating the weights call.
  mm <- model.matrix(delete.response(terms(mod)), data = as.data.frame(m))
  pred <- as.numeric(mm %*% coef(mod))
  m[[target]] - pred
}
m[, whiff_resid := fit_resid("whiff_pct")]
m[, xwoba_resid := fit_resid("xwobacon")]

# ---- Classify & report -----------------------------------------------------

m[, archetype := ifelse(axis_diff <= AXIS_DIFF_MAX, "parachute (same-spin)",
                        "turnover (diff-spin)")]

parachute <- m[velo_diff > VELO_DIFF_MIN & axis_diff <= AXIS_DIFF_MAX]
setorder(parachute, -whiff_pct)

# League-average changeup benchmark (all qualified CH).
league <- m[, .(
  ch_n = sum(ch_n),
  velo_diff = round(weighted.mean(velo_diff, ch_n), 1),
  axis_diff = round(weighted.mean(axis_diff, ch_n), 1),
  whiff_pct = round(weighted.mean(whiff_pct, ch_n, na.rm = TRUE), 1),
  chase_pct = round(weighted.mean(chase_pct, ch_n, na.rm = TRUE), 1),
  csw_pct   = round(weighted.mean(csw_pct, ch_n, na.rm = TRUE), 1),
  xwobacon  = round(weighted.mean(xwobacon, ch_n, na.rm = TRUE), 3),
  gb_pct    = round(weighted.mean(gb_pct, ch_n, na.rm = TRUE), 1),
  rv100     = round(weighted.mean(rv100, ch_n, na.rm = TRUE), 2)
)]

round_cols <- c("ch_velo","velo_diff","axis_diff","spin_diff","ivb_drop","rel_diff",
                "ext_diff","arm_angle","whiff_pct","chase_pct","csw_pct","gb_pct",
                "avg_ev","rv100","whiff_resid","xwoba_resid","ch_ivb","ch_hb")
parachute[, (round_cols) := lapply(.SD, round, 2), .SDcols = round_cols]
parachute[, ch_spin := round(ch_spin, 0)]
parachute[, xwobacon := round(xwobacon, 3)]

cat("=== League-average changeup (all", nrow(m), "qualified pitchers) ===\n")
print(league)

cat(sprintf("\n=== Same-spin PARACHUTE changeups (velo diff > %d, spin-axis diff <= %d deg) ===\n",
            VELO_DIFF_MIN, AXIS_DIFF_MAX))
print(parachute[, .(player_name, p_throws, ch_n,
                    velo_diff, axis_diff, ivb_drop, rel_diff, arm_angle,
                    whiff_pct, chase_pct, csw_pct, xwobacon, gb_pct, rv100,
                    whiff_resid, xwoba_resid)])

data.table::fwrite(parachute, out_csv)
cat(sprintf("\nWrote %s\n", out_csv))
