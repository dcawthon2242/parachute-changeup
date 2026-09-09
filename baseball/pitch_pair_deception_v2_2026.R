#!/usr/bin/env Rscript

# Pitch-pair deception, v2. Three upgrades over v1:
#   (1) ANCHOR on the four-seam (fallback SI, then FC) instead of most-thrown
#       fastball -- so cutter-primary arsenals (e.g. Bibee) compare the changeup
#       to the riding fastball it actually mirrors.
#   (2) MECHANISM test: split "matched-spin" secondaries by spin-rate kill and
#       velo drop to see whether same-axis + big-kill is the real whiff engine.
#   (3) MISS-OVER-THE-TOP model: compute vertical approach angle (VAA) and a
#       location+geometry-aware whiff model; "whiff over expected" is whiffs
#       beyond what location/approach/stuff predict, and we tag whether those
#       whiffs happen low in the zone (hitter swinging over the top).

suppressPackageStartupMessages({ library(data.table) })

in_csv  <- file.path("data", "statcast_2026", "statcast_2026_all.csv")
out_dir <- file.path("data", "statcast_2026")

REGULAR_ONLY <- TRUE
MIN_PRIMARY  <- 100
MIN_SECOND   <- 40
AXIS_MATCH   <- 25
AXIS_MIRROR  <- 155
KILL_BIG     <- 500   # rpm; "big" spin-rate kill vs the fastball

dt <- data.table::fread(in_csv, showProgress = FALSE, select = c(
  "pitch_type","game_type","pitcher","player_name","p_throws",
  "release_speed","release_spin_rate","spin_axis","pfx_x","pfx_z",
  "release_pos_x","release_pos_z","release_extension",
  "description","zone","plate_x","plate_z","sz_top","sz_bot",
  "vy0","vz0","ay","az",
  "estimated_woba_using_speedangle","delta_run_exp",
  "attack_angle","swing_path_tilt"))

if (REGULAR_ONLY) dt <- dt[game_type == "R"]
dt <- dt[!is.na(pfx_x) & !is.na(pfx_z) & pitch_type != ""]

# ---- Pitch-level derived fields -------------------------------------------

yf <- 17/12
dt[, t_plate := (-vy0 - sqrt(vy0^2 - 2*ay*(50 - yf))) / ay]
dt[, vaa := -atan((vz0 + az*t_plate) / (vy0 + ay*t_plate)) * 180/pi]
dt[, zone_rel := (plate_z - sz_bot) / (sz_top - sz_bot)]   # 0=bottom,1=top,<0 below

swing_desc <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip",
                "hit_into_play","foul_bunt","missed_bunt","bunt_foul_tip")
whiff_desc <- c("swinging_strike","swinging_strike_blocked","foul_tip","missed_bunt")
dt[, is_swing  := description %in% swing_desc]
dt[, is_whiff  := description %in% whiff_desc]
dt[, is_called := description == "called_strike"]
dt[, out_zone  := !is.na(zone) & zone >= 11]

cmean <- function(a){ a<-a[!is.na(a)]; if(!length(a)) return(NA_real_)
  r<-a*pi/180; ((atan2(mean(sin(r)),mean(cos(r)))*180/pi)+360)%%360 }
circd <- function(a,b){ d<-abs(a-b)%%360; pmin(d,360-d) }

# ---- (3) Location+geometry whiff model (pitch level, on swings) ------------
# Expected whiff from WHERE and HOW the pitch arrives (location, approach angle,
# velo, movement, pitch type) -- NOT the spin-axis relationship. Whiffs beyond
# this = deception. We also record whether whiffs are low (over-the-top misses).

sw <- dt[is_swing == TRUE & is.finite(vaa) & is.finite(zone_rel) &
         !is.na(plate_x) & !is.na(plate_z) & !is.na(release_speed) &
         !is.na(pfx_z) & !is.na(pfx_x) & pitch_type != ""]
sw[, pt := factor(pitch_type)]
whmod <- glm(is_whiff ~ pt + release_speed + pfx_z + pfx_x + plate_x +
               plate_z + I(plate_z^2) + vaa,
             data = sw, family = binomial(), model = FALSE)
sw[, pred_whiff := predict(whmod, type = "response")]

swagg <- sw[, .(
  swings2   = .N,
  whiffs2   = sum(is_whiff),
  exp_whiff = sum(pred_whiff),
  whiff_zone_rel = mean(zone_rel[is_whiff], na.rm = TRUE),      # low = over the top
  underbarrel    = mean(zone_rel[is_whiff] < 0.34, na.rm = TRUE) # frac of whiffs low
), by = .(pitcher, pitch_type)]
swagg[, whiff_over_exp := (whiffs2 - exp_whiff) / swings2 * 100]

# ---- Per pitcher x pitch type ---------------------------------------------

agg <- dt[, .(
  n = .N,
  velo = mean(release_speed, na.rm=TRUE), spin = mean(release_spin_rate, na.rm=TRUE),
  axis = cmean(spin_axis), ivb = mean(pfx_z)*12, hb = mean(pfx_x)*12,
  relx = mean(release_pos_x, na.rm=TRUE), relz = mean(release_pos_z, na.rm=TRUE),
  ext = mean(release_extension, na.rm=TRUE), vaa = mean(vaa, na.rm=TRUE),
  swings = sum(is_swing), whiffs = sum(is_whiff), called = sum(is_called),
  oz = sum(out_zone), ozsw = sum(out_zone & is_swing),
  xwobacon = mean(estimated_woba_using_speedangle[!is.na(estimated_woba_using_speedangle)]),
  wh_attack = mean(attack_angle[is_whiff], na.rm=TRUE),
  wh_tilt   = mean(swing_path_tilt[is_whiff], na.rm=TRUE)
), by = .(pitcher, player_name, p_throws, pitch_type)]
agg <- merge(agg, swagg, by = c("pitcher","pitch_type"), all.x = TRUE)

# ---- (1) FOUR-SEAM anchor (fallback SI, then FC) --------------------------

pick_anchor <- function(sub) {
  for (ft in c("FF","SI","FC")) {
    r <- sub[pitch_type == ft & n >= MIN_PRIMARY]
    if (nrow(r)) return(r[which.max(n)])
  }
  NULL
}
anchor <- agg[, { a <- pick_anchor(.SD); if (is.null(a)) NULL else a },
              by = pitcher, .SDcols = names(agg)]
anchor <- anchor[, .(pitcher, fb_type = pitch_type, fb_velo = velo, fb_spin = spin,
                     fb_axis = axis, fb_ivb = ivb, fb_hb = hb, fb_relx = relx,
                     fb_relz = relz, fb_ext = ext, fb_vaa = vaa)]

# ---- Secondary pitches vs their four-seam ----------------------------------

pairs <- merge(agg[n >= MIN_SECOND], anchor, by = "pitcher")
pairs <- pairs[pitch_type != fb_type]

pairs[, axis_diff := circd(axis, fb_axis)]
pairs[, spin_diff := fb_spin - spin]      # + = secondary spins slower (kill)
pairs[, velo_diff := fb_velo - velo]
pairs[, vert_sep  := fb_ivb - ivb]
pairs[, horz_sep  := hb - fb_hb]
pairs[, vaa_diff  := vaa - fb_vaa]        # - = steeper descent than FB
pairs[, rel_diff  := sqrt((relx-fb_relx)^2 + (relz-fb_relz)^2)*12]
pairs[, ext_diff  := fb_ext - ext]
pairs[, whiff_pct := whiffs / swings * 100]
pairs[, chase_pct := ozsw / oz * 100]
pairs[, csw_pct   := (called + whiffs) / n * 100]

# ---- Stuff-model whiff residual (as v1) -----------------------------------

model_types <- c("FF","SI","FC","SL","ST","CU","KC","CH","FS","SV","CS")
mp <- pairs[pitch_type %in% model_types & is.finite(whiff_pct) & swings >= 15 &
            is.finite(rel_diff) & is.finite(ext_diff)]
mp[, pt := factor(pitch_type)]
stuff <- lm(whiff_pct ~ pt + velo + ivb + hb + spin + velo_diff +
              I(abs(vert_sep)) + I(abs(horz_sep)) + rel_diff + ext_diff,
            data = mp, weights = n)
mp[, whiff_resid := residuals(stuff)]
pairs <- merge(pairs, mp[, .(pitcher, pitch_type, whiff_resid)],
               by = c("pitcher","pitch_type"), all.x = TRUE)

# ---- (2) MECHANISM: matched-axis x spin-kill ------------------------------

valid <- pairs[is.finite(whiff_resid)]
valid[, grp := fifelse(axis_diff <= AXIS_MATCH & spin_diff >= KILL_BIG, "matched + big kill",
               fifelse(axis_diff <= AXIS_MATCH & spin_diff <  KILL_BIG, "matched + small kill",
               fifelse(axis_diff >= AXIS_MIRROR, "mirror (opposite)", "middle")))]

mech <- valid[, .(
  n = .N,
  mean_axis_diff = round(weighted.mean(axis_diff, n),1),
  mean_spin_kill = round(weighted.mean(spin_diff, n),0),
  mean_velo_drop = round(weighted.mean(velo_diff, n),1),
  whiff_pct = round(weighted.mean(whiff_pct, n, na.rm=TRUE),1),
  whiff_vs_stuff = round(weighted.mean(whiff_resid, n),2),
  whiff_over_exp = round(weighted.mean(whiff_over_exp, n, na.rm=TRUE),2),
  whiff_zone_rel = round(weighted.mean(whiff_zone_rel, n, na.rm=TRUE),2)
), by = grp][order(-whiff_vs_stuff)]

cat("=== (2) MECHANISM: whiff overperformance by spin relationship x spin-kill ===\n")
cat("    (whiff_zone_rel: lower = whiffs happen lower in zone = swinging over the top)\n")
print(mech)

# regression: does spin-kill add whiffs among matched-axis pitches?
mm <- valid[axis_diff <= AXIS_MATCH & is.finite(whiff_resid)]
fit <- lm(whiff_resid ~ spin_diff + velo_diff + axis_diff, data = mm, weights = mm$n)
cat("\nAmong matched-axis (<=25 deg) secondaries, whiff_resid ~ spin_diff + velo_diff + axis_diff:\n")
print(round(summary(fit)$coefficients, 4))

# ---- (1) Leaderboards with FF anchor --------------------------------------

lb <- function(d) d[, .(player_name, p_throws, pitch=pitch_type, sec_n=n, fb=fb_type,
  Dvelo=round(velo_diff,1), axis=round(axis_diff,1), Dspin=round(spin_diff),
  vertSep=round(vert_sep,1), whiff=round(whiff_pct,1),
  wOverExp=round(whiff_over_exp,1), wVsStuff=round(whiff_resid,1))]

cat("\n=== (1) Top MATCHED-SPIN look-alikes (FF-anchored, axis<=25), by whiff-over-expected ===\n")
print(lb(valid[axis_diff <= AXIS_MATCH & velo_diff > 6 & is.finite(whiff_over_exp)][order(-whiff_over_exp)][1:20]))

cat("\n=== (1) Top MIRROR pairs (axis>=155), by whiff-over-expected ===\n")
print(lb(valid[axis_diff >= AXIS_MIRROR & is.finite(whiff_over_exp)][order(-whiff_over_exp)][1:15]))

# ---- (3) Over-the-top signature: examples ---------------------------------

examples <- c("Cease, Dylan","Ribalta, Orlando","Dion, Will","Cantillo, Joey",
              "Ragans, Cole","Vesia, Alex","Boyd, Matthew","Bibee, Tanner",
              "Sanchez, Cristopher","Sánchez, Cristopher","Skubal, Tarik","Brown, Hunter")
ex <- pairs[player_name %in% examples & pitch_type == "CH"]
setorder(ex, -whiff_over_exp)
cat("\n=== (3) Changeup examples (FF-anchored; Luzardo removed, Vesia added) ===\n")
cat("    whiff_zone_rel low & underbarrel high = misses are over-the-top\n")
print(ex[, .(player_name, ch_n=n, axis=round(axis_diff,1), Dspin=round(spin_diff),
  Dvelo=round(velo_diff,1), vaa_diff=round(vaa_diff,1), whiff=round(whiff_pct,1),
  wOverExp=round(whiff_over_exp,1), zoneRel=round(whiff_zone_rel,2),
  underBrl=round(underbarrel,2))])

# league CH baseline for whiff_zone_rel
chall <- pairs[pitch_type=="CH" & is.finite(whiff_zone_rel)]
cat(sprintf("\nLeague CH whiff_zone_rel (baseline): %.2f | underbarrel: %.2f\n",
            weighted.mean(chall$whiff_zone_rel, chall$n), weighted.mean(chall$underbarrel, chall$n)))

fwrite(pairs, file.path(out_dir, "pitch_pair_deception_v2_2026.csv"))
cat(sprintf("\nWrote %s\n", file.path(out_dir, "pitch_pair_deception_v2_2026.csv")))
