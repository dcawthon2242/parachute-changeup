#!/usr/bin/env Rscript

# Does the "matched look + kill -> over-the-top whiffs" rule hold for pitch-type
# pairs OTHER than FF->CH?  We anchor every secondary on the HARDER pitch it most
# resembles (min spin-axis diff among harder pitches), so a slider gets tested
# against a cutter, a sweeper against a slider, a split against a sinker, etc.
# Then we ask, within each anchor->secondary TYPE pairing: among matched-look
# instances, does a big velo/spin kill produce whiff overperformance and misses
# low in the zone (swinging over the top)?

suppressPackageStartupMessages({ library(data.table) })

in_csv  <- file.path("data", "statcast_2026", "statcast_2026_all.csv")
out_dir <- file.path("data", "statcast_2026")

MIN_ANCHOR <- 60
MIN_SECOND <- 40
AXIS_MATCH <- 30      # "same look" ceiling: spin-direction diff (deg)
EFF_MATCH  <- 0.20    # "same look" ceiling: spin-efficiency diff
VELO_GAP   <- 3       # anchor must be at least this much harder
KILL_VELO  <- 7       # mph; "big" velo kill

dt <- data.table::fread(in_csv, showProgress = FALSE, select = c(
  "pitch_type","game_type","pitcher","player_name","p_throws",
  "release_speed","release_spin_rate","spin_axis","pfx_x","pfx_z",
  "description","zone","plate_x","plate_z","sz_top","sz_bot",
  "vx0","vy0","vz0","ax","ay","az"))

dt <- dt[game_type == "R" & !is.na(pfx_x) & !is.na(pfx_z) & pitch_type != ""]

g <- 32.174; yf <- 17/12
dt[, t_plate := (-vy0 - sqrt(vy0^2 - 2*ay*(50 - yf))) / ay]
dt[, vaa := -atan((vz0 + az*t_plate) / (vy0 + ay*t_plate)) * 180/pi]
dt[, zone_rel := (plate_z - sz_bot) / (sz_top - sz_bot)]

# ---- spin efficiency (Nathan Magnus decomposition) -------------------------
# Magnus accel = non-gravity accel minus the drag (velocity-parallel) part.
# efficiency ~ |a_magnus| / (v * total_spin); normalized so the 99th pct = 1.
dt[, `:=`(vxm = vx0 + ax*t_plate/2, vym = vy0 + ay*t_plate/2, vzm = vz0 + az*t_plate/2)]
dt[, vmag := sqrt(vxm^2 + vym^2 + vzm^2)]
dt[, dotp := (ax*vxm + ay*vym + (az+g)*vzm) / vmag^2]
dt[, amag := sqrt((ax-dotp*vxm)^2 + (ay-dotp*vym)^2 + ((az+g)-dotp*vzm)^2)]
dt[, eff_raw := fifelse(!is.na(release_spin_rate) & release_spin_rate > 0,
                        amag/(vmag*release_spin_rate), NA_real_)]
kq <- quantile(dt$eff_raw, 0.99, na.rm = TRUE)
dt[, spin_eff := pmin(eff_raw/kq, 1.05)]

swing_desc <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip",
                "hit_into_play","foul_bunt","missed_bunt","bunt_foul_tip")
whiff_desc <- c("swinging_strike","swinging_strike_blocked","foul_tip","missed_bunt")
dt[, is_swing := description %in% swing_desc]
dt[, is_whiff := description %in% whiff_desc]

cmean <- function(a){ a<-a[!is.na(a)]; if(!length(a)) return(NA_real_)
  r<-a*pi/180; ((atan2(mean(sin(r)),mean(cos(r)))*180/pi)+360)%%360 }
circd <- function(a,b){ d<-abs(a-b)%%360; pmin(d,360-d) }

# ---- location+geometry whiff model (expected whiff, over-the-top signature) --
sw <- dt[is_swing == TRUE & is.finite(vaa) & is.finite(zone_rel) &
         !is.na(plate_x) & !is.na(plate_z) & !is.na(release_speed) &
         !is.na(pfx_z) & !is.na(pfx_x)]
sw[, pt := factor(pitch_type)]
whmod <- glm(is_whiff ~ pt + release_speed + pfx_z + pfx_x + plate_x +
               plate_z + I(plate_z^2) + vaa, data = sw, family = binomial())
sw[, pred_whiff := predict(whmod, type = "response")]
swagg <- sw[, .(swings2=.N, whiffs2=sum(is_whiff), exp_whiff=sum(pred_whiff),
                whiff_zone_rel=mean(zone_rel[is_whiff], na.rm=TRUE),
                underbarrel=mean(zone_rel[is_whiff] < 0.34, na.rm=TRUE)),
            by=.(pitcher, pitch_type)]
swagg[, whiff_over_exp := (whiffs2 - exp_whiff)/swings2*100]

# ---- per pitcher x pitch type ----------------------------------------------
agg <- dt[, .(
  n=.N, velo=mean(release_speed,na.rm=TRUE), spin=mean(release_spin_rate,na.rm=TRUE),
  axis=cmean(spin_axis), eff=mean(spin_eff,na.rm=TRUE),
  ivb=mean(pfx_z)*12, hb=mean(pfx_x)*12, vaa=mean(vaa,na.rm=TRUE),
  swings=sum(is_swing), whiffs=sum(is_whiff)
), by=.(pitcher, player_name, pitch_type)]
agg <- merge(agg, swagg, by=c("pitcher","pitch_type"), all.x=TRUE)
agg[, whiff_pct := whiffs/swings*100]
agg <- agg[swings >= 15 & is.finite(whiff_pct)]

# ---- clean per-pitch stuff residual (no pair terms) ------------------------
model_types <- c("FF","SI","FC","SL","ST","CU","KC","CH","FS","SV","CS")
ms <- agg[pitch_type %in% model_types & n >= MIN_SECOND]
ms[, pt := factor(pitch_type)]
stuff <- lm(whiff_pct ~ pt + velo + ivb + hb + spin, data=ms, weights=n)
ms[, whiff_resid := residuals(stuff)]
agg <- merge(agg, ms[, .(pitcher, pitch_type, whiff_resid)],
             by=c("pitcher","pitch_type"), all.x=TRUE)

# ---- all harder->softer pairs; anchor = most-similar-looking harder pitch ---
A <- agg[n >= MIN_ANCHOR, .(pitcher, a_type=pitch_type, a_velo=velo, a_spin=spin,
                            a_axis=axis, a_eff=eff, a_ivb=ivb, a_vaa=vaa)]
S <- agg[n >= MIN_SECOND & is.finite(whiff_resid),
         .(pitcher, player_name, s_type=pitch_type, s_n=n, s_velo=velo, s_spin=spin,
           s_axis=axis, s_eff=eff, s_ivb=ivb, s_vaa=vaa, whiff_pct, whiff_resid,
           whiff_over_exp, whiff_zone_rel, underbarrel)]

pp <- merge(S, A, by="pitcher", allow.cartesian=TRUE)
pp <- pp[a_type != s_type]
pp[, axis_diff := circd(a_axis, s_axis)]
pp[, eff_diff  := abs(a_eff - s_eff)]
pp[, velo_diff := a_velo - s_velo]
pp[, spin_diff := a_spin - s_spin]
pp[, vert_sep  := a_ivb - s_ivb]
pp[, vaa_diff  := s_vaa - a_vaa]
pp <- pp[velo_diff >= VELO_GAP]                     # anchor genuinely harder
# LOOK-ALIKE gate: same spin DIRECTION and same spin EFFICIENCY (so they truly
# move alike out of the hand -- rejects gyro-cutter vs efficient-sinker artifacts)
pp <- pp[axis_diff <= AXIS_MATCH & eff_diff <= EFF_MATCH]
# among qualifying harder look-alikes, keep the closest by combined axis+eff score
pp[, look_score := axis_diff/AXIS_MATCH + eff_diff/EFF_MATCH]
setorder(pp, pitcher, s_type, look_score)
pp <- pp[, .SD[1], by=.(pitcher, s_type)]
pp[, pair := paste0(a_type, "->", s_type)]
pp[, big_kill := velo_diff >= KILL_VELO]

# ---- (A) rule test by pitch-type pairing -----------------------------------
byp <- pp[, .(
  n=.N,
  velo_drop=round(mean(velo_diff),1), spin_kill=round(mean(spin_diff)),
  axis=round(mean(axis_diff),1), effD=round(mean(eff_diff),2),
  n_big=sum(big_kill),
  wResid_big=round(mean(whiff_resid[big_kill]),2),
  zone_big=round(mean(whiff_zone_rel[big_kill], na.rm=TRUE),2),
  wResid_small=round(mean(whiff_resid[!big_kill]),2)
), by=pair][n_big >= 5][order(-wResid_big)]

cat("=== (A) 'Matched look + velo kill' rule, by pitch-type pairing ===\n")
cat("    anchor = harder look-alike pitch | wResid = whiff over stuff-model\n")
cat("    zone_big<0 = big-kill whiffs land low = hitters swinging over the top\n")
cat("    (follows the rule when wResid_big > 0 and > wResid_small)\n\n")
print(byp)

# ---- (B) faceted mechanism: matched + big kill vs small kill, by secondary --
fac <- pp[, .(
  n_big=sum(big_kill), wResid_big=round(mean(whiff_resid[big_kill]),2),
  zone_big=round(mean(whiff_zone_rel[big_kill], na.rm=TRUE),2),
  n_small=sum(!big_kill), wResid_small=round(mean(whiff_resid[!big_kill]),2)
), by=.(secondary=s_type)][n_big >= 5][order(-wResid_big)]
cat("\n=== (B) Big-kill vs small-kill whiff overperformance, by SECONDARY type ===\n")
print(fac)

# ---- (C) overall: within matched-look pairs, does kill drive whiffs? --------
fit <- lm(whiff_resid ~ velo_diff + spin_diff + axis_diff, data=pp, weights=pp$s_n)
cat("\n=== (C) All matched-look pairs: whiff_resid ~ velo_diff + spin_diff + axis_diff ===\n")
print(round(summary(fit)$coefficients, 4))
cat(sprintf("\nBig-kill matched pairs: mean whiff_resid = %.2f (n=%d)\n",
            mean(pp[big_kill==TRUE]$whiff_resid), nrow(pp[big_kill==TRUE])))
cat(sprintf("Small-kill matched pairs: mean whiff_resid = %.2f (n=%d)\n",
            mean(pp[big_kill==FALSE]$whiff_resid), nrow(pp[big_kill==FALSE])))

# ---- top individual non-CH examples that follow the rule -------------------
cat("\n=== Top individual big-kill matched-look pairs (secondary != CH) ===\n")
top <- pp[big_kill==TRUE & s_type!="CH"][order(-whiff_resid)][1:20]
print(top[, .(player_name, pair, s_n, Dvelo=round(velo_diff,1), axis=round(axis_diff,1),
  effD=round(eff_diff,2), Dspin=round(spin_diff), whiff=round(whiff_pct,1),
  wResid=round(whiff_resid,1), zoneRel=round(whiff_zone_rel,2))])

fwrite(pp, file.path(out_dir, "pitch_pair_all_types_2026.csv"))
cat(sprintf("\nWrote %s (%d matched-look pairs)\n",
            file.path(out_dir, "pitch_pair_all_types_2026.csv"), nrow(pp)))
