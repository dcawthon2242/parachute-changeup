#!/usr/bin/env Rscript

# 1) Are the missing fouls 2-strike fouls?
# 2) Does a high velo-gap BREAKING ball help xRV (not just offspeed)?
# 3) Do hitters chase these pitches less?

suppressPackageStartupMessages({ library(data.table); library(splines) })
set.seed(7)
options(width = 215)

P <- fread("data/statcast_2026/velo_gap_tunnel_pairs.csv")
P[, vg3 := cut(velo_gap, quantile(velo_gap, 0:3/3), include.lowest=TRUE,
               labels=c("small gap","mid gap","big gap"))]
BR <- c("SL","ST","CU","KC","SV","CS")
OS <- c("CH","FS","FO")
P[, family := fifelse(pitch_type %in% BR, "breaking",
               fifelse(pitch_type %in% OS, "offspeed", "other"))]

cols <- c("game_year","game_type","pitcher","pitch_type","description","balls","strikes",
          "plate_x","plate_z","sz_top","sz_bot","estimated_woba_using_speedangle","delta_run_exp")
dt <- rbindlist(lapply(c(2025, 2026), function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress=FALSE, select=cols)))
setnames(dt, "estimated_woba_using_speedangle", "xw")
dt <- dt[game_type=="R" & pitch_type!="" & balls<=3 & strikes<=2 & is.finite(delta_run_exp) &
         is.finite(plate_x) & is.finite(plate_z)]
dt[, `:=`(bat_rv=delta_run_exp, in_zone=abs(plate_x)<=0.83 & plate_z>=sz_bot & plate_z<=sz_top)]
WH <- c("swinging_strike","swinging_strike_blocked","foul_tip")
SW <- c(WH,"foul","hit_into_play")
dt[, swung := description %in% SW]
dt[, chase := swung & !in_zone]
dt[, twoK  := strikes==2]
mb <- lm(bat_rv ~ ns(xw,5)+factor(paste0(balls,"-",strikes)),
         data=dt[description=="hit_into_play" & is.finite(xw)])
dt[, xrv_p := bat_rv]
dt[description=="hit_into_play" & is.finite(xw),
   xrv_p := predict(mb, newdata=dt[description=="hit_into_play" & is.finite(xw)])]
dt <- dt[!(description=="hit_into_play" & !is.finite(xw))]

sec <- merge(dt, P[, .(pitcher, pitch_type, vg3, family, velo_gap)], by=c("pitcher","pitch_type"))
cat(sprintf("Secondary pitches: %s\n\n", format(nrow(sec), big.mark=",")))

# =============================================================================
# 1. Fouls split by 2-strike vs not
# =============================================================================
cat("############ 1. Foul rates and foul xRV contribution, by count ############\n\n")
cat("  Contribution = pitcher runs per 100 pitches of the secondary from that event.\n")
cat("  A 2-strike foul is near-zero RV (at-bat continues). A 0-/1-strike foul is a cheap strike.\n\n")

F <- sec[, .(
  n=.N,
  foul0  = 100*mean(description=="foul" & twoK==FALSE),
  foul2  = 100*mean(description=="foul" & twoK==TRUE),
  c0     = -100*mean(fifelse(description=="foul" & twoK==FALSE, xrv_p, 0)),
  c2     = -100*mean(fifelse(description=="foul" & twoK==TRUE,  xrv_p, 0)),
  # what happens instead with 2 strikes
  whiff2 = 100*mean(description %in% WH & twoK==TRUE),
  bip2   = 100*mean(description=="hit_into_play" & twoK==TRUE),
  ball2  = 100*mean(description=="ball" & twoK==TRUE),
  w2c    = -100*mean(fifelse(description %in% WH & twoK==TRUE, xrv_p, 0)),
  b2c    = -100*mean(fifelse(description=="hit_into_play" & twoK==TRUE, xrv_p, 0))
), by=.(pitcher, pitch_type, vg3, family)]

for (v in c("foul0","foul2","c0","c2","whiff2","bip2","ball2","w2c","b2c"))
  F[, paste0(v,"_a") := residuals(lm(get(v) ~ factor(pitch_type)+n))]

cat("  Type-adjusted rates (per 100 pitches of the secondary):\n")
print(F[, .(pairs=.N,
            foul_0_1K=round(mean(foul0_a),2), foul_2K=round(mean(foul2_a),2),
            contrib_0_1K=round(mean(c0_a),3), contrib_2K=round(mean(c2_a),3),
            whiff_2K=round(mean(whiff2_a),2), bip_2K=round(mean(bip2_a),2)),
        by=vg3][order(vg3)], row.names=FALSE)

cat("\n  Big vs small, type-adjusted:\n")
tt <- function(v, lab) {
  t <- t.test(F[vg3=="big gap"][[v]], F[vg3=="small gap"][[v]])
  cat(sprintf("     %-16s  %+0.3f vs %+0.3f   diff %+0.3f   p = %.3f\n",
              lab, mean(F[vg3=="big gap"][[v]]), mean(F[vg3=="small gap"][[v]]),
              mean(F[vg3=="big gap"][[v]])-mean(F[vg3=="small gap"][[v]]), t$p.value))
}
tt("foul0_a", "foul rate, 0-1K")
tt("foul2_a", "foul rate, 2K")
tt("c0_a",    "foul xRV, 0-1K")
tt("c2_a",    "foul xRV, 2K")
tt("whiff2_a","whiff rate, 2K")
tt("bip2_a",  "BIP rate, 2K")
tt("w2c_a",   "whiff xRV, 2K")
tt("b2c_a",   "BIP xRV, 2K")

cat("\n  Raw 2-strike mix (what the 2K swing becomes):\n")
print(sec[twoK==TRUE & swung==TRUE, .(
  swings=.N,
  foul=round(100*mean(description=="foul"),1),
  whiff=round(100*mean(description %in% WH),1),
  bip=round(100*mean(description=="hit_into_play"),1)
), by=vg3][order(vg3)], row.names=FALSE)

# =============================================================================
# 2. High velo-gap breaking balls vs offspeed
# =============================================================================
cat("\n############ 2. xRV of a high-gap breaking ball vs high-gap offspeed ############\n\n")
PR <- sec[, .(
  n=.N, xrv=-100*mean(xrv_p),
  whiff=100*mean(description %in% WH),
  chase=100*mean(chase),
  zone=100*mean(in_zone),
  swing=100*mean(swung),
  cs=100*mean(description=="called_strike")
), by=.(pitcher, pitch_type, vg3, family, velo_gap)]

# residualize within family so we compare big vs small sliders, not sliders vs changeups
PR[, xrv_a    := residuals(lm(xrv ~ factor(pitch_type)+n)),    by=family]
PR[, whiff_a  := residuals(lm(whiff ~ factor(pitch_type)+n)),  by=family]
PR[, chase_a  := residuals(lm(chase ~ factor(pitch_type)+n)),  by=family]
PR[, zone_a   := residuals(lm(zone ~ factor(pitch_type)+n)),   by=family]

cat("  Type-adjusted within family (big-gap slider vs small-gap slider, etc.):\n")
print(PR[family %in% c("breaking","offspeed"),
         .(pairs=.N, velo_gap=round(mean(velo_gap),1),
           xrv=round(mean(xrv),3), xrv_a=round(mean(xrv_a),3),
           whiff_a=round(mean(whiff_a),2), chase_a=round(mean(chase_a),2),
           zone_a=round(mean(zone_a),2)),
         by=.(family, vg3)][order(family, vg3)], row.names=FALSE)

cat("\n  Continuous: xRV residual on velocity gap, by family.\n")
cat("  Coefficient is per 1 mph of gap, type held fixed.\n")
for (fam in c("breaking","offspeed")) {
  d <- PR[family==fam]
  m <- lm(xrv ~ velo_gap + factor(pitch_type) + n, data=d)
  s <- summary(m)$coefficients["velo_gap",]
  cat(sprintf("     %-10s  β = %+.4f / mph   t = %+.2f   p = %.3f   R2 = %.3f   n = %d\n",
              fam, s[1], s[3], s[4], summary(m)$r.squared, nrow(d)))
}

cat("\n  Same, on raw xRV (no type FE) — does a bigger-gap curve just beat a smaller-gap slider?\n")
for (fam in c("breaking","offspeed")) {
  d <- PR[family==fam]
  m <- lm(xrv ~ velo_gap + n, data=d)
  s <- summary(m)$coefficients["velo_gap",]
  cat(sprintf("     %-10s  β = %+.4f / mph   p = %.3f   n = %d\n", fam, s[1], s[4], nrow(d)))
}

cat("\n  By breaking-ball type, big vs small (cells with 15+ pairs):\n")
BT <- PR[family=="breaking", .(
  pairs=.N, gap=round(mean(velo_gap),1),
  xrv=round(mean(xrv),3), xrv_a=round(mean(xrv_a),3),
  whiff_a=round(mean(whiff_a),2), chase_a=round(mean(chase_a),2)
), by=.(pitch_type, vg3)]
print(BT[order(pitch_type, vg3)], row.names=FALSE)

# terciles within breaking only, so "big" is big among breakers
cat("\n  Breaking balls only, re-terciled on their own gap distribution:\n")
BRP <- PR[family=="breaking"]
BRP[, vgB := cut(velo_gap, quantile(velo_gap, 0:3/3), include.lowest=TRUE,
                 labels=c("small (among BR)","mid","big (among BR)"))]
print(BRP[, .(pairs=.N, gap=round(mean(velo_gap),1),
              xrv=round(mean(xrv),3),
              xrv_a=round(mean(residuals(lm(xrv ~ factor(pitch_type)+n))),3),
              chase=round(mean(chase),1),
              whiff=round(mean(whiff),1)), by=vgB][order(vgB)], row.names=FALSE)

# =============================================================================
# 3. Chase
# =============================================================================
cat("\n############ 3. Do hitters chase these less? ############\n\n")
cat("  Chase = swing at a pitch outside the zone, per 100 pitches of the secondary.\n")
cat("  Also: chase given the pitch is already outside the zone (the real chase rate).\n\n")

CH <- sec[, .(
  n=.N,
  chase=100*mean(chase),
  chase_ooz=100*mean(swung[!in_zone]),
  ooz=100*mean(!in_zone),
  swing_z=100*mean(swung[in_zone])
), by=.(pitcher, pitch_type, vg3, family)]
for (v in c("chase","chase_ooz","ooz","swing_z"))
  CH[, paste0(v,"_a") := residuals(lm(get(v) ~ factor(pitch_type)+n))]

print(CH[, .(pairs=.N,
             chase=round(mean(chase),2), chase_a=round(mean(chase_a),2),
             chase_ooz=round(mean(chase_ooz),2), chase_ooz_a=round(mean(chase_ooz_a),2),
             swing_in_zone_a=round(mean(swing_z_a),2)),
         by=vg3][order(vg3)], row.names=FALSE)
cat("\n  Big vs small:\n")
for (v in c("chase_a","chase_ooz_a","swing_z_a")) {
  t <- t.test(CH[vg3=="big gap"][[v]], CH[vg3=="small gap"][[v]])
  cat(sprintf("     %-16s  %+0.3f vs %+0.3f   diff %+0.3f   p = %.3f\n",
              v, mean(CH[vg3=="big gap"][[v]]), mean(CH[vg3=="small gap"][[v]]),
              mean(CH[vg3=="big gap"][[v]])-mean(CH[vg3=="small gap"][[v]]), t$p.value))
}
cat("\n  Chase-outside-zone residual, by family:\n")
print(CH[family %in% c("breaking","offspeed"),
         .(pairs=.N, chase_ooz=round(mean(chase_ooz),1),
           chase_ooz_a=round(mean(chase_ooz_a),2),
           swing_z_a=round(mean(swing_z_a),2)),
         by=.(family, vg3)][order(family, vg3)], row.names=FALSE)
