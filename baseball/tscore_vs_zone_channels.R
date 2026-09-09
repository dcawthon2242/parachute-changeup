#!/usr/bin/env Rscript
# In-zone channel split for high vs low tscore pitches, plus type-adjusted
# within-pitcher zone slopes.

suppressPackageStartupMessages({ library(data.table); library(splines) })
options(width = 215)

Q <- fread("data/statcast_2026/tscore_vs_zone_pitches.csv")
A <- fread("data/statcast_2026/tscore_vs_zone_pitchers.csv")

cols <- c("game_year","game_type","pitcher","pitch_type","description","balls","strikes",
          "plate_x","plate_z","sz_top","sz_bot","estimated_woba_using_speedangle",
          "delta_run_exp")
dt <- rbindlist(lapply(c(2025, 2026), function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress=FALSE, select=cols)))
setnames(dt, "estimated_woba_using_speedangle", "xw")
dt <- dt[game_type=="R" & pitch_type!="" & balls<=3 & strikes<=2 & is.finite(delta_run_exp)]
dt[, in_zone := is.finite(plate_x) & abs(plate_x)<=0.83 &
                 is.finite(plate_z) & plate_z>=sz_bot & plate_z<=sz_top]
dt[, bat_rv := delta_run_exp]
dt[, cnt := paste0(balls,"-",strikes)]
WH <- c("swinging_strike","swinging_strike_blocked","foul_tip")
mb <- lm(bat_rv ~ ns(xw,5)+factor(cnt), data=dt[description=="hit_into_play" & is.finite(xw)])
dt[, xrv_p := bat_rv]
dt[description=="hit_into_play" & is.finite(xw),
   xrv_p := predict(mb, newdata=dt[description=="hit_into_play" & is.finite(xw)])]
dt <- dt[!(description=="hit_into_play" & !is.finite(xw))]
dt[, chan := fifelse(description %in% WH, "whiff",
             fifelse(description=="hit_into_play", "bip",
             fifelse(description=="foul", "foul",
             fifelse(description=="called_strike", "cs", "other"))))]

key <- unique(Q[, .(pitcher, pitch_type, ts3)])
d2 <- merge(dt, key, by=c("pitcher","pitch_type"))

cat("############ In-zone channels (pitcher runs / 100 pitches of that type) ############\n\n")
# contribution: -100 * mean(xrv_p * I(in_zone & chan))  over ALL pitches of the type
CH <- d2[, .(
  n=.N,
  z_cs    = -100*mean(fifelse(in_zone & chan=="cs",    xrv_p, 0)),
  z_whiff = -100*mean(fifelse(in_zone & chan=="whiff", xrv_p, 0)),
  z_foul  = -100*mean(fifelse(in_zone & chan=="foul",  xrv_p, 0)),
  z_bip   = -100*mean(fifelse(in_zone & chan=="bip",   xrv_p, 0)),
  ooz_cs  = -100*mean(fifelse(!in_zone & chan=="cs",   xrv_p, 0)),
  ooz_whiff=-100*mean(fifelse(!in_zone & chan=="whiff",xrv_p, 0)),
  ooz_foul= -100*mean(fifelse(!in_zone & chan=="foul", xrv_p, 0)),
  ooz_bip = -100*mean(fifelse(!in_zone & chan=="bip",  xrv_p, 0)),
  ooz_ball= -100*mean(fifelse(!in_zone & chan=="other",xrv_p, 0)),
  total   = -100*mean(xrv_p)
), by=.(pitcher, pitch_type, ts3)]

SUM <- CH[, lapply(.SD, mean),
          .SDcols=c("z_cs","z_whiff","z_foul","z_bip","ooz_whiff","ooz_foul","ooz_bip","ooz_ball","total"),
          by=ts3]
print(SUM[, lapply(.SD, function(x) if(is.numeric(x)) round(x,3) else x)][order(ts3)],
      row.names=FALSE)

cat("\n  High minus low, in-zone channels:\n")
h <- CH[ts3=="high tscore pitch"]; l <- CH[ts3=="low tscore pitch"]
for (v in c("z_cs","z_whiff","z_foul","z_bip","ooz_whiff","ooz_bip","ooz_ball","total")) {
  t <- t.test(h[[v]], l[[v]])
  cat(sprintf("     %-10s  %+6.3f vs %+6.3f   diff %+6.3f   p = %.3f\n",
              v, mean(h[[v]]), mean(l[[v]]), mean(h[[v]])-mean(l[[v]]), t$p.value))
}

# rate when the pitch IS in the zone
cat("\n############ Per in-zone pitch (not contribution) ############\n\n")
Z <- d2[in_zone==TRUE, .(
  n=.N,
  xrv=-100*mean(xrv_p),
  cs=100*mean(chan=="cs"),
  whiff=100*mean(chan=="whiff"),
  foul=100*mean(chan=="foul"),
  bip=100*mean(chan=="bip"),
  swing=100*mean(chan %in% c("whiff","foul","bip"))
), by=.(pitcher, pitch_type, ts3)]
print(Z[, .(cells=.N,
            xrv=round(mean(xrv),3),
            cs=round(mean(cs),1),
            whiff=round(mean(whiff),1),
            foul=round(mean(foul),1),
            bip=round(mean(bip),1),
            swing=round(mean(swing),1)), by=ts3][order(ts3)], row.names=FALSE)

# type-adjusted within-pitcher: book xRV ~ best zone residual + rest zone residual
cat("\n############ Type-adjusted within-pitcher zone slopes ############\n\n")
Q[, zone_a := residuals(lm(zone ~ factor(pitch_type)))]
BEST <- Q[, .SD[which.max(tscore_a)], by=pitcher]
REST <- Q[!BEST, on=.(pitcher, pitch_type)]
REST_A <- REST[, .(rest_zone=weighted.mean(zone, n),
                   rest_zone_a=weighted.mean(zone_a, n),
                   rest_xrv=weighted.mean(xrv, n)), by=pitcher]
W <- merge(BEST[, .(pitcher, best_type=pitch_type, best_zone=zone,
                    best_zone_a=zone_a, best_xrv=xrv, tscore, tscore_a)],
           REST_A, by="pitcher")
W <- merge(W, A[, .(pitcher, pit_xrv=xrv, pit_tscore=tscore)], by="pitcher")
W <- W[is.finite(rest_zone)]
cat("  Book xRV ~ type-adj timing-pitch zone + rest zone + tscore\n")
print(round(summary(lm(pit_xrv ~ best_zone_a + rest_zone + pit_tscore, data=W))$coefficients, 4))
cat("\n  Timing-pitch xRV ~ type-adj own zone  (does THIS pitch get more from zone?)\n")
print(round(summary(lm(best_xrv ~ best_zone_a + factor(best_type), data=W))$coefficients[1:2,,drop=FALSE], 4))

# slope of pitch xRV ~ zone, high vs low, OFF only (the interesting group)
cat("\n  OFF only, xRV ~ zone + type, by tscore tercile:\n")
for (lab in c("low tscore pitch","mid tscore pitch","high tscore pitch")) {
  d <- Q[pgrp != "FB" & ts3==lab]
  if (nrow(d) < 30) next
  m <- lm(xrv ~ zone + factor(pitch_type), data=d)
  cf <- summary(m)$coefficients
  cat(sprintf("     %-20s  n=%d  β = %+.4f  p = %.3f\n",
              lab, nrow(d), cf["zone","Estimate"], cf["zone","Pr(>|t|)"]))
}
cat("\n  OFF interaction xRV ~ zone * hi + type\n")
Qo <- Q[pgrp != "FB"]
Qo[, hi := as.integer(ts3=="high tscore pitch")]
print(round(summary(lm(xrv ~ zone * hi + factor(pitch_type), data=Qo))$coefficients[
  c("zone","hi","zone:hi"),], 4))
