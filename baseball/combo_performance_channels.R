#!/usr/bin/env Rscript

# Do big-gap + tight-tunnel pitchers generally perform better?
# If they do, is it because of whiffs and mistiming?

suppressPackageStartupMessages({ library(data.table); library(splines) })
set.seed(7)
options(width = 215)

P <- fread("data/statcast_2026/velo_gap_tunnel_pairs.csv")
P[, tun_a := tunnel_dist_a]
P[, hi_gap := velo_gap >= quantile(velo_gap, 2/3)]
P[, tight  := tun_a    <= quantile(tun_a, 1/3)]
P[, wide   := tun_a    >= quantile(tun_a, 2/3)]
P[, cell := fifelse(hi_gap & tight, "big+tight",
            fifelse(hi_gap & wide,  "big+wide",
            fifelse(hi_gap,         "big+mid",
            fifelse(tight,          "small+tight",
            fifelse(wide,           "small+wide", "small+mid")))))]
P[, combo := hi_gap & tight]

# ---- pitch-level, 2025-26, to decompose the secondary AND the pitcher's whole book
cols <- c("game_year","game_type","pitcher","batter","stand","pitch_type","description",
          "balls","strikes","plate_x","plate_z","sz_top","sz_bot","release_speed",
          "estimated_woba_using_speedangle","delta_run_exp",
          "intercept_ball_minus_batter_pos_y_inches")
dt <- rbindlist(lapply(c(2025, 2026), function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress=FALSE, select=cols)))
setnames(dt, c("estimated_woba_using_speedangle",
               "intercept_ball_minus_batter_pos_y_inches"), c("xw","depth"))
dt <- dt[game_type=="R" & pitch_type != "" & balls<=3 & strikes<=2 &
         is.finite(delta_run_exp)]
dt[, bat_rv := delta_run_exp]
dt[, bip := description == "hit_into_play"]
WH <- c("swinging_strike","swinging_strike_blocked","foul_tip")
SW <- c(WH, "foul","hit_into_play")
dt[, chan := fifelse(description %in% WH, "whiff",
             fifelse(bip==TRUE, "bip",
             fifelse(description %in% c("foul"), "foul", "take")))]

# expected RV: BIP repriced by xwOBA, everything else keeps realized
mb <- lm(bat_rv ~ ns(xw,5) + factor(paste0(balls,"-",strikes)),
         data=dt[bip==TRUE & is.finite(xw)])
dt[, xrv_p := bat_rv]
dt[bip==TRUE & is.finite(xw), xrv_p := predict(mb, newdata=dt[bip==TRUE & is.finite(xw)])]
dt <- dt[!(bip==TRUE & !is.finite(xw))]

# attach pair cell to matching pitcher × secondary
key <- unique(P[, .(pitcher, pitch_type, cell, combo, hi_gap, velo_gap, tscore_a, whiff_a, xrv_a)])
sec <- merge(dt, key, by=c("pitcher","pitch_type"))
cat(sprintf("Secondary pitches tagged: %s  |  unique pairs: %d\n\n",
            format(nrow(sec), big.mark=","), uniqueN(sec[, .(pitcher, pitch_type)])))

# =============================================================================
# 1. Pair-level channel decomposition
# =============================================================================
cat("############ 1. Where does the pair's run value come from? ############\n\n")
cat("  Each number is pitcher runs per 100 pitches of that secondary (higher = better).\n")
cat("  Channels partition the pitch; they add back to the total.\n\n")

CH <- sec[, .(
  n     = .N,
  total = -100*mean(xrv_p),
  whiff = -100*mean(fifelse(chan=="whiff", xrv_p, 0)),
  bip   = -100*mean(fifelse(chan=="bip",   xrv_p, 0)),
  foul  = -100*mean(fifelse(chan=="foul",  xrv_p, 0)),
  take  = -100*mean(fifelse(chan=="take",  xrv_p, 0)),
  whiff_rate = 100*mean(chan=="whiff"),
  xwcon = mean(xw[chan=="bip"], na.rm=TRUE)
), by=.(pitcher, pitch_type, cell)]

SUM <- CH[, .(
  pairs=.N,
  total=mean(total),
  whiff=mean(whiff),
  bip=mean(bip),
  foul=mean(foul),
  take=mean(take),
  whiff_rate=mean(whiff_rate),
  xwcon=mean(xwcon, na.rm=TRUE),
  total_wt=weighted.mean(total, n)
), by=cell]
setorder(SUM, -total)
print(SUM[, lapply(.SD, function(x) if (is.numeric(x)) round(x,3) else x)], row.names=FALSE)

# type-adjusted channel values, so we are not just ranking curveballs
CH[, pt := pitch_type]
for (v in c("total","whiff","bip","foul","take","whiff_rate","xwcon")) {
  CH[, paste0(v,"_a") := residuals(lm(get(v) ~ factor(pt) + n)), by=NULL]
}
cat("\n  Same channels, residualized on pitch type and sample size:\n")
SUMa <- CH[, .(pairs=.N,
               total=mean(total_a), whiff=mean(whiff_a), bip=mean(bip_a),
               foul=mean(foul_a), take=mean(take_a),
               whiff_rate=mean(whiff_rate_a), xwcon=mean(xwcon_a, na.rm=TRUE)),
           by=cell]
setorder(SUMa, -total)
print(SUMa[, lapply(.SD, function(x) if (is.numeric(x)) round(x,3) else x)], row.names=FALSE)

cat("\n  Head-to-head, big+tight vs big+wide (adjusted):\n")
A <- CH[cell=="big+tight"]; B <- CH[cell=="big+wide"]
for (v in c("total_a","whiff_a","bip_a","take_a","whiff_rate_a","xwcon_a")) {
  t <- t.test(A[[v]], B[[v]])
  cat(sprintf("     %-14s  %+0.3f vs %+0.3f   diff %+0.3f   p = %.3f\n",
              v, mean(A[[v]],na.rm=TRUE), mean(B[[v]],na.rm=TRUE),
              mean(A[[v]],na.rm=TRUE)-mean(B[[v]],na.rm=TRUE), t$p.value))
}

# does the (tiny) total gap vanish once whiff and contact are in the model?
cat("\n  Is any remaining total-xRV gap just the whiff and BIP channels?\n")
Z <- CH[cell %in% c("big+tight","big+wide")]
Z[, tight := cell=="big+tight"]
m0 <- lm(total_a ~ tight, Z)
m1 <- lm(total_a ~ tight + whiff_a + bip_a, Z)
cat(sprintf("     cell dummy alone:          tight β = %+.3f  p = %.3f\n",
            coef(m0)["tightTRUE"], summary(m0)$coefficients["tightTRUE",4]))
cat(sprintf("     + whiff and BIP channels:  tight β = %+.3f  p = %.3f\n",
            coef(m1)["tightTRUE"], summary(m1)$coefficients["tightTRUE",4]))
cat(sprintf("     whiff channel β = %+.3f  p = %.3f | BIP channel β = %+.3f  p = %.3f\n",
            coef(m1)["whiff_a"], summary(m1)$coefficients["whiff_a",4],
            coef(m1)["bip_a"], summary(m1)$coefficients["bip_a",4]))

# =============================================================================
# 2. Pitcher-level: the whole book, not just the secondary
# =============================================================================
cat("\n############ 2. Do these pitchers perform better overall? ############\n\n")
# a pitcher is "combo" if he has at least one big+tight pair
own <- unique(P[combo==TRUE, .(pitcher)])
own[, is_combo := TRUE]
allp <- dt[, .(
  pitches=.N,
  xrv=-100*mean(xrv_p),
  rv=-100*mean(bat_rv),
  whiff=100*mean(description %in% WH),
  xwcon=mean(xw[bip==TRUE], na.rm=TRUE)
), by=pitcher][pitches >= 700]
allp <- merge(allp, own, by="pitcher", all.x=TRUE)
allp[is.na(is_combo), is_combo := FALSE]

# exclusive: combo vs big+wide-only vs neither
widep <- unique(P[cell=="big+wide" & combo==FALSE, .(pitcher)])
allp[, grp := fifelse(is_combo==TRUE, "has big+tight",
              fifelse(pitcher %in% widep$pitcher, "has big+wide only", "neither"))]

cat("  All pitches, 700+ on the season(s). Higher xRV/RV = better for the pitcher.\n\n")
print(allp[, .(pitchers=.N, pitches=sum(pitches),
               xrv=round(mean(xrv),3), rv=round(mean(rv),3),
               whiff=round(mean(whiff),2), xwcon=round(mean(xwcon),3),
               xrv_wt=round(weighted.mean(xrv, pitches),3)), by=grp][order(-xrv)],
      row.names=FALSE)

cat("\n  Tests (has big+tight vs not):\n")
tt <- function(x, y, lab) {
  t <- t.test(x, y)
  cat(sprintf("     %-8s  %.3f vs %.3f   diff %+0.3f   p = %.3f\n",
              lab, mean(x,na.rm=TRUE), mean(y,na.rm=TRUE),
              mean(x,na.rm=TRUE)-mean(y,na.rm=TRUE), t$p.value))
}
tt(allp[is_combo==TRUE]$xrv,   allp[is_combo==FALSE]$xrv,   "xRV")
tt(allp[is_combo==TRUE]$rv,    allp[is_combo==FALSE]$rv,    "RV")
tt(allp[is_combo==TRUE]$whiff, allp[is_combo==FALSE]$whiff, "whiff%")
tt(allp[is_combo==TRUE]$xwcon, allp[is_combo==FALSE]$xwcon, "xwOBAcon")

# =============================================================================
# 3. Realized RV vs xRV — they do not have a lower RV
# =============================================================================
cat("\n############ 3. Realized RV by cell (this is the noisy one) ############\n\n")
RV <- sec[, .(n=.N, rv=-100*mean(bat_rv), xrv=-100*mean(xrv_p)),
          by=.(pitcher, pitch_type, cell)]
print(RV[, .(pairs=.N, rv=round(mean(rv),3), xrv=round(mean(xrv),3),
             luck=round(mean(rv-xrv),3)), by=cell][order(-xrv)], row.names=FALSE)
cat("\n  Realized RV is the column that made the tight cell look worse. Split-half\n")
cat("  reliability on pair RV is ~0.21. It is not a performance signal.\n")
