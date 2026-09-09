#!/usr/bin/env Rscript

# Does a high-tscore arsenal throw more strikes?
# And does zone% on high-tscore pitches buy more (or less) success than zone%
# on everything else?

suppressPackageStartupMessages({ library(data.table); library(splines) })
set.seed(41)
options(width = 215)

cols <- c("game_year","game_type","player_name","pitcher","batter","stand","pitch_type",
          "description","bb_type","events","balls","strikes","plate_x","plate_z",
          "sz_top","sz_bot","release_speed","launch_speed","launch_speed_angle",
          "estimated_woba_using_speedangle","delta_run_exp",
          "intercept_ball_minus_batter_pos_y_inches")

dt <- rbindlist(lapply(c(2025, 2026), function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress = FALSE, select = cols)))
setnames(dt, c("estimated_woba_using_speedangle",
               "intercept_ball_minus_batter_pos_y_inches"), c("xw","depth"))
dt <- dt[game_type == "R" & pitch_type != "" & balls <= 3 & strikes <= 2 &
         is.finite(delta_run_exp)]

FB <- c("FF","SI","FC"); BR <- c("SL","ST","CU","KC","SV","CS"); OS <- c("CH","FS","FO")
dt[, pgrp := fifelse(pitch_type %in% FB,"FB", fifelse(pitch_type %in% BR,"BR",
             fifelse(pitch_type %in% OS,"OS",NA_character_)))]
dt <- dt[!is.na(pgrp)]
dt[, pg2 := fifelse(pgrp=="FB","FB","OFF")]
dt[, in_zone := is.finite(plate_x) & abs(plate_x) <= 0.83 &
                 is.finite(plate_z) & plate_z >= sz_bot & plate_z <= sz_top]
dt[, bat_rv := delta_run_exp]
dt[, cnt := paste0(balls,"-",strikes)]

WH <- c("swinging_strike","swinging_strike_blocked","foul_tip")
SW <- c(WH, "foul", "hit_into_play")
CS <- "called_strike"
dt[, swung := description %in% SW]
dt[, strike := description %in% c(WH, CS, "foul")]

# xRV: BIP repriced by xwOBA + count; takes/whiffs keep realized
mb <- lm(bat_rv ~ ns(xw,5) + factor(cnt),
         data=dt[description=="hit_into_play" & is.finite(xw)])
dt[, xrv_p := bat_rv]
dt[description=="hit_into_play" & is.finite(xw),
   xrv_p := predict(mb, newdata=dt[description=="hit_into_play" & is.finite(xw)])]
dt <- dt[!(description=="hit_into_play" & !is.finite(xw))]

# ---- tscore on swings (same construction as timing_directional_score.R)
SWING <- SW
sw <- dt[description %in% SWING & is.finite(depth) & is.finite(plate_x) &
         is.finite(plate_z) & is.finite(release_speed)]
sw[, px_bat := fifelse(stand=="R", -plate_x, plate_x)]
sw[, pz_rel := (plate_z - sz_bot)/pmax(sz_top - sz_bot, 0.1)]
fitd <- lm(depth ~ ns(px_bat,5)*pgrp + ns(pz_rel,5) + ns(release_speed,4) + stand, data=sw)
sw[, r1 := residuals(fitd)]
sw[, nb := .N, by=batter]; sw <- sw[nb >= 200]
sw[, tdev := r1 - mean(r1), by=batter]
surf <- lm(bat_rv ~ ns(tdev, 6)*pgrp + factor(cnt), data=sw)
sw[, trv := predict(surf, newdata=sw)]
sw[, trv_b := trv - mean(trv), by=batter]
sw[, trv_ex := trv_b - mean(trv_b), by=pg2]

# pitcher tscore
PIT <- sw[, .(
  swings = .N,
  tscore = -100*mean(trv_ex),
  fb_rate = -100*mean(trv_ex[pg2=="FB"]),
  off_rate = -100*mean(trv_ex[pg2=="OFF"]),
  fb_part = -100*sum(trv_ex[pg2=="FB"])/.N,
  off_part = -100*sum(trv_ex[pg2=="OFF"])/.N
), by=pitcher]

# pitcher book, all pitches (not just swings)
BK <- dt[, .(
  pitches = .N,
  name = player_name[1],
  xrv = -100*mean(xrv_p),
  rv = -100*mean(bat_rv),
  zone = 100*mean(in_zone, na.rm=TRUE),
  strike = 100*mean(strike),
  csw = 100*mean(description %in% c(WH, CS)),
  cs = 100*mean(description==CS),
  ball = 100*mean(description=="ball"),
  bb = 100*mean(events=="walk", na.rm=TRUE),
  k = 100*mean(events %in% c("strikeout","strikeout_double_play"), na.rm=TRUE),
  chase_ooz = 100*mean(swung[!in_zone]),
  swing_z = 100*mean(swung[in_zone]),
  fb_sh = 100*mean(pg2=="FB"),
  zone_fb = 100*mean(in_zone[pg2=="FB"], na.rm=TRUE),
  zone_off = 100*mean(in_zone[pg2=="OFF"], na.rm=TRUE),
  xrv_z = -100*mean(xrv_p[in_zone]),
  xrv_ooz = -100*mean(xrv_p[!in_zone])
), by=pitcher]

A <- merge(PIT, BK, by="pitcher")
A <- A[pitches >= 700 & swings >= 300]
A[, ts3 := cut(tscore, quantile(tscore, c(0,1/3,2/3,1)),
               labels=c("low tscore","mid tscore","high tscore"), include.lowest=TRUE)]

cat(sprintf("Pitchers (700+ pitches, 300+ scored swings): %d\n\n", nrow(A)))

# =============================================================================
# 1. Arsenal tscore vs strike throwing
# =============================================================================
cat("############ 1. High-tscore arsenals and strike throwing ############\n\n")
print(A[, .(
  n=.N,
  tscore=round(mean(tscore),2),
  zone=round(mean(zone),1),
  strike=round(mean(strike),1),
  csw=round(mean(csw),1),
  cs=round(mean(cs),1),
  ball=round(mean(ball),1),
  bb=round(mean(bb),2),
  k=round(mean(k),2),
  chase=round(mean(chase_ooz),1),
  swing_z=round(mean(swing_z),1),
  zone_fb=round(mean(zone_fb),1),
  zone_off=round(mean(zone_off),1),
  xrv=round(mean(xrv),3)
), by=ts3][order(ts3)], row.names=FALSE)

cat("\n  Correlations with tscore (pitcher level):\n")
for (v in c("zone","strike","csw","cs","ball","bb","k","chase_ooz","swing_z",
            "zone_fb","zone_off","xrv","fb_sh")) {
  r <- cor(A$tscore, A[[v]], use="complete.obs")
  p <- cor.test(A$tscore, A[[v]])$p.value
  cat(sprintf("     %-10s r = %+.3f   p = %.3f\n", v, r, p))
}

cat("\n  Partial: zone ~ tscore + FB share + FB tscore rate + off tscore rate\n")
print(round(summary(lm(zone ~ tscore + fb_sh, data=A))$coefficients, 4))
cat("\n  zone_off ~ off_rate + fb_sh\n")
print(round(summary(lm(zone_off ~ off_rate + fb_sh, data=A))$coefficients, 4))
cat("\n  zone_fb ~ fb_rate + fb_sh\n")
print(round(summary(lm(zone_fb ~ fb_rate + fb_sh, data=A))$coefficients, 4))

# Welch high vs low
tt <- function(v, lab) {
  a <- A[ts3=="high tscore"][[v]]; b <- A[ts3=="low tscore"][[v]]
  t <- t.test(a, b)
  cat(sprintf("     %-10s high %+6.2f vs low %+6.2f   diff %+6.2f   p = %.3f\n",
              lab, mean(a,na.rm=TRUE), mean(b,na.rm=TRUE),
              mean(a,na.rm=TRUE)-mean(b,na.rm=TRUE), t$p.value))
}
cat("\n  High vs low tscore tercile:\n")
for (v in c("zone","strike","csw","cs","ball","bb","k","chase_ooz","swing_z",
            "zone_fb","zone_off","xrv")) tt(v, v)

# =============================================================================
# 2. Pitcher x pitch_type: is zone more valuable on high-tscore pitches?
# =============================================================================
cat("\n############ 2. Zone% value on high-tscore pitches vs other pitches ############\n\n")

PT <- sw[, .(
  swings=.N,
  tscore=-100*mean(trv_ex)
), by=.(pitcher, pitch_type, pgrp, pg2)]

BKPT <- dt[, .(
  n=.N,
  zone=100*mean(in_zone, na.rm=TRUE),
  strike=100*mean(strike),
  xrv=-100*mean(xrv_p),
  xrv_z=-100*mean(xrv_p[in_zone]),
  xrv_ooz=-100*mean(xrv_p[!in_zone]),
  n_z=sum(in_zone, na.rm=TRUE),
  n_ooz=sum(!in_zone, na.rm=TRUE),
  cs=100*mean(description==CS),
  whiff=100*mean(description %in% WH),
  chase_ooz=100*mean(swung[!in_zone]),
  swing_z=100*mean(swung[in_zone])
), by=.(pitcher, pitch_type)]

Q <- merge(PT, BKPT, by=c("pitcher","pitch_type"))
Q <- Q[n >= 150 & swings >= 60]
Q <- merge(Q, A[, .(pitcher, pit_tscore=tscore, pit_xrv=xrv, pit_zone=zone)], by="pitcher")

# residualize by pitch type so curves vs sliders don't drive it
Q[, tscore_a := residuals(lm(tscore ~ factor(pitch_type)))]
Q[, zone_a   := residuals(lm(zone ~ factor(pitch_type)))]
Q[, xrv_a    := residuals(lm(xrv ~ factor(pitch_type)))]
Q[, hi_ts := tscore_a >= quantile(tscore_a, 2/3)]
Q[, lo_ts := tscore_a <= quantile(tscore_a, 1/3)]
Q[, ts3 := fifelse(hi_ts, "high tscore pitch",
            fifelse(lo_ts, "low tscore pitch", "mid tscore pitch"))]
Q[, ts3 := factor(ts3, levels=c("low tscore pitch","mid tscore pitch","high tscore pitch"))]

cat(sprintf("  Pitcher x pitch_type cells: %d (150+ pitches, 60+ scored swings)\n\n", nrow(Q)))
print(Q[, .(
  cells=.N,
  tscore=round(mean(tscore),2),
  zone=round(mean(zone),1),
  zone_a=round(mean(zone_a),2),
  strike=round(mean(strike),1),
  xrv=round(mean(xrv),3),
  xrv_z=round(mean(xrv_z),3),
  xrv_ooz=round(mean(xrv_ooz),3),
  z_minus_ooz=round(mean(xrv_z - xrv_ooz),3),
  cs=round(mean(cs),1),
  whiff=round(mean(whiff),1),
  chase=round(mean(chase_ooz),1),
  swing_z=round(mean(swing_z),1)
), by=ts3][order(ts3)], row.names=FALSE)

cat("\n  Type-adjusted cor(tscore_a, zone_a) = ",
    sprintf("%+.3f  p = %.3f\n",
            cor(Q$tscore_a, Q$zone_a),
            cor.test(Q$tscore_a, Q$zone_a)$p.value))
cat(sprintf("  Type-adjusted cor(tscore_a, xrv_a)  = %+.3f\n", cor(Q$tscore_a, Q$xrv_a)))
cat(sprintf("  Type-adjusted cor(zone_a, xrv_a)    = %+.3f\n", cor(Q$zone_a, Q$xrv_a)))

cat("\n  xRV ~ zone  (type FE), split by tscore tercile:\n")
for (lab in levels(Q$ts3)) {
  d <- Q[ts3==lab]
  m <- lm(xrv ~ zone + factor(pitch_type), data=d)
  cf <- summary(m)$coefficients
  cat(sprintf("     %-20s  β(zone) = %+.4f per pt   p = %.3f   n = %d\n",
              lab, cf["zone","Estimate"], cf["zone","Pr(>|t|)"], nrow(d)))
}

cat("\n  Interaction: xRV ~ zone * high-tscore + pitch type\n")
Q[, hi := as.integer(hi_ts)]
m_int <- lm(xrv ~ zone * hi + factor(pitch_type), data=Q)
print(round(summary(m_int)$coefficients[c("zone","hi","zone:hi"),,drop=FALSE], 4))

cat("\n  Same, type-adjusted residuals: xrv_a ~ zone_a * hi\n")
print(round(summary(lm(xrv_a ~ zone_a * hi, data=Q))$coefficients, 4))

cat("\n  Does zone buy MORE xRV on high-tscore pitches? Interaction p above.\n")
cat("  Positive zone:hi means each extra zone point is worth more when the pitch is a timing pitch.\n")

# =============================================================================
# 3. Within-pitcher: their high-tscore pitch vs their other pitches
# =============================================================================
cat("\n############ 3. Within the same pitcher: timing pitch vs the rest ############\n\n")

# each pitcher's highest-tscore pitch type (min 150) vs the rest of their book
BEST <- Q[, .SD[which.max(tscore_a)], by=pitcher]
REST <- Q[!BEST, on=.(pitcher, pitch_type)]
REST_A <- REST[, .(
  rest_n=sum(n),
  rest_zone=weighted.mean(zone, n),
  rest_xrv=weighted.mean(xrv, n),
  rest_xrv_z=weighted.mean(xrv_z, n),
  rest_xrv_ooz=weighted.mean(xrv_ooz, n),
  rest_tscore=weighted.mean(tscore, n)
), by=pitcher]
W <- merge(BEST[, .(pitcher, pit_tscore, pit_xrv, pit_zone,
                    best_type=pitch_type, best_ts=tscore, best_ts_a=tscore_a,
                    best_zone=zone, best_xrv=xrv, best_xrv_z=xrv_z,
                    best_xrv_ooz=xrv_ooz, best_n=n)],
           REST_A, by="pitcher")
W <- W[is.finite(rest_zone) & rest_n >= 200]
W[, zone_gap := best_zone - rest_zone]
W[, xrv_gap  := best_xrv - rest_xrv]
W[, zval_best := best_xrv_z - best_xrv_ooz]
W[, zval_rest := rest_xrv_z - rest_xrv_ooz]

cat(sprintf("  Pitchers with a scored timing pitch and 200+ other pitches: %d\n\n", nrow(W)))
cat("  Mean zone% / xRV of the pitcher's highest-tscore pitch vs the rest of his book:\n")
cat(sprintf("     timing pitch   zone %.1f   xRV %+.3f   in-zone xRV %+.3f   ooz xRV %+.3f   (in-out %+0.3f)\n",
            mean(W$best_zone), mean(W$best_xrv), mean(W$best_xrv_z),
            mean(W$best_xrv_ooz), mean(W$zval_best)))
cat(sprintf("     rest of book   zone %.1f   xRV %+.3f   in-zone xRV %+.3f   ooz xRV %+.3f   (in-out %+0.3f)\n",
            mean(W$rest_zone), mean(W$rest_xrv), mean(W$rest_xrv_z),
            mean(W$rest_xrv_ooz), mean(W$zval_rest)))
cat(sprintf("     paired t on zone%%:     p = %.3f\n", t.test(W$best_zone, W$rest_zone, paired=TRUE)$p.value))
cat(sprintf("     paired t on (in-ooz) xRV premium: p = %.3f\n",
            t.test(W$zval_best, W$zval_rest, paired=TRUE)$p.value))

# among high-arsenal-tscore pitchers, is the timing pitch's zone even more decisive?
W[, hi_ars := pit_tscore >= quantile(pit_tscore, 2/3)]
cat("\n  In-zone minus out-of-zone xRV premium, by arsenal tscore:\n")
print(W[, .(
  n=.N,
  timing_premium=round(mean(zval_best),3),
  rest_premium=round(mean(zval_rest),3),
  extra=round(mean(zval_best - zval_rest),3)
), by=.(hi_ars)][order(-hi_ars)], row.names=FALSE)

# does zone on the timing pitch predict book xRV more than zone on the rest?
cat("\n  Book xRV ~ timing-pitch zone + rest-of-book zone + arsenal tscore\n")
print(round(summary(lm(pit_xrv ~ best_zone + rest_zone + pit_tscore, data=W))$coefficients, 4))

cat("\n  Standardized (1 SD) for comparison:\n")
Ws <- copy(W)
for (v in c("best_zone","rest_zone","pit_tscore","best_xrv","rest_xrv"))
  Ws[[paste0(v,"_z")]] <- as.numeric(scale(Ws[[v]]))
print(round(summary(lm(pit_xrv ~ best_zone_z + rest_zone_z + pit_tscore_z, data=Ws))$coefficients, 4))

# =============================================================================
# 4. Where does the run value of high-tscore pitches actually live?
# =============================================================================
cat("\n############ 4. Share of xRV coming from in-zone vs out-of-zone ############\n\n")
cat("  Channel: pitcher xRV contribution = -100 * mean(xrv_p * I(location)).\n")
cat("  These add back to total xRV. Higher = more pitcher value from that location.\n\n")

# pitch-level contribution, then average within tscore tercile of the pitch type
CH <- dt[, .(
  n=.N,
  total=-100*mean(xrv_p),
  from_z=-100*mean(fifelse(in_zone, xrv_p, 0)),
  from_ooz=-100*mean(fifelse(!in_zone, xrv_p, 0)),
  zone=100*mean(in_zone)
), by=.(pitcher, pitch_type)]
CH <- merge(CH, Q[, .(pitcher, pitch_type, tscore, tscore_a, ts3, xrv)],
            by=c("pitcher","pitch_type"))

print(CH[, .(
  cells=.N,
  zone=round(mean(zone),1),
  xrv=round(mean(total),3),
  from_in_zone=round(mean(from_z),3),
  from_out=round(mean(from_ooz),3),
  share_from_zone=round(100*mean(from_z)/mean(total),1)
), by=ts3][order(ts3)], row.names=FALSE)

# in-zone vs ooz RATE (already have xrv_z / xrv_ooz) — this is the per-pitch value
cat("\n  Per-pitch xRV when the pitch is in the zone vs out, by tscore tercile:\n")
print(Q[, .(
  cells=.N,
  in_zone_xrv=round(mean(xrv_z),3),
  ooz_xrv=round(mean(xrv_ooz),3),
  premium=round(mean(xrv_z - xrv_ooz),3)
), by=ts3][order(ts3)], row.names=FALSE)

hi <- Q[hi_ts==TRUE]; lo <- Q[lo_ts==TRUE]
cat(sprintf("\n  High vs low tscore pitches, in-zone xRV:  %+.3f vs %+.3f  p = %.3f\n",
            mean(hi$xrv_z), mean(lo$xrv_z), t.test(hi$xrv_z, lo$xrv_z)$p.value))
cat(sprintf("  High vs low tscore pitches, ooz xRV:      %+.3f vs %+.3f  p = %.3f\n",
            mean(hi$xrv_ooz), mean(lo$xrv_ooz), t.test(hi$xrv_ooz, lo$xrv_ooz)$p.value))
cat(sprintf("  High vs low, in-zone premium (in-ooz):    %+.3f vs %+.3f  p = %.3f\n",
            mean(hi$xrv_z - hi$xrv_ooz), mean(lo$xrv_z - lo$xrv_ooz),
            t.test(hi$xrv_z - hi$xrv_ooz, lo$xrv_z - lo$xrv_ooz)$p.value))

# FB vs OFF separately — timing direction differs
cat("\n  Same premium, fastballs only / secondaries only:\n")
for (g in c("FB","OFF")) {
  d <- Q[pg2==g]
  cat(sprintf("  %s\n", g))
  print(d[, .(cells=.N,
              in_z=round(mean(xrv_z),3), ooz=round(mean(xrv_ooz),3),
              prem=round(mean(xrv_z-xrv_ooz),3),
              zone=round(mean(zone),1)), by=ts3][order(ts3)], row.names=FALSE)
}

# =============================================================================
# 5. Year-ahead: 2025 tscore -> 2026 zone?  2025 zone on timing pitch -> 2026 xRV?
# =============================================================================
cat("\n############ 5. Year-ahead (does this stick?) ############\n\n")

sw[, yr := game_year]
PIT_Y <- sw[, .(swings=.N, tscore=-100*mean(trv_ex)), by=.(pitcher, yr)]
BK_Y <- dt[, .(
  pitches=.N,
  zone=100*mean(in_zone),
  xrv=-100*mean(xrv_p),
  strike=100*mean(strike)
), by=.(pitcher, game_year)]
setnames(BK_Y, "game_year", "yr")
Y <- merge(PIT_Y, BK_Y, by=c("pitcher","yr"))
Y25 <- Y[yr==2025 & pitches>=400 & swings>=200]
Y26 <- Y[yr==2026 & pitches>=400 & swings>=200]
YY <- merge(Y25[, .(pitcher, ts25=tscore, zone25=zone, xrv25=xrv, st25=strike)],
            Y26[, .(pitcher, ts26=tscore, zone26=zone, xrv26=xrv, st26=strike)],
            by="pitcher")
cat(sprintf("  Pitchers in both seasons (400+ each): %d\n", nrow(YY)))
cat(sprintf("  cor(2025 tscore, 2026 zone)    = %+.3f\n", cor(YY$ts25, YY$zone26)))
cat(sprintf("  cor(2025 tscore, 2026 strike)  = %+.3f\n", cor(YY$ts25, YY$st26)))
cat(sprintf("  cor(2025 tscore, 2026 xRV)     = %+.3f\n", cor(YY$ts25, YY$xrv26)))
cat(sprintf("  cor(2025 zone,   2026 xRV)     = %+.3f\n", cor(YY$zone25, YY$xrv26)))
cat("\n  2026 xRV ~ 2025 xRV + 2025 tscore + 2025 zone\n")
print(round(summary(lm(xrv26 ~ xrv25 + ts25 + zone25, data=YY))$coefficients, 4))

# save a slim table for a canvas if needed
fwrite(A[, .(pitcher, name, pitches, swings, tscore, fb_rate, off_rate, ts3,
             zone, strike, csw, cs, ball, bb, k, chase_ooz, swing_z,
             zone_fb, zone_off, xrv)],
       "data/statcast_2026/tscore_vs_zone_pitchers.csv")
fwrite(Q[, .(pitcher, pitch_type, pgrp, n, swings, tscore, tscore_a, ts3,
             zone, zone_a, xrv, xrv_a, xrv_z, xrv_ooz, strike, cs, whiff)],
       "data/statcast_2026/tscore_vs_zone_pitches.csv")
cat("\nWrote pitcher and pitch-type tables.\n")
