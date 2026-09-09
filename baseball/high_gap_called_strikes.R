#!/usr/bin/env Rscript

# High velocity-gap secondaries: called-strike ability, and where they underperform.

suppressPackageStartupMessages({ library(data.table); library(splines) })
set.seed(7)
options(width = 215)

P <- fread("data/statcast_2026/velo_gap_tunnel_pairs.csv")
P[, vg3 := cut(velo_gap, quantile(velo_gap, 0:3/3), include.lowest=TRUE,
               labels=c("small gap","mid gap","big gap"))]

cols <- c("game_year","game_type","pitcher","batter","stand","p_throws","pitch_type",
          "description","balls","strikes","plate_x","plate_z","sz_top","sz_bot",
          "release_speed","pfx_x","pfx_z","estimated_woba_using_speedangle",
          "delta_run_exp","zone")
dt <- rbindlist(lapply(c(2025, 2026), function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress=FALSE, select=cols)))
setnames(dt, "estimated_woba_using_speedangle", "xw")
dt <- dt[game_type=="R" & pitch_type!="" & balls<=3 & strikes<=2 & is.finite(delta_run_exp) &
         is.finite(plate_x) & is.finite(plate_z)]
dt[, `:=`(bat_rv=delta_run_exp,
          px_bat=fifelse(stand=="R", -plate_x, plate_x),
          pz_rel=(plate_z-sz_bot)/pmax(sz_top-sz_bot, 0.1),
          cnt=paste0(balls,"-",strikes),
          zc=(sz_top+sz_bot)/2)]
# shadow / heart / chase geometry (Statcast-ish)
dt[, in_heart := abs(px_bat) <= 0.56 & pz_rel >= 0.25 & pz_rel <= 0.75]
dt[, in_zone  := abs(plate_x) <= 0.83 & plate_z >= sz_bot & plate_z <= sz_top]
dt[, in_shadow:= !in_zone & abs(plate_x) <= 1.17 &
                 plate_z >= sz_bot-0.28 & plate_z <= sz_top+0.28]
WH <- c("swinging_strike","swinging_strike_blocked","foul_tip")
SW <- c(WH,"foul","hit_into_play")
CS <- "called_strike"
dt[, swung := description %in% SW]
dt[, called := description == CS]
dt[, ball   := description == "ball"]
dt[, chase  := swung & !in_zone]
dt[, take   := !swung]

mb <- lm(bat_rv ~ ns(xw,5)+factor(cnt), data=dt[description=="hit_into_play" & is.finite(xw)])
dt[, xrv_p := bat_rv]
dt[description=="hit_into_play" & is.finite(xw),
   xrv_p := predict(mb, newdata=dt[description=="hit_into_play" & is.finite(xw)])]
dt <- dt[!(description=="hit_into_play" & !is.finite(xw))]

key <- P[, .(pitcher, pitch_type, vg3, velo_gap, tunnel_dist, cell_tun=fifelse(
  velo_gap>=quantile(P$velo_gap,2/3) & tunnel_dist_a<=quantile(P$tunnel_dist_a,1/3), "big+tight",
  fifelse(velo_gap>=quantile(P$velo_gap,2/3) & tunnel_dist_a>=quantile(P$tunnel_dist_a,2/3),
          "big+wide", "other")))]
sec <- merge(dt, key, by=c("pitcher","pitch_type"))
cat(sprintf("Secondary pitches: %s\n\n", format(nrow(sec), big.mark=",")))

# =============================================================================
# 1. Called-strike profile by velo-gap tercile
# =============================================================================
cat("############ 1. Called-strike ability by velocity-gap tercile ############\n\n")
cat("  Pair-level rates, then type-adjusted. Higher xRV = better for pitcher.\n\n")

PR <- sec[, .(
  n=.N,
  xrv=-100*mean(xrv_p),
  zone=100*mean(in_zone),
  heart=100*mean(in_heart),
  shadow=100*mean(in_shadow),
  cs=100*mean(called),
  cs_zone=100*mean(called[in_zone]),          # called-strike rate when in zone
  take_zone=100*mean(take[in_zone]),
  cs_shadow=100*mean(called[in_shadow]),
  chase=100*mean(chase),
  swing=100*mean(swung),
  swing_zone=100*mean(swung[in_zone]),
  whiff=100*mean(description %in% WH),
  ball=100*mean(ball),
  take_xrv=-100*mean(fifelse(take, xrv_p, 0)),
  swing_xrv=-100*mean(fifelse(swung, xrv_p, 0)),
  cs_xrv=-100*mean(fifelse(called, xrv_p, 0)),
  ball_xrv=-100*mean(fifelse(ball, xrv_p, 0))
), by=.(pitcher, pitch_type, vg3)]

# type-adjust
for (v in setdiff(names(PR), c("pitcher","pitch_type","vg3","n"))) {
  PR[, paste0(v,"_a") := residuals(lm(get(v) ~ factor(pitch_type) + n))]
}

show <- function(cols, title) {
  cat(sprintf("  -- %s --\n", title))
  print(PR[, c(.(pairs=.N), lapply(.SD, function(x) round(mean(x),2))),
           by=vg3, .SDcols=cols][order(vg3)], row.names=FALSE)
}
show(c("xrv","zone","heart","shadow","cs","cs_zone","cs_shadow","chase",
       "swing","swing_zone","whiff","ball"), "Raw rates")
show(c("xrv_a","zone_a","heart_a","cs_a","cs_zone_a","cs_shadow_a","chase_a",
       "swing_a","swing_zone_a","whiff_a","ball_a","take_xrv_a","swing_xrv_a",
       "cs_xrv_a","ball_xrv_a"), "Type-adjusted")

cat("\n  Big vs small gap, type-adjusted t-tests:\n")
A <- PR[vg3=="big gap"]; B <- PR[vg3=="small gap"]
for (v in c("xrv_a","zone_a","heart_a","cs_a","cs_zone_a","cs_shadow_a","chase_a",
            "swing_a","swing_zone_a","whiff_a","ball_a","take_xrv_a","swing_xrv_a",
            "cs_xrv_a","ball_xrv_a")) {
  t <- t.test(A[[v]], B[[v]])
  cat(sprintf("     %-14s  %+0.3f vs %+0.3f   diff %+0.3f   p = %.3f\n",
              v, mean(A[[v]]), mean(B[[v]]), mean(A[[v]])-mean(B[[v]]), t$p.value))
}

# =============================================================================
# 2. Called-strike rate in the shadow, the actual command test
# =============================================================================
cat("\n############ 2. Shadow-zone called strikes — the steal ############\n\n")
SH <- sec[in_shadow==TRUE, .(
  n=.N,
  cs=100*mean(called),
  swing=100*mean(swung),
  xrv=-100*mean(xrv_p)
), by=.(pitcher, pitch_type, vg3, pitch_type)]
# wait, duplicate pitch_type in by
SH <- sec[in_shadow==TRUE, .(
  n=.N, pt=pitch_type[1],
  cs=100*mean(called), swing=100*mean(swung), xrv=-100*mean(xrv_p)
), by=.(pitcher, pitch_type, vg3)]
SH <- SH[n >= 20]
SH[, cs_a := residuals(lm(cs ~ factor(pitch_type)))]
cat("  Shadow pitches only, 20+ per pair:\n")
print(SH[, .(pairs=.N, cs=round(mean(cs),1), cs_a=round(mean(cs_a),2),
             swing=round(mean(swing),1), xrv=round(mean(xrv),2)), by=vg3][order(vg3)],
      row.names=FALSE)
cat(sprintf("  Big vs small shadow CS residual: p = %.3f\n",
            t.test(SH[vg3=="big gap"]$cs_a, SH[vg3=="small gap"]$cs_a)$p.value))

# =============================================================================
# 3. Where the take channel loses: count and location
# =============================================================================
cat("\n############ 3. Take-channel xRV by count ############\n\n")
sec[, kgrp := fifelse(strikes==2, "2-strike",
               fifelse(balls>=2 & strikes==0, "behind 2-0/3-0",
               fifelse(balls>strikes, "behind", "even/ahead")))]
TK <- sec[take==TRUE, .(n=.N, xrv=-100*mean(xrv_p), cs=100*mean(called),
                        zone=100*mean(in_zone)),
          by=.(pitcher, pitch_type, vg3, kgrp)]
print(TK[, .(pairs=.N, xrv=round(mean(xrv),2), cs=round(mean(cs),1),
             zone=round(mean(zone),1)), by=.(vg3, kgrp)][order(kgrp, vg3)],
      row.names=FALSE)

cat("\n  Type-adjusted take xRV by count group:\n")
TK[, xrv_a := residuals(lm(xrv ~ factor(pitch_type)+factor(kgrp)))]
print(dcast(TK[, .(v=round(mean(xrv_a),3), k=.N), by=.(vg3, kgrp)],
            kgrp ~ vg3, value.var="v"), row.names=FALSE)

# location bins on takes
sec[, loc := fifelse(in_heart, "heart",
              fifelse(in_zone, "zone (not heart)",
              fifelse(in_shadow, "shadow", "chase/waste")))]
cat("\n  Type-adjusted take xRV by location:\n")
LOC <- sec[take==TRUE, .(n=.N, xrv=-100*mean(xrv_p)),
           by=.(pitcher, pitch_type, vg3, loc)]
LOC[, xrv_a := residuals(lm(xrv ~ factor(pitch_type)+factor(loc)))]
print(dcast(LOC[, .(v=round(mean(xrv_a),3)), by=.(vg3, loc)],
            loc ~ vg3, value.var="v"), row.names=FALSE)
cat("\n  Raw called-strike % by location:\n")
print(dcast(sec[take==TRUE, .(cs=round(100*mean(called),1)), by=.(vg3, loc)],
            loc ~ vg3, value.var="cs"), row.names=FALSE)

cat("\n  Even/ahead takes only — why is take xRV worse for big-gap pairs?\n")
EA <- sec[take==TRUE & kgrp=="even/ahead", .(
  n=.N, cs=100*mean(called), ball=100*mean(ball), zone=100*mean(in_zone),
  shadow=100*mean(in_shadow), xrv=-100*mean(xrv_p)
), by=.(pitcher, pitch_type, vg3)]
for (v in c("cs","ball","zone","shadow","xrv"))
  EA[, paste0(v,"_a") := residuals(lm(get(v) ~ factor(pitch_type)+n))]
print(EA[, .(pairs=.N, cs_a=round(mean(cs_a),2), ball_a=round(mean(ball_a),2),
             zone_a=round(mean(zone_a),2), xrv_a=round(mean(xrv_a),3)),
         by=vg3][order(vg3)], row.names=FALSE)
cat(sprintf("  Even/ahead take xRV residual, big vs small: p = %.3f\n",
            t.test(EA[vg3=="big gap"]$xrv_a, EA[vg3=="small gap"]$xrv_a)$p.value))

# =============================================================================
# 4. Full underperformance map: every channel × gap
# =============================================================================
cat("\n############ 4. Full underperformance map (type-adjusted xRV contribution) ############\n\n")
sec[, chan := fifelse(called, "called strike",
               fifelse(ball, "ball",
               fifelse(description %in% WH, "whiff",
               fifelse(description=="foul", "foul",
               fifelse(description=="hit_into_play", "bip", "other")))))]
chans <- c("called strike","ball","whiff","foul","bip")
MAP <- rbindlist(lapply(chans, function(ch) {
  d <- sec[, .(n=.N,
               contrib=-100*mean(fifelse(chan==ch, xrv_p, 0)),
               rate=100*mean(chan==ch)),
           by=.(pitcher, pitch_type, vg3)]
  d[, `:=`(contrib_a=residuals(lm(contrib ~ factor(pitch_type)+n)),
           rate_a=residuals(lm(rate ~ factor(pitch_type)+n)),
           channel=ch)]
  d
}))
print(MAP[, .(pairs=.N,
              contrib=round(mean(contrib),3),
              contrib_a=round(mean(contrib_a),3),
              rate=round(mean(rate),2),
              rate_a=round(mean(rate_a),2)),
          by=.(channel, vg3)][order(channel, vg3)], row.names=FALSE)

cat("\n  Big minus small, type-adjusted contribution (negative = underperforms):\n")
CMP <- MAP[, {
  a <- contrib_a[vg3=="big gap"]; b <- contrib_a[vg3=="small gap"]
  t <- t.test(a, b)
  .(diff=round(mean(a)-mean(b),3), p=round(t$p.value,3),
    big=round(mean(a),3), small=round(mean(b),3))
}, by=channel]
setorder(CMP, diff)
print(CMP, row.names=FALSE)

# =============================================================================
# 5. Same map by pitch type, so we can see if it's curves
# =============================================================================
cat("\n############ 5. Called-strike residual by secondary type ############\n\n")
print(PR[, .(pairs=.N, cs=round(mean(cs),1), cs_a=round(mean(cs_a),2),
             zone_a=round(mean(zone_a),2), chase_a=round(mean(chase_a),2),
             take_xrv_a=round(mean(take_xrv_a),3),
             xrv_a=round(mean(xrv_a),3)),
         by=.(pitch_type, vg3)][order(pitch_type, vg3)], row.names=FALSE)

# =============================================================================
# 6. Fastball side: does a big gap also cost the primary called strikes?
# =============================================================================
cat("\n############ 6. The primary fastball of high-gap pitchers ############\n\n")
# classify the pitcher by his largest pair gap so each pitcher appears once
fbkey <- P[, .SD[which.max(velo_gap)], by=pitcher][, .(pitcher, vg3, fb_type)]
fb <- merge(dt, fbkey, by="pitcher")[pitch_type==fb_type]
FB <- fb[, .(n=.N,
             xrv=-100*mean(xrv_p),
             cs=100*mean(called),
             zone=100*mean(in_zone),
             chase=100*mean(chase),
             whiff=100*mean(description %in% WH),
             take_xrv=-100*mean(fifelse(take, xrv_p, 0))),
         by=.(pitcher, vg3, fb_type)]
FB <- FB[n >= 200]
for (v in c("xrv","cs","zone","chase","whiff","take_xrv"))
  FB[, paste0(v,"_a") := residuals(lm(get(v) ~ factor(fb_type)+n))]
print(FB[, .(pitchers=.N,
             xrv_a=round(mean(xrv_a),3), cs_a=round(mean(cs_a),2),
             zone_a=round(mean(zone_a),2), chase_a=round(mean(chase_a),2),
             whiff_a=round(mean(whiff_a),2), take_xrv_a=round(mean(take_xrv_a),3)),
         by=vg3][order(vg3)], row.names=FALSE)
cat(sprintf("  Fastball CS residual, big vs small: p = %.3f\n",
            t.test(FB[vg3=="big gap"]$cs_a, FB[vg3=="small gap"]$cs_a)$p.value))
cat(sprintf("  Fastball take-xRV residual, big vs small: p = %.3f\n",
            t.test(FB[vg3=="big gap"]$take_xrv_a, FB[vg3=="small gap"]$take_xrv_a)$p.value))
