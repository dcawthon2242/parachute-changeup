#!/usr/bin/env Rscript

# Does contact depth explain which pitch types manage contact?
#
# Motivating question: the cutter is conventionally described as a contact-management
# pitch, but the directional timing work found it is the fastball hitters get FURTHEST
# out front on — supposedly the worst place for a fastball to be met. Either the timing
# story is incomplete or the cutter's reputation is.
#
# Statcast tracks BOTH intercept axes:
#   ..._y_inches = depth  (fore/aft; how far in front of the batter contact happens)
#   ..._x_inches = lateral (how far from the batter's body, i.e. barrel vs. jammed)
# The timing work only used depth. A cutter's job is described as breaking bats and
# producing end-of-bat/handle contact, which lives on the LATERAL axis. This script
# checks both.

suppressPackageStartupMessages({ library(data.table); library(splines) })

cols <- c("game_year","game_type","player_name","pitcher","batter","stand","p_throws",
          "pitch_type","description","bb_type","events","balls","strikes",
          "plate_x","plate_z","sz_top","sz_bot","release_speed",
          "launch_speed","launch_angle","launch_speed_angle",
          "estimated_woba_using_speedangle","woba_value","woba_denom","delta_run_exp",
          "intercept_ball_minus_batter_pos_y_inches",
          "intercept_ball_minus_batter_pos_x_inches")

dt <- rbindlist(lapply(c(2025, 2026), function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress = FALSE, select = cols)))
setnames(dt, c("intercept_ball_minus_batter_pos_y_inches",
               "intercept_ball_minus_batter_pos_x_inches"), c("depth","lat_raw"))
dt <- dt[game_type == "R" & pitch_type != "" & is.finite(delta_run_exp)]

SWING <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip","hit_into_play")
sw <- dt[description %in% SWING & is.finite(depth) & is.finite(plate_x) &
         is.finite(plate_z) & is.finite(release_speed)]
FB<-c("FF","SI","FC"); BR<-c("SL","ST","CU","KC","SV","CS"); OS<-c("CH","FS","FO")
sw[, pgrp := fifelse(pitch_type %in% FB,"FB", fifelse(pitch_type %in% BR,"BR",
             fifelse(pitch_type %in% OS,"OS",NA_character_)))]
sw <- sw[!is.na(pgrp)]
sw[, pg2 := fifelse(pgrp=="FB","FB","OFF")]
sw[, px_bat := fifelse(stand=="R", -plate_x, plate_x)]   # positive = inside to the batter
sw[, pz_rel := (plate_z - sz_bot)/pmax(sz_top - sz_bot, 0.1)]

# ---- orient the lateral axis empirically -------------------------------------
cat("############ 0. Orienting the lateral intercept axis ############\n")
print(sw[, .(n=.N, mean_lat_raw=round(mean(lat_raw, na.rm=TRUE),2)), by=stand],
      row.names=FALSE)
cat("Contact happens out away from the body, so the sign that yields a positive mean\n")
cat("is 'away from the batter'. Flipping RHH to put both handednesses on one scale.\n")
sw[, lat := fifelse(stand=="R", -lat_raw, lat_raw)]
if (mean(sw$lat, na.rm=TRUE) < 0) sw[, lat := -lat]
cat(sprintf("After flip: mean lateral reach = %.2f in (positive = farther from the body,\n", 
            mean(sw$lat, na.rm=TRUE)))
cat("            i.e. out toward the barrel/end; negative = jammed in on the hands)\n\n")

# ---- the two deviation axes --------------------------------------------------
sw <- sw[is.finite(lat)]
fitd <- lm(depth ~ ns(px_bat,5)*pgrp + ns(pz_rel,5) + ns(release_speed,4) + stand, data=sw)
sw[, r1 := residuals(fitd)]
fitl <- lm(lat ~ ns(px_bat,5)*pgrp + ns(pz_rel,5) + ns(release_speed,4) + stand, data=sw)
sw[, l1 := residuals(fitl)]

sw[, nb := .N, by=batter]; sw <- sw[nb >= 200]
sw[, tdev := r1 - mean(r1), by=batter]      # + = out front, - = deep
sw[, ldev := l1 - mean(l1), by=batter]      # + = met farther out on the bat, - = jammed
sw[, whiff := description %in% c("swinging_strike","swinging_strike_blocked","foul_tip")]
sw[, isbip := description == "hit_into_play" & bb_type != "" & is.finite(launch_speed)]
sw[, pit_rv := -delta_run_exp]

bp <- sw[isbip == TRUE & is.finite(estimated_woba_using_speedangle)]
bp[, `:=`(xw = estimated_woba_using_speedangle,
          brl = is.finite(launch_speed_angle) & launch_speed_angle == 6,
          hard = launch_speed >= 95,
          pop = bb_type == "popup",
          gb  = bb_type == "ground_ball")]
LGXW <- mean(bp$xw)
cat(sprintf("Balls in play with tracked contact point: %s | league xwOBAcon = %.3f\n\n",
            format(nrow(bp), big.mark=","), LGXW))

# =============================================================================
# 1. Pitch-type table: where each pitch is met, and how much damage it allows
# =============================================================================
cat("############ 1. Mean contact depth vs. xwOBAcon, by pitch type ############\n")
PT <- merge(
  bp[, .(bip=.N, tdev=mean(tdev), ldev=mean(ldev), xw=mean(xw),
         wcon=sum(woba_value)/pmax(sum(woba_denom),1),
         ev=mean(launch_speed), la=mean(launch_angle, na.rm=TRUE),
         brl=100*mean(brl), hard=100*mean(hard), pop=100*mean(pop), gb=100*mean(gb)),
     by=.(pitch_type, pg2)],
  sw[, .(sw=.N, whiff=100*mean(whiff), rv=mean(pit_rv)), by=pitch_type],
  by="pitch_type")
PT <- PT[sw >= 6000]
# signed toward the pitcher's target for that group: deep on FB, out front on OFF
PT[, toward := fifelse(pg2=="FB", -tdev, tdev)]
PT[, xw_vs_lg := xw - LGXW]
setorder(PT, pg2, -sw)
print(PT[, .(pitch_type, pg2, bip, tdev=round(tdev,2), toward=round(toward,2),
             ldev=round(ldev,2), xw=round(xw,3), xw_vs_lg=round(xw_vs_lg,3),
             wcon=round(wcon,3), ev=round(ev,1), brl=round(brl,1), hard=round(hard,1),
             pop=round(pop,1), gb=round(gb,1), whiff=round(whiff,1), rv=round(rv,4))],
      row.names=FALSE)

cat("\n  Across pitch types, does depth placement track contact quality?\n")
cat(sprintf("     cor(mean tdev, xwOBAcon)   all types      = %+.3f  (n=%d types)\n",
    cor(PT$tdev, PT$xw), nrow(PT)))
cat(sprintf("     cor(toward-target, xwOBAcon) all types    = %+.3f\n", cor(PT$toward, PT$xw)))
cat(sprintf("     cor(mean tdev, xwOBAcon)   fastballs only = %+.3f  (n=%d)\n",
    PT[pg2=="FB", cor(tdev, xw)], PT[pg2=="FB", .N]))
cat(sprintf("     cor(mean tdev, xwOBAcon)   offspeed only  = %+.3f  (n=%d)\n",
    PT[pg2=="OFF", cor(tdev, xw)], PT[pg2=="OFF", .N]))
cat(sprintf("     cor(mean ldev, xwOBAcon)   all types      = %+.3f\n", cor(PT$ldev, PT$xw)))

# =============================================================================
# 2. Within pitch type, across pitchers
# =============================================================================
cat("\n############ 2. Within a pitch type, does deeper contact mean less damage? ############\n")
cat("Pitcher x pitch-type cells, min 60 balls in play. Correlations are computed inside\n")
cat("each pitch type, so they ask: among cutters, does the cutter that is met deeper\n")
cat("allow a lower xwOBAcon?\n\n")
PP <- merge(
  bp[, .(bip=.N, tdev=mean(tdev), ldev=mean(ldev), absl=mean(abs(ldev)),
         xw=mean(xw), brl=100*mean(brl), hard=100*mean(hard)),
     by=.(pitcher, pitch_type, pg2)],
  sw[, .(nsw=.N, whiff=100*mean(whiff), rv=mean(pit_rv)), by=.(pitcher, pitch_type)],
  by=c("pitcher","pitch_type"))
PP <- PP[bip >= 60]
W <- PP[, .(cells=.N,
            r_tdev_xw = cor(tdev, xw), r_tdev_brl = cor(tdev, brl),
            r_tdev_hard = cor(tdev, hard),
            r_ldev_xw = cor(ldev, xw), r_absl_xw = cor(absl, xw)),
        by=.(pitch_type, pg2)][cells >= 25]
setorder(W, pg2, -cells)
print(W[, .(pitch_type, pg2, cells, r_tdev_xw=round(r_tdev_xw,3),
            r_tdev_brl=round(r_tdev_brl,3), r_tdev_hard=round(r_tdev_hard,3),
            r_ldev_xw=round(r_ldev_xw,3), r_absl_xw=round(r_absl_xw,3))],
      row.names=FALSE)
cat("\n  Expected sign on r_tdev_xw: NEGATIVE for offspeed (out front is good for the\n")
cat("  pitcher), POSITIVE for fastballs (out front is bad, so deeper = lower xwOBAcon\n")
cat("  means the correlation should be positive).\n")

# =============================================================================
# 3. The cutter
# =============================================================================
cat("\n############ 3. Is the cutter a contact-management pitch? ############\n")
THREE <- c("FF","SI","FC")
cat("\n  -- Raw comparison against the other fastballs --\n")
print(PT[pitch_type %in% THREE,
  .(pitch_type, bip, tdev=round(tdev,2), ldev=round(ldev,2), xw=round(xw,3),
    wcon=round(wcon,3), ev=round(ev,1), la=round(la,1), brl=round(brl,1),
    hard=round(hard,1), pop=round(pop,1), gb=round(gb,1), whiff=round(whiff,1),
    rv=round(rv,4))], row.names=FALSE)

cat("\n  -- xwOBAcon by contact depth, within each fastball type --\n")
bp[, tb := cut(tdev, c(-Inf,-12,-6,-3,3,6,12,Inf),
   labels=c("deep >12","deep 6-12","deep 3-6","on time","front 3-6","front 6-12","front >12"))]
D <- bp[pitch_type %in% THREE, .(bip=.N, xw=mean(xw), brl=100*mean(brl)),
        by=.(pitch_type, tb)]
print(dcast(D, tb ~ pitch_type, value.var="xw")[, lapply(.SD, function(x)
  if (is.numeric(x)) round(x,3) else x)], row.names=FALSE)
cat("\n  Share of each pitch's contact in each depth bin (%):\n")
D[, sh := 100*bip/sum(bip), by=pitch_type]
print(dcast(D, tb ~ pitch_type, value.var="sh")[, lapply(.SD, function(x)
  if (is.numeric(x)) round(x,1) else x)], row.names=FALSE)

cat("\n  -- Depth-matched comparison: reweight FC and SI to the FF depth mix --\n")
mix <- D[pitch_type=="FF", .(tb, w = bip/sum(bip))]
adj <- merge(D, mix, by="tb")[, .(xw_adj = sum(w*xw)/sum(w),
                                  brl_adj = sum(w*brl)/sum(w)), by=pitch_type]
print(merge(PT[pitch_type %in% THREE, .(pitch_type, xw_raw=round(xw,3), brl_raw=round(brl,1))],
            adj[, .(pitch_type, xw_adj=round(xw_adj,3), brl_adj=round(brl_adj,1))],
            by="pitch_type"), row.names=FALSE)
cat("  If the cutter's edge survives depth matching, depth is not where it comes from.\n")

# =============================================================================
# 4. The lateral axis
# =============================================================================
cat("\n############ 4. The lateral axis: jam vs. barrel ############\n")
bp[, lb := cut(ldev, c(-Inf,-4,-2,2,4,Inf),
   labels=c("jammed >4in","jammed 2-4in","centered","reached 2-4in","reached >4in"))]
cat("\n  xwOBAcon, exit velo and barrel rate by lateral contact point (all pitches):\n")
L <- bp[, .(bip=.N, share=NA_real_, xw=mean(xw), ev=mean(launch_speed),
            brl=100*mean(brl), hard=100*mean(hard), pop=100*mean(pop)), by=lb]
L[, share := 100*bip/sum(bip)]
setorder(L, lb)
print(L[, .(lb, bip, share=round(share,1), xw=round(xw,3), ev=round(ev,1),
            brl=round(brl,1), hard=round(hard,1), pop=round(pop,1))], row.names=FALSE)

cat("\n  Which pitch types push contact off-center laterally?\n")
LT <- bp[, .(bip=.N, ldev=mean(ldev), absl=mean(abs(ldev)),
             jam=100*mean(ldev <= -2), reach=100*mean(ldev >= 2), xw=mean(xw)),
         by=.(pitch_type, pg2)][bip >= 2000]
setorder(LT, -absl)
print(LT[, .(pitch_type, pg2, bip, ldev=round(ldev,2), absl=round(absl,2),
             jam=round(jam,1), reach=round(reach,1), xw=round(xw,3))], row.names=FALSE)
cat(sprintf("\n     cor(mean |ldev|, xwOBAcon) across pitch types = %+.3f\n",
            cor(LT$absl, LT$xw)))
cat(sprintf("     cor(jam rate,    xwOBAcon) across pitch types = %+.3f\n",
            cor(LT$jam, LT$xw)))

cat("\n  -- Two-axis view for the three fastballs: xwOBAcon by depth x lateral --\n")
for (p in THREE) {
  cat(sprintf("\n     %s\n", p))
  M <- bp[pitch_type==p, .(n=.N, xw=mean(xw)), by=.(tb, lb)][n >= 120]
  print(dcast(M, tb ~ lb, value.var="xw")[, lapply(.SD, function(x)
    if (is.numeric(x)) round(x,3) else x)], row.names=FALSE)
}

# =============================================================================
# 5. What actually separates the three fastballs
# =============================================================================
cat("\n############ 5. Decomposing the fastball gap ############\n")
cat("Additive contribution to xwOBAcon of depth mix and lateral mix, using league\n")
cat("cell means. 'residual' is everything else (velocity, shape, location, sequencing).\n\n")
cellT <- bp[, .(mT = mean(xw)), by=tb]
cellL <- bp[, .(mL = mean(xw)), by=lb]
bp2 <- merge(merge(bp, cellT, by="tb"), cellL, by="lb")
Z <- bp2[pitch_type %in% THREE, .(
  bip = .N,
  actual   = mean(xw),
  from_dep = mean(mT) - LGXW,
  from_lat = mean(mL) - LGXW
), by=pitch_type]
Z[, total_vs_lg := actual - LGXW]
Z[, residual := total_vs_lg - from_dep - from_lat]
print(Z[, .(pitch_type, bip, actual=round(actual,3), total_vs_lg=round(total_vs_lg,4),
            from_dep=round(from_dep,4), from_lat=round(from_lat,4),
            residual=round(residual,4))], row.names=FALSE)

cat("\n  -- Best cutters in 2026 (min 60 BIP) --\n")
NAME <- unique(sw[, .(pitcher, player_name)], by="pitcher")
FCP <- merge(PP[pitch_type=="FC"], NAME, by="pitcher")
print(FCP[order(xw)][1:12, .(player_name, bip, nsw, tdev=round(tdev,1),
  ldev=round(ldev,1), xw=round(xw,3), brl=round(brl,1), hard=round(hard,1),
  whiff=round(whiff,1), rv=round(rv,4))], row.names=FALSE)
cat(sprintf("\n     Among %d qualified cutters: cor(tdev, xwOBAcon) = %+.3f, cor(ldev, xwOBAcon) = %+.3f\n",
    nrow(FCP), cor(FCP$tdev, FCP$xw), cor(FCP$ldev, FCP$xw)))

# =============================================================================
# 6. All swings vs. contact only — the ordering of the fastballs flips
# =============================================================================
cat("\n############ 6. Mean tdev on ALL swings vs. on CONTACT only ############\n")
cat("The directional-timing work measured depth over every swing, including whiffs and\n")
cat("fouls. Restricting to balls in play changes the fastball ordering, so the two\n")
cat("numbers are not interchangeable.\n\n")
cmp <- merge(
  sw[pitch_type %in% c("FF","SI","FC","SL","CH","CU","ST","FS"),
     .(all_sw=.N, tdev_all=mean(tdev)), by=pitch_type],
  bp[pitch_type %in% c("FF","SI","FC","SL","CH","CU","ST","FS"),
     .(bip=.N, tdev_bip=mean(tdev), xw=mean(xw)), by=pitch_type], by="pitch_type")
cmp <- merge(cmp, sw[, .(tdev_whiff = mean(tdev[whiff==TRUE]),
                         tdev_foul  = mean(tdev[description=="foul"])), by=pitch_type],
             by="pitch_type")
setorder(cmp, -all_sw)
print(cmp[, .(pitch_type, all_sw, bip, tdev_all=round(tdev_all,2),
              tdev_whiff=round(tdev_whiff,2), tdev_foul=round(tdev_foul,2),
              tdev_bip=round(tdev_bip,2),
              shift=round(tdev_bip - tdev_all,2), xw=round(xw,3))], row.names=FALSE)
cat("\n  cor(tdev on all swings, xwOBAcon) across these types = ")
cat(sprintf("%+.3f\n", cor(cmp$tdev_all, cmp$xw)))
cat("  cor(tdev on contact only, xwOBAcon)                  = ")
cat(sprintf("%+.3f\n", cor(cmp$tdev_bip, cmp$xw)))
cat("  Fastballs only, contact:  ")
cat(sprintf("%+.3f | all swings: %+.3f\n",
    cmp[pitch_type %in% c("FF","SI","FC"), cor(tdev_bip, xw)],
    cmp[pitch_type %in% c("FF","SI","FC"), cor(tdev_all, xw)]))

# =============================================================================
# 7. Platoon split — a cutter is two different pitches
# =============================================================================
cat("\n############ 7. Platoon split on the three fastballs ############\n")
bp[, plat := fifelse(stand == p_throws, "same-handed", "opposite")]
PL <- bp[pitch_type %in% THREE, .(bip=.N, tdev=mean(tdev), ldev=mean(ldev),
         jam=100*mean(ldev <= -2), reach=100*mean(ldev >= 2),
         xw=mean(xw), ev=mean(launch_speed), brl=100*mean(brl),
         pop=100*mean(pop), gb=100*mean(gb)), by=.(pitch_type, plat)]
setorder(PL, pitch_type, plat)
print(PL[, .(pitch_type, plat, bip, tdev=round(tdev,2), ldev=round(ldev,2),
             jam=round(jam,1), reach=round(reach,1), xw=round(xw,3),
             ev=round(ev,1), brl=round(brl,1), pop=round(pop,1), gb=round(gb,1))],
      row.names=FALSE)
cat("\n  A cutter breaks toward the pitcher's glove side: in on an opposite-handed hitter,\n")
cat("  away from a same-handed one. If the jamming story is right it should show up as a\n")
cat("  negative ldev against opposite-handed hitters specifically.\n")

# =============================================================================
# 8. Why barrel rate and xwOBAcon disagree at the out-front end
# =============================================================================
cat("\n############ 8. Batted-ball mix along the depth axis ############\n")
cat("Out-front contact raises barrel rate, but xwOBAcon peaks and turns back. The mix\n")
cat("explains it: out front adds pop-ups as fast as it adds barrels.\n\n")
MX <- bp[, .(bip=.N, xw=mean(xw), brl=100*mean(brl), pop=100*mean(pop),
             fb=100*mean(bb_type=="fly_ball"), gb=100*mean(gb)), by=tb]
setorder(MX, tb)
print(MX[, .(tb, bip, xw=round(xw,3), brl=round(brl,1), pop=round(pop,1),
             fb=round(fb,1), gb=round(gb,1))], row.names=FALSE)
cat("\n  Fastballs only:\n")
MF <- bp[pg2=="FB", .(bip=.N, xw=mean(xw), brl=100*mean(brl), pop=100*mean(pop)), by=tb]
setorder(MF, tb)
print(MF[, .(tb, bip, xw=round(xw,3), brl=round(brl,1), pop=round(pop,1))], row.names=FALSE)

cat("\n  Cross-pitch-type correlations are built on 9 points and are fragile:\n")
cat(sprintf("     all 9 types            r = %+.3f\n", cor(PT$tdev, PT$xw)))
cat(sprintf("     3 fastballs            r = %+.3f  (n=3; only says they rank-order alike)\n",
    PT[pg2=="FB", cor(tdev, xw)]))
cat(sprintf("     6 offspeed types       r = %+.3f\n", PT[pg2=="OFF", cor(tdev, xw)]))
cat(sprintf("     5 offspeed, no KC      r = %+.3f  <- the offspeed signal is entirely KC\n",
    PT[pg2=="OFF" & pitch_type != "KC", cor(tdev, xw)]))
