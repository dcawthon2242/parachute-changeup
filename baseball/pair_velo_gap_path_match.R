#!/usr/bin/env Rscript

# Does a fastball/offspeed PAIR that keeps a big velocity gap while arriving at its
# location the same way produce the most timing disruption?
#
# The hypothesis: same apparent route to the plate, different arrival time. A pitch
# pair with matched path-to-location ratio looks like it is getting to its spot by
# the same mechanism, so the only thing the hitter has to resolve is timing.
#
# Unit of analysis is a PAIR: a pitcher's primary fastball crossed with each
# secondary he throws at least 5% of the time.
#   velo_gap = primary fastball velo - secondary velo
#   p2l_gap  = |path-to-location ratio of secondary - that of the primary fastball|
#              (small = the two pitches reach their locations the same way)
#
# Outcomes are measured on the SECONDARY, since that is the pitch the pairing is
# supposed to set up:
#   off_tscore  runs saved per 100 swings from getting hitters out front (modeled)
#   off_rv      actual pitcher run value per 100 pitches (the practical test)
#   off_xw      xwOBA on contact
#   off_whiff   whiff rate per swing

suppressPackageStartupMessages({ library(data.table); library(splines) })
set.seed(7)
options(width = 210)

cols <- c("game_year","game_type","player_name","pitcher","batter","stand","p_throws",
          "pitch_type","description","balls","strikes","plate_x","plate_z",
          "sz_top","sz_bot","release_speed","release_pos_x","release_pos_z",
          "release_extension","pfx_x","pfx_z","estimated_woba_using_speedangle",
          "delta_run_exp","game_pk","intercept_ball_minus_batter_pos_y_inches")

dt <- fread("data/statcast_2026/statcast_2026_all.csv", showProgress=FALSE, select=cols)
setnames(dt, c("intercept_ball_minus_batter_pos_y_inches",
               "estimated_woba_using_speedangle"), c("depth","xw"))
dt <- dt[game_type=="R" & pitch_type != "" & balls<=3 & strikes<=2 &
         is.finite(delta_run_exp) & is.finite(plate_x) & is.finite(plate_z) &
         is.finite(release_speed)]
FB<-c("FF","SI","FC"); BR<-c("SL","ST","CU","KC","SV","CS"); OS<-c("CH","FS","FO")
dt[, pgrp := fifelse(pitch_type %in% FB,"FB", fifelse(pitch_type %in% BR,"BR",
             fifelse(pitch_type %in% OS,"OS",NA_character_)))]
dt <- dt[!is.na(pgrp)]
dt[, pg2 := fifelse(pgrp=="FB","FB","OFF")]
dt[, `:=`(px_bat = fifelse(stand=="R", -plate_x, plate_x),
          pz_rel = (plate_z - sz_bot)/pmax(sz_top - sz_bot, 0.1),
          zc     = (sz_top + sz_bot)/2,
          cnt    = paste0(balls,"-",strikes),
          bat_rv = delta_run_exp)]
dt[, d_act := sqrt(plate_x^2 + (plate_z - zc)^2)]
dt[, d_aim := sqrt((plate_x - pfx_x)^2 + (plate_z - pfx_z - zc)^2)]
cat(sprintf("2026 pitches: %s\n\n", format(nrow(dt), big.mark=",")))

# ---- timing score machinery, identical to timing_directional_score.R ---------
SWING <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip","hit_into_play")
sw <- dt[description %in% SWING & is.finite(depth)]
fitd <- lm(depth ~ ns(px_bat,5)*pgrp + ns(pz_rel,5) + ns(release_speed,4) + stand, data=sw)
sw[, r1 := residuals(fitd)]
sw[, nb := .N, by=batter]; sw <- sw[nb >= 200]
sw[, tdev := r1 - mean(r1), by=batter]
surf <- lm(bat_rv ~ ns(tdev,6)*pgrp + factor(cnt), data=sw)
sw[, trv := predict(surf, newdata=sw)]
sw[, trv_b := trv - mean(trv), by=batter]
sw[, trv_ex := trv_b - mean(trv_b), by=pg2]

# =============================================================================
# 1. Pitch-type level profile per pitcher
# =============================================================================
PT <- dt[, .(n=.N, velo=mean(release_speed), p2l=mean(d_act)/mean(d_aim),
             mx=mean(pfx_x,na.rm=TRUE), mz=mean(pfx_z,na.rm=TRUE),
             rx=mean(release_pos_x,na.rm=TRUE), rz=mean(release_pos_z,na.rm=TRUE),
             ext=mean(release_extension,na.rm=TRUE),
             rv=-100*mean(bat_rv)), by=.(pitcher, pitch_type, pg2)]
PT[, use := 100*n/sum(n), by=pitcher]

SWPT <- sw[, .(swings=.N, tscore=-100*mean(trv_ex), tdev=mean(tdev),
               whiff=100*mean(description %in%
                 c("swinging_strike","swinging_strike_blocked","foul_tip"))),
           by=.(pitcher, pitch_type)]
BPT <- dt[description=="hit_into_play" & is.finite(xw),
          .(bip=.N, xw=mean(xw)), by=.(pitcher, pitch_type)]
PT <- merge(merge(PT, SWPT, by=c("pitcher","pitch_type"), all.x=TRUE),
            BPT, by=c("pitcher","pitch_type"), all.x=TRUE)

# primary fastball = most-thrown fastball, must be a real pitch
PRIM <- PT[pg2=="FB" & n >= 100][order(pitcher, -n)][, .SD[1], by=pitcher]
setnames(PRIM, c("pitch_type","velo","p2l","mx","mz","rx","rz","ext","tscore","tdev","n","use"),
         c("fb_type","fb_velo","fb_p2l","fb_mx","fb_mz","fb_rx","fb_rz","fb_ext",
           "fb_tscore","fb_tdev","fb_n","fb_use"))
PRIM <- PRIM[, .(pitcher, fb_type, fb_velo, fb_p2l, fb_mx, fb_mz, fb_rx, fb_rz,
                 fb_ext, fb_tscore, fb_tdev, fb_n, fb_use)]

SEC <- PT[pg2=="OFF" & use >= 5 & n >= 60 & swings >= 60]
P <- merge(SEC, PRIM, by="pitcher")
P[, `:=`(
  velo_gap = fb_velo - velo,
  p2l_gap  = abs(p2l - fb_p2l),
  p2l_sign = p2l - fb_p2l,
  mov_gap  = sqrt((mx - fb_mx)^2 + (mz - fb_mz)^2),
  rel_gap  = 12*sqrt((rx - fb_rx)^2 + (rz - fb_rz)^2),
  ext_gap  = 12*abs(ext - fb_ext)
)]
P <- P[is.finite(velo_gap) & is.finite(p2l_gap) & is.finite(tscore)]
setnames(P, c("tscore","rv","xw","whiff","tdev"),
         c("off_tscore","off_rv","off_xw","off_whiff","off_tdev"))
NAME <- unique(dt[, .(pitcher, player_name)], by="pitcher")
P <- merge(P, NAME, by="pitcher")

cat(sprintf("############ 1. Pair panel: %d fastball/secondary pairs, %d pitchers ############\n\n",
            nrow(P), uniqueN(P$pitcher)))
cat("  Distribution of the two pairing variables:\n")
print(rbindlist(list(
  data.table(v="Velocity gap (mph)", mean=mean(P$velo_gap), sd=sd(P$velo_gap),
             p10=quantile(P$velo_gap,.1), p90=quantile(P$velo_gap,.9)),
  data.table(v="Path-to-location gap", mean=mean(P$p2l_gap), sd=sd(P$p2l_gap),
             p10=quantile(P$p2l_gap,.1), p90=quantile(P$p2l_gap,.9)),
  data.table(v="Movement gap (ft)", mean=mean(P$mov_gap), sd=sd(P$mov_gap),
             p10=quantile(P$mov_gap,.1), p90=quantile(P$mov_gap,.9)),
  data.table(v="Release gap (in)", mean=mean(P$rel_gap), sd=sd(P$rel_gap),
             p10=quantile(P$rel_gap,.1), p90=quantile(P$rel_gap,.9))
))[, lapply(.SD, function(x) if (is.numeric(x)) round(x,3) else x)], row.names=FALSE)
cat("\n  Are the two pairing variables entangled?\n")
cat(sprintf("     cor(velo gap, path gap)  = %+.3f\n", cor(P$velo_gap, P$p2l_gap)))
cat(sprintf("     cor(velo gap, mov gap)   = %+.3f\n", cor(P$velo_gap, P$mov_gap)))
cat(sprintf("     cor(path gap, mov gap)   = %+.3f\n", cor(P$p2l_gap, P$mov_gap)))

cat("\n  Pairs by secondary pitch type:\n")
print(P[, .(pairs=.N, velo_gap=round(mean(velo_gap),1), p2l_gap=round(mean(p2l_gap),3),
            mov_gap=round(mean(mov_gap),2), off_tscore=round(mean(off_tscore),3),
            off_rv=round(mean(off_rv),2)), by=pitch_type][order(-pairs)], row.names=FALSE)

# =============================================================================
# 2. The 2x2: big gap + matched path vs everything else
# =============================================================================
cat("\n############ 2. The 2x2 ############\n\n")
P[, vg_hi := velo_gap >= median(velo_gap)]
P[, pg_lo := p2l_gap <= median(p2l_gap)]
P[, cell := paste0(fifelse(vg_hi,"big velo gap","small velo gap"), " + ",
                   fifelse(pg_lo,"matched path","mismatched path"))]
CELL <- P[, .(pairs=.N,
              velo_gap=mean(velo_gap), p2l_gap=mean(p2l_gap),
              off_tscore=mean(off_tscore), off_tdev=mean(off_tdev),
              off_rv=mean(off_rv), off_xw=mean(off_xw, na.rm=TRUE),
              off_whiff=mean(off_whiff)), by=cell]
setorder(CELL, -off_tscore)
print(CELL[, lapply(.SD, function(x) if (is.numeric(x)) round(x,3) else x)], row.names=FALSE)

cat("\n  Same cells, but the practical column only (pitcher run value per 100 pitches):\n")
base_rv <- P[, mean(off_rv)]
print(CELL[, .(cell, pairs, off_rv=round(off_rv,3),
               vs_avg=round(off_rv - base_rv,3))], row.names=FALSE)

# The raw cells are confounded: changeups have small velo gaps AND matched paths,
# curveballs have big gaps AND mismatched paths. Strip pitch type and usage out.
cat("\n  Cells again, now with secondary pitch type and usage residualized out.\n")
cat("  This is the comparison that matters - the raw cells above are mostly\n")
cat("  changeups versus curveballs.\n")
P[, `:=`(ts_adj = residuals(lm(off_tscore ~ factor(pitch_type) + use)),
         rv_adj = residuals(lm(off_rv ~ factor(pitch_type) + use)),
         wh_adj = residuals(lm(off_whiff ~ factor(pitch_type) + use)))]
ADJ <- P[, .(pairs=.N, ts_adj=mean(ts_adj), rv_adj=mean(rv_adj),
             wh_adj=mean(wh_adj),
             rv_adj_wt = weighted.mean(rv_adj, n)), by=cell]
setorder(ADJ, -ts_adj)
print(ADJ[, lapply(.SD, function(x) if (is.numeric(x)) round(x,3) else x)], row.names=FALSE)

cat("\n  3x3 grid on off_tscore (rows = velo gap tercile, cols = path gap tercile):\n")
q3 <- function(x) cut(x, quantile(x, 0:3/3), include.lowest=TRUE, labels=c("low","mid","high"))
P[, vg3 := q3(velo_gap)][, pg3 := q3(p2l_gap)]
print(dcast(P[, .(v=round(mean(off_tscore),3), k=.N), by=.(vg3, pg3)],
            vg3 ~ pg3, value.var="v"), row.names=FALSE)
cat("\n  ... and cell counts:\n")
print(dcast(P[, .(k=.N), by=.(vg3, pg3)], vg3 ~ pg3, value.var="k"), row.names=FALSE)
cat("\n  ... and the same grid on actual run value per 100 pitches:\n")
print(dcast(P[, .(v=round(mean(off_rv),2)), by=.(vg3, pg3)],
            vg3 ~ pg3, value.var="v"), row.names=FALSE)

# =============================================================================
# 3. Is the interaction real?
# =============================================================================
cat("\n############ 3. Regression: does path matching amplify the velocity gap? ############\n")
Z <- copy(P)
for (v in c("velo_gap","p2l_gap","mov_gap","rel_gap","ext_gap","fb_velo","use","p2l_sign"))
  Z[[v]] <- as.numeric(scale(Z[[v]]))
# sign flipped so "match" is high = more similar path, easier to read the interaction
Z[, p2l_match := -p2l_gap]

runmod <- function(y, lab, wt = FALSE) {
  f0 <- as.formula(paste(y, "~ velo_gap + p2l_match + factor(pitch_type)"))
  f1 <- as.formula(paste(y, "~ velo_gap * p2l_match + factor(pitch_type)"))
  f2 <- as.formula(paste(y, "~ velo_gap * p2l_match + mov_gap + rel_gap + ext_gap +",
                         "fb_velo + use + factor(pitch_type)"))
  w <- if (wt) Z$n else rep(1, nrow(Z))
  m0<-lm(f0,Z,weights=w); m1<-lm(f1,Z,weights=w); m2<-lm(f2,Z,weights=w)
  cat(sprintf("\n  -- %s --\n", lab))
  cat(sprintf("     main effects only            R2 = %.4f\n", summary(m0)$r.squared))
  cat(sprintf("     + interaction                R2 = %.4f  (LRT p = %.4f)\n",
              summary(m1)$r.squared, anova(m0,m1)$`Pr(>F)`[2]))
  cat(sprintf("     + movement/release controls  R2 = %.4f\n", summary(m2)$r.squared))
  co <- as.data.table(summary(m2)$coefficients, keep.rownames="term")
  setnames(co, c("term","beta","se","t","p"))
  co <- co[!grepl("factor\\(pitch_type\\)|Intercept", term)]
  print(co[, .(term, beta=round(beta,4), t=round(t,2),
               p=formatC(p, format="f", digits=4),
               sig=fifelse(p<0.001,"***",fifelse(p<0.01,"**",
                    fifelse(p<0.05,"*",fifelse(p<0.10,".","")))))], row.names=FALSE)
}
cat("\n  Standardized predictors; p2l_match = -(path gap), so positive means\n")
cat("  the two pitches arrive at their locations more similarly.\n")
runmod("off_tscore", "Modeled timing score on the secondary")
runmod("off_rv",     "Actual pitcher run value per 100 pitches (the practical test)")
runmod("off_rv",     "Same, weighted by pitch count (run value is noisy on small samples)",
       wt = TRUE)
runmod("off_whiff",  "Whiff rate")
runmod("off_xw",     "xwOBA on contact")

# Which of these outcome columns can actually adjudicate anything at pair sample
# sizes? Split each pitcher-pitch's own pitches in half and correlate.
cat("\n############ 3b. Which outcome columns are even measurable here? ############\n\n")
key <- unique(P[, .(pitcher, pitch_type)])
h <- merge(dt, key, by=c("pitcher","pitch_type"))
h[, half := sample(c(1,2), .N, replace=TRUE), by=.(pitcher, pitch_type)]
hs <- merge(sw, key, by=c("pitcher","pitch_type"))
hs[, half := sample(c(1,2), .N, replace=TRUE), by=.(pitcher, pitch_type)]

rel <- function(d, expr, minn, lab) {
  e <- d[, .(v = eval(expr), k = .N), by=.(pitcher, pitch_type, half)]
  w <- dcast(e, pitcher + pitch_type ~ half, value.var=c("v","k"))
  w <- w[k_1 >= minn & k_2 >= minn & is.finite(v_1) & is.finite(v_2)]
  r <- cor(w$v_1, w$v_2)
  data.table(outcome=lab, pairs=nrow(w), split_half=round(r,3),
             reliability=round(2*r/(1+r), 3))
}
REL <- rbindlist(list(
  rel(hs, quote(-100*mean(trv_ex)), 30, "Modeled timing score (off_tscore)"),
  rel(hs, quote(mean(tdev)), 30, "Mean timing deviation (off_tdev)"),
  rel(hs, quote(100*mean(description %in%
      c("swinging_strike","swinging_strike_blocked","foul_tip"))), 30, "Whiff rate"),
  rel(h[description=="hit_into_play" & is.finite(xw)], quote(mean(xw)), 15, "xwOBA on contact"),
  rel(h, quote(-100*mean(bat_rv)), 30, "Actual run value per 100")
))
setorder(REL, -reliability)
print(REL, row.names=FALSE)
cat("\n  Reliability is the Spearman-Brown projection to the full sample. Anything\n")
cat("  under ~0.4 cannot separate these cells no matter what the means look like.\n")

# =============================================================================
# 4. Marginal effect of velo gap at each level of path matching
# =============================================================================
cat("\n############ 4. Slope of the velocity gap, by how well the paths match ############\n\n")
P[, pgQ := cut(p2l_gap, quantile(p2l_gap, 0:4/4), include.lowest=TRUE,
               labels=c("Q1 closest match","Q2","Q3","Q4 most mismatched"))]
SL <- P[, {
  m1 <- lm(off_tscore ~ velo_gap); m2 <- lm(off_rv ~ velo_gap)
  .(pairs=.N, mean_p2l_gap=round(mean(p2l_gap),3),
    slope_tscore=round(coef(m1)[2],4), p_tscore=round(summary(m1)$coefficients[2,4],4),
    slope_rv=round(coef(m2)[2],4), p_rv=round(summary(m2)$coefficients[2,4],4))
}, by=pgQ][order(pgQ)]
print(SL, row.names=FALSE)
cat("\n  Slopes are per 1 mph of velocity gap.\n")

# =============================================================================
# 5. Does the pairing help the pitcher overall, not just that one pitch?
# =============================================================================
cat("\n############ 5. Pitcher level: best pair in the arsenal ############\n\n")
P[, `:=`(vg_z = as.numeric(scale(velo_gap)), pg_z = as.numeric(scale(p2l_gap)))]
PIT <- P[, .(n_pairs = .N,
             best_combo = max(vg_z - pg_z),   # high gap and closely matched path
             max_vg = max(velo_gap),
             min_pg = min(p2l_gap)), by=pitcher]
# a pitcher "has the combo" if any single pair is top-third velo gap AND top-third match
vg_cut <- quantile(P$velo_gap, 2/3); pg_cut <- quantile(P$p2l_gap, 1/3)
HAS <- P[, .(combo = any(velo_gap >= vg_cut & p2l_gap <= pg_cut)), by=pitcher]
TOT <- sw[, .(swings=.N, tscore=-100*mean(trv_ex),
              off_rate=-100*mean(trv_ex[pg2=="OFF"])), by=pitcher][swings >= 250]
CMP <- merge(merge(HAS, TOT, by="pitcher"), PIT, by="pitcher")
print(CMP[, .(pitchers=.N, tscore=round(mean(tscore),3),
              off_rate=round(mean(off_rate),3)), by=combo], row.names=FALSE)
cat(sprintf("\n  Welch t-test on tscore: p = %.4f\n",
    t.test(tscore ~ combo, data=CMP)$p.value))
cat(sprintf("  cor(best combo index, pitcher tscore) = %+.3f (n = %d)\n",
    cor(CMP$best_combo, CMP$tscore), nrow(CMP)))

# =============================================================================
# 6. Leaderboards
# =============================================================================
cat("\n############ 6. Pairs with big gap AND matched path ############\n\n")
SHOW <- c("player_name","fb_type","pitch_type","use","swings","velo_gap","p2l_gap",
          "mov_gap","rel_gap","off_tdev","off_tscore","off_rv","off_xw","off_whiff")
rnd <- function(d) d[, lapply(.SD, function(x) if (is.numeric(x)) round(x,3) else x)]
Q <- P[velo_gap >= vg_cut & p2l_gap <= pg_cut]
cat(sprintf("  %d pairs qualify (top third velo gap, closest third path match).\n", nrow(Q)))
cat("  Best 15 by run value:\n")
print(rnd(Q[order(-off_rv)][1:15, ..SHOW]), row.names=FALSE)
cat("\n  Worst 8 of the qualifying pairs, to show the spread within the cell:\n")
print(rnd(Q[order(off_rv)][1:8, ..SHOW]), row.names=FALSE)

fwrite(P, "data/statcast_2026/pair_velo_path_2026.csv")
cat("\nWrote data/statcast_2026/pair_velo_path_2026.csv\n")
