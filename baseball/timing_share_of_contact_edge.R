#!/usr/bin/env Rscript

# How much of each contact-manager type's edge is actually TIMING?
#
# A pitcher's hard-hit rate can be written as a weighted average of league hard-hit rates
# within timing bins, plus whatever he does to contact inside those bins:
#
#   hardhit_actual = SUM_b w_b * league_hard_b        <- "timing mix" component
#                  + (hardhit_actual - that)          <- "everything else" component
#
# The first term is what his placement of hitters on the timing axis alone implies. The
# second is stuff, location, and bat-speed suppression that survives after timing is
# accounted for. Splitting the four types on this says whether weak-contact pitchers get
# there by upsetting timing or by something timing does not see.
#
# Also reports how much of the timing component comes from the extreme tail (|tdev| >= 12),
# since Section 1 of the prior script showed displacement only pays past that point.

suppressPackageStartupMessages({ library(data.table); library(splines) })
set.seed(17)

cols <- c("game_year","game_type","player_name","pitcher","batter","stand","pitch_type",
          "description","bb_type","balls","strikes","plate_x","plate_z","sz_top","sz_bot",
          "release_speed","launch_speed","launch_angle","launch_speed_angle","delta_run_exp",
          "estimated_woba_using_speedangle","intercept_ball_minus_batter_pos_y_inches")

dt <- rbindlist(lapply(c(2025, 2026), function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress = FALSE, select = cols)))
setnames(dt, "intercept_ball_minus_batter_pos_y_inches", "depth")
dt <- dt[game_type == "R" & pitch_type != "" & balls <= 3 & strikes <= 2 &
         is.finite(delta_run_exp)]

SWING <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip","hit_into_play")
sw <- dt[description %in% SWING & is.finite(depth) & is.finite(plate_x) &
         is.finite(plate_z) & is.finite(release_speed)]
FB<-c("FF","SI","FC"); BR<-c("SL","ST","CU","KC","SV","CS"); OS<-c("CH","FS","FO")
sw[, pgrp := fifelse(pitch_type %in% FB,"FB", fifelse(pitch_type %in% BR,"BR",
             fifelse(pitch_type %in% OS,"OS",NA_character_)))]
sw <- sw[!is.na(pgrp)]
sw[, px_bat := fifelse(stand=="R", -plate_x, plate_x)]
sw[, pz_rel := (plate_z - sz_bot)/pmax(sz_top - sz_bot, 0.1)]
fit <- lm(depth ~ ns(px_bat,5)*pgrp + ns(pz_rel,5) + ns(release_speed,4) + stand, data=sw)
sw[, r1 := residuals(fit)]
sw[, nb := .N, by=batter]; sw <- sw[nb >= 200]
sw[, tdev := r1 - mean(r1), by=batter]
sw[, adev := abs(tdev)]
sw[, b_adev := mean(adev), by=batter]
sw[, excess := adev - b_adev]

bp <- sw[description == "hit_into_play" & bb_type != "" & is.finite(launch_speed) &
         is.finite(launch_angle) & is.finite(launch_speed_angle)]
bp[, hard   := launch_speed >= 95]
bp[, barrel := launch_speed_angle == 6]
bp[, tb := cut(tdev, c(-Inf,-15,-12,-9,-6,-3,0,3,6,9,12,15,Inf))]
bp[, extreme := adev >= 12]

# ---- league rates within timing bin ---------------------------------------
LG <- bp[, .(lg_hard = mean(hard), lg_brl = mean(barrel), lg_n = .N), by=tb]
bp <- merge(bp, LG, by="tb")
cat(sprintf("Balls in play: %s | league hard-hit %.1f%% | barrel %.1f%%\n\n",
  format(nrow(bp), big.mark=","), 100*mean(bp$hard), 100*mean(bp$barrel)))

# ---- what mean tdev, spread and disruption actually are -------------------
cat("############ 0. The three timing statistics, on the same swings ############\n")
cat(sprintf("  tdev is signed inches of contact depth vs the hitter's own norm.\n"))
cat(sprintf("    league mean tdev        = %+.3f in  (zero by construction, demeaned per hitter)\n",
            mean(sw$tdev)))
cat(sprintf("    league SD of tdev       = %6.2f in  (a single swing's typical displacement)\n",
            sd(sw$tdev)))
cat(sprintf("    league mean |tdev|      = %6.2f in\n", mean(sw$adev)))
PS <- sw[game_year==2026, .(swings=.N, mean_tdev=mean(tdev), sd_tdev=sd(tdev),
                            mean_abs=mean(adev), disrupt=mean(excess)), by=pitcher][swings>=250]
cat(sprintf("\n  Across %d pitchers (2026, min 250 swings), the BETWEEN-pitcher spread is much smaller:\n",
            nrow(PS)))
cat(sprintf("    mean tdev:  SD %.2f in, range %+.2f to %+.2f\n",
            sd(PS$mean_tdev), min(PS$mean_tdev), max(PS$mean_tdev)))
cat(sprintf("    disruption: SD %.2f in, range %+.2f to %+.2f\n",
            sd(PS$disrupt), min(PS$disrupt), max(PS$disrupt)))
cat(sprintf("    SD of tdev: SD %.2f in, range %.2f to %.2f\n",
            sd(PS$sd_tdev), min(PS$sd_tdev), max(PS$sd_tdev)))
cat("\n  So a pitcher's mean tdev is a small shift in the CENTER of a wide distribution,\n")
cat("  not a per-swing displacement and not a standard deviation.\n")
cat(sprintf("  cor(mean tdev, disruption) = %+.3f | cor(mean tdev, SD of tdev) = %+.3f\n",
  cor(PS$mean_tdev, PS$disrupt), cor(PS$mean_tdev, PS$sd_tdev)))

# ---- decomposition --------------------------------------------------------
cat("\n############ 1. Timing-mix vs everything-else decomposition ############\n")
P <- bp[game_year==2026, .(
  bip      = .N,
  hardhit  = 100*mean(hard),
  barrelpc = 100*mean(barrel),
  conv     = 100*sum(barrel)/pmax(sum(hard),1),
  t_hard   = 100*mean(lg_hard),      # timing mix alone
  t_brl    = 100*mean(lg_brl),
  extpct   = 100*mean(extreme),
  mean_tdev= mean(tdev),
  la       = mean(launch_angle),
  gb       = 100*mean(bb_type=="ground_ball")
), by=pitcher][bip >= 120]
# Two versions of disruption: on balls in play only, and on ALL swings (the version the
# leaderboards and the archetype figure use). The all-swings one is the relevant comparison.
D <- merge(bp[game_year==2026, .(disrupt_bip = mean(excess)), by=pitcher],
           sw[game_year==2026, .(swings = .N, disrupt = mean(excess)), by=pitcher],
           by="pitcher")
P <- merge(P, D, by="pitcher")
P[, `:=`(e_hard = hardhit - t_hard, e_brl = barrelpc - t_brl)]

cat(sprintf("  %d pitchers, 2026, min 120 balls in play.\n", nrow(P)))
cat(sprintf("  Timing-mix component of hard-hit%%: SD across pitchers = %.2f pts\n", sd(P$t_hard)))
cat(sprintf("  Everything-else component:         SD across pitchers = %.2f pts\n", sd(P$e_hard)))
cat(sprintf("  -> timing explains %.0f%% of the variance in hard-hit%% allowed\n",
  100*var(P$t_hard)/var(P$hardhit)))
cat(sprintf("  Timing-mix component of barrel%%:   SD = %.2f pts | else SD = %.2f pts -> timing %.0f%%\n",
  sd(P$t_brl), sd(P$e_brl), 100*var(P$t_brl)/var(P$barrelpc)))

hmed <- median(P$hardhit); cmed <- median(P$conv)
P[, ctype := fifelse(hardhit <= hmed & conv <= cmed, "Weak contact",
            fifelse(hardhit >  hmed & conv <= cmed, "Angle manager",
            fifelse(hardhit <= hmed & conv >  cmed, "Velocity suppressor", "Vulnerable")))]

cat("\n  By type -- hard-hit%% split into the part timing implies and the rest:\n")
cat("  (t_hard = league hard-hit rate averaged over HIS timing mix; e_hard = actual minus that)\n\n")
print(P[, .(pitchers=.N,
            hardhit=round(mean(hardhit),1),
            timing_implied=round(mean(t_hard),1),
            everything_else=round(mean(e_hard),2),
            disrupt=round(mean(disrupt),3),
            ext_pct=round(mean(extpct),1),
            mean_tdev=round(mean(mean_tdev),2)),
        by=ctype][order(hardhit)], row.names=FALSE)

cat("\n  By type -- barrel%% split the same way:\n\n")
print(P[, .(pitchers=.N,
            barrel=round(mean(barrelpc),1),
            timing_implied=round(mean(t_brl),1),
            everything_else=round(mean(e_brl),2),
            mean_tdev=round(mean(mean_tdev),2),
            la=round(mean(la),1),
            gb=round(mean(gb),1)),
        by=ctype][order(barrel)], row.names=FALSE)

# ---- is the disruption difference between types real? ---------------------
cat("\n############ 2. Is the type difference in disruption statistically real? ############\n")
for (v in c("disrupt","disrupt_bip","mean_tdev","extpct")) {
  a <- aov(as.formula(paste(v, "~ ctype")), data=P); s <- summary(a)[[1]]
  cat(sprintf("  ANOVA of %-12s by type: F = %5.2f, p = %8.3g, eta2 = %.3f | type means %s\n",
    v, s$`F value`[1], s$`Pr(>F)`[1], s$`Sum Sq`[1]/sum(s$`Sum Sq`),
    paste(sprintf("%+.2f", P[, mean(get(v)), by=ctype][order(ctype)]$V1), collapse=" ")))
}
cat(sprintf("  (type order: %s)\n", paste(sort(unique(P$ctype)), collapse=", ")))

cat("\n  Weak contact vs Vulnerable, the two diagonal corners:\n")
for (v in c("disrupt","mean_tdev","extpct","hardhit","barrelpc","t_hard","e_hard")) {
  x <- P[ctype=="Weak contact"][[v]]; y <- P[ctype=="Vulnerable"][[v]]
  tt <- t.test(x, y)
  cat(sprintf("     %-10s weak %+8.3f | vuln %+8.3f | diff %+7.3f | p = %.3g\n",
    v, mean(x), mean(y), mean(x)-mean(y), tt$p.value))
}
cat("\n  Weak contact vs Velocity suppressor (same hard-hit tier, opposite conversion):\n")
for (v in c("disrupt","mean_tdev","extpct","hardhit","conv","t_hard","e_hard")) {
  x <- P[ctype=="Weak contact"][[v]]; y <- P[ctype=="Velocity suppressor"][[v]]
  tt <- t.test(x, y)
  cat(sprintf("     %-10s weak %+8.3f | velo %+8.3f | diff %+7.3f | p = %.3g\n",
    v, mean(x), mean(y), mean(x)-mean(y), tt$p.value))
}

# ---- what is the disruption gap worth? -----------------------------------
cat("\n############ 3. What is the disruption gap actually worth? ############\n")
mh <- lm(hardhit ~ disrupt, data=P); mb <- lm(barrelpc ~ disrupt, data=P)
cat(sprintf("  hard-hit%% per inch of disruption: %+.2f pts (SE %.2f)\n",
  coef(mh)[2], summary(mh)$coefficients[2,2]))
cat(sprintf("  barrel%%   per inch of disruption: %+.2f pts (SE %.2f)\n",
  coef(mb)[2], summary(mb)$coefficients[2,2]))
gap <- mean(P[ctype=="Weak contact"]$disrupt) - mean(P[ctype=="Vulnerable"]$disrupt)
cat(sprintf("\n  Weak-contact minus vulnerable disruption gap = %.3f in\n", gap))
cat(sprintf("  At the fitted slope that is worth %.2f pts of hard-hit%% and %.2f pts of barrel%%,\n",
  abs(gap*coef(mh)[2]), abs(gap*coef(mb)[2])))
cat(sprintf("  against an ACTUAL gap of %.1f pts hard-hit and %.1f pts barrel.\n",
  mean(P[ctype=="Weak contact"]$hardhit) - mean(P[ctype=="Vulnerable"]$hardhit),
  mean(P[ctype=="Weak contact"]$barrelpc) - mean(P[ctype=="Vulnerable"]$barrelpc)))
cat(sprintf("  -> disruption accounts for %.0f%% of the hard-hit gap and %.0f%% of the barrel gap.\n",
  100*abs(gap*coef(mh)[2])/abs(mean(P[ctype=="Weak contact"]$hardhit) -
                               mean(P[ctype=="Vulnerable"]$hardhit)),
  100*abs(gap*coef(mb)[2])/abs(mean(P[ctype=="Weak contact"]$barrelpc) -
                               mean(P[ctype=="Vulnerable"]$barrelpc))))

cat("\n############ 4. Do weak-contact pitchers suppress WELL-TIMED contact too? ############\n")
cat("If their edge were really timing, they should look ordinary on swings the hitter timed well.\n\n")
W <- merge(bp[game_year==2026 & adev < 6, .(bip_wt=.N, hard_wt=100*mean(hard),
             brl_wt=100*mean(barrel)), by=pitcher],
           bp[game_year==2026 & adev >= 12, .(bip_ex=.N, hard_ex=100*mean(hard),
             brl_ex=100*mean(barrel)), by=pitcher], by="pitcher")
W <- merge(P[, .(pitcher, ctype, hardhit, barrelpc, disrupt)], W, by="pitcher")
W <- W[bip_wt >= 40]
cat(sprintf("  %d pitchers with >=40 well-timed balls in play (|tdev| < 6 in)\n\n", nrow(W)))
print(W[, .(pitchers=.N,
            hard_all=round(mean(hardhit),1),
            hard_welltimed=round(mean(hard_wt),1),
            brl_welltimed=round(mean(brl_wt),1),
            hard_extreme=round(mean(hard_ex, na.rm=TRUE),1)),
       by=ctype][order(hard_all)], row.names=FALSE)
cat(sprintf("\n  cor(hard-hit%% overall, hard-hit%% on well-timed contact) = %+.3f\n",
            cor(W$hardhit, W$hard_wt)))
cat(sprintf("  cor(disruption, hard-hit%% on well-timed contact)        = %+.3f\n",
            cor(W$disrupt, W$hard_wt)))
cat("\n  Weak contact vs Vulnerable on well-timed contact only:\n")
tt <- t.test(W[ctype=="Weak contact"]$hard_wt, W[ctype=="Vulnerable"]$hard_wt)
cat(sprintf("     %.1f%% vs %.1f%% (diff %+.1f pts, p = %.3g)\n",
  mean(W[ctype=="Weak contact"]$hard_wt), mean(W[ctype=="Vulnerable"]$hard_wt),
  mean(W[ctype=="Weak contact"]$hard_wt) - mean(W[ctype=="Vulnerable"]$hard_wt), tt$p.value))
