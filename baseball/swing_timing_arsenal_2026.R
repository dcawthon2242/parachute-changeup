#!/usr/bin/env Rscript

# 2026 pitcher x pitch type: who owns the best run value in each count state, and do the
# leaders get there through the SOFT CONTACT objective or the SWING-AND-MISS objective?
#
# Objectives, as established by swing_timing_by_count.R:
#   SOFT CONTACT   = drive the hitter LATE (tdev <= -6in), minimise xwOBAcon  -> pays at 0 strikes
#   SWING-AND-MISS = drive the hitter EARLY, maximise whiff rate              -> pays at 2 strikes
#
# Two measurement problems are handled explicitly rather than papered over:
#   (1) Run value per swing at pitcher x pitch-type x count grain is mostly sampling noise.
#       Section 0 measures how much true spread survives, and leaderboards are shrunk.
#   (2) xwOBAcon and whiff rate are mechanically components of run value in the same sample,
#       so a same-sample correlation is a decomposition, not a forecast. Section 3b splits
#       each arsenal in half and predicts one half's run value from the other half's objectives.
#
# The timing axis is built on 2025+2026 pooled, then subset to 2026.

suppressPackageStartupMessages({ library(data.table); library(splines) })
set.seed(11)

cols <- c("game_year","game_type","player_name","pitcher","batter","stand","pitch_type",
          "description","bb_type","balls","strikes","plate_x","plate_z","sz_top","sz_bot",
          "release_speed","launch_speed","estimated_woba_using_speedangle","delta_run_exp",
          "intercept_ball_minus_batter_pos_y_inches")

dt <- rbindlist(lapply(c(2025, 2026), function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress = FALSE, select = cols)))
setnames(dt, "intercept_ball_minus_batter_pos_y_inches", "depth")
dt <- dt[game_type == "R" & pitch_type != "" & balls <= 3 & strikes <= 2 & is.finite(delta_run_exp)]
dt[, pit_rv := -delta_run_exp]
dt[, sgrp := paste0(strikes, fifelse(strikes == 1, " strike", " strikes"))]

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
sw[, whiff := description %in% c("swinging_strike","swinging_strike_blocked","foul_tip")]
sw[, bipf  := description == "hit_into_play" & bb_type != "" & is.finite(launch_speed)]

s26 <- sw[game_year == 2026]
p26 <- dt[game_year == 2026]
s26[, half := sample(rep_len(1:2, .N)), by=.(pitcher, pitch_type, sgrp)]
cat(sprintf("2026: %s pitches | %s swings on the timing axis | %d pitchers\n\n",
  format(nrow(p26), big.mark=","), format(nrow(s26), big.mark=","), uniqueN(s26$pitcher)))

# ---- aggregation ----------------------------------------------------------
agg <- function(swd, pd, by) {
  a <- swd[, .(
    swings=.N, rv_swing=mean(pit_rv), se=sd(pit_rv)/sqrt(.N),
    whiff=100*mean(whiff), se_whiff=100*sqrt(mean(whiff)*(1-mean(whiff))/.N),
    bip=sum(bipf),
    xwobacon=mean(estimated_woba_using_speedangle[bipf], na.rm=TRUE),
    se_xw=sd(estimated_woba_using_speedangle[bipf], na.rm=TRUE)/sqrt(pmax(sum(bipf),1)),
    tdev=mean(tdev), late6=100*mean(tdev <= -6), late12=100*mean(tdev <= -12),
    early12=100*mean(tdev >= 12)
  ), by = by]
  b <- pd[, .(pitches=.N, rv100=100*mean(pit_rv)), by = by]
  merge(a, b, by = by)
}
BYS <- c("pitcher","player_name","pitch_type","sgrp")
Araw <- agg(s26, p26, BYS)
Oraw <- agg(s26, p26, c("pitcher","player_name","pitch_type"))

# =============================================================================
# 0. Signal-to-noise: is there real arm-level spread to rank at all?
# =============================================================================
audit1 <- function(x, val, sev, lab) {
  x <- x[is.finite(get(sev)) & is.finite(get(val))]
  obs <- var(x[[val]]); noi <- mean(x[[sev]]^2); tru <- obs - noi
  cat(sprintf("    %-12s observed SD %7.4f | noise SD %7.4f | true SD %-12s | reliability %s\n",
    lab, sqrt(obs), sqrt(noi),
    if (tru > 0) sprintf("%.4f", sqrt(tru)) else "~0",
    if (tru > 0) sprintf("%.2f", tru/obs) else "0.00"))
}
cat("############ 0. How much of the arm-to-arm spread is real, not sampling noise? ############\n")
cat("  -- all counts pooled, min 120 swings --\n")
x <- Oraw[swings >= 120]
audit1(x, "rv_swing", "se", "run value"); audit1(x, "whiff", "se_whiff", "whiff %")
audit1(x, "xwobacon", "se_xw", "xwOBAcon")
for (g in c("0 strikes","1 strike","2 strikes")) {
  cat(sprintf("  -- %s, min 80 swings --\n", g))
  x <- Araw[sgrp == g & swings >= 80]
  audit1(x, "rv_swing", "se", "run value"); audit1(x, "whiff", "se_whiff", "whiff %")
  audit1(x, "xwobacon", "se_xw", "xwOBAcon")
}

# ---- empirical-Bayes shrinkage within each qualifying population -----------
shrink <- function(d, key = NULL) {
  d <- copy(d); d[, .g := if (is.null(key)) "all" else as.character(get(key))]
  d[!is.finite(se), se := NA_real_]
  d[, mu := mean(rv_swing[is.finite(se)], na.rm=TRUE), by=.g]
  d[, tau2 := pmax(var(rv_swing[is.finite(se)], na.rm=TRUE) -
                   mean(se[is.finite(se)]^2, na.rm=TRUE), 1e-8), by=.g]
  d[, rel := fifelse(is.finite(se), tau2/(tau2 + se^2), 0)]
  d[, rv_adj := mu + rel*(rv_swing - mu)]
  d[, `:=`(mu=NULL, tau2=NULL, .g=NULL)][]
}
O   <- shrink(Oraw[swings >= 120])
A   <- shrink(Araw[swings >= 80], "sgrp")
A40 <- shrink(Araw[swings >= 40], "sgrp")

lab_type <- function(t) fifelse(t <= -2, "late-inducing",
                       fifelse(t >=  2, "early-inducing", "neutral"))
for (d in list(O, A, A40)) d[, otype := lab_type(tdev)]

fmt <- function(d, cs) { d <- copy(d)
  for (v in c("rv_adj","rv_swing","se")) d[[v]] <- round(d[[v]], 4)
  for (v in c("xwobacon")) d[[v]] <- round(d[[v]], 3)
  for (v in c("rel")) d[[v]] <- round(d[[v]], 2)
  for (v in c("rv100","whiff","tdev","late6","early12")) d[[v]] <- round(d[[v]], 1)
  d[, ..cs] }
SHOW <- c("player_name","pitch_type","swings","rv_adj","rv_swing","rel","rv100",
          "whiff","xwobacon","tdev","late6","early12","otype")

# =============================================================================
# 1 / 2. Leaderboards
# =============================================================================
cat("\n############ 1. Best 2026 pitches, all counts, shrunk run value per swing (min 120 swings) ############\n")
cat(sprintf("  %d qualify | median reliability %.2f\n\n", nrow(O), median(O$rel)))
print(fmt(O[order(-rv_adj)][1:25], SHOW))
cat("\n  -- bottom 8 --\n"); print(fmt(O[order(rv_adj)][1:8], SHOW))

for (g in c("0 strikes","1 strike","2 strikes")) {
  d <- A[sgrp == g][order(-rv_adj)]
  cat(sprintf("\n############ 2. %s: top 15, shrunk run value per swing (min 80 swings) ############\n", g))
  cat(sprintf("   optimal objective per the count analysis: %s\n", switch(g,
    "0 strikes"="SOFT CONTACT (drive hitters late)",
    "1 strike" ="either -- statistical dead heat",
    "2 strikes"="SWING-AND-MISS (drive hitters early)")))
  cat(sprintf("   %d qualify | median reliability %.2f\n\n", nrow(d), median(d$rel)))
  print(fmt(d[1:15], SHOW))
}

cat("\n############ 2b. Top 18 arsenals overall: shrunk run value by count state ############\n")
top <- O[swings >= 200][order(-rv_adj)][1:18,
  .(pitcher, player_name, pitch_type, rv_all=rv_adj, whiff_all=whiff,
    xw_all=xwobacon, tdev_all=tdev, late6_all=late6, early12_all=early12)]
sp <- dcast(A40, pitcher + pitch_type ~ sgrp, value.var=c("rv_adj","whiff","xwobacon"))
setnames(sp, gsub(" ", "", names(sp)))
m <- merge(top, sp, by=c("pitcher","pitch_type"), all.x=TRUE)[order(-rv_all)]
print(m[, .(player_name, pitch_type,
  rv_all=round(rv_all,4), whiff=round(whiff_all,1), xwobacon=round(xw_all,3),
  tdev=round(tdev_all,1), late6=round(late6_all,1), early12=round(early12_all,1),
  rv_0str=round(rv_adj_0strikes,4), rv_1str=round(rv_adj_1strike,4),
  rv_2str=round(rv_adj_2strikes,4),
  wh_0str=round(whiff_0strikes,1), wh_2str=round(whiff_2strikes,1),
  xw_0str=round(xwobacon_0strikes,3), xw_2str=round(xwobacon_2strikes,3))])

# =============================================================================
# 3a. Where the run value came from (same-sample decomposition, NOT a forecast)
# =============================================================================
cat("\n############ 3a. Same-sample decomposition: which channel carries the run value? ############\n")
cat("  NOTE: xwOBAcon and whiff are mechanically inside run value here, so these are\n")
cat("  attributions, not predictions. The forecast test is 3b.\n\n")
for (g in c("0 strikes","1 strike","2 strikes")) {
  d <- A40[sgrp == g & swings >= 60 & is.finite(xwobacon)]
  f <- lm(scale(rv_swing) ~ scale(xwobacon) + scale(whiff) + scale(late6), data=d)
  co <- summary(f)$coefficients
  cat(sprintf("  %-10s n=%3d R2=%.3f | xwOBAcon %+.3f (p=%.1g) | whiff %+.3f (p=%.1g) | late6%% %+.3f (p=%.1g)\n",
    g, nrow(d), summary(f)$r.squared, co[2,1], co[2,4], co[3,1], co[3,4], co[4,1], co[4,4]))
}

# =============================================================================
# 3b. Split-half forecast: do the objectives in one half predict RV in the other?
# =============================================================================
cat("\n############ 3b. Split-half forecast (objectives in half A -> run value in half B) ############\n")
h <- function(hh, by) s26[half == hh, .(
  n=.N, rv=mean(pit_rv), whiff=100*mean(whiff),
  xw=mean(estimated_woba_using_speedangle[bipf], na.rm=TRUE),
  late6=100*mean(tdev <= -6), early12=100*mean(tdev >= 12)), by=by]
pc <- function(x,y){ok<-is.finite(x)&is.finite(y); if(sum(ok)<25) return(c(NA,NA,sum(ok)))
  ct<-cor.test(x[ok],y[ok]); c(ct$estimate, ct$p.value, sum(ok))}
sh_report <- function(by, minn, lab) {
  A1 <- h(1, by); B1 <- h(2, by)
  mm <- merge(A1[n >= minn], B1[n >= minn], by=by, suffixes=c("_a","_b"))
  if (!nrow(mm)) return(invisible())
  cat(sprintf("\n  -- %s (min %d swings per half, n=%d) --\n", lab, minn, nrow(mm)))
  for (v in c("xw","whiff","late6","early12")) {
    r1 <- pc(mm[[paste0(v,"_a")]], mm$rv_b)
    r2 <- pc(mm[[paste0(v,"_a")]], mm[[paste0(v,"_b")]])
    cat(sprintf("     %-8s -> other-half run value r=%+.3f (p=%.2g) | own split-half stability r=%+.3f\n",
      v, r1[1], r1[2], r2[1]))
  }
}
sh_report(c("pitcher","pitch_type"), 80, "all counts pooled")
for (g in c("0 strikes","1 strike","2 strikes")) {
  A1 <- h(1, c("pitcher","pitch_type","sgrp"))[sgrp==g]
  B1 <- h(2, c("pitcher","pitch_type","sgrp"))[sgrp==g]
  mm <- merge(A1[n>=40], B1[n>=40], by=c("pitcher","pitch_type","sgrp"), suffixes=c("_a","_b"))
  cat(sprintf("\n  -- %s (min 40 swings per half, n=%d) --\n", g, nrow(mm)))
  for (v in c("xw","whiff","late6","early12")) {
    r1 <- pc(mm[[paste0(v,"_a")]], mm$rv_b)
    r2 <- pc(mm[[paste0(v,"_a")]], mm[[paste0(v,"_b")]])
    cat(sprintf("     %-8s -> other-half run value r=%+.3f (p=%.2g) | own split-half stability r=%+.3f\n",
      v, r1[1], r1[2], r2[1]))
  }
}

# =============================================================================
# 4. Archetypes
# =============================================================================
cat("\n############ 4. Late-inducing vs early-inducing arsenals, by count state ############\n")
cat("(classified on mean tdev: <=-2in late-inducing, >=+2in early-inducing; min 60 swings)\n")
for (g in c("0 strikes","1 strike","2 strikes")) {
  d <- A40[sgrp == g & swings >= 60]
  cat(sprintf("\n  -- %s --\n", g))
  print(d[, .(arsenals=.N, swings=sum(swings), tdev=round(mean(tdev),1),
    whiff=round(mean(whiff),1), xwobacon=round(mean(xwobacon,na.rm=TRUE),3),
    rv_swing=round(mean(rv_swing),4), rv100=round(mean(rv100),2)),
    by=otype][order(-rv_swing)])
}

cat("\n############ 5. Arsenals top-quartile on BOTH objectives (all counts, min 120 swings) ############\n")
lc <- quantile(O$late6, 0.75); wc <- quantile(O$whiff, 0.75)
both <- O[late6 >= lc & whiff >= wc][order(-rv_adj)]
cat(sprintf("  thresholds: late6 >= %.1f%%, whiff >= %.1f%% | %d of %d arsenals clear both\n\n",
            lc, wc, nrow(both), nrow(O)))
print(fmt(both[1:min(20,nrow(both))], SHOW))

fwrite(O[order(-rv_adj)], file.path("data","statcast_2026","timing_arsenal_2026_all.csv"))
fwrite(A40[order(sgrp,-rv_adj)], file.path("data","statcast_2026","timing_arsenal_2026_bycount.csv"))
cat("\nWrote data/statcast_2026/timing_arsenal_2026_all.csv and _bycount.csv\n")
