#!/usr/bin/env Rscript

# Does the velo-gap benefit on sliders (and breaking balls) only exist for
# hard throwers? 79 mph sliders as the test case.

suppressPackageStartupMessages({ library(data.table); library(splines) })
set.seed(7)
options(width = 210)

P <- fread("data/statcast_2026/velo_gap_tunnel_pairs.csv")

cols <- c("game_year","game_type","pitcher","pitch_type","description","balls","strikes",
          "release_speed","estimated_woba_using_speedangle","delta_run_exp")
dt <- rbindlist(lapply(c(2025, 2026), function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress=FALSE, select=cols)))
setnames(dt, "estimated_woba_using_speedangle", "xw")
dt <- dt[game_type=="R" & pitch_type!="" & balls<=3 & strikes<=2 & is.finite(delta_run_exp)]
dt[, bat_rv := delta_run_exp]
mb <- lm(bat_rv ~ ns(xw,5)+factor(paste0(balls,"-",strikes)),
         data=dt[description=="hit_into_play" & is.finite(xw)])
dt[, xrv_p := bat_rv]
dt[description=="hit_into_play" & is.finite(xw),
   xrv_p := predict(mb, newdata=dt[description=="hit_into_play" & is.finite(xw)])]
dt <- dt[!(description=="hit_into_play" & !is.finite(xw))]

# pair-level slider / sweeper velocity (the secondary itself)
SEC <- dt[, .(n=.N, sl_velo=mean(release_speed),
              xrv=-100*mean(xrv_p), rv=-100*mean(bat_rv)),
          by=.(pitcher, pitch_type)]
P <- merge(P, SEC[, .(pitcher, pitch_type, sl_velo, n_sec=n, xrv2=xrv)],
           by=c("pitcher","pitch_type"))
# prefer the already-computed xrv if present
if (!"xrv" %in% names(P) || all(is.na(P$xrv))) P[, xrv := xrv2]
P[, xrv := fifelse(is.finite(xrv), xrv, xrv2)]

BR <- c("SL","ST")
S <- P[pitch_type %in% BR]
cat(sprintf("SL/ST pairs: %d  (SL %d, ST %d)\n", nrow(S), nrow(S[pitch_type=="SL"]),
            nrow(S[pitch_type=="ST"])))
cat(sprintf("  Slider velo: mean %.1f  p10 %.1f  p50 %.1f  p90 %.1f\n",
            mean(S[pitch_type=="SL"]$sl_velo),
            quantile(S[pitch_type=="SL"]$sl_velo, .1),
            quantile(S[pitch_type=="SL"]$sl_velo, .5),
            quantile(S[pitch_type=="SL"]$sl_velo, .9)))
cat(sprintf("  Sweeper velo: mean %.1f  p10 %.1f  p50 %.1f  p90 %.1f\n\n",
            mean(S[pitch_type=="ST"]$sl_velo),
            quantile(S[pitch_type=="ST"]$sl_velo, .1),
            quantile(S[pitch_type=="ST"]$sl_velo, .5),
            quantile(S[pitch_type=="ST"]$sl_velo, .9)))

# =============================================================================
# 1. Raw: do slow sliders perform poorly?
# =============================================================================
cat("############ 1. Slider / sweeper xRV by the pitch's own velocity ############\n\n")
band <- function(x) cut(x, c(0, 79, 82, 85, 88, 100),
                        labels=c("<79","79-82","82-85","85-88","88+"),
                        right=FALSE)
S[, vband := band(sl_velo)]
print(S[, .(pairs=.N,
            sl_velo=round(mean(sl_velo),1),
            fb_velo=round(mean(fb_velo),1),
            velo_gap=round(mean(velo_gap),1),
            xrv=round(mean(xrv),3),
            xrv_wt=round(weighted.mean(xrv, n_sec),3),
            whiff=round(mean(whiff),1)),
        by=.(pitch_type, vband)][order(pitch_type, vband)], row.names=FALSE)

cat("\n  SL only, continuous: xRV on slider velo (no gap in the model yet):\n")
m <- lm(xrv ~ sl_velo + n_sec, data=S[pitch_type=="SL"])
s <- summary(m)$coefficients["sl_velo",]
cat(sprintf("     β = %+.4f / mph   p = %.3f   R2 = %.3f   n = %d\n",
            s[1], s[4], summary(m)$r.squared, nobs(m)))

# =============================================================================
# 2. Does the gap slope depend on how hard the slider is?
# =============================================================================
cat("\n############ 2. Gap slope inside each slider-velocity band ############\n\n")
cat("  SL only. Coefficient is xRV per 1 mph of fastball-minus-slider gap.\n\n")
SL <- S[pitch_type=="SL"]
for (b in levels(SL$vband)) {
  d <- SL[vband==b]
  if (nrow(d) < 20) {
    cat(sprintf("     %-6s  n = %2d   skipped\n", b, nrow(d)))
    next
  }
  m <- lm(xrv ~ velo_gap + n_sec, data=d)
  s <- summary(m)$coefficients["velo_gap",]
  cat(sprintf("     %-6s  n = %3d   mean SL %.1f mph   mean gap %.1f   β = %+.4f / mph   p = %.3f   mean xRV = %+.3f\n",
              b, nrow(d), mean(d$sl_velo), mean(d$velo_gap), s[1], s[4], mean(d$xrv)))
}

cat("\n  Interaction test, SL: xRV ~ velo_gap * sl_velo\n")
mI <- lm(xrv ~ velo_gap * sl_velo + n_sec, data=SL)
print(round(summary(mI)$coefficients, 4))
cat(sprintf("  R2 = %.3f\n", summary(mI)$r.squared))

cat("\n  Same interaction, ST:\n")
ST <- S[pitch_type=="ST"]
mS <- lm(xrv ~ velo_gap * sl_velo + n_sec, data=ST)
print(round(summary(mS)$coefficients, 4))

cat("\n  SL+ST pooled, pitch-type FE:\n")
mP <- lm(xrv ~ velo_gap * sl_velo + factor(pitch_type) + n_sec, data=S)
print(round(summary(mP)$coefficients, 4))

# =============================================================================
# 3. Hard thrower vs slow thrower — defined on the FASTBALL
# =============================================================================
cat("\n############ 3. Is it the fastball that's required to be hard? ############\n\n")
SL[, fbQ := cut(fb_velo, quantile(fb_velo, 0:3/3), include.lowest=TRUE,
                labels=c("softest FB third","mid FB","hardest FB third"))]
cat("  Gap slope inside each fastball-velocity tercile (SL only):\n")
for (b in levels(SL$fbQ)) {
  d <- SL[fbQ==b]
  m <- lm(xrv ~ velo_gap + n_sec, data=d)
  s <- summary(m)$coefficients["velo_gap",]
  cat(sprintf("     %-18s  n = %3d   FB %.1f   SL %.1f   gap %.1f   β = %+.4f   p = %.3f   xRV = %+.3f\n",
              b, nrow(d), mean(d$fb_velo), mean(d$sl_velo), mean(d$velo_gap),
              s[1], s[4], mean(d$xrv)))
}
mF <- lm(xrv ~ velo_gap * fb_velo + n_sec, data=SL)
cat("\n  Interaction velo_gap × FB velo (SL):\n")
print(round(summary(mF)$coefficients, 4))

# =============================================================================
# 4. The 79 mph slice, directly
# =============================================================================
cat("\n############ 4. Sliders at 77-81 mph vs 85+ ############\n\n")
SL[, grp := fifelse(sl_velo < 81, "SL under 81",
             fifelse(sl_velo >= 85, "SL 85+", "SL 81-85"))]
print(SL[, .(pairs=.N, sl_velo=round(mean(sl_velo),1), fb_velo=round(mean(fb_velo),1),
             gap=round(mean(velo_gap),1),
             xrv=round(mean(xrv),3),
             xrv_wt=round(weighted.mean(xrv, n_sec),3),
             r_gap=round(cor(velo_gap, xrv),3)),
         by=grp][order(grp)], row.names=FALSE)

cat("\n  Gap slope in each slice:\n")
for (g in c("SL under 81","SL 81-85","SL 85+")) {
  d <- SL[grp==g]
  m <- lm(xrv ~ velo_gap + n_sec, data=d)
  s <- summary(m)$coefficients["velo_gap",]
  cat(sprintf("     %-14s  n = %3d   β = %+.4f / mph   p = %.3f   10th-90th gap %.1f-%.1f\n",
              g, nrow(d), s[1], s[4], quantile(d$velo_gap,.1), quantile(d$velo_gap,.9)))
}

cat("\n  Among sliders under 81: does a bigger gap rescue them?\n")
slow <- SL[sl_velo < 81]
slow[, g3 := cut(velo_gap, quantile(velo_gap, 0:3/3), include.lowest=TRUE,
                 labels=c("small gap","mid","big gap"))]
print(slow[, .(pairs=.N, sl_velo=round(mean(sl_velo),1), gap=round(mean(velo_gap),1),
               fb_velo=round(mean(fb_velo),1), xrv=round(mean(xrv),3),
               xrv_wt=round(weighted.mean(xrv, n_sec),3)),
           by=g3][order(g3)], row.names=FALSE)

cat("\n  Among sliders 85+: same terciles of gap\n")
hard <- SL[sl_velo >= 85]
hard[, g3 := cut(velo_gap, quantile(velo_gap, 0:3/3), include.lowest=TRUE,
                 labels=c("small gap","mid","big gap"))]
print(hard[, .(pairs=.N, sl_velo=round(mean(sl_velo),1), gap=round(mean(velo_gap),1),
               fb_velo=round(mean(fb_velo),1), xrv=round(mean(xrv),3),
               xrv_wt=round(weighted.mean(xrv, n_sec),3)),
           by=g3][order(g3)], row.names=FALSE)

# =============================================================================
# 5. Names: slow sliders with a big gap vs slow sliders with a small gap
# =============================================================================
cat("\n############ 5. Slow sliders, biggest and smallest gaps ############\n\n")
SHOW <- c("player_name","fb_type","sl_velo","fb_velo","velo_gap","n_sec","xrv","whiff")
cat("  SL under 81, largest 8 gaps:\n")
print(slow[order(-velo_gap)][1:8, .(player_name, fb_type,
  sl=round(sl_velo,1), fb=round(fb_velo,1), gap=round(velo_gap,1),
  n=n_sec, xrv=round(xrv,2), whiff=round(whiff,1))], row.names=FALSE)
cat("\n  SL under 81, smallest 8 gaps:\n")
print(slow[order(velo_gap)][1:8, .(player_name, fb_type,
  sl=round(sl_velo,1), fb=round(fb_velo,1), gap=round(velo_gap,1),
  n=n_sec, xrv=round(xrv,2), whiff=round(whiff,1))], row.names=FALSE)

# =============================================================================
# 6. Partial: gap holding slider velo fixed, and slider velo holding gap fixed
# =============================================================================
cat("\n############ 6. Which number actually matters once both are in the model? ############\n\n")
mBoth <- lm(xrv ~ velo_gap + sl_velo + fb_velo + n_sec, data=SL)
cat("  SL: xRV ~ gap + slider velo + FB velo\n")
print(round(summary(mBoth)$coefficients, 4))
cat(sprintf("  R2 = %.3f   adj = %.3f   n = %d\n",
            summary(mBoth)$r.squared, summary(mBoth)$adj.r.squared, nobs(mBoth)))
cat(sprintf("  VIF-ish: cor(gap, SL velo) = %+.3f   cor(gap, FB velo) = %+.3f   cor(SL, FB) = %+.3f\n",
            cor(SL$velo_gap, SL$sl_velo), cor(SL$velo_gap, SL$fb_velo),
            cor(SL$sl_velo, SL$fb_velo)))
