#!/usr/bin/env Rscript

# The user's question:
#   98 FB + 90 SL vs 98 FB + 86 SL, same movement.
#   Why does the 90 look better even though the 86 has the bigger velo gap?
# And: do high-gap sliders lose run value by not drawing chase?

suppressPackageStartupMessages({ library(data.table); library(splines) })
set.seed(7)
options(width = 215)

cols <- c("game_year","game_type","pitcher","player_name","p_throws","stand",
          "pitch_type","description","balls","strikes",
          "release_speed","pfx_x","pfx_z","release_spin_rate","release_extension",
          "plate_x","plate_z","sz_top","sz_bot",
          "launch_speed","launch_speed_angle",
          "estimated_woba_using_speedangle","delta_run_exp")
dt <- rbindlist(lapply(c(2025, 2026), function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress=FALSE, select=cols)))
setnames(dt, "estimated_woba_using_speedangle", "xw")
dt <- dt[game_type=="R" & pitch_type!="" & balls<=3 & strikes<=2 &
         is.finite(delta_run_exp)]
dt[, bat_rv := delta_run_exp]
dt[, in_zone := is.finite(plate_x) & abs(plate_x)<=0.83 &
                is.finite(plate_z) & plate_z>=sz_bot & plate_z<=sz_top]
WH <- c("swinging_strike","swinging_strike_blocked","foul_tip")
SW <- c(WH,"foul","hit_into_play")
dt[, swung := description %in% SW]
dt[, cnt := paste0(balls,"-",strikes)]
mb <- lm(bat_rv ~ ns(xw,5)+factor(cnt),
         data=dt[description=="hit_into_play" & is.finite(xw)])
dt[, xrv_p := bat_rv]
dt[description=="hit_into_play" & is.finite(xw),
   xrv_p := predict(mb, newdata=dt[description=="hit_into_play" & is.finite(xw)])]
dt <- dt[!(description=="hit_into_play" & !is.finite(xw))]

# glove-side horizontal (positive = sweep away from RHH for RHP, etc.)
dt[, pfx_glove := fifelse(p_throws=="R", -pfx_x, pfx_x)]

FB <- dt[pitch_type %in% c("FF","SI","FC") & is.finite(release_speed)]
# primary FB = most-thrown of FF/SI/FC
PRI <- FB[, .N, by=.(pitcher, pitch_type)][order(-N), .SD[1], by=pitcher]
FBm <- FB[, .(fb_n=.N, fb_velo=mean(release_speed)), by=pitcher]
FBm <- merge(FBm, PRI[, .(pitcher, fb_type=pitch_type)], by="pitcher")
# FB velo of the primary type only
FBp <- FB[PRI, on=.(pitcher, pitch_type),
          .(fb_pri=mean(release_speed), fb_pri_n=.N), by=pitcher]
FBm <- merge(FBm, FBp, by="pitcher")

sl <- dt[pitch_type=="SL"]
# pitcher-level slider
SL <- sl[, .(
  n=.N,
  name=player_name[1],
  sl_velo=mean(release_speed),
  hmov=mean(pfx_glove, na.rm=TRUE),
  vmov=mean(pfx_z, na.rm=TRUE),
  spin=mean(release_spin_rate, na.rm=TRUE),
  ext=mean(release_extension, na.rm=TRUE),
  xrv=-100*mean(xrv_p),
  rv=-100*mean(bat_rv),
  zone=100*mean(in_zone, na.rm=TRUE),
  swing=100*mean(swung),
  swing_z=100*mean(swung[in_zone]),
  chase=100*mean(swung[!in_zone]),
  whiff=100*mean(description %in% WH),
  wswing=100*mean((description %in% WH)[swung]),
  cs=100*mean(description=="called_strike"),
  ball=100*mean(description=="ball"),
  hard=100*mean(launch_speed[description=="hit_into_play"]>=95, na.rm=TRUE),
  brl=100*mean(launch_speed_angle[description=="hit_into_play"]==6, na.rm=TRUE),
  xwcon=mean(xw[description=="hit_into_play"], na.rm=TRUE),
  # xRV channels (add to total)
  ch_whiff=-100*mean(fifelse(description %in% WH, xrv_p, 0)),
  ch_bip=-100*mean(fifelse(description=="hit_into_play", xrv_p, 0)),
  ch_foul=-100*mean(fifelse(description=="foul", xrv_p, 0)),
  ch_cs=-100*mean(fifelse(description=="called_strike", xrv_p, 0)),
  ch_ball=-100*mean(fifelse(description=="ball", xrv_p, 0)),
  ch_chase=-100*mean(fifelse(!in_zone & swung, xrv_p, 0)),
  ch_ooz_take=-100*mean(fifelse(!in_zone & !swung, xrv_p, 0)),
  ch_iz_swing=-100*mean(fifelse(in_zone & swung, xrv_p, 0)),
  ch_iz_take=-100*mean(fifelse(in_zone & !swung, xrv_p, 0))
), by=pitcher]
SL <- merge(SL, FBm, by="pitcher")
SL[, gap := fb_pri - sl_velo]
SL <- SL[n>=80 & fb_pri_n>=80 & is.finite(gap) & is.finite(hmov)]
cat(sprintf("SL pitchers, 80+ SL and 80+ primary FB: %d\n", nrow(SL)))
cat(sprintf("  SL velo mean %.1f  p10 %.1f  p50 %.1f  p90 %.1f\n",
            mean(SL$sl_velo), quantile(SL$sl_velo,.1),
            quantile(SL$sl_velo,.5), quantile(SL$sl_velo,.9)))
cat(sprintf("  FB pri  mean %.1f  p10 %.1f  p50 %.1f  p90 %.1f\n\n",
            mean(SL$fb_pri), quantile(SL$fb_pri,.1),
            quantile(SL$fb_pri,.5), quantile(SL$fb_pri,.9)))

tt <- function(a, b, lab) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a)<8 || length(b)<8) { cat(sprintf("     %-12s n too small\n", lab)); return() }
  t <- t.test(a, b)
  cat(sprintf("     %-12s  A %+7.3f vs B %+7.3f   diff %+7.3f   p = %.3f   n=%d/%d\n",
              lab, mean(a), mean(b), mean(a)-mean(b), t$p.value, length(a), length(b)))
}

# =============================================================================
# 1. The unmatched view: harder slider looks better
# =============================================================================
cat("############ 1. Unmatched: SL velo bands (this is the stuff-model view) ############\n\n")
SL[, vband := fifelse(sl_velo>=88.5 & sl_velo<91.5, "89-91",
               fifelse(sl_velo>=85.5 & sl_velo<87.5, "85-87",
               fifelse(sl_velo>=88, "88+",
               fifelse(sl_velo>=85, "85-88", "<85"))))]
print(SL[, .(
  n=.N, sl=round(mean(sl_velo),1), fb=round(mean(fb_pri),1), gap=round(mean(gap),1),
  hmov=round(mean(hmov),2), vmov=round(mean(vmov),2),
  xrv=round(mean(xrv),3), whiff=round(mean(whiff),1), chase=round(mean(chase),1),
  zone=round(mean(zone),1), hard=round(mean(hard),1), xw=round(mean(xwcon),3)
), by=vband][order(vband)], row.names=FALSE)

cat("\n  Continuous, unmatched:\n")
print(round(summary(lm(xrv ~ sl_velo, data=SL))$coefficients, 4))
cat(sprintf("  cor(SL velo, FB velo) = %+.3f   cor(SL velo, gap) = %+.3f\n",
            cor(SL$sl_velo, SL$fb_pri), cor(SL$sl_velo, SL$gap)))
cat(sprintf("  cor(SL velo, xRV) = %+.3f   cor(gap, xRV) = %+.3f   cor(FB, xRV) = %+.3f\n",
            cor(SL$sl_velo, SL$xrv), cor(SL$gap, SL$xrv), cor(SL$fb_pri, SL$xrv)))

# =============================================================================
# 2. THE experiment: same FB velo, 90 SL vs 86 SL
# =============================================================================
cat("\n############ 2. Hold FB at 96–100: 89–91 SL vs 85–87 SL ############\n\n")
HARD <- SL[fb_pri>=96 & fb_pri<=100]
A <- HARD[sl_velo>=89 & sl_velo<91.5]   # ~90
B <- HARD[sl_velo>=85 & sl_velo<87.5]   # ~86
cat(sprintf("  Hard-FB pool: %d pitchers.  ~90 SL: %d    ~86 SL: %d\n",
            nrow(HARD), nrow(A), nrow(B)))
cat(sprintf("  Mean FB: 90-group %.1f    86-group %.1f\n",
            mean(A$fb_pri), mean(B$fb_pri)))
cat(sprintf("  Mean gap: 90-group %.1f    86-group %.1f\n\n",
            mean(A$gap), mean(B$gap)))

cat("  Outcomes (90 minus 86; positive = 90 is better for pitcher):\n")
for (v in c("xrv","whiff","wswing","chase","swing_z","zone","cs","ball",
            "hard","brl","xwcon","hmov","vmov")) tt(A[[v]], B[[v]], v)

cat("\n  xRV channels (90 minus 86):\n")
for (v in c("ch_whiff","ch_bip","ch_foul","ch_cs","ch_ball",
            "ch_chase","ch_ooz_take","ch_iz_swing","ch_iz_take")) tt(A[[v]], B[[v]], v)

# slightly wider for power
cat("\n  Wider bands, still FB 96–100: 88–92 vs 84–88\n")
Aw <- HARD[sl_velo>=88 & sl_velo<92]
Bw <- HARD[sl_velo>=84 & sl_velo<88]
cat(sprintf("  n = %d vs %d\n", nrow(Aw), nrow(Bw)))
for (v in c("xrv","whiff","chase","hard","xwcon","hmov","vmov",
            "ch_chase","ch_ooz_take","ch_iz_swing","ch_whiff","ch_bip"))
  tt(Aw[[v]], Bw[[v]], v)

# =============================================================================
# 3. Movement-matched: nearest neighbor on (hmov, vmov) within hard-FB pool
# =============================================================================
cat("\n############ 3. Same FB, movement-matched 90 vs 86 ############\n\n")
# standardize movement in the hard-FB 90/86 pool
POOL <- rbind(A[, tag:="90"], B[, tag:="86"])
POOL[, hz := as.numeric(scale(hmov))]
POOL[, vz := as.numeric(scale(vmov))]
# for each 90, nearest 86 on movement
nn <- rbindlist(lapply(seq_len(nrow(A)), function(i) {
  a <- POOL[tag=="90"][i]
  cand <- POOL[tag=="86"]
  d <- (cand$hz - a$hz)^2 + (cand$vz - a$vz)^2
  j <- which.min(d)
  data.table(p90=a$pitcher, p86=cand$pitcher[j], dist=sqrt(d[j]),
             h90=a$hmov, h86=cand$hmov[j], v90=a$vmov, v86=cand$vmov[j])
}))
# unique 86s (greedy: keep closest assignment per 86)
setorder(nn, dist)
nn <- nn[!duplicated(p86)]
cat(sprintf("  Matched pairs (unique 86s, nearest 90): %d\n", nrow(nn)))
cat(sprintf("  Mean |hmov| 90/86: %.2f / %.2f    vmov: %.2f / %.2f    match dist: %.2f SD\n\n",
            mean(nn$h90), mean(nn$h86), mean(nn$v90), mean(nn$v86), mean(nn$dist)))
M90 <- SL[pitcher %in% nn$p90]
M86 <- SL[pitcher %in% nn$p86]
cat("  Movement-matched outcomes (90 minus 86):\n")
for (v in c("xrv","whiff","wswing","chase","swing_z","zone","cs","ball",
            "hard","xwcon","sl_velo","gap","fb_pri","hmov","vmov",
            "ch_chase","ch_ooz_take","ch_iz_swing","ch_whiff","ch_bip","ch_cs","ch_ball"))
  tt(M90[[v]], M86[[v]], v)

# =============================================================================
# 4. Continuous, FB-held: does SL velo still help?
# =============================================================================
cat("\n############ 4. Continuous models ############\n\n")
cat("  All SL (unmatched):\n")
print(round(summary(lm(xrv ~ sl_velo + fb_pri + hmov + vmov, data=SL))$coefficients, 4))
cat("\n  FB 96–100 only:\n")
print(round(summary(lm(xrv ~ sl_velo + fb_pri + hmov + vmov, data=HARD))$coefficients, 4))
cat("\n  Same, replace sl_velo with gap (gap = fb - sl, so sign flips):\n")
print(round(summary(lm(xrv ~ gap + fb_pri + hmov + vmov, data=HARD))$coefficients, 4))
cat("\n  Chase ~ sl_velo | FB + movement, FB 96–100:\n")
print(round(summary(lm(chase ~ sl_velo + fb_pri + hmov + vmov, data=HARD))$coefficients, 4))
cat("\n  In-zone swing ~ sl_velo | FB + movement:\n")
print(round(summary(lm(swing_z ~ sl_velo + fb_pri + hmov + vmov, data=HARD))$coefficients, 4))
cat("\n  xwOBAcon ~ sl_velo | FB + movement:\n")
print(round(summary(lm(xwcon ~ sl_velo + fb_pri + hmov + vmov, data=HARD))$coefficients, 4))
cat("\n  Whiff ~ sl_velo | FB + movement:\n")
print(round(summary(lm(whiff ~ sl_velo + fb_pri + hmov + vmov, data=HARD))$coefficients, 4))

# =============================================================================
# 5. Gap terciles among hard-FB sliders: chase as an RV channel
# =============================================================================
cat("\n############ 5. Among 96–100 FB, gap terciles (not SL-velo bands) ############\n\n")
HARD[, g3 := cut(gap, quantile(gap, 0:3/3), include.lowest=TRUE,
                 labels=c("small gap (~90 SL)","mid","big gap (~86 SL)"))]
print(HARD[, .(
  n=.N, sl=round(mean(sl_velo),1), fb=round(mean(fb_pri),1), gap=round(mean(gap),1),
  hmov=round(mean(hmov),2), vmov=round(mean(vmov),2),
  xrv=round(mean(xrv),3),
  chase=round(mean(chase),1), swing_z=round(mean(swing_z),1),
  zone=round(mean(zone),1), ball=round(mean(ball),1), cs=round(mean(cs),1),
  whiff=round(mean(whiff),1), hard=round(mean(hard),1), xw=round(mean(xwcon),3),
  ch_chase=round(mean(ch_chase),3), ch_ooz_take=round(mean(ch_ooz_take),3),
  ch_iz_swing=round(mean(ch_iz_swing),3), ch_iz_take=round(mean(ch_iz_take),3)
), by=g3][order(g3)], row.names=FALSE)
hi <- HARD[g3=="big gap (~86 SL)"]; lo <- HARD[g3=="small gap (~90 SL)"]
cat("\n  Big-gap minus small-gap (positive = bigger gap is better):\n")
for (v in c("xrv","chase","swing_z","ball","cs","whiff","hard","xwcon",
            "ch_chase","ch_ooz_take","ch_iz_swing","ch_iz_take","ch_whiff","ch_bip","ch_cs","ch_ball"))
  tt(hi[[v]], lo[[v]], v)

# =============================================================================
# 6. Is the chase loss big enough to offset the contact win?
# =============================================================================
cat("\n############ 6. Channel sum: does chase give back the contact win? ############\n\n")
cat("  Big-gap minus small-gap, FB 96–100 sliders:\n")
chs <- c("ch_chase","ch_ooz_take","ch_iz_swing","ch_iz_take","ch_whiff","ch_bip","ch_foul","ch_cs","ch_ball")
# note some overlap: ch_chase is ooz swings (whiff+foul+bip ooz); ch_whiff is all whiffs
# Use a partition:
#   ooz swing, ooz take, iz swing, iz take  — these four add to total xRV
cat("  Partition (adds to xRV):\n")
for (v in c("ch_iz_swing","ch_iz_take","ch_chase","ch_ooz_take","xrv")) {
  d <- mean(hi[[v]],na.rm=TRUE) - mean(lo[[v]],na.rm=TRUE)
  cat(sprintf("     %-14s  %+7.3f\n", v, d))
}

# =============================================================================
# 7. Names: the actual 98/90 and 98/86 pitchers
# =============================================================================
cat("\n############ 7. Who is in the 96–100 FB, 89–91 vs 85–87 cells ############\n\n")
cat("  90 mph SL group, by xRV:\n")
print(A[order(-xrv), .(name, sl=round(sl_velo,1), fb=round(fb_pri,1), gap=round(gap,1),
                       hmov=round(hmov,2), vmov=round(vmov,2),
                       xrv=round(xrv,2), chase=round(chase,1), whiff=round(whiff,1),
                       hard=round(hard,1))][1:12], row.names=FALSE)
cat("\n  86 mph SL group, by xRV:\n")
print(B[order(-xrv), .(name, sl=round(sl_velo,1), fb=round(fb_pri,1), gap=round(gap,1),
                       hmov=round(hmov,2), vmov=round(vmov,2),
                       xrv=round(xrv,2), chase=round(chase,1), whiff=round(whiff,1),
                       hard=round(hard,1))][1:12], row.names=FALSE)

# save for canvas
fwrite(SL, "data/statcast_2026/hard_slider_vs_gap.csv")
cat("\nWrote pitcher slider table.\n")
