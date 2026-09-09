#!/usr/bin/env Rscript

# What kind of pitcher disrupts hitter timing?
#
# Outcome is the directional timing score built in timing_directional_score.R:
#   tscore = runs saved per 100 swings from putting hitters DEEP on fastballs and
#            OUT FRONT on breaking/offspeed, net of each hitter's own norms.
#   fb_rate / off_rate = the two halves, each per 100 swings of that group.
#   disrupt = the older direction-agnostic version (mean excess |mistiming|, inches).
#
# Candidate traits, grouped by the hypothesis they test:
#   ARSENAL SHAPE  entropy of pitch mix, effective pitch count, fastball share
#   SEPARATION     velocity spread, fastball-to-offspeed gap, movement spread
#   DECEPTION      release-point consistency across pitch types, arm angle, extension
#   PATH           path-to-location ratio (does movement create the location, or does aim?)
#                  vertical approach angle flatness
#   QUALITY        Location+ and Stuff+ analogs fit here from run value

suppressPackageStartupMessages({ library(data.table); library(splines) })
set.seed(7)
options(width = 210)

cols <- c("game_year","game_type","player_name","pitcher","batter","stand","p_throws",
          "pitch_type","description","balls","strikes","plate_x","plate_z",
          "sz_top","sz_bot","release_speed","release_pos_x","release_pos_z",
          "release_extension","release_spin_rate","pfx_x","pfx_z","spin_axis",
          "arm_angle","vx0","vy0","vz0","ax","ay","az","delta_run_exp",
          "game_pk","inning","intercept_ball_minus_batter_pos_y_inches")

dt <- rbindlist(lapply(c(2025, 2026), function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress = FALSE, select = cols)))
setnames(dt, "intercept_ball_minus_batter_pos_y_inches", "depth")
dt <- dt[game_type == "R" & pitch_type != "" & balls <= 3 & strikes <= 2 &
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
cat(sprintf("Pitches: %s (2025 %s, 2026 %s)\n\n", format(nrow(dt), big.mark=","),
  format(nrow(dt[game_year==2025]), big.mark=","),
  format(nrow(dt[game_year==2026]), big.mark=",")))

# ---- vertical approach angle -------------------------------------------------
dt[, vyf := -sqrt(pmax(vy0^2 - 2*ay*(50 - 17/12), 0))]
dt[, tf  := (vyf - vy0)/ay]
dt[, `:=`(vzf = vz0 + az*tf, vxf = vx0 + ax*tf)]
dt[, vaa := atan2(vzf, sqrt(vxf^2 + vyf^2)) * 180/pi]

# ---- path-to-location ratio --------------------------------------------------
# Where would this pitch have crossed with no spin-induced movement? pfx_x/pfx_z are
# exactly that displacement (gravity already removed), so subtracting them gives the
# "aim point". If movement carries the ball AWAY from the middle of the zone the ratio
# exceeds 1: the pitch's path is what creates its location. Below 1 means the pitcher
# aimed further off-center and movement pulled the ball back toward the heart.
dt[, d_act := sqrt(plate_x^2 + (plate_z - zc)^2)]
dt[, d_aim := sqrt((plate_x - pfx_x)^2 + (plate_z - pfx_z - zc)^2)]

# =============================================================================
# 1. Location+ and Stuff+ analogs
# =============================================================================
cat("############ 1. Fitting Location+ and Stuff+ analogs ############\n")
dt[, swung := description %in% c("swinging_strike","swinging_strike_blocked",
                                 "foul","foul_tip","hit_into_play")]
samp <- dt[sample(.N, min(350000, .N))]
cnt_only <- summary(lm(bat_rv ~ factor(cnt), data=samp))$r.squared
mLoc <- lm(bat_rv ~ ns(px_bat,6)*ns(pz_rel,6) + factor(cnt) + stand, data=samp)

# Stuff has almost no purchase on the run value of a TAKEN pitch — that channel is
# location deciding ball vs. strike. Fitting the physical-characteristic model on
# swings only gives stuff a fair test.
ok <- dt[, is.finite(release_spin_rate) & is.finite(release_extension) &
           is.finite(pfx_x) & is.finite(pfx_z)]
sfit <- dt[ok & swung == TRUE][sample(.N, min(250000, .N))]
mStf <- lm(bat_rv ~ ns(release_speed,5)*pg2 + ns(pfx_x,5)*ns(pfx_z,5) +
             ns(release_extension,4) + ns(release_spin_rate,4) + factor(cnt) +
             stand + p_throws, data=sfit)
mStf_all <- lm(bat_rv ~ ns(release_speed,5)*pg2 + ns(pfx_x,5)*ns(pfx_z,5) +
             ns(release_extension,4) + ns(release_spin_rate,4) + factor(cnt) +
             stand + p_throws,
             data=samp[is.finite(release_spin_rate) & is.finite(release_extension) &
                       is.finite(pfx_x) & is.finite(pfx_z)])
cat(sprintf("  Count fixed effects alone            R2 = %.4f\n", cnt_only))
cat(sprintf("  Location + count (all pitches)       R2 = %.4f\n", summary(mLoc)$r.squared))
cat(sprintf("  Stuff + count (all pitches)          R2 = %.4f  <- near zero, unusable\n",
            summary(mStf_all)$r.squared))
cat(sprintf("  Stuff + count (swings only)          R2 = %.4f  <- used below\n",
            summary(mStf)$r.squared))
cat("  A taken pitch's run value is decided by location, which is why the all-pitch\n")
cat("  stuff model is empty. Restricting to swings gives stuff a channel to work in.\n")
dt[, loc_rv := predict(mLoc, newdata=dt)]
dt[, stf_rv := NA_real_]
dt[ok, stf_rv := predict(mStf, newdata=dt[ok])]
cat("\n  Both are rescaled across pitchers to a + scale: 100 = league average,\n")
cat("  1 SD = 10 points, higher = better for the pitcher.\n\n")

# =============================================================================
# 2. Pitcher traits
# =============================================================================
traits <- function(d) {
  # pitch-type level first, for the arsenal-shape measures
  a <- d[, .(n=.N, velo=mean(release_speed), mx=mean(pfx_x, na.rm=TRUE),
             mz=mean(pfx_z, na.rm=TRUE), rx=mean(release_pos_x, na.rm=TRUE),
             rz=mean(release_pos_z, na.rm=TRUE)), by=.(pitcher, pitch_type, pg2)]
  a[, w := n/sum(n), by=pitcher]
  shape <- a[, {
    ww <- w[w >= 0.02]; ww <- ww/sum(ww)
    used <- .SD[w >= 0.05]
    # usage-weighted spread of each characteristic across the arsenal
    wm <- function(x, wt) sum(wt*x)/sum(wt)
    wsd <- function(x, wt) sqrt(pmax(sum(wt*(x - wm(x,wt))^2)/sum(wt), 0))
    # mean pairwise separation in movement space, usage-weighted
    ms <- if (.N > 1) {
      g <- expand.grid(i=seq_len(.N), j=seq_len(.N))
      g <- g[g$i < g$j, ]
      sum(w[g$i]*w[g$j]*sqrt((mx[g$i]-mx[g$j])^2 + (mz[g$i]-mz[g$j])^2)) /
        sum(w[g$i]*w[g$j])
    } else 0
    .(entropy   = -sum(ww*log(ww)),
      eff_pitch = exp(-sum(ww*log(ww))),
      n_used    = nrow(used),
      velo_sd   = wsd(velo, w),
      velo_gap  = if (nrow(used) > 1) max(used$velo) - min(used$velo) else 0,
      mov_sep   = ms,
      rel_sd    = sqrt(wsd(rx, w)^2 + wsd(rz, w)^2) * 12,   # inches
      fb_velo   = if (any(pg2=="FB")) wm(velo[pg2=="FB"], w[pg2=="FB"]) else NA_real_,
      off_velo  = if (any(pg2=="OFF")) wm(velo[pg2=="OFF"], w[pg2=="OFF"]) else NA_real_)
  }, by=pitcher]
  shape[, fb_off_gap := fb_velo - off_velo]

  # pitch level for everything else
  base <- d[, .(
    pitches   = .N,
    fbsh      = 100*mean(pg2=="FB"),
    ext       = mean(release_extension, na.rm=TRUE),
    armang    = mean(arm_angle, na.rm=TRUE),
    spin      = mean(release_spin_rate, na.rm=TRUE),
    vaa_fb    = mean(vaa[pg2=="FB"], na.rm=TRUE),
    path_ratio= mean(d_act, na.rm=TRUE)/mean(d_aim, na.rm=TRUE),
    zone_dist = mean(d_act, na.rm=TRUE),
    # role, so the trait set is not just separating starters from relievers
    ppa       = .N/uniqueN(game_pk),
    loc_raw   = -mean(loc_rv, na.rm=TRUE),
    stuff_raw = -mean(stf_rv, na.rm=TRUE)
  ), by=pitcher]
  out <- merge(shape, base, by="pitcher")
  # rescale command and stuff to the familiar + convention within this population
  out[, loc_plus   := 100 + 10*as.numeric(scale(loc_raw))]
  out[, stuff_plus := 100 + 10*as.numeric(scale(stuff_raw))]
  out[]
}

# ---- the outcome: directional timing score -----------------------------------
SWING <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip","hit_into_play")
sw <- dt[description %in% SWING & is.finite(depth)]
fitd <- lm(depth ~ ns(px_bat,5)*pgrp + ns(pz_rel,5) + ns(release_speed,4) + stand, data=sw)
sw[, r1 := residuals(fitd)]
sw[, nb := .N, by=batter]; sw <- sw[nb >= 200]
sw[, tdev := r1 - mean(r1), by=batter]
sw[, adev := abs(tdev)]
sw[, b_adev := mean(adev), by=batter]
surf <- lm(bat_rv ~ ns(tdev, 6)*pgrp + factor(cnt), data=sw)
sw[, trv := predict(surf, newdata=sw)]
sw[, trv_b := trv - mean(trv), by=batter]
sw[, trv_ex := trv_b - mean(trv_b), by=pg2]

outcome <- function(d) d[, .(
  swings   = .N,
  tscore   = -100*mean(trv_ex),
  fb_rate  = -100*mean(trv_ex[pg2=="FB"]),
  off_rate = -100*mean(trv_ex[pg2=="OFF"]),
  tdev_fb  = mean(tdev[pg2=="FB"]),
  tdev_off = mean(tdev[pg2=="OFF"]),
  disrupt  = mean(adev - b_adev),
  whiffpct = 100*mean(description %in%
              c("swinging_strike","swinging_strike_blocked","foul_tip"))
), by=pitcher]

build <- function(yr, minp = 700, mins = 250) {
  X <- traits(dt[game_year==yr][, .N, by=pitcher][N >= minp][
        dt[game_year==yr], on="pitcher", nomatch=0])
  Y <- outcome(sw[game_year==yr])[swings >= mins]
  merge(X, Y, by="pitcher")
}
A <- build(2026)
NAME <- unique(dt[, .(pitcher, player_name)], by="pitcher")
A <- merge(A, NAME, by="pitcher")
cat(sprintf("############ 2. Trait panel: %d pitchers, 2026 (700+ pitches, 250+ swings) ############\n\n",
            nrow(A)))

TR <- c("entropy","eff_pitch","n_used","velo_sd","velo_gap","fb_off_gap","mov_sep",
        "rel_sd","armang","ext","fb_velo","spin","vaa_fb","path_ratio","zone_dist",
        "fbsh","loc_plus","stuff_plus","ppa")
LAB <- c(entropy="Pitch-mix entropy", eff_pitch="Effective pitch count",
  n_used="Pitch types used 5%+", velo_sd="Velocity spread (SD, mph)",
  velo_gap="Fastest minus slowest (mph)", fb_off_gap="Fastball minus offspeed (mph)",
  mov_sep="Movement separation (ft)", rel_sd="Release scatter across arsenal (in)",
  armang="Arm angle (deg)", ext="Extension (ft)", fb_velo="Fastball velocity (mph)",
  spin="Spin rate (rpm)", vaa_fb="Fastball VAA (deg, less negative = flatter)",
  path_ratio="Path-to-location ratio", zone_dist="Distance from zone center (ft)",
  fbsh="Fastball usage %", loc_plus="Location+", stuff_plus="Stuff+",
  ppa="Pitches per appearance (role)")

cat("  Trait distribution:\n")
print(rbindlist(lapply(TR, function(v) data.table(trait=LAB[[v]],
  mean=round(mean(A[[v]], na.rm=TRUE),3), sd=round(sd(A[[v]], na.rm=TRUE),3),
  p10=round(quantile(A[[v]], .1, na.rm=TRUE),3),
  p90=round(quantile(A[[v]], .9, na.rm=TRUE),3)))), row.names=FALSE)

# =============================================================================
# 3. Univariate: which traits move the timing score?
# =============================================================================
cat("\n############ 3. Univariate correlations with the timing score ############\n\n")
U <- rbindlist(lapply(TR, function(v) {
  f <- function(y) cor(A[[v]], A[[y]], use="complete.obs")
  data.table(trait=LAB[[v]], tscore=f("tscore"), fb_rate=f("fb_rate"),
             off_rate=f("off_rate"), tdev_fb=f("tdev_fb"), tdev_off=f("tdev_off"),
             disrupt=f("disrupt"), whiff=f("whiffpct"))
}))
setorder(U, -tscore)
print(U[, lapply(.SD, function(x) if (is.numeric(x)) round(x,3) else x)], row.names=FALSE)
cat("\n  Reminder on signs: tdev_fb should be NEGATIVE (deep) and tdev_off POSITIVE\n")
cat("  (out front) for a good disruptor, so a trait helping the fastball half shows a\n")
cat("  negative tdev_fb correlation and a positive fb_rate correlation.\n")

# =============================================================================
# 4. Multivariate
# =============================================================================
cat("\n############ 4. Multivariate: which traits survive together? ############\n")
Z <- copy(A)
for (v in TR) Z[[v]] <- as.numeric(scale(Z[[v]]))
KEEP <- c("entropy","velo_gap","fb_off_gap","mov_sep","rel_sd","armang","ext",
          "fb_velo","spin","vaa_fb","path_ratio","fbsh","loc_plus","stuff_plus","ppa")

# collinearity check — several of these traits are near-duplicates by construction
vifs <- sapply(KEEP, function(v)
  1/(1 - summary(lm(as.formula(paste(v, "~", paste(setdiff(KEEP, v), collapse=" + "))),
                    data=Z))$r.squared))
cat("\n  Variance inflation factors (>5 means the coefficient is unstable):\n")
print(data.table(trait = LAB[names(vifs)], vif = round(vifs, 2))[order(-vif)],
      row.names=FALSE)

show_fit <- function(y, lab, d = Z) {
  m <- lm(as.formula(paste(y, "~", paste(KEEP, collapse=" + "))), data=d)
  s <- summary(m)
  co <- as.data.table(s$coefficients, keep.rownames="trait")[trait != "(Intercept)"]
  setnames(co, c("trait","beta","se","tval","p"))
  co[, label := LAB[trait]]
  co <- co[order(-abs(beta))]
  cat(sprintf("\n  -- %s | R2 = %.3f, adj = %.3f, n = %d --\n", lab, s$r.squared,
              s$adj.r.squared, nobs(m)))
  print(co[, .(label, beta=round(beta,3), t=round(tval,2),
               p=formatC(p, format="f", digits=3),
               sig=fifelse(p<0.001,"***",fifelse(p<0.01,"**",
                    fifelse(p<0.05,"*",fifelse(p<0.10,".","")))))], row.names=FALSE)
}
cat("\n  Standardized betas: change in the outcome, in SDs, per 1 SD of the trait.\n")
show_fit("tscore",   "Overall timing score")
show_fit("fb_rate",  "Fastball half (getting hitters deep on fastballs)")
show_fit("off_rate", "Offspeed half (getting hitters out front on breaking/offspeed)")
show_fit("disrupt",  "Magnitude-only disruption, for contrast")

cat("\n  ---- Robustness ----\n")
ZS <- Z[A$ppa >= 50]
cat(sprintf("  Starters only (50+ pitches per appearance), n = %d:", nrow(ZS)))
show_fit("tscore", "Overall timing score, starters only", ZS)

ZO <- Z[A$armang > 0]
cat(sprintf("\n  Excluding submarine/sidearm outliers (arm angle > 0 deg), n = %d:", nrow(ZO)))
show_fit("tscore", "Overall timing score, no submariners", ZO)

cat("\n  Rank (Spearman) correlations, immune to single outliers like Tyler Rogers:\n")
RK <- rbindlist(lapply(TR, function(v) data.table(trait=LAB[[v]],
  pearson = cor(A[[v]], A$tscore, use="complete.obs"),
  spearman = cor(A[[v]], A$tscore, method="spearman", use="complete.obs"))))
print(RK[order(-abs(spearman))][, lapply(.SD, function(x)
  if (is.numeric(x)) round(x,3) else x)], row.names=FALSE)

# =============================================================================
# 5. Out-of-sample: do 2025 traits forecast 2026 timing?
# =============================================================================
cat("\n############ 5. Do 2025 traits forecast 2026 timing? ############\n")
B <- build(2025)
M <- merge(B[, c("pitcher", TR), with=FALSE],
           outcome(sw[game_year==2026])[swings >= 250,
             .(pitcher, tscore26=tscore, fb26=fb_rate, off26=off_rate)], by="pitcher")
M <- merge(M, B[, .(pitcher, tscore25=tscore)], by="pitcher")
Zm <- copy(M); for (v in TR) Zm[[v]] <- as.numeric(scale(Zm[[v]]))
cat(sprintf("\n  n = %d pitchers present in both seasons\n\n", nrow(M)))
r2 <- function(f, d) summary(lm(as.formula(f), data=d))$r.squared
cat(sprintf("     2026 tscore from 2025 tscore alone        R2 = %.3f\n",
            r2("tscore26 ~ tscore25", Zm)))
cat(sprintf("     2026 tscore from 2025 traits alone        R2 = %.3f\n",
            r2(paste("tscore26 ~", paste(KEEP, collapse=" + ")), Zm)))
cat(sprintf("     2026 tscore from both                     R2 = %.3f\n",
            r2(paste("tscore26 ~ tscore25 +", paste(KEEP, collapse=" + ")), Zm)))
cat("\n  Univariate trait correlations with NEXT season's timing score:\n")
UO <- rbindlist(lapply(TR, function(v) data.table(trait=LAB[[v]],
  next_tscore = cor(Zm[[v]], Zm$tscore26, use="complete.obs"),
  next_fb     = cor(Zm[[v]], Zm$fb26, use="complete.obs"),
  next_off    = cor(Zm[[v]], Zm$off26, use="complete.obs"))))
setorder(UO, -next_tscore)
print(UO[, lapply(.SD, function(x) if (is.numeric(x)) round(x,3) else x)], row.names=FALSE)

# =============================================================================
# 6. What the best disruptors look like
# =============================================================================
cat("\n############ 6. Trait profile of the top and bottom of the leaderboard ############\n")
A[, band := fifelse(tscore >= quantile(tscore,.9), "top 10%",
            fifelse(tscore >= quantile(tscore,.7), "next 20%",
            fifelse(tscore <= quantile(tscore,.1), "bottom 10%", "middle")))]
PROF <- A[, c(.(pitchers=.N), lapply(.SD, function(x) round(mean(x, na.rm=TRUE),2))),
          by=band, .SDcols=c("tscore","fb_rate","off_rate", TR)]
setorder(PROF, -tscore)
print(PROF[, .(band, pitchers, tscore, fb_rate, off_rate, entropy, velo_gap,
               fb_off_gap, mov_sep, rel_sd, armang, ext, fb_velo, vaa_fb,
               path_ratio, fbsh, loc_plus, stuff_plus)], row.names=FALSE)

cat("\n  Ten best disruptors and their traits:\n")
print(A[order(-tscore)][1:10, .(player_name, swings, tscore=round(tscore,2),
  fb_rate=round(fb_rate,2), off_rate=round(off_rate,2), entropy=round(entropy,2),
  velo_gap=round(velo_gap,1), rel_sd=round(rel_sd,1), armang=round(armang,0),
  ext=round(ext,1), fb_velo=round(fb_velo,1), vaa=round(vaa_fb,1),
  path=round(path_ratio,2), locp=round(loc_plus,1), stfp=round(stuff_plus,1))],
  row.names=FALSE)
cat("\n  Ten worst:\n")
print(A[order(tscore)][1:10, .(player_name, swings, tscore=round(tscore,2),
  fb_rate=round(fb_rate,2), off_rate=round(off_rate,2), entropy=round(entropy,2),
  velo_gap=round(velo_gap,1), rel_sd=round(rel_sd,1), armang=round(armang,0),
  ext=round(ext,1), fb_velo=round(fb_velo,1), vaa=round(vaa_fb,1),
  path=round(path_ratio,2), locp=round(loc_plus,1), stfp=round(stuff_plus,1))],
  row.names=FALSE)

fwrite(A, file.path("data","statcast_2026","timing_disruptor_traits_2026.csv"))
cat("\nWrote data/statcast_2026/timing_disruptor_traits_2026.csv\n")
