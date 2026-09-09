#!/usr/bin/env Rscript

# PARACHUTE CHANGEUP ANALYSIS ON NCAA D1 TRACKMAN, 2023-2025.
#
# Mirrors the MLB pipeline: a four-seam anchor per pitcher-season, an out-of-fold whiff model on
# tjStuff+ v3.0 features plus location, and residuals read as overperformance. Three differences
# forced by the data, each flagged where it bites:
#
#   spin axis      TrackMan's is within ~5 degrees of a deterministic function of the break for
#                  changeups (script 01), so the axis gap here is close to a movement-direction
#                  gap. Movement is already a model feature, which makes this a harder test than
#                  the MLB one rather than an easier one.
#   efficiency     inferred from spin rate and Magnus magnitude, calibrated on MLB measured
#                  active spin (script 02: r = +.92 on four-seamers, +.75 on changeups). It runs
#                  about four points higher on NCAA than MLB, so thresholds are set by PERCENTILE
#                  within this population rather than copied across as absolute numbers.
#   arm slot       TrackMan reports none. A release-point proxy is built and its accuracy is
#                  measured against MLB arm angle before it is used for anything.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(11); options(width = 205)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")
CACHE <- file.path(MDIR, "ncaa_whiff_resid.rds"); MINSW <- 40L

D <- readRDS(file.path(MDIR, "ncaa_d1_pitches.rds"))
SH <- readRDS(file.path(MDIR, "ncaa_pitcher_season_shapes.rds"))

## ---- arm-slot proxy: build on MLB, measure it, then apply -------------------------------------
cat("=== how good is a release-point arm-slot proxy? (validated on MLB) ===\n")
M <- readRDS(file.path(MDIR, "parachute_ff.rds"))
MS <- M[is.finite(arm_angle), .(arm = mean(arm_angle), rx = mean(abs(release_pos_x)),
        rz = mean(release_pos_z), ext = mean(release_extension, na.rm=TRUE),
        velo = mean(release_speed)), by = .(pitcher, season)]
# The geometric proxy: height of the release above a nominal shoulder, against its distance to
# the side. No listed heights here, so the shoulder is a constant - which is exactly the
# approximation that limits this.
MS[, proxy := atan2(rz - 4.7, rx) * 180/pi]
cat(sprintf("  geometric proxy vs measured arm angle: r = %+.3f (R2 = %.2f) over %d seasons\n",
            cor(MS$proxy, MS$arm), cor(MS$proxy, MS$arm)^2, nrow(MS)))
fitp <- lm(arm ~ proxy + rz + rx + ext, MS)
cat(sprintf("  with a fitted correction:              r = %+.3f (R2 = %.2f)\n",
            cor(fitted(fitp), MS$arm), summary(fitp)$r.squared))
cat("  Good enough to rank slots, not to place a hard 44-degree gate. Used descriptively below.\n")

## ---- four-seam anchor and the changeup's gaps to it -------------------------------------------
D[, L := PitcherThrows == "Left"]
FF <- D[pt == "Four-Seam" & is.finite(RelSpeed),
        .(nff = .N, ff_speed = mean(RelSpeed), ff_ax = mean(ax, na.rm=TRUE),
          ff_az = mean(az, na.rm=TRUE), ff_axis = mean(SpinAxis, na.rm=TRUE)),
        by = .(PitcherId, season)][nff >= 30]
C <- D[pt == "Changeup"]
C <- merge(C, FF, by = c("PitcherId","season"))
C <- C[is.finite(RelSpeed) & is.finite(ax) & is.finite(az) & is.finite(SpinAxis) &
       is.finite(px) & is.finite(pz) & is.finite(SpinRate) & is.finite(Extension) &
       is.finite(VertApprAngle) & is.finite(HorzApprAngle)]
C[, `:=`(speed_diff = RelSpeed - ff_speed, ax_diff = ax - ff_ax, az_diff = az - ff_az,
         axis_diff = pmin(abs(SpinAxis - ff_axis), 360 - abs(SpinAxis - ff_axis)))]
# tjStuff+ v3.0 mirrors every x-dimension quantity for left-handers so both hands share a scale.
C[, `:=`(tj_x0 = fifelse(L, -RelSide, RelSide), tj_ax = fifelse(L, -ax, ax),
         tj_ax_diff = fifelse(L, -ax_diff, ax_diff), tj_px = fifelse(L, -px, px),
         tj_haa = fifelse(L, -HorzApprAngle, HorzApprAngle),
         tj_axis = fifelse(L, (360 - SpinAxis) %% 360, SpinAxis),
         same_hand = as.integer((PitcherThrows == "Left") == (BatterSide == "Left")))]
cat(sprintf("\nD1 changeups with a four-seam anchor: %s from %d pitcher-seasons\n",
            format(nrow(C), big.mark=","), uniqueN(paste(C$PitcherId, C$season))))

TJ  <- c("RelSpeed","SpinRate","Extension","tj_ax","az","tj_x0","RelHeight","tj_axis",
         "speed_diff","tj_ax_diff","az_diff")
LOC <- c("tj_px","pz","VertApprAngle","tj_haa","same_hand","Balls","Strikes")
SW  <- C[swing == TRUE]
cat(sprintf("changeup swings: %s   whiff rate %.3f\n", format(nrow(SW), big.mark=","), mean(SW$whiff)))

if (!file.exists(CACHE)) {
  fit <- function(FEAT, tag) {
    K <- 4; fold <- sample(rep(1:K, length.out = nrow(SW))); p <- rep(NA_real_, nrow(SW))
    for (f in 1:K) {
      tr <- SW[fold != f]; n <- nrow(tr); vi <- sample(n, floor(.12*n))
      dtr <- lgb.Dataset(as.matrix(tr[-vi, ..FEAT]), label = tr$whiff[-vi])
      dva <- lgb.Dataset.create.valid(dtr, as.matrix(tr[vi, ..FEAT]), label = tr$whiff[vi])
      m <- lgb.train(params = list(objective = "binary", metric = "binary_logloss",
                     learning_rate = .06, num_leaves = 31, min_data_in_leaf = 300,
                     feature_fraction = .8, bagging_fraction = .8, bagging_freq = 1),
                     data = dtr, nrounds = 1500, valids = list(v = dva),
                     early_stopping_rounds = 50, verbose = -1)
      p[fold == f] <- predict(m, as.matrix(SW[fold == f, ..FEAT]))
    }
    cat(sprintf("  %-28s %2d features  R2=%.4f\n", tag, length(FEAT), 1-var(SW$whiff-p)/var(SW$whiff)))
    SW$whiff - p
  }
  cat("\n=== out-of-fold whiff models on D1 changeups ===\n")
  SW[, r_tj  := fit(TJ, "tjStuff+ v3.0")]
  SW[, r_all := fit(c(TJ, LOC), "tjStuff+ v3.0 + location")]
  saveRDS(SW[, .(PitcherId, Pitcher, season, axis_diff, speed_diff, RelHeight, RelSide,
                 r_tj, r_all)], CACHE)
} else cat("\n(using cached fits)\n")
R <- readRDS(CACHE)

## ---- pitcher-season table ---------------------------------------------------------------------
S <- R[, .(nsw = .N, axis = mean(axis_diff), velo_sep = -mean(speed_diff), rz = mean(RelHeight),
           rx = mean(abs(RelSide)), w = 100*mean(r_all), w_tj = 100*mean(r_tj),
           name = Pitcher[1]), by = .(PitcherId, season)]
S[, slot := atan2(rz - 4.7, rx) * 180/pi]
S <- merge(S, SH[pt == "Changeup", .(PitcherId, season, eff_ch = eff, n_ch = n)],
           by = c("PitcherId","season"))
S <- merge(S, SH[pt == "Four-Seam", .(PitcherId, season, eff_ff = eff)],
           by = c("PitcherId","season"))
S <- S[nsw >= MINSW]
S[, eff_gap := eff_ch - eff_ff]
cat(sprintf("\npitcher-seasons with %d+ changeup swings: %d (%d pitchers)\n",
            MINSW, nrow(S), uniqueN(S$PitcherId)))
cat(sprintf("axis gap: mean %.1f, median %.1f | inferred eff CH %.3f FF %.3f | slot %.1f\n",
            mean(S$axis), median(S$axis), mean(S$eff_ch), mean(S$eff_ff), mean(S$slot)))

# Percentile thresholds, since the inferred efficiency sits on a different scale to MLB's
# measured figure and copying absolute cut points across would move the bin size.
q <- function(x, p) as.numeric(quantile(x, p))
EFF85 <- q(pmin(S$eff_ch, S$eff_ff), .35)   # matches the share of MLB seasons clearing .85 on both
SLOTQ <- q(S$slot, .70)
cat(sprintf("thresholds: efficiency floor %.3f (35th pct of the lower of the two), slot %.1f (70th pct)\n",
            EFF85, SLOTQ))

S[, `:=`(matched   = axis <= 10,
         eff_floor = pmin(eff_ch, eff_ff) >= EFF85,
         eff_close = abs(eff_gap) <= .10,
         hi_slot   = slot >= SLOTQ)]
S[, `:=`(core = matched & eff_close & hi_slot,
         corex = matched & eff_close & hi_slot & eff_floor)]

cat("\n=== does any cue predict whiff overperformance in D1? ===\n")
cues <- list(c("Axis gap <= 10 deg","matched"), c("Inferred efficiency floor","eff_floor"),
             c("|Efficiency gap| <= .10","eff_close"), c("High slot (70th pct)","hi_slot"),
             c("Core: matched + close + slot","core"), c("Core + efficiency floor","corex"))
print(rbindlist(lapply(cues, function(x) { i <- S[[x[2]]]
  t <- t.test(S$w[i], S$w[!i]); tt <- t.test(S$w_tj[i], S$w_tj[!i])
  data.table(cue = x[1], n = sum(i), whiff_vs_full = round(diff(rev(t$estimate)),2),
             p = round(t$p.value,4), whiff_vs_tjonly = round(diff(rev(tt$estimate)),2),
             p_tj = round(tt$p.value,4)) })), row.names = FALSE)

cat("\n=== velocity separation vs whiff residual: the MLB result, retested ===\n")
vc <- function(i, lab) { D <- S[i]; if (nrow(D) < 8) return(NULL); ct <- cor.test(D$velo_sep, D$w)
  data.table(group = lab, n = nrow(D), r = round(ct$estimate,3),
             ci = sprintf("[%+.2f, %+.2f]", ct$conf.int[1], ct$conf.int[2]),
             slope = round(coef(lm(w ~ velo_sep, D))[2],3), p = round(ct$p.value,4),
             spearman = round(cor(D$velo_sep, D$w, method="spearman"),3)) }
print(rbindlist(list(vc(S$core, "Core bin"), vc(S$corex, "Core + eff floor"),
                     vc(S$matched, "axis <= 10 only"), vc(!S$core, "everything else"),
                     vc(rep(TRUE, nrow(S)), "all D1 seasons"))), row.names = FALSE)
for (b in c("core","corex")) { S[, bb := get(b)]
  d <- summary(lm(w ~ velo_sep*bb, S, weights = nsw))$coefficients["velo_sep:bbTRUE",]
  cat(sprintf("  interaction, %-6s %+.2f pp per mph, p = %.4f\n", b, d[1], d[4])) }

fwrite(S, file.path(AST, "ncaa_parachute_seasons.csv"))
cat("\nwrote ncaa_parachute_seasons.csv\n")
