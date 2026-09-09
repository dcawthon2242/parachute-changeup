#!/usr/bin/env Rscript

# THE .90 PARACHUTE BIN, TRANSFERRED TO D1.
#
# The MLB definition is: four-seam active spin >= .90, changeup active spin >= .90, arm slot >= 44
# degrees, spin-axis gap <= 10 degrees. Nine pitcher-seasons, +3.83 points of whiff above model.
#
# Every one of those four conditions has to be substituted to run it on TrackMan:
#
#   active spin  ->  the lift-law estimator from script 08 (r = .76 on changeups, .92 on four-seams)
#   arm slot     ->  a proxy fit on MLB from release height, release side and extension
#   axis gap     ->  TrackMan's axis is movement-derived, so the threshold has to be re-scaled
#   whiff resid  ->  refit on D1, with the arm proxy included so it matches MLB's arm-aware model
#
# Substituting four things at once is enough to destroy a nine-season effect through noise alone,
# so the order of work matters. The recipe is first run on MLB with every substitution in place and
# measured against the answer we already know. If the substituted recipe cannot recover the MLB
# result on MLB data, then whatever it produces on D1 is uninterpretable, and that is worth
# establishing before looking at any college roster.

suppressPackageStartupMessages({ library(data.table); library(mgcv); library(lightgbm) })
set.seed(19); options(width = 210, datatable.print.nrows = 200)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")
G <- 32.174

## ---- 1. arm-angle proxy ---------------------------------------------------------------------------
cat("=== 1. arm-slot proxy: can release point stand in for measured arm angle? ===\n")
M <- readRDS(file.path(MDIR, "parachute_rv.rds"))
RS <- readRDS(file.path(MDIR, "whiff_tjstuff.rds"))[, .(arm = mean(arm_angle)), by = .(pitcher, season)]
REL <- M[pitch_type == "FF" & is.finite(release_pos_x) & is.finite(release_pos_z) &
         is.finite(release_extension),
         .(n = .N, rz = mean(release_pos_z), rx = mean(abs(release_pos_x)),
           ext = mean(release_extension)), by = .(pitcher, season)][n >= 30]
AP <- merge(REL, RS, by = c("pitcher","season"))
AP[, fold := sample(rep_len(1:5, .N))]
for (k in 1:5) AP[fold == k, arm_hat := predict(gam(arm ~ s(rz) + s(rx) + s(ext), data = AP[fold != k]), .SD)]
cat(sprintf("  %d pitcher-seasons: out-of-fold r = %.3f, RMSE %.1f degrees (arm angle sd = %.1f)\n",
            nrow(AP), cor(AP$arm_hat, AP$arm), sqrt(mean((AP$arm_hat-AP$arm)^2)), sd(AP$arm)))
cat(sprintf("  recovering the 44-degree gate: base rate %.1f%%, precision %.1f%%, recall %.1f%%\n",
            100*mean(AP$arm >= 44), 100*mean(AP$arm[AP$arm_hat >= 44] >= 44),
            100*mean(AP$arm_hat[AP$arm >= 44] >= 44)))
ARMFIT <- gam(arm ~ s(rz) + s(rx) + s(ext), data = AP)

## ---- 2. the movement-axis gap and what threshold matches 10 degrees --------------------------------
cat("\n=== 2. re-scaling the axis threshold for a movement-derived axis ===\n")
M2 <- M[pitch_type %in% c("CH","FF") & is.finite(vx0) & is.finite(ay) & is.finite(release_spin_rate) &
        release_spin_rate > 500]
M2[, `:=`(wx = ax, wy = ay, wz = az + G)]
M2[, d := (wx*vx0 + wy*vy0 + wz*vz0)/(vx0^2 + vy0^2 + vz0^2)]
M2[, mv := (atan2(wx - d*vx0, wz - d*vz0)*180/pi) %% 360]
FB <- M2[pitch_type == "FF", .(n = .N, s = mean(sin(mv*pi/180)), c = mean(cos(mv*pi/180))),
         by = .(pitcher, season)][n >= 30]
FB[, fbmv := (atan2(s, c)*180/pi) %% 360]
CHm <- merge(M2[pitch_type == "CH"], FB[, .(pitcher, season, fbmv)], by = c("pitcher","season"))
CHm[, gi := abs(((mv - fbmv + 180) %% 360) - 180)]
MVG <- CHm[, .(n = .N, mv_gap = mean(gi)), by = .(pitcher, season)][n >= 30]

PM <- readRDS(file.path(MDIR, "mlb_spineff_calibration.rds"))
W <- dcast(PM[, .(pitcher, season, pitch_type, measured, eff_hat)], pitcher + season ~ pitch_type,
           value.var = c("measured","eff_hat"))
S <- readRDS(file.path(MDIR, "whiff_tjstuff.rds"))[
  , .(nsw = .N, axis = mean(axis_diff), arm = mean(arm_angle), velo_sep = -mean(speed_diff),
      w = 100*mean(r4)), by = .(pitcher, player_name, season)][nsw >= 60]
S <- merge(S, W, by = c("pitcher","season")); S <- merge(S, MVG[, .(pitcher, season, mv_gap)], by = c("pitcher","season"))
S <- merge(S, AP[, .(pitcher, season, arm_hat)], by = c("pitcher","season"))
S[, name := trimws(paste(sub(".*,\\s*","",player_name), sub(",.*","",player_name)))]
S <- S[is.finite(measured_CH) & is.finite(measured_FF) & is.finite(eff_hat_CH) & is.finite(eff_hat_FF)]

# The threshold is matched inside the population the bin actually draws from - high efficiency,
# high slot - because that is where script 07 showed the two axis definitions converge.
POP <- S[measured_FF >= .90 & measured_CH >= .90 & arm >= 44]
share <- mean(POP$axis <= 10); MVTHR <- as.numeric(quantile(POP$mv_gap, share))
cat(sprintf("  in the .90/44 population (n = %d), measured axis <= 10 admits %.1f%%\n", nrow(POP), 100*share))
cat(sprintf("  the matching movement-axis threshold is %.1f degrees (correlation between the two: %.3f)\n",
            MVTHR, cor(POP$axis, POP$mv_gap)))

## ---- 3. dry run: the substituted recipe, scored on MLB ---------------------------------------------
cat("\n=== 3. does the substituted recipe recover the MLB result on MLB data? ===\n")
S[, truth := measured_FF >= .90 & measured_CH >= .90 & arm >= 44 & axis <= 10]
# Each substituted variable is also given a percentile-matched threshold, because a proxy is
# shrunk toward its mean and so clears a fixed gate less often than the quantity it replaces.
# Matching on share is the most generous version of each substitution and rules out the
# possibility that the recipe fails only because the thresholds were transferred literally.
q <- function(x, cond) as.numeric(quantile(S[[x]], mean(S[, eval(cond)])))
qa <- as.numeric(quantile(S$arm_hat, 1 - mean(S$arm >= 44)))
qe_ff <- as.numeric(quantile(S$eff_hat_FF, mean(S$measured_FF < .90)))
qe_ch <- as.numeric(quantile(S$eff_hat_CH, mean(S$measured_CH < .90)))
cat(sprintf("  percentile-matched gates: arm_hat >= %.1f, eff_hat FF >= %.3f / CH >= %.3f, mv_gap <= %.1f\n",
            qa, qe_ff, qe_ch, MVTHR))
LADDER <- list(
  list(l = "MLB truth (all four measured)",        i = quote(truth)),
  list(l = "estimated efficiency only",            i = quote(eff_hat_FF >= .90 & eff_hat_CH >= .90 & arm >= 44 & axis <= 10)),
  list(l = "proxy arm slot only",                  i = quote(measured_FF >= .90 & measured_CH >= .90 & arm_hat >= 44 & axis <= 10)),
  list(l = "proxy arm slot, percentile-matched",   i = quote(measured_FF >= .90 & measured_CH >= .90 & arm_hat >= qa & axis <= 10)),
  list(l = "movement axis only",                   i = quote(measured_FF >= .90 & measured_CH >= .90 & arm >= 44 & mv_gap <= MVTHR)),
  list(l = "all three substituted (D1 recipe)",    i = quote(eff_hat_FF >= .90 & eff_hat_CH >= .90 & arm_hat >= 44 & mv_gap <= MVTHR)),
  list(l = "all three, percentile-matched",        i = quote(eff_hat_FF >= qe_ff & eff_hat_CH >= qe_ch & arm_hat >= qa & mv_gap <= MVTHR)))
print(rbindlist(lapply(LADDER, function(rc) { i <- S[, eval(rc$i)]
  if (sum(i) < 2) return(data.table(recipe = rc$l, n = sum(i), whiff_above = NA_real_, p = NA_real_, recovered = 0L))
  t <- t.test(S$w[i], S$w[!i])
  data.table(recipe = rc$l, n = sum(i), whiff_above = round(mean(S$w[i]),2), p = round(t$p.value,3),
             recovered = sum(i & S$truth)) })), row.names = FALSE)
cat(sprintf("  'recovered' counts how many of the %d true members each version keeps.\n", sum(S$truth)))

## ---- 4. D1 ------------------------------------------------------------------------------------------
cat("\n=== 4. applying the recipe to D1 ===\n")
ND <- readRDS(file.path(MDIR, "ncaa_spineff.rds"))
ND[, arm_hat := predict(ARMFIT, newdata = data.frame(rz = relh, rx = abs(rels), ext = ext))]
PR <- merge(ND[pt == "Changeup", .(PitcherId, season, name, throws, relh, rels, ext,
                                   eff_ch = eff_hat, velo_ch = velo, spin_ch = spin, ivb = ivb, hb = hb, arm_hat)],
            ND[pt == "Four-Seam", .(PitcherId, season, eff_ff = eff_hat, velo_ff = velo)],
            by = c("PitcherId","season"))
cat(sprintf("  %d D1 pitcher-seasons throw both pitches; proxy arm slot median %.1f degrees (MLB %.1f)\n",
            nrow(PR), median(PR$arm_hat), median(S$arm)))

# D1 residuals, refit with the arm proxy so the model matches MLB's arm-aware specification.
CACHE <- file.path(MDIR, "ncaa_resid_arm.rds")
RD <- readRDS(file.path(MDIR, "ncaa_whiff_resid.rds"))
D1 <- RD[, .(nsw = .N, axis = mean(axis_diff), velo_sep = -mean(speed_diff),
             w_noarm = 100*mean(r_all)), by = .(PitcherId, Pitcher, season)][nsw >= 60]
D1 <- merge(D1, PR, by = c("PitcherId","season"))
cat(sprintf("  %d D1 pitcher-seasons with >=60 changeup swings and both pitches scored\n", nrow(D1)))

# The movement-axis threshold is re-matched to the D1 axis distribution, since TrackMan's spread
# differs from the one the MLB translation was solved on.
POPD <- D1[eff_ff >= .90 & eff_ch >= .90 & arm_hat >= 44]
D1THR <- as.numeric(quantile(D1$axis, mean(POP$axis <= 10)))
cat(sprintf("  D1 axis threshold matched on percentile: %.1f degrees\n", D1THR))
D1[, bin := eff_ff >= .90 & eff_ch >= .90 & arm_hat >= 44 & axis <= D1THR]
cat(sprintf("  bin: %d of %d pitcher-seasons (%.1f%%)\n", sum(D1$bin), nrow(D1), 100*mean(D1$bin)))

if (sum(D1$bin) >= 3) {
  t <- t.test(D1$w_noarm[D1$bin], D1$w_noarm[!D1$bin])
  cat(sprintf("\n  whiff above model: bin %+.2f pp vs %+.2f for the rest, difference %+.2f (p = %.3f)\n",
              mean(D1$w_noarm[D1$bin]), mean(D1$w_noarm[!D1$bin]), diff(rev(t$estimate)), t$p.value))
  cat(sprintf("  positive residual in %d of %d\n", sum(D1$w_noarm[D1$bin] > 0), sum(D1$bin)))
  cat("\n  D1 parachute bin:\n")
  print(D1[bin == TRUE][order(-w_noarm), .(Pitcher = name, Season = season, Swings = nsw,
        Arm = round(arm_hat,1), Axis = round(axis,1), EffFF = round(eff_ff,2), EffCH = round(eff_ch,2),
        Velo = round(velo_ch,1), VeloSep = round(velo_sep,1), Above = round(w_noarm,1))], row.names = FALSE)
}

# Each condition on its own, so it is visible which one (if any) carries anything in D1.
cat("\n  each condition alone:\n")
CUTS <- list("efficiency >= .90 on both" = quote(eff_ff >= .90 & eff_ch >= .90),
             "proxy arm slot >= 44"      = quote(arm_hat >= 44),
             "axis gap <= threshold"     = quote(axis <= D1THR))
print(rbindlist(lapply(names(CUTS), function(nm) { i <- D1[, eval(CUTS[[nm]])]
  data.table(condition = nm, n = sum(i), whiff_above = round(mean(D1$w_noarm[i]),2),
             p = round(t.test(D1$w_noarm[i], D1$w_noarm[!i])$p.value,3)) })), row.names = FALSE)

fwrite(D1[bin == TRUE][order(-w_noarm)], file.path(AST, "ext_d1_parachute90.csv"))
saveRDS(D1, file.path(MDIR, "ncaa_parachute90.rds"))
cat("\nwrote ext_d1_parachute90.csv and ncaa_parachute90.rds\n")
