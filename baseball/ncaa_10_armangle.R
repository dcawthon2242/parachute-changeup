#!/usr/bin/env Rscript

# Rebuild the D1 arm-slot estimate now that listed heights exist, and re-test the parachute gate.
#
# The arm condition is what broke the college transfer. Fit on release height, release side and
# extension alone the proxy reached r = .765 against Statcast's measured arm angle, but at the
# 44-degree gate the bin actually uses it ran 58 percent precision and 42 percent recall - and a
# four-condition conjunction cannot survive a component that misses more than half its members.
#
# Height is the missing input because arm angle is measured from the release point relative to the
# SHOULDER, and shoulder position scales with stature: two pitchers releasing at the same absolute
# height have different arm angles if they are different sizes. The estimator is calibrated on MLB,
# where measured arm angle exists as truth, and then applied to D1, where it does not.
#
# Note what this does and does not fix. It repairs one of the two substituted conditions. The other
# - the spin axis gap - is substituted with a movement-derived axis on the TrackMan side, and that
# recovered none of the true MLB bin members on its own. So the honest expectation is a better arm
# gate, not a working parachute bin.

suppressPackageStartupMessages({ library(data.table); library(mgcv) })
set.seed(23); options(width = 205)
MDIR <- "data/statcast_model"; DIR <- "data/ncaa_rosters"

## ---- calibrate on MLB, where arm angle is measured -------------------------------------------
H  <- fread("data/pitcher_heights.csv")
M  <- readRDS(file.path(MDIR, "parachute_rv.rds"))
RS <- readRDS(file.path(MDIR, "whiff_tjstuff.rds"))[, .(arm = mean(arm_angle)), by = .(pitcher, season)]
REL <- M[pitch_type == "FF" & is.finite(release_pos_x) & is.finite(release_pos_z) &
         is.finite(release_extension),
         .(n = .N, rz = mean(release_pos_z), rx = mean(abs(release_pos_x)),
           ext = mean(release_extension)), by = .(pitcher, season)][n >= 50]
A <- merge(merge(REL, RS, by = c("pitcher","season")), H[, .(pitcher, ht = ht_in)], by = "pitcher")[is.finite(ht)]
A[, `:=`(sh_z = 0.70*ht/12, sh_x = 0.115*ht/12)]
A[, geo := atan2(rz - sh_z, pmax(rx - sh_x, .01))*180/pi]

FORM <- arm ~ s(geo) + s(rz) + s(rx) + s(ext) + s(ht)
A[, fold := sample(rep_len(1:5, .N))]
for (k in 1:5) A[fold == k, hat := predict(gam(FORM, data = A[fold != k]), .SD)]
cat(sprintf("MLB calibration on %d pitcher-seasons: out-of-fold r = %.3f, RMSE %.2f degrees\n",
            nrow(A), cor(A$hat, A$arm), sqrt(mean((A$hat-A$arm)^2))))
g <- A$hat >= 44; t <- A$arm >= 44
cat(sprintf("  44-degree gate: precision %.1f%%, recall %.1f%% (release point alone was 58/41)\n",
            100*mean(t[g]), 100*mean(g[t])))
FIT <- gam(FORM, data = A)

## ---- apply to D1 -----------------------------------------------------------------------------
HT <- fread(file.path(DIR, "D1PitcherHeights.csv"))
HT <- HT[, .(ht = median(PitcherHeight)), by = .(PitcherId, season = year)]
ND <- as.data.table(readRDS(file.path(MDIR, "ncaa_spineff.rds")))
setnames(ND, old = intersect(names(ND), c("relh","rels")), new = c("rz","rx"), skip_absent = TRUE)
D <- merge(ND[pt == "Four-Seam"], HT, by.x = c("PitcherId","season"), by.y = c("PitcherId","season"))
cat(sprintf("\nD1: %d four-seam pitcher-seasons, %d (%.1f%%) now carry a listed height\n",
            nrow(ND[pt == "Four-Seam"]), nrow(D), 100*nrow(D)/nrow(ND[pt == "Four-Seam"]))) 

D[, `:=`(rx = abs(rx), sh_z = 0.70*ht/12, sh_x = 0.115*ht/12)]
D[, geo := atan2(rz - sh_z, pmax(rx - sh_x, .01))*180/pi]
D[, arm_hat := as.numeric(predict(FIT, D))]
cat(sprintf("estimated D1 arm slot: median %.1f, iqr %.1f-%.1f, share >= 44 deg: %.1f%%\n",
            median(D$arm_hat), quantile(D$arm_hat,.25), quantile(D$arm_hat,.75),
            100*mean(D$arm_hat >= 44)))
cat(sprintf("  MLB for comparison:  median %.1f, share >= 44 deg: %.1f%%\n",
            median(A$arm), 100*mean(A$arm >= 44)))
saveRDS(D[, .(PitcherId, season, ht, arm_hat, geo)], file.path(MDIR, "ncaa_armangle.rds"))

## ---- does the parachute gate work now? -------------------------------------------------------
PR <- as.data.table(readRDS(file.path(MDIR, "ncaa_spineff_pairs.rds")))
PR <- merge(PR, D[, .(PitcherId, season, arm_hat, ht)], by = c("PitcherId","season"))
# The residual file is per pitch, so it has to be collapsed to the pitcher-season the bin is
# defined on. r_tj is the residual against the tjStuff-style feature set, which is the like-for-like
# comparison with the MLB work; it is averaged and expressed in percentage points.
W <- readRDS(file.path(MDIR, "ncaa_whiff_resid.rds"))
if (!is.data.table(W)) setDT(W)
W <- W[!is.na(PitcherId), .(w_above = 100*mean(r_tj, na.rm = TRUE), n_ch = .N),
       by = .(PitcherId, season)][n_ch >= 60]
PR <- merge(PR, W[, .(PitcherId, season, w_above)], by = c("PitcherId","season"), all.x = TRUE)
rc <- "w_above"
cat(sprintf("\n%d of %d pairs have a whiff-above-model value (60+ changeups)\n",
            sum(is.finite(PR$w_above)), nrow(PR)))

show <- function(lab, idx) {
  x <- PR[[rc]][idx & is.finite(PR[[rc]])]; y <- PR[[rc]][!idx & is.finite(PR[[rc]])]
  if (length(x) < 3) { cat(sprintf("  %-34s n = %d (too few)\n", lab, length(x))); return(invisible()) }
  tt <- t.test(x, y)
  cat(sprintf("  %-34s n = %4d  whiff above model %+.2f pp  (rest %+.2f)  p = %.3f\n",
              lab, length(x), mean(x), mean(y), tt$p.value))
}

# The efficiency floor is swept rather than fixed. Lowering it trades purity for sample size: the
# MLB calibration put roughly 95 percent purity against a true .85 floor at an ESTIMATE of .90, so
# .85 admits a meaningful share of pitches that are not actually high efficiency. If the effect is
# real it should strengthen as the floor rises and the bin gets cleaner; if it is noise it will
# wander. That pattern across the sweep is more informative than any single threshold.
for (TH in c(.85, .88, .90)) {
  PR[, bin := eff_ch >= TH & eff_ff >= TH & arm_hat >= 44]
  cat(sprintf("\n=== efficiency floor %.2f, height-informed arm slot ===\n", TH))
  show(sprintf("efficiency only (%.2f/%.2f)", TH, TH), PR[, eff_ch >= TH & eff_ff >= TH])
  show("+ arm slot >= 44",              PR$bin)
  show("+ arm slot and axis gap <= 10", PR[, bin & axis_gap <= 10])
  B <- PR[bin == TRUE & is.finite(get(rc))]
  if (nrow(B) >= 10) {
    ct <- cor.test(B$velo_sep, B[[rc]])
    cat(sprintf("  velo separation inside the bin: r = %+.3f (p = %.3f, n = %d)\n",
                ct$estimate, ct$p.value, nrow(B)))
  } else cat("  velo separation: too few pairs to test\n")
}
