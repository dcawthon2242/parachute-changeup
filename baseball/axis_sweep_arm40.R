#!/usr/bin/env Rscript

# ARM SLOT HELD AT 40, THE AXIS THRESHOLD SWEPT.
#
# 40 degrees sits just under the D1 top third and just under MLB's, so it is roughly comparable
# across the two leagues without being as punishing as 44 was in the compressed D1 distribution.
# With slot fixed, the axis threshold is the only free parameter left, and sweeping it trades
# purity for sample in the one dimension the hypothesis actually cares about.
#
# What to look for: the effect should DECAY as the threshold loosens if matched spin is what
# matters, because each step admits pitches whose axes are further from the fastball's. A flat
# profile would mean the axis condition is not selecting on anything relevant, and a profile that
# only appears at the tightest cut is more likely to be the tail of a noisy distribution than a
# mechanism.

suppressPackageStartupMessages(library(data.table)); options(width = 205)
MDIR <- "data/statcast_model"; GATE <- 40L; EFF <- .85; ARMMIN <- 40

ARM <- as.data.table(readRDS(file.path(MDIR, "ncaa_armangle.rds")))
PR  <- as.data.table(readRDS(file.path(MDIR, "ncaa_spineff_pairs.rds")))
PR  <- merge(PR, ARM[, .(PitcherId, season, arm_hat)], by = c("PitcherId","season"))
setnames(PR, "axis_gap", "axis", skip_absent = TRUE)
CR  <- readRDS(file.path(MDIR, "ncaa_chase_resid.rds")); setDT(CR)
D <- merge(PR, CR[!is.na(PitcherId), .(y = 100*mean(r_all, na.rm = TRUE), n = .N),
                  by = .(PitcherId, season)], by = c("PitcherId","season"))[n >= GATE]
D[, `:=`(arm = arm_hat, id = PitcherId)]

F  <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(F); F <- F[pitch_type == "CH"]
AS <- readRDS(file.path(MDIR, "active_spin_long.rds")); setDT(AS)
R  <- readRDS(file.path(MDIR, "mlb_chase_resid.rds")); setDT(R)
S <- F[, .(axis = mean(axis_diff, na.rm = TRUE), arm = mean(arm_angle, na.rm = TRUE)),
       by = .(pitcher, season)]
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, eff_ch = active_spin)], by = c("pitcher","season"))
S <- merge(S, AS[pitch_type == "FF", .(pitcher, season, eff_ff = active_spin)], by = c("pitcher","season"))
S <- merge(S, R[, .(y = 100*mean(r_all), n = .N), by = .(pitcher, season)], by = c("pitcher","season"))
S <- S[n >= GATE & is.finite(arm)][, id := pitcher]

sweep <- function(X, lg) {
  P <- X[eff_ch >= EFF & eff_ff >= EFF & arm >= ARMMIN]
  cat(sprintf("\n=== %s: eff >= %.2f both + arm slot >= %d, axis gap swept ===\n", lg, EFF, ARMMIN))
  cat(sprintf("  that base pool is %d seasons from %d arms; its axis gaps run median %.1f, quartiles %.1f / %.1f\n",
              nrow(P), uniqueN(P$id), median(P$axis), quantile(P$axis,.25), quantile(P$axis,.75)))
  rbindlist(lapply(c(8, 10, 12, 15, 20, 25, 30, 999), function(a) {
    i <- X[, eff_ch >= EFF & eff_ff >= EFF & arm >= ARMMIN & axis <= a]
    x <- X$y[i]; y <- X$y[!i]
    if (length(x) < 3) return(NULL)
    t <- t.test(x, y)
    A <- data.table(id = X$id[i], v = x)[, .(m = mean(v)), by = id]
    data.table(league = lg, axis_max = if (a == 999) NA_integer_ else a,
               n = length(x), arms = nrow(A), bin = round(mean(x),2),
               diff = round(mean(x)-mean(y),2), p = round(t$p.value,4),
               se = round(sd(x)/sqrt(length(x)),2),
               se_by_arm = if (nrow(A) > 1) round(sd(A$m)/sqrt(nrow(A)),2) else NA_real_) }))
}
print(sweep(D, "D1"), row.names = FALSE)
print(sweep(S, "MLB"), row.names = FALSE)

# A threshold sweep hides the shape. Binning the same pitchers into disjoint axis bands shows
# whether chase actually falls off with axis distance, or whether one band is carrying everything.
cat("\n=== disjoint axis bands within eff >= .85 + arm slot >= 40 ===\n")
for (lg in c("D1","MLB")) { X <- if (lg == "D1") D else S
  P <- X[eff_ch >= EFF & eff_ff >= EFF & arm >= ARMMIN]
  B <- P[, .(n = .N, arms = uniqueN(id), chase_above = round(mean(y),2),
             se = round(sd(y)/sqrt(.N),2)),
         by = .(band = cut(axis, c(-1,10,20,30,45,360),
                           c("0-10","10-20","20-30","30-45","45+")))][order(band)]
  cat(sprintf("  %s:\n", lg)); print(B, row.names = FALSE)
  ct <- cor.test(-P$axis, P$y)
  cat(sprintf("    axis match vs chase above model across the whole base pool: r = %+.3f (p = %.3f, n = %d)\n",
              ct$estimate, ct$p.value, nrow(P))) }
