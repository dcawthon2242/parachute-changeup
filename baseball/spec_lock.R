#!/usr/bin/env Rscript

# THE FROZEN SPECIFICATION.
#
# Everything downstream reads its bin definition from this file and nothing else. The point is to
# make further threshold searching impossible by construction: the numbers below were arrived at by
# sweeping, that sweeping is already spent, and any additional tuning would be fitting noise that
# has already been fitted once.
#
# THE LOCKED CELL
#   changeup active spin  >= .85
#   four-seam active spin >= .85
#   spin axis gap         <= 10 degrees
#   arm slot              >= each league's own 67th percentile
#   qualifying            >= 40 changeup swings in the pitcher-season
#
# Arm slot is a percentile rather than a fixed angle because the two distributions are not on the
# same scale. D1's is a regression estimate and therefore compressed - its standard deviation is 68
# percent of MLB's measured spread - so a fixed 44 degrees selects the top 17 percent of D1 against
# the top 27 percent of MLB. A percentile at least asks the same question of both.
#
# The residual is season-level velocity-adjusted in both leagues. D1 needed it (the raw season-level
# residual correlated with mean velocity separation at r = +.062, p = .014, despite the pitch-level
# correlation being +.0004 - an aggregation effect, not model error). MLB does not need it, but it
# is applied there too so the two sides are treated identically.

# bit64 is required, not optional. The NCAA PitcherId is stored as integer64, and as.character on an
# integer64 without bit64 attached reinterprets the underlying bits as a double and returns strings
# like "3.37562447470698e-318". That mapping is a bijection, so merges inside a single script still
# work and nothing looks wrong - but a second script that happens to have bit64 loaded produces the
# real digits instead, and the two conventions silently fail to join.
suppressPackageStartupMessages({ library(data.table); library(bit64) }); options(width = 200)
set.seed(11)
MDIR <- "data/statcast_model"
SPEC <- list(eff_min = .85, axis_max = 10, arm_pctile = 2/3, swing_gate = 40L,
             require_ff_primary = TRUE)

## ---- MLB ---------------------------------------------------------------------------------------
AS <- readRDS(file.path(MDIR, "active_spin_long.rds")); setDT(AS)
F  <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(F); FC <- F[pitch_type == "CH"]
RW <- readRDS(file.path(MDIR, "mlb_whiff_locaware.rds")); setDT(RW)
S <- FC[, .(axis = mean(axis_diff, na.rm = TRUE), arm = mean(arm_angle, na.rm = TRUE),
            vs = -mean(speed_diff, na.rm = TRUE)), by = .(pitcher, season)]
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, ec = active_spin)], by = c("pitcher","season"))
S <- merge(S, AS[pitch_type == "FF", .(pitcher, season, ef = active_spin)], by = c("pitcher","season"))
S <- merge(S, RW[, .(y_raw = 100*mean(r_all), nsw = .N), by = .(pitcher, season)],
           by = c("pitcher","season"))
S <- S[nsw >= SPEC$swing_gate & is.finite(arm)]
S[, `:=`(league = "MLB", id = as.character(pitcher))]

## ---- D1 ----------------------------------------------------------------------------------------
ARM <- as.data.table(readRDS(file.path(MDIR, "ncaa_armangle.rds")))
PR  <- as.data.table(readRDS(file.path(MDIR, "ncaa_spineff_pairs.rds")))
PR  <- merge(PR, ARM[, .(PitcherId, season, arm_hat)], by = c("PitcherId","season"))
setnames(PR, "axis_gap", "axis", skip_absent = TRUE)
PR[, PitcherId := as.character(PitcherId)]
DR <- readRDS(file.path(MDIR, "ncaa_whiff_resid.rds")); setDT(DR)
DR[, PitcherId := as.character(PitcherId)]
DS <- DR[, .(y_raw = 100*mean(r_all, na.rm = TRUE), nsw = .N, vs = -mean(speed_diff, na.rm = TRUE)),
         by = .(PitcherId, season)][nsw >= SPEC$swing_gate]
D <- merge(PR[, .(PitcherId, season, axis, arm = arm_hat, ec = eff_ch, ef = eff_ff, name)],
           DS, by = c("PitcherId","season"))
D[, `:=`(league = "D1", id = PitcherId)]

## ---- assemble, attach usage, filter ------------------------------------------------------------
COLS <- c("league","id","season","y_raw","nsw","vs","axis","arm","ec","ef")
P <- rbind(S[, ..COLS], D[, ..COLS])
n_pre <- P[, .N, by = league]

# Usage is required, not optional. A blank Savant cell means the pitcher did not throw that pitch
# and is treated as zero. A missing arsenal *row* is dropped. Imputing would quietly put
# sinker-primary arms back into the pool.
AR <- rbindlist(lapply(2020:2026, function(y) {
  x <- fread(sprintf("data/savant/arsenal_%d.csv", y), showProgress = FALSE)
  setnames(x, names(x)[2], "pitcher"); x[, season := y]
  x[, .(id = as.character(pitcher), season,
        ff_use = suppressWarnings(as.numeric(n_ff)),
        si_use = suppressWarnings(as.numeric(n_si)))] }), use.names = TRUE)
AR[is.na(ff_use), ff_use := 0][is.na(si_use), si_use := 0]
UD <- file.path(MDIR, "ncaa_usage.rds")
if (!file.exists(UD)) stop("ncaa_usage.rds is missing; build it from mech_01_usage.R before locking")
UN <- readRDS(UD); setDT(UN)
UN[, id := as.character(PitcherId)]
U <- rbind(AR[, .(league = "MLB", id, season, ff_use, si_use)],
           UN[, .(league = "D1", id, season, ff_use, si_use)])
P <- merge(P, U, by = c("league","id","season"), all.x = TRUE)
n_miss <- P[is.na(ff_use) | is.na(si_use), .N]
if (n_miss > 0) cat(sprintf("dropping %d pitcher-seasons with missing usage (fail closed)\n", n_miss))
P <- P[is.finite(ff_use) & is.finite(si_use)]
if (isTRUE(SPEC$require_ff_primary)) {
  n_si <- P[si_use > ff_use, .N]
  P <- P[ff_use >= si_use]
  cat(sprintf("four-seam-primary filter: dropped %d sinker-primary seasons\n", n_si))
}
cat("\nattrition, before vs after the population filter:\n")
print(merge(n_pre[, .(league, before = N)], P[, .(after = .N), by = league], by = "league"),
      row.names = FALSE)

# Season-level velocity adjustment, fit within league so one league's slope cannot contaminate the
# other's residual. Arm percentile is cut AFTER the filter so the top third is the
# four-seam-primary population, not the old mixed one.
P[, y := resid(lm(y_raw ~ vs)), by = league]
P[, arm_thr := as.numeric(quantile(arm, SPEC$arm_pctile)), by = league]
P[, bin := ec >= SPEC$eff_min & ef >= SPEC$eff_min & axis <= SPEC$axis_max & arm >= arm_thr]
saveRDS(list(spec = SPEC, data = P), file.path(MDIR, "locked_spec.rds"))

cat("\n=== the locked specification ===\n")
cat(sprintf("  active spin >= %.2f on both pitches | axis gap <= %d deg | arm slot >= league p%.0f | %d+ swings\n",
            SPEC$eff_min, SPEC$axis_max, 100*SPEC$arm_pctile, SPEC$swing_gate))
cat(sprintf("  population: four-seam usage >= sinker usage (missing usage dropped)\n"))
print(P[, .(pool = .N, arms = uniqueN(id), arm_threshold = round(arm_thr[1],1),
            bin_seasons = sum(bin), bin_arms = uniqueN(id[bin])), by = league], row.names = FALSE)

## ---- per-league effect and the pooled estimate --------------------------------------------------
est <- P[, { t <- t.test(y[bin], y[!bin]); d <- as.numeric(diff(rev(t$estimate)))
             .(n = sum(bin), arms = uniqueN(id[bin]), diff = d, se = as.numeric(d/t$statistic),
               p = t$p.value) }, by = league]
cat("\n=== per league, on the velocity-adjusted residual ===\n")
print(est[, .(league, n, arms, diff = round(diff,2), se = round(se,2), p = round(p,4))],
      row.names = FALSE)

w <- 1/est$se^2; m <- sum(est$diff*w)/sum(w); se <- sqrt(1/sum(w)); z <- m/se
Q <- sum(w*(est$diff - m)^2)
cat(sprintf("\n=== pooled, inverse-variance ===\n  %+.2f +/- %.2f | z = %.2f | p = %.4f | 95%% CI %+.2f to %+.2f\n",
            m, se, z, 2*pnorm(-abs(z)), m - 1.96*se, m + 1.96*se))
cat(sprintf("  heterogeneity Q = %.2f on 1 df (p = %.3f)%s\n", Q, 1 - pchisq(Q, 1),
            if (1 - pchisq(Q,1) > .1) " - the leagues are statistically compatible" else " - the leagues disagree"))

# Arm-clustered SE is the number that decides. One pitcher with four qualifying seasons is one
# observation, not four. Season-level p is kept as the secondary figure.
cat("\n=== pooled, one row per arm (the deciding number) ===\n")
A <- P[, .(y = mean(y), bin = any(bin), n = .N), by = .(league, id)]
estA <- A[, { t <- t.test(y[bin], y[!bin]); d <- as.numeric(diff(rev(t$estimate)))
              .(n = sum(bin), diff = d, se = as.numeric(d/t$statistic), p = t$p.value) },
          by = league]
print(estA[, .(league, n, diff = round(diff,2), se = round(se,2), p = round(p,4))], row.names = FALSE)
wA <- 1/estA$se^2; mA <- sum(estA$diff*wA)/sum(wA); seA <- sqrt(1/sum(wA)); zA <- mA/seA
cat(sprintf("  pooled arms: %+.2f +/- %.2f | z = %.2f | p = %.4f | 95%% CI %+.2f to %+.2f\n",
            mA, seA, zA, 2*pnorm(-abs(zA)), mA - 1.96*seA, mA + 1.96*seA))

## ---- permutation: how often does a bin this size look this good in both leagues at once? --------
# Shuffling the residual within league preserves every bin size and every marginal distribution and
# destroys only the association. Requiring BOTH leagues to clear the observed effect is the joint
# null the pooled estimate is really being tested against.
perm <- replicate(4000, {
  Pp <- copy(P)[, y := sample(y), by = league]
  e <- Pp[, .(d = mean(y[bin]) - mean(y[!bin])), by = league]
  ep <- merge(e, est[, .(league, obs = diff)], by = "league")
  c(both = all(ep$d >= ep$obs), pooled = sum(ep$d*w)/sum(w)) })
cat(sprintf("\n=== permutation, 4000 shuffles within league ===\n  both leagues at or above observed: p = %.4f\n",
            mean(perm["both",] > 0)))
cat(sprintf("  pooled estimate at or above observed:    p = %.4f (null mean %+.3f, sd %.3f)\n",
            mean(perm["pooled",] >= m), mean(perm["pooled",]), sd(perm["pooled",])))
cat(sprintf("\nwrote %s\n", file.path(MDIR, "locked_spec.rds")))
