#!/usr/bin/env Rscript

# WOULD 2022 AND 2026 ACTUALLY CHANGE THE ANSWER?
#
# The bin at arm slot 40 sits at p = 0.19, which is the kind of number that invites "we just need
# more data". Sometimes that is true and sometimes the gap is far larger than two extra seasons can
# close, and the difference is worth knowing BEFORE paying for the data and the scraping that comes
# with it.
#
# The projection has three parts. How many bin seasons each year of tracking currently yields, which
# sets the per-season rate. What sample size the observed effect would need in order to clear .05,
# which is fixed by the effect size and the residual spread. And how those two numbers compare.
#
# One thing the arithmetic below cannot capture, and it cuts against the optimistic reading: the
# +1.42 being extrapolated was itself chosen after sweeping thresholds. Effects selected that way
# shrink when they are re-estimated on new data, so treating the current point estimate as the truth
# to power against is the most generous assumption available, not a neutral one.

suppressPackageStartupMessages(library(data.table)); options(width = 200)
MDIR <- "data/statcast_model"; GATE <- 40L; EFF <- .85; ARMMIN <- 40; AXIS <- 10

ARM <- as.data.table(readRDS(file.path(MDIR, "ncaa_armangle.rds")))
PR  <- as.data.table(readRDS(file.path(MDIR, "ncaa_spineff_pairs.rds")))
PR  <- merge(PR, ARM[, .(PitcherId, season, arm_hat)], by = c("PitcherId","season"))
setnames(PR, "axis_gap", "axis", skip_absent = TRUE)
CR  <- readRDS(file.path(MDIR, "ncaa_chase_resid.rds")); setDT(CR)
D <- merge(PR, CR[!is.na(PitcherId), .(y = 100*mean(r_all, na.rm = TRUE), n = .N),
                  by = .(PitcherId, season)], by = c("PitcherId","season"))[n >= GATE]
D[, `:=`(arm = arm_hat, inbin = eff_ch >= EFF & eff_ff >= EFF & arm_hat >= ARMMIN & axis <= AXIS)]

cat(sprintf("=== the bin at arm slot >= %d, axis gap <= %d, eff >= %.2f both ===\n", ARMMIN, AXIS, EFF))
B <- D[inbin == TRUE]; O <- D[inbin == FALSE]
tt <- t.test(B$y, O$y)
cat(sprintf("  n = %d seasons from %d arms | bin %+.2f, rest %+.2f | diff %+.2f, se %.2f, p = %.3f\n",
            nrow(B), uniqueN(B$PitcherId), mean(B$y), mean(O$y),
            mean(B$y)-mean(O$y), sd(B$y)/sqrt(nrow(B)), tt$p.value))
cat(sprintf("  95%% CI on the difference: %+.2f to %+.2f\n", tt$conf.int[1], tt$conf.int[2]))

cat("\n=== per-season yield, which sets the rate new years would arrive at ===\n")
print(D[, .(pool = .N, bin = sum(inbin), bin_arms = uniqueN(PitcherId[inbin]),
            bin_effect = round(mean(y[inbin]) - mean(y[!inbin]), 2)), by = season][order(season)],
      row.names = FALSE)
per_season <- nrow(B)/uniqueN(D$season)
cat(sprintf("  average %.1f bin seasons per year of tracking data\n", per_season))

## ---- what sample would .05 require? -------------------------------------------------------------
# The comparison group is two orders of magnitude larger, so its contribution to the standard error
# is negligible and the requirement is set almost entirely by the bin's own size.
sdb <- sd(B$y)
need <- function(eff, target_p = .05) {
  f <- function(n) abs(eff)/ (sdb/sqrt(n)) - qt(1 - target_p/2, n - 1)
  n <- 4; while (f(n) < 0 && n < 1e5) n <- n + 1; n }
cat(sprintf("\n=== sample required, bin residual sd = %.2f ===\n", sdb))
print(rbindlist(lapply(c(1.42, 1.20, 1.00, 0.80, 0.50), function(e)
  data.table(assumed_effect = e, n_for_p05 = need(e),
             years_of_data = round(need(e)/per_season, 1)))), row.names = FALSE)

## ---- what two more seasons actually buys ---------------------------------------------------------
cat("\n=== projecting 2022 and 2026 onto the current rate ===\n")
print(rbindlist(lapply(c(1.42, 1.20, 1.00, 0.80), function(e) {
  n5 <- round(per_season*5)
  se5 <- sdb/sqrt(n5); t5 <- e/se5
  data.table(assumed_effect = e, years = 5, projected_n = n5,
             projected_se = round(se5,2), projected_t = round(t5,2),
             projected_p = round(2*pt(-abs(t5), n5-1),3),
             power_at_05 = round(pnorm(abs(e)/se5 - 1.96),2)) })), row.names = FALSE)

cat("\n=== and if every remaining D1 season through 2030 were added ===\n")
print(rbindlist(lapply(c(9, 12), function(yr) {
  n <- round(per_season*yr); se <- sdb/sqrt(n)
  data.table(years = yr, projected_n = n, se = round(se,2),
             p_if_effect_is_142 = round(2*pt(-abs(1.42/se), n-1),4),
             p_if_effect_is_070 = round(2*pt(-abs(0.70/se), n-1),4)) })), row.names = FALSE)
