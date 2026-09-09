#!/usr/bin/env Rscript

# WHAT DOES THE D1 ARM-SLOT DISTRIBUTION ACTUALLY LOOK LIKE?
#
# One caution has to come first, because it changes how a percentile threshold should be read. D1
# arm slot is not measured. It is the output of a height-informed proxy calibrated on MLB, and a
# regression prediction is always narrower than the thing it predicts - the proxy explains part of
# the variance and the rest shows up as shrinkage toward the mean. So the D1 spread below is
# compressed relative to true arm angles, and a cut at the estimated 67th percentile admits a
# different set than a cut at the true 67th percentile would.
#
# That makes the percentile framing better than the fixed-degree one regardless. A fixed 44 degrees
# means something different in a compressed distribution than in a real one; a percentile is at
# least self-consistent within whichever distribution it is applied to.

suppressPackageStartupMessages(library(data.table)); options(width = 200)
MDIR <- "data/statcast_model"; GATE <- 40L

ARM <- as.data.table(readRDS(file.path(MDIR, "ncaa_armangle.rds")))
PR  <- as.data.table(readRDS(file.path(MDIR, "ncaa_spineff_pairs.rds")))
PR  <- merge(PR, ARM[, .(PitcherId, season, arm_hat)], by = c("PitcherId","season"))
setnames(PR, "axis_gap", "axis", skip_absent = TRUE)
CR  <- readRDS(file.path(MDIR, "ncaa_chase_resid.rds")); setDT(CR)
D <- merge(PR, CR[!is.na(PitcherId), .(y = 100*mean(r_all, na.rm = TRUE), n = .N),
                  by = .(PitcherId, season)], by = c("PitcherId","season"))[n >= GATE][, arm := arm_hat]

F  <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(F); F <- F[pitch_type == "CH"]
AS <- readRDS(file.path(MDIR, "active_spin_long.rds")); setDT(AS)
R  <- readRDS(file.path(MDIR, "mlb_chase_resid.rds")); setDT(R)
S <- F[, .(axis = mean(axis_diff, na.rm = TRUE), arm = mean(arm_angle, na.rm = TRUE)),
       by = .(pitcher, season)]
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, eff_ch = active_spin)], by = c("pitcher","season"))
S <- merge(S, AS[pitch_type == "FF", .(pitcher, season, eff_ff = active_spin)], by = c("pitcher","season"))
S <- merge(S, R[, .(y = 100*mean(r_all), n = .N), by = .(pitcher, season)], by = c("pitcher","season"))
S <- S[n >= GATE & is.finite(arm)]

cat("=== arm slot distributions, gate-40 pools ===\n")
desc <- function(v, lab) data.table(pool = lab, n = length(v), mean = round(mean(v),1),
  sd = round(sd(v),1), min = round(min(v),1), p10 = round(quantile(v,.10),1),
  p25 = round(quantile(v,.25),1), median = round(median(v),1), p67 = round(quantile(v,.667),1),
  p75 = round(quantile(v,.75),1), p90 = round(quantile(v,.90),1), max = round(max(v),1))
print(rbind(desc(D$arm, "D1 (estimated)"), desc(S$arm, "MLB (measured)")), row.names = FALSE)
cat(sprintf("\n  the proxy is compressed as expected: D1 sd %.1f vs MLB measured sd %.1f (%.0f%% as wide)\n",
            sd(D$arm), sd(S$arm), 100*sd(D$arm)/sd(S$arm)))
cat(sprintf("  a fixed 44-degree cut takes the top %.0f%% of D1 but the top %.0f%% of MLB\n",
            100*mean(D$arm >= 44), 100*mean(S$arm >= 44)))

cat("\n=== D1 deciles ===\n")
print(data.table(decile = paste0("p", seq(10,90,10)),
                 arm_slot = round(as.numeric(quantile(D$arm, seq(.1,.9,.1))),1)), row.names = FALSE)

## ---- the bin at percentile cuts rather than fixed degrees --------------------------------------
cat("\n=== eff >= .85 both + axis gap <= 10, slot cut by PERCENTILE of each league's own pool ===\n")
run <- function(X, lg) rbindlist(lapply(c(0, .50, .60, .667, .75, .80), function(q) {
  thr <- as.numeric(quantile(X$arm, q))
  i <- X[, eff_ch >= .85 & eff_ff >= .85 & axis <= 10 & arm >= thr]
  x <- X$y[i]; y <- X$y[!i]
  if (length(x) < 3) return(NULL)
  t <- t.test(x, y)
  id <- if (lg == "D1") X$PitcherId[i] else X$pitcher[i]
  data.table(league = lg, slot_pctile = sprintf("top %.0f%%", 100*(1-q)),
             threshold_deg = round(thr,1), n = length(x), arms = uniqueN(id),
             bin = round(mean(x),2), diff = round(mean(x)-mean(y),2),
             p = round(t$p.value,4), se = round(sd(x)/sqrt(length(x)),2)) }))
print(rbind(run(D, "D1"), run(S, "MLB")), row.names = FALSE)
