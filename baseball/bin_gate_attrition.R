#!/usr/bin/env Rscript

# WHAT IS ACTUALLY LIMITING THE BIN TO FOURTEEN SEASONS - THE GATE, OR THE CONDITIONS?
#
# Two different things can starve a bin, and they call for opposite responses. If the residual gate
# is the binding constraint then loosening it recovers seasons cheaply, at the cost of noisier
# per-season estimates. If the bin conditions are what is cutting, then no gate setting helps and
# the only route to more evidence is more seasons of tracking data.
#
# The attrition is therefore walked twice: the gate swept with the conditions held fixed, and the
# conditions applied one at a time with the gate held at 40.

suppressPackageStartupMessages(library(data.table)); options(width = 200)
MDIR <- "data/statcast_model"

ARM <- as.data.table(readRDS(file.path(MDIR, "ncaa_armangle.rds")))
PR  <- as.data.table(readRDS(file.path(MDIR, "ncaa_spineff_pairs.rds")))
PR  <- merge(PR, ARM[, .(PitcherId, season, arm_hat)], by = c("PitcherId","season"))
setnames(PR, "axis_gap", "axis", skip_absent = TRUE)
CR  <- readRDS(file.path(MDIR, "ncaa_chase_resid.rds")); setDT(CR)
CNT <- CR[!is.na(PitcherId), .(y = 100*mean(r_all, na.rm = TRUE), nooz = .N),
          by = .(PitcherId, season)]
D0 <- merge(PR, CNT, by = c("PitcherId","season"))
cat(sprintf("D1 pitcher-seasons with a screened CH/FF pair AND any out-of-zone changeups: %d\n",
            nrow(D0)))
cat(sprintf("out-of-zone changeup counts: median %d, quartiles %d / %d, max %d\n",
            median(D0$nooz), quantile(D0$nooz,.25), quantile(D0$nooz,.75), max(D0$nooz)))

cond <- function(X) list(
  full = X[, eff_ch >= .85 & eff_ff >= .85 & arm_hat >= 44 & axis <= 10],
  noax = X[, eff_ch >= .85 & eff_ff >= .85 & arm_hat >= 44],
  noarm = X[, eff_ch >= .85 & eff_ff >= .85 & axis <= 10])

cat("\n=== gate swept, bin conditions held fixed ===\n")
print(rbindlist(lapply(c(0, 10, 15, 25, 40, 60, 100, 150), function(g) {
  X <- D0[nooz >= g]; C <- cond(X)
  t <- if (sum(C$full) >= 3 && sum(!C$full) >= 3) t.test(X$y[C$full], X$y[!C$full]) else NULL
  data.table(gate = g, pool = nrow(X), full_bin = sum(C$full),
             arms = uniqueN(X$PitcherId[C$full]),
             median_ooz_in_bin = if (sum(C$full)) median(X$nooz[C$full]) else NA_integer_,
             effect = if (!is.null(t)) round(diff(rev(t$estimate)),2) else NA_real_,
             p = if (!is.null(t)) round(t$p.value,3) else NA_real_) })), row.names = FALSE)

cat("\n=== conditions applied one at a time, gate held at 40 ===\n")
X <- D0[nooz >= 40]
step <- function(l, i) data.table(step = l, seasons = sum(i),
                                  pct_of_pool = round(100*mean(i),1),
                                  arms = uniqueN(X$PitcherId[i]))
print(rbind(step("pool (gate 40)",        rep(TRUE, nrow(X))),
            step("+ CH efficiency >= .85", X[, eff_ch >= .85]),
            step("+ FF efficiency >= .85", X[, eff_ch >= .85 & eff_ff >= .85]),
            step("+ arm slot >= 44",       X[, eff_ch >= .85 & eff_ff >= .85 & arm_hat >= 44]),
            step("+ axis gap <= 10",       X[, eff_ch >= .85 & eff_ff >= .85 & arm_hat >= 44 & axis <= 10])),
      row.names = FALSE)

# Which single condition is the scarce one? Each is applied alone so the orderings above do not
# disguise a condition that only looks cheap because an earlier one already removed its casualties.
cat("\n=== each condition alone, gate 40 ===\n")
print(rbind(step("CH eff >= .85 alone",  X[, eff_ch >= .85]),
            step("FF eff >= .85 alone",  X[, eff_ff >= .85]),
            step("arm slot >= 44 alone", X[, arm_hat >= 44]),
            step("axis gap <= 10 alone", X[, axis <= 10])), row.names = FALSE)
cat(sprintf("\naxis gap: median %.1f deg, only %.1f%% of the gate-40 pool sits at or under 10\n",
            median(X$axis), 100*mean(X$axis <= 10)))
