#!/usr/bin/env Rscript

# The .90/.90 + 44-degree bin, held fixed, tested across every changeup-swing gate.
#
# D1 changeup usage is thin - the median pitcher-season has 14 changeup swings - so the swing gate
# does most of the work in deciding how many pitcher-seasons survive to be tested. That creates a
# genuine tension rather than an obvious choice. A loose gate buys sample size but each pitcher's
# whiff residual is noisier, so the bin mean gets a bigger standard error per member; a tight gate
# buys precise residuals but leaves cells too small to test at all, which is what reduced the
# axis-gap row to six pitcher-seasons.
#
# Sweeping it separates those two effects from a real signal. If the bin genuinely overperforms, the
# estimate should stay on one side of zero as the gate moves and tighten as residuals get cleaner.
# If it is noise, the sign will wander and significance will track nothing but sample size.

suppressPackageStartupMessages({ library(data.table); library(mgcv) })
set.seed(23); options(width = 205)
MDIR <- "data/statcast_model"

EFF <- { a <- commandArgs(TRUE); if (length(a)) as.numeric(a[1]) else 0.90 }
ARM <- as.data.table(readRDS(file.path(MDIR, "ncaa_armangle.rds")))
PR  <- as.data.table(readRDS(file.path(MDIR, "ncaa_spineff_pairs.rds")))
PR  <- merge(PR, ARM[, .(PitcherId, season, arm_hat)], by = c("PitcherId","season"))
W   <- readRDS(file.path(MDIR, "ncaa_whiff_resid.rds")); setDT(W)

cat(sprintf("bin held at: CH and FF efficiency >= %.2f, estimated arm slot >= 44, axis gap <= 10\n", EFF))
cat(sprintf("%d pitcher-seasons have both an efficiency pair and an arm estimate\n\n", nrow(PR)))

res <- rbindlist(lapply(c(20, 30, 40, 50, 60, 80), function(g) {
  S <- W[!is.na(PitcherId), .(w = 100*mean(r_tj, na.rm = TRUE), nsw = .N),
         by = .(PitcherId, season)][nsw >= g]
  D <- merge(PR, S, by = c("PitcherId","season"))
  D[, bin  := eff_ch >= EFF & eff_ff >= EFF & arm_hat >= 44]
  D[, bina := bin & axis_gap <= 10]
  # The arm condition costs about three quarters of the members, so the axis bin is also shown
  # without it. If the axis gap is what matters, that version is the better-powered test of it.
  D[, binx := eff_ch >= EFF & eff_ff >= EFF & axis_gap <= 10]
  row <- function(lab, i) {
    x <- D$w[i]; y <- D$w[!i]
    if (length(x) < 3 || length(y) < 3)
      return(data.table(gate = g, cut = lab, n = length(x),
                        bin_w = if (length(x)) round(mean(x),2) else NA_real_,
                        rest_w = NA_real_, diff = NA_real_, p = NA_real_,
                        se = if (length(x) > 1) round(sd(x)/sqrt(length(x)),2) else NA_real_))
    t <- t.test(x, y)
    data.table(gate = g, cut = lab, n = length(x), bin_w = round(mean(x),2),
               rest_w = round(mean(y),2), diff = round(mean(x)-mean(y),2),
               p = round(t$p.value,3), se = round(sd(x)/sqrt(length(x)),2))
  }
  rbind(data.table(gate = g, cut = "qualifying pool", n = nrow(D), bin_w = round(mean(D$w),2),
                   rest_w = NA_real_, diff = NA_real_, p = NA_real_, se = NA_real_),
        row("eff + arm 44", D$bin),
        row("eff + arm 44 + axis<=10", D$bina),
        row("eff + axis<=10 (no arm)", D$binx))
}))

cat("=== the fixed bin across swing gates ===\n")
print(res[cut != "qualifying pool"], row.names = FALSE)
cat("\npool size at each gate, for reference:\n")
print(res[cut == "qualifying pool", .(gate, pitcher_seasons = n, mean_resid = bin_w)],
      row.names = FALSE)

cat("\n=== reading the sweep ===\n")
M <- res[cut == "eff + arm 44 + axis<=10" & is.finite(diff)]
cat(sprintf("  sign of the bin effect across gates: %s\n",
            paste(ifelse(M$diff > 0, "+", "-"), collapse = " ")))
cat(sprintf("  range %+.2f to %+.2f pp, smallest p = %.3f\n",
            min(M$diff), max(M$diff), min(M$p)))
cat(sprintf("  bin standard error runs %.2f to %.2f pp, so an effect under about %.1f pp\n",
            min(M$se), max(M$se), 2*max(M$se)))
cat("  could not be detected at any of these gates regardless of whether it exists.\n")

# What size of effect could this design actually find? Reporting a null without its detectable
# floor invites reading "no effect" where the honest statement is "no effect this large".
cat("\n=== minimum detectable effect (80% power, two-sided) ===\n")
for (g in c(20, 40, 60)) {
  r <- res[gate == g & cut == "eff + arm 44 + axis<=10"]
  if (!is.finite(r$se)) next
  cat(sprintf("  gate %2d: n = %3d, MDE = %.2f pp of whiff above model\n", g, r$n, 2.8*r$se))
}
