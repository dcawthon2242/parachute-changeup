#!/usr/bin/env Rscript

# VELOCITY SEPARATION AND THE WHIFF RESIDUAL.
#
# speed_diff is a feature in both whiff models, so across the full population the residual should be
# flat in velocity by construction. That flat baseline is the control: any slope that appears inside
# a bin is a conditional effect, meaning velocity pays off differently for those pitchers than the
# model's global fit expects.
#
# Three questions, in the order they should be asked. Whether the baseline really is flat, which
# validates the whole exercise. Whether velocity correlates with the residual inside the bin.
# Whether gating on velocity improves the one cell that survived - axis gap under 10 with a
# top-third arm slot - or merely shrinks it.
#
# On the third: a gate that "improves" a 15-season cell by cutting it to 8 is not an improvement,
# it is a smaller sample with a larger point estimate. The standard error is reported alongside so
# that trade is visible rather than hidden behind a p-value.

suppressPackageStartupMessages(library(data.table)); options(width = 200)
MDIR <- "data/statcast_model"; GATE <- 40L; EFF <- .85

AS <- readRDS(file.path(MDIR, "active_spin_long.rds")); setDT(AS)
F  <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(F); F <- F[pitch_type == "CH"]
RW <- readRDS(file.path(MDIR, "mlb_whiff_locaware.rds")); setDT(RW)
S <- F[, .(axis = mean(axis_diff, na.rm = TRUE), arm = mean(arm_angle, na.rm = TRUE),
           vs = -mean(speed_diff, na.rm = TRUE)), by = .(pitcher, season)]
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, ec = active_spin)], by = c("pitcher","season"))
S <- merge(S, AS[pitch_type == "FF", .(pitcher, season, ef = active_spin)], by = c("pitcher","season"))
S <- merge(S, RW[, .(y = 100*mean(r_all), n = .N), by = .(pitcher, season)], by = c("pitcher","season"))
S <- S[n >= GATE & is.finite(arm), .(id = pitcher, y, axis, arm, ec, ef, vs)]

ARM <- as.data.table(readRDS(file.path(MDIR, "ncaa_armangle.rds")))
PR  <- as.data.table(readRDS(file.path(MDIR, "ncaa_spineff_pairs.rds")))
PR  <- merge(PR, ARM[, .(PitcherId, season, arm_hat)], by = c("PitcherId","season"))
setnames(PR, "axis_gap", "axis", skip_absent = TRUE)
DW  <- readRDS(file.path(MDIR, "ncaa_whiff_resid.rds")); setDT(DW)
D <- merge(PR, DW[!is.na(PitcherId), .(y = 100*mean(r_all, na.rm = TRUE), n = .N),
                  by = .(PitcherId, season)], by = c("PitcherId","season"))[n >= GATE]
D <- D[, .(id = PitcherId, y, axis, arm = arm_hat, ec = eff_ch, ef = eff_ff, vs = velo_sep)]

cat("=== 1. is the baseline flat, as it must be if speed_diff is properly absorbed? ===\n")
for (lg in c("D1","MLB")) { X <- if (lg == "D1") D else S; ct <- cor.test(X$vs, X$y)
  cat(sprintf("  %-3s whole pool (n = %d): r = %+.3f (p = %.3f), slope %+.3f whiff pts per mph\n",
              lg, nrow(X), ct$estimate, ct$p.value, coef(lm(y ~ vs, X))["vs"])) }

cat("\n=== 2. velocity against whiff above model, nested inside each cell ===\n")
print(rbindlist(lapply(c("D1","MLB"), function(lg) { X <- if (lg == "D1") D else S
  thr <- as.numeric(quantile(X$arm, 2/3))
  rbindlist(lapply(list(
      list("eff .85 both",              quote(ec >= EFF & ef >= EFF)),
      list("+ axis <= 10",              quote(ec >= EFF & ef >= EFF & axis <= 10)),
      list("+ axis <= 10, top third",   quote(ec >= EFF & ef >= EFF & axis <= 10 & arm >= thr)),
      list("+ axis <= 15, top third",   quote(ec >= EFF & ef >= EFF & axis <= 15 & arm >= thr))),
    function(r) { A <- X[eval(r[[2]])]
      if (nrow(A) < 5) return(NULL); ct <- cor.test(A$vs, A$y)
      data.table(league = lg, cell = r[[1]], n = nrow(A), r = round(ct$estimate,3),
                 slope = round(coef(lm(y ~ vs, A))["vs"],3), p = round(ct$p.value,3)) })) })),
  row.names = FALSE)

cat("\n=== 3. does a velocity gate improve the surviving cell? ===\n")
print(rbindlist(lapply(c("D1","MLB"), function(lg) { X <- if (lg == "D1") D else S
  thr <- as.numeric(quantile(X$arm, 2/3))
  base <- X[, ec >= EFF & ef >= EFF & axis <= 10 & arm >= thr]
  rbindlist(lapply(list(list("no velo gate", NULL), list("velo sep >= 7", 7), list("velo sep >= 8", 8),
                        list("velo sep >= 9", 9), list("velo sep <= 8", -8), list("velo sep <= 7", -7)),
    function(g) {
      i <- if (is.null(g[[2]])) base else if (g[[2]] > 0) base & X$vs >= g[[2]] else base & X$vs <= -g[[2]]
      if (sum(i) < 4) return(data.table(league = lg, gate = g[[1]], n = sum(i), arms = uniqueN(X$id[i]),
                                        diff = NA_real_, se = NA_real_, p = NA_real_))
      t <- t.test(X$y[i], X$y[!i])
      data.table(league = lg, gate = g[[1]], n = sum(i), arms = uniqueN(X$id[i]),
                 diff = round(mean(X$y[i]) - mean(X$y[!i]), 2),
                 se = round(sd(X$y[i])/sqrt(sum(i)), 2), p = round(t$p.value, 3)) })) })),
  row.names = FALSE)
