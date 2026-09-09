#!/usr/bin/env Rscript

# THE FALSIFICATION BATTERY.
#
# Three ways the +1.82 could be an artefact of the procedure rather than a property of the pitches,
# each with a prediction that does not depend on believing the mechanism.
#
# F1  Negative control. Rerun the identical machinery on spin RATE gap instead of spin AXIS gap.
#     Rate has no mimicry story - a hitter cannot see revolutions per minute - so if "efficiency
#     floors plus a high slot plus any narrow gap" produces a positive number, the finding is about
#     the construction and not about axis.
# F2  Out of sample. Everything so far used the same MLB seasons that the thresholds were found on.
#     Fit the arm-slot percentile and the pool on 2020-2023 and evaluate on 2024-2026.
# F3  Within pitcher. Some arms have seasons on both sides of the bin. If the effect is a property
#     of the pitch rather than the pitcher, those arms should whiff better in the seasons they meet
#     the definition than in the seasons they do not.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 200); set.seed(23)
MDIR <- "data/statcast_model"
LOCK <- readRDS(file.path(MDIR, "locked_spec.rds"))
SPEC <- LOCK$spec; L <- LOCK$data; setDT(L); L[, id := as.character(id)]
stopifnot(is.numeric(SPEC$arm_pctile), length(SPEC$arm_pctile) == 1)
gap <- function(A, B) c(nrow(A), uniqueN(A$id), mean(A$y) - mean(B$y),
                        sqrt(var(A$y)/nrow(A) + var(B$y)/nrow(B)))
show <- function(lab, e) cat(sprintf("  %-28s %3d seasons / %3d arms   %+6.2f +/- %.2f   p %.3f\n",
  lab, e[1], e[2], e[3], e[4], 2*pnorm(-abs(e[3]/e[4]))))

## ---- F1: spin rate gap as a negative control ---------------------------------------------------
SP <- rbindlist(lapply(2020:2026, function(y) {
  x <- fread(sprintf("data/savant/spin_%d.csv", y), showProgress = FALSE)
  setnames(x, names(x)[2], "pitcher"); x[, season := y]
  x[, .(id = as.character(pitcher), season, spin_ch = suppressWarnings(as.numeric(ch_avg_spin)),
        spin_ff = suppressWarnings(as.numeric(ff_avg_spin)))] }), use.names = TRUE)
PR <- as.data.table(readRDS(file.path(MDIR, "ncaa_spineff_pairs.rds")))
PR[, id := as.character(PitcherId)]
RG <- rbind(SP[, .(league = "MLB", id, season, spin_ch, spin_ff)],
            PR[, .(league = "D1", id, season, spin_ch, spin_ff)])
G <- merge(L, RG, by = c("league","id","season"))
G <- G[is.finite(spin_ch) & is.finite(spin_ff)]
G[, rate_gap := abs(spin_ch - spin_ff)]
# The axis cut keeps roughly the tightest sixth of the base pool, so the rate cut is set to the same
# quantile. Matching the selectivity is what makes it a control rather than a different test.
G[, base := ec >= .85 & ef >= .85 & arm >= arm_thr]
q <- G[base == TRUE, mean(axis <= 10)]
G[, rate_thr := as.numeric(quantile(rate_gap[base == TRUE], q, na.rm = TRUE)), by = league]
cat(sprintf("=== F1: negative control on spin RATE gap (cut at the same %.0f%% selectivity as axis) ===\n",
            100*q))
for (lgv in c("MLB","D1")) {
  Z <- G[league == lgv]; B <- Z[base == FALSE]
  cat(sprintf("\n%s   rate-gap threshold %.0f rpm   (axis threshold 10 deg)\n", lgv, Z$rate_thr[1]))
  show("axis gap <= 10 (the claim)", gap(Z[base == TRUE & axis <= 10], B))
  show("rate gap tight (the control)", gap(Z[base == TRUE & rate_gap <= rate_thr], B))
}
est <- G[, { Z <- .SD; B <- Z[base == FALSE]
  a <- gap(Z[base == TRUE & axis <= 10], B); r <- gap(Z[base == TRUE & rate_gap <= rate_thr], B)
  .(a_d = a[3], a_se = a[4], r_d = r[3], r_se = r[4]) }, by = league]
pw <- function(d, s) { w <- 1/s^2; m <- sum(d*w)/sum(w); se <- sqrt(1/sum(w))
  sprintf("%+.2f +/- %.2f (p %.3f)", m, se, 2*pnorm(-abs(m/se))) }
cat(sprintf("\n  pooled, axis    : %s\n  pooled, rate    : %s\n",
            pw(est$a_d, est$a_se), pw(est$r_d, est$r_se)))
cat("  the control must be near zero. If it is not, the bin is selecting on something other than\n")
cat("  axis - most likely pitchers whose changeup and fastball are alike in every way at once.\n")

## ---- F2: out of sample -------------------------------------------------------------------------
# The arm-slot percentile is recomputed inside the training era only, so nothing about the later
# seasons enters the definition. MLB alone, since D1 has only three seasons.
cat("\n=== F2: MLB fit on 2020-2023, evaluated on 2024-2026 ===\n")
M <- L[league == "MLB"]
for (era in list(c(2020,2023), c(2024,2026))) {
  E <- M[season >= era[1] & season <= era[2]]
  thr <- as.numeric(quantile(M[season <= 2023]$arm, SPEC$arm_pctile))
  E[, b := ec >= .85 & ef >= .85 & axis <= 10 & arm >= thr]
  show(sprintf("%d-%d (slot >= %.1f from train)", era[1], era[2], thr),
       gap(E[b == TRUE], E[b == FALSE]))
}

## ---- F3: within pitcher ------------------------------------------------------------------------
# The comparison that no amount of "these arms are simply good" can explain, because the arm is held
# fixed. Only pitchers observed on both sides contribute.
cat("\n=== F3: arms with seasons both inside and outside the bin ===\n")
for (lgv in c("MLB","D1")) {
  Z <- L[league == lgv]
  sw <- Z[, .(nin = sum(bin), nout = sum(!bin)), by = id][nin > 0 & nout > 0]
  W <- Z[id %in% sw$id]
  d <- W[, .(din = mean(y[bin]) - mean(y[!bin])), by = id]
  t <- if (nrow(d) >= 3) t.test(d$din) else NULL
  cat(sprintf("  %-4s %2d arms qualify (%d in-bin seasons, %d out)  within-pitcher %+6.2f +/- %.2f  p %.3f\n",
              lgv, nrow(d), sum(W$bin), sum(!W$bin), mean(d$din),
              if (is.null(t)) NA_real_ else sd(d$din)/sqrt(nrow(d)),
              if (is.null(t)) NA_real_ else t$p.value))
  if (lgv == "MLB") {
    FC <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(FC)
    NM <- unique(FC[pitch_type == "CH", .(id = as.character(pitcher), player_name)])
    print(merge(d, NM, by = "id")[order(-din), .(player_name, within_diff = round(din,2))],
          row.names = FALSE)
  }
}
