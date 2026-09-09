#!/usr/bin/env Rscript

# ADDING A VELOCITY GATE - AND KEEPING HONEST ABOUT WHAT THAT COSTS.
#
# Searching a grid for a cell where two leagues agree will always find one. With a dozen cells per
# league and noise on both sides, agreement in sign is roughly a coin flip per cell, so several
# matches are the EXPECTED outcome under a pure null. The grid below is therefore reported together
# with the number of cells tested and a permutation estimate of how often this much agreement shows
# up when the residuals are shuffled.
#
# There is one genuine prior to lean on, which separates this from pure fishing. The velocity slope
# inside the bin came back negative in both leagues independently - D1 at r = -.19 and MLB at -.33,
# neither individually significant but agreeing in sign and surviving a jackknife. That predicts the
# effect should concentrate at LOW velocity separation. A hit in the low band is therefore a
# confirmation of a stated direction; a hit anywhere else is a grid artifact and should be read that
# way even if its p-value looks better.

suppressPackageStartupMessages(library(data.table)); options(width = 210)
set.seed(4); MDIR <- "data/statcast_model"; GATE <- 40L; EFF <- .85

ARM <- as.data.table(readRDS(file.path(MDIR, "ncaa_armangle.rds")))
PR  <- as.data.table(readRDS(file.path(MDIR, "ncaa_spineff_pairs.rds")))
PR  <- merge(PR, ARM[, .(PitcherId, season, arm_hat)], by = c("PitcherId","season"))
setnames(PR, "axis_gap", "axis", skip_absent = TRUE)
CR  <- readRDS(file.path(MDIR, "ncaa_chase_resid.rds")); setDT(CR)
D <- merge(PR, CR[!is.na(PitcherId), .(y = 100*mean(r_all, na.rm = TRUE), n = .N),
                  by = .(PitcherId, season)], by = c("PitcherId","season"))[n >= GATE]
D <- D[, .(id = PitcherId, season, y, axis, arm = arm_hat, ec = eff_ch, ef = eff_ff, vs = velo_sep)]

F  <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(F); F <- F[pitch_type == "CH"]
AS <- readRDS(file.path(MDIR, "active_spin_long.rds")); setDT(AS)
R  <- readRDS(file.path(MDIR, "mlb_chase_resid.rds")); setDT(R)
S <- F[, .(axis = mean(axis_diff, na.rm = TRUE), arm = mean(arm_angle, na.rm = TRUE),
           vs = -mean(speed_diff, na.rm = TRUE)), by = .(pitcher, season)]
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, ec = active_spin)], by = c("pitcher","season"))
S <- merge(S, AS[pitch_type == "FF", .(pitcher, season, ef = active_spin)], by = c("pitcher","season"))
S <- merge(S, R[, .(y = 100*mean(r_all), n = .N), by = .(pitcher, season)], by = c("pitcher","season"))
S <- S[n >= GATE & is.finite(arm), .(id = pitcher, season, y, axis, arm, ec, ef, vs)]

# Velocity separation is banded by each league's own tertiles, since D1 separates less than MLB.
for (X in list(D, S)) X[, vband := cut(vs, quantile(vs, c(0,1/3,2/3,1)), c("low","mid","high"),
                                       include.lowest = TRUE)]
D[, thr := as.numeric(quantile(arm, 2/3))]; S[, thr := as.numeric(quantile(arm, 2/3))]
cat(sprintf("velo separation tertile cuts - D1: %.1f / %.1f mph, MLB: %.1f / %.1f mph\n",
            quantile(D$vs,1/3), quantile(D$vs,2/3), quantile(S$vs,1/3), quantile(S$vs,2/3)))

cell <- function(X, ax, slot, vb) {
  i <- X[, ec >= EFF & ef >= EFF & axis <= ax & vband == vb & (if (slot) arm >= thr else TRUE)]
  x <- X$y[i]; y <- X$y[!i]
  if (length(x) < 4) return(list(n = length(x), d = NA_real_, p = NA_real_))
  list(n = length(x), d = mean(x)-mean(y), p = t.test(x, y)$p.value, arms = uniqueN(X$id[i]))
}
GRID <- CJ(ax = c(10, 15, 20), slot = c(FALSE, TRUE), vb = c("low","mid","high"), sorted = FALSE)
OUT <- rbindlist(lapply(seq_len(nrow(GRID)), function(k) {
  g <- GRID[k]; a <- cell(D, g$ax, g$slot, g$vb); b <- cell(S, g$ax, g$slot, g$vb)
  data.table(axis = g$ax, slot = ifelse(g$slot, "top third", "any"), velo = g$vb,
             D1_n = a$n, D1_diff = round(a$d,2), D1_p = round(a$p,3),
             MLB_n = b$n, MLB_diff = round(b$d,2), MLB_p = round(b$p,3),
             agree = !is.na(a$d) && !is.na(b$d) && sign(a$d) == sign(b$d) && a$d > 0) }))
cat("\n=== the grid: both leagues, same cell definitions ===\n")
print(OUT, row.names = FALSE)

ok <- OUT[!is.na(D1_diff) & !is.na(MLB_diff)]
cat(sprintf("\ncells evaluated: %d | both positive: %d | both positive AND both p<.10: %d\n",
            nrow(ok), sum(ok$agree), ok[agree == TRUE & D1_p < .10 & MLB_p < .10, .N]))

## ---- how much agreement does noise produce? -----------------------------------------------------
# The residuals are shuffled within each league, destroying any real association while preserving
# the bin sizes and the marginal spread, and the grid is rebuilt. The distribution of "both
# positive" counts under that null is what the observed count has to beat.
perm <- replicate(600, {
  Dp <- copy(D)[, y := sample(y)]; Sp <- copy(S)[, y := sample(y)]
  sum(sapply(seq_len(nrow(GRID)), function(k) {
    g <- GRID[k]; a <- cell(Dp, g$ax, g$slot, g$vb); b <- cell(Sp, g$ax, g$slot, g$vb)
    !is.na(a$d) && !is.na(b$d) && sign(a$d) == sign(b$d) && a$d > 0 })) })
cat(sprintf("under 600 shuffles: expected %.1f agreeing cells (95%% of draws fall at or below %d)\n",
            mean(perm), quantile(perm, .95)))
cat(sprintf("observed %d -> permutation p = %.3f\n", sum(ok$agree), mean(perm >= sum(ok$agree))))

## ---- the pre-registered direction, tested once ---------------------------------------------------
cat("\n=== the one prediction that was stated in advance: low velo separation ===\n")
for (lg in c("D1","MLB")) { X <- if (lg == "D1") D else S
  i <- X[, ec >= EFF & ef >= EFF & axis <= 10 & vband == "low"]
  x <- X$y[i]; y <- X$y[!i]; t <- t.test(x, y)
  cat(sprintf("  %-3s axis<=10 + low velo band: n = %2d from %2d arms | %+.2f vs %+.2f | diff %+.2f, p = %.3f\n",
              lg, sum(i), uniqueN(X$id[i]), mean(x), mean(y), mean(x)-mean(y), t$p.value)) }
