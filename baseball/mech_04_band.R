#!/usr/bin/env Rscript

# THE 30-45 DEGREE BAND.
#
# Holding the efficiency floors and the arm-slot percentile fixed and cutting on axis gap alone, MLB
# reads larger at 30-45 degrees than at 0-10. That is a problem for the whole framing: a changeup
# whose spin axis sits 40 degrees off the fastball's is not a look-alike, so if that band is where
# the effect really lives, "matched axis" is a mislabel for whatever is going on.
#
# Four things could produce it, and they are separable:
#   1. it is a handful of arms - check the distinct-pitcher count and collapse repeats
#   2. the anchor is wrong - a pitcher whose real fastball is a sinker is being scored against a
#      four-seamer he barely throws, so his "gap" is measured off a pitch the hitter rarely sees
#   3. it is noise - check whether D1 agrees, and what the band looks like under permutation
#   4. it is real, and axis similarity is not the operative variable

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 200); set.seed(19)
MDIR <- "data/statcast_model"
L <- readRDS(file.path(MDIR, "locked_spec.rds"))$data; setDT(L); L[, id := as.character(id)]
U <- readRDS(file.path(MDIR, "mech_usage.rds")); setDT(U); U[, id := as.character(id)]

# Everything except the axis cut, so the bands are comparable to the locked cell in all other ways.
L[, base := ec >= .85 & ef >= .85 & arm >= arm_thr]
BR <- c(0, 10, 20, 30, 45, 360)
L[, band := cut(axis, BR, right = TRUE, include.lowest = TRUE,
                labels = c("0-10","10-20","20-30","30-45","45+"))]

eff <- function(D, sel) { A <- D[sel]; B <- D[base == FALSE]
  if (nrow(A) < 3) return(c(NA, NA, NA, NA))
  c(nrow(A), uniqueN(A$id), mean(A$y) - mean(B$y), sqrt(var(A$y)/nrow(A) + var(B$y)/nrow(B))) }

cat("=== axis bands within the efficiency + arm-slot base, versus everything outside the base ===\n")
for (lgv in c("MLB","D1")) {
  D <- L[league == lgv]
  cat(sprintf("\n%s (base pool %d seasons, comparison pool %d)\n", lgv, sum(D$base), sum(!D$base)))
  for (b in levels(L$band)) { e <- eff(D, D$base & D$band == b)
    if (!is.na(e[1])) cat(sprintf("  axis %-6s  %3d seasons / %3d arms   %+6.2f +/- %.2f\n",
                                  b, e[1], e[2], e[3], e[4])) }
}

## ---- 1. is the band a handful of arms? ---------------------------------------------------------
cat("\n=== who is in the MLB 30-45 band ===\n")
B <- L[league == "MLB" & base & band == "30-45"][order(-y)]
FC <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(FC)
NM <- unique(FC[pitch_type == "CH", .(id = as.character(pitcher), player_name)])
B <- merge(B, NM, by = "id", all.x = TRUE)[order(-y)]
print(B[, .(player_name, season, axis = round(axis,1), arm = round(arm,1), ec, ef,
            nsw, y = round(y,2))], row.names = FALSE)
cat(sprintf("\n  %d seasons from %d arms; the top contributor appears %d times\n",
            nrow(B), uniqueN(B$id), max(B[, .N, by = id]$N)))

# One row per arm, averaging his seasons, so a pitcher with four good years counts once. If the band
# survives this it is not one man repeated.
cat("\n=== every MLB band, collapsed to one row per pitcher ===\n")
D <- L[league == "MLB"]
CO <- D[, .(y = mean(y), base = base[1], band = band[1], nb = .N), by = id]
out <- D[base == FALSE][, .(y = mean(y)), by = id]
for (b in levels(L$band)) { A <- CO[base == TRUE & band == b]
  if (nrow(A) >= 3) cat(sprintf("  axis %-6s  %3d arms  %+6.2f +/- %.2f\n", b, nrow(A),
    mean(A$y) - mean(out$y), sqrt(var(A$y)/nrow(A) + var(out$y)/nrow(out)))) }

## ---- 2. is the anchor wrong? -------------------------------------------------------------------
# The axis gap is measured against the four-seamer. For a sinker-primary pitcher that is the wrong
# reference: the hitter's expectation is set by the sinker, so a large gap to the four-seamer says
# little about whether the changeup looks like the pitch he actually sees.
cat("\n=== four-seam versus sinker usage inside each MLB band ===\n")
# After the four-seam-primary lock, usage is already on L. Re-merging it would duplicate columns.
if (all(c("ff_use","si_use") %in% names(L))) {
  M <- L[league == "MLB" & base == TRUE]
} else {
  M <- merge(L[league == "MLB" & base == TRUE], U[league == "MLB", .(id, season, ff_use, si_use)],
             by = c("id","season"))
}
print(M[, .(seasons = .N, ff_use = round(mean(ff_use),1), si_use = round(mean(si_use),1),
            sinker_primary = sum(si_use > ff_use)), by = band][order(band)], row.names = FALSE)
cat("\n  the band re-measured after dropping sinker-primary arms:\n")
MB <- M[si_use <= ff_use]; ob <- L[league == "MLB" & base == FALSE]
for (b in levels(L$band)) { A <- MB[band == b]
  if (nrow(A) >= 3) cat(sprintf("    axis %-6s  %3d seasons / %2d arms  %+6.2f +/- %.2f\n", b,
    nrow(A), uniqueN(A$id), mean(A$y) - mean(ob$y),
    sqrt(var(A$y)/nrow(A) + var(ob$y)/nrow(ob)))) }

## ---- 3. is it noise? ---------------------------------------------------------------------------
# Five bands were examined, so the largest of them is expected to look good even under the null.
# The permutation asks how often ANY band reaches +2.37 when the residual is shuffled - the right
# comparison for a number that was found by looking across bands.
cat("\n=== permutation: how often does the best of five bands reach the observed 30-45 value? ===\n")
D <- L[league == "MLB"]; obs <- eff(D, D$base & D$band == "30-45")[3]
best <- replicate(4000, { Dp <- copy(D)[, y := sample(y)]
  max(sapply(levels(L$band), function(b) { A <- Dp[base & band == b]
    if (nrow(A) < 3) -Inf else mean(A$y) - mean(Dp[base == FALSE]$y) })) })
cat(sprintf("  observed 30-45 = %+.2f | best-of-five null mean %+.2f, sd %.2f | p = %.4f\n",
            obs, mean(best), sd(best), mean(best >= obs)))
cat(sprintf("  the same shuffle test for the locked 0-10 cell: p = %.4f\n",
            mean(replicate(2000, { Dp <- copy(D)[, y := sample(y)]
              A <- Dp[base & band == "0-10"]; mean(A$y) - mean(Dp[base == FALSE]$y) }) >=
                 eff(D, D$base & D$band == "0-10")[3])))

## ---- 4. does D1 see it? ------------------------------------------------------------------------
cat("\n=== the two leagues side by side, band for band ===\n")
for (b in levels(L$band)) {
  r <- sapply(c("MLB","D1"), function(lg) { D <- L[league == lg]
    e <- eff(D, D$base & D$band == b); if (is.na(e[1])) NA_real_ else e[3] })
  cat(sprintf("  axis %-6s  MLB %+6.2f   D1 %+6.2f%s\n", b, r[1], r[2],
              if (!is.na(r[1]) && !is.na(r[2]) && sign(r[1]) == sign(r[2])) "   (agree)" else ""))
}
