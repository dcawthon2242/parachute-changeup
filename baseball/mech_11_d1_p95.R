#!/usr/bin/env Rscript

# D1 AXIS CONDITION THAT ACTUALLY TRANSFERS.
#
# TrackMan's SpinAxis is movement-derived. At high inferred efficiency the movement axis and
# Hawk-Eye's measured axis converge, so a D1 test at .95 on both pitches is the first one where
# "axis gap" means approximately the same thing as the MLB construct. A fixed 10-degree cut still
# does not transfer: the movement gap runs wider, and MLB's own .95 pool puts only 4.5 percent of
# seasons at or under 10 degrees. The D1 cut is therefore the matching percentile of D1's .95 pool,
# not 10 degrees copied across.
#
# This is not a new search. The floors and the matching rule were specified before looking at the
# D1 residual at this cut. Arm slot is reported as a robustness check against the locked
# specification, not as a gate that can be retuned.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 200)
MDIR <- "data/statcast_model"
L <- readRDS(file.path(MDIR, "locked_spec.rds"))$data; setDT(L)

MLB <- L[league == "MLB" & ec >= .95 & ef >= .95]
D1  <- L[league == "D1"  & ec >= .95 & ef >= .95]
share <- mean(MLB$axis <= 10)
cut <- as.numeric(quantile(D1$axis, share, na.rm = TRUE))
cat("=== transfer rule, declared from the MLB .95 pool ===\n")
cat(sprintf("  MLB seasons with both active spins >= .95: %d\n", nrow(MLB)))
cat(sprintf("  of those, axis <= 10 degrees: %d (%.1f%%)\n", sum(MLB$axis <= 10), 100*share))
cat(sprintf("  D1 seasons with both inferred efficiencies >= .95: %d\n", nrow(D1)))
cat(sprintf("  matching percentile cut on D1 axis: %.2f degrees\n", cut))

D1[, in_p := axis <= cut]
cat(sprintf("  D1 seasons inside the matched cut: %d from %d arms\n",
            sum(D1$in_p), uniqueN(D1[in_p == TRUE]$id)))

gap <- function(A, B, lab) {
  if (nrow(A) < 3) { cat(sprintf("  %-40s  n=%d  (too small)\n", lab, nrow(A))); return(invisible(NULL)) }
  e <- mean(A$y) - mean(B$y)
  se <- sqrt(var(A$y)/nrow(A) + var(B$y)/nrow(B))
  cat(sprintf("  %-40s  %2d seasons / %2d arms   %+6.2f +/- %.2f   p %.3f\n",
              lab, nrow(A), uniqueN(A$id), e, se, 2*pnorm(-abs(e/se))))
}

# Comparison group is the rest of the locked D1 pool, not the rest of the .95 pool, so this is
# the same contrast the locked cell uses.
OUT <- L[league == "D1" & !(ec >= .95 & ef >= .95 & axis <= cut)]
cat("\n=== D1 .95 + percentile-matched axis, vs the rest of the locked D1 pool ===\n")
gap(D1[in_p == TRUE], OUT, ".95 + matched axis (the test)")
gap(D1[in_p == TRUE & arm >= arm_thr], OUT, ".95 + matched axis + top-third slot")
gap(D1[axis <= 10], L[league == "D1" & axis > 10], "fixed 10 deg inside the .95 pool (old cut)")

# MLB at the same .95 floor, fixed 10 degrees, for the side-by-side the transfer is about.
MOUT <- L[league == "MLB" & !(ec >= .95 & ef >= .95 & axis <= 10)]
cat("\n=== MLB at the same .95 floor, axis <= 10 (not a new MLB search) ===\n")
gap(MLB[axis <= 10], MOUT, "MLB .95 + axis <= 10")
gap(MLB[axis <= 10 & arm >= arm_thr], MOUT, "MLB .95 + axis <= 10 + top-third")

cat("\n=== D1 seasons in the transferred cell ===\n")
print(D1[in_p == TRUE][order(-y), .(id, season, axis = round(axis,1), arm = round(arm,1),
      ec = round(ec,3), ef = round(ef,3), nsw, y = round(y,2), top_third = arm >= arm_thr)],
      row.names = FALSE)

saveRDS(list(share = share, cut = cut, d1 = D1, mlb = MLB),
        file.path(MDIR, "d1_p95_transfer.rds"))
cat("\nwrote d1_p95_transfer.rds\n")
