#!/usr/bin/env Rscript

# MULTIPLICITY AUDIT OF THE SEARCHES THAT PRODUCED PUBLISHED CELLS.
#
# This is a correction of old searches, not a new one. The outcome is the four-seam-primary
# residual already locked in spec_lock.R. Each family is the grid that was actually examined;
# the permutation asks how often the best cell of that family looks as good as the observed
# best when the residual is shuffled within league.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 200); set.seed(37)
MDIR <- "data/statcast_model"
L <- readRDS(file.path(MDIR, "locked_spec.rds"))$data; setDT(L)
L[, id := as.character(id)]
NPERM <- 2000L

cell_eff <- function(D, sel) {
  A <- D[sel]; B <- D[!sel]
  if (nrow(A) < 4 || nrow(B) < 4) return(NA_real_)
  mean(A$y) - mean(B$y)
}

audit <- function(name, cells, D) {
  obs <- sapply(cells, function(fn) cell_eff(D, fn(D)))
  names(obs) <- names(cells)
  best <- max(obs, na.rm = TRUE)
  who <- names(obs)[which.max(replace(obs, is.na(obs), -Inf))]
  cat(sprintf("\n=== %s (%d cells) ===\n", name, length(cells)))
  for (nm in names(obs)) {
    if (is.na(obs[nm])) cat(sprintf("  %-28s  (too small)\n", nm))
    else cat(sprintf("  %-28s  %+6.2f%s\n", nm, obs[nm], if (nm == who) "   <- best" else ""))
  }
  null_best <- replicate(NPERM, {
    Dp <- copy(D)[, y := sample(y)]
    max(sapply(cells, function(fn) cell_eff(Dp, fn(Dp))), na.rm = TRUE)
  })
  cat(sprintf("  observed best = %+.2f (%s) | null mean of best-of-%d = %+.2f (sd %.2f) | p = %.4f\n",
              best, who, length(cells), mean(null_best), sd(null_best),
              mean(null_best >= best)))
  invisible(list(family = name, best = who, est = best, p = mean(null_best >= best)))
}

## ---- 1. axis bands, efficiency floors + top-third slot held fixed ------------------------------
# Five bands, the search that produced the 30-45 anomaly and the locked 0-10 cell.
BR <- c(0, 10, 20, 30, 45, 360)
L[, band := cut(axis, BR, right = TRUE, include.lowest = TRUE,
                labels = c("0-10","10-20","20-30","30-45","45+"))]
base <- function(D) D$ec >= .85 & D$ef >= .85 & D$arm >= D$arm_thr
bands <- setNames(lapply(levels(L$band), function(b) {
  force(b); function(D) base(D) & D$band == b
}), paste0("axis ", levels(L$band)))

## ---- 2. arm-slot x axis-gap grid that selected the locked cell ---------------------------------
# Efficiency floors held at .85; axis at 10/15/20; slot either unrestricted or the league top third.
# That is the six-cell search the locked specification came out of.
grid2 <- list()
for (ax in c(10, 15, 20)) for (slot in c(FALSE, TRUE)) {
  lab <- sprintf("axis<=%d, %s", ax, if (slot) "top third" else "any slot")
  grid2[[lab]] <- local({
    ax <- ax; slot <- slot
    function(D) D$ec >= .85 & D$ef >= .85 & D$axis <= ax & (if (slot) D$arm >= D$arm_thr else TRUE)
  })
}

## ---- 3. velocity-separation gates --------------------------------------------------------------
# The 18-cell grid from velo_gate_grid.R: three axis cuts, slot on/off, velo tertiles.
L[, vband := cut(vs, quantile(vs, c(0, 1/3, 2/3, 1)), c("low","mid","high"),
                 include.lowest = TRUE), by = league]
grid3 <- list()
for (ax in c(10, 15, 20)) for (slot in c(FALSE, TRUE)) for (vb in c("low","mid","high")) {
  lab <- sprintf("axis<=%d, %s, velo %s", ax, if (slot) "top third" else "any slot", vb)
  grid3[[lab]] <- local({
    ax <- ax; slot <- slot; vb <- vb
    function(D) D$ec >= .85 & D$ef >= .85 & D$axis <= ax & D$vband == vb &
      (if (slot) D$arm >= D$arm_thr else TRUE)
  })
}

cat("Outcome is the locked four-seam-primary residual. Permutations are within league, then the\n")
cat("two league best-of-k p-values are reported separately so one noisy league cannot carry both.\n")

out <- list()
for (lgv in c("MLB","D1")) {
  D <- L[league == lgv]
  cat(sprintf("\n######## %s  (%d seasons) ########\n", lgv, nrow(D)))
  out[[paste(lgv,"bands")]] <- audit("axis bands", bands, D)
  out[[paste(lgv,"armaxis")]] <- audit("arm x axis grid", grid2, D)
  out[[paste(lgv,"velo")]] <- audit("velo-separation grid", grid3, D)
}

# The locked cell itself, tested once, for comparison with the best-of-k numbers.
cat("\n=== locked cell, tested once (not a search) ===\n")
for (lgv in c("MLB","D1")) {
  D <- L[league == lgv]
  e <- cell_eff(D, D$bin)
  null <- replicate(NPERM, { Dp <- copy(D)[, y := sample(y)]; cell_eff(Dp, Dp$bin) })
  cat(sprintf("  %-4s  %+6.2f | one-cell permutation p = %.4f\n", lgv, e, mean(null >= e)))
}
saveRDS(out, file.path(MDIR, "multiplicity_audit.rds"))
cat("\nwrote multiplicity_audit.rds\n")
