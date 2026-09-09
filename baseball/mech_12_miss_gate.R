#!/usr/bin/env Rscript

# MISS DISTANCE ON THE LOCKED MLB RECIPE, AND WHAT A LOOSER SWING GATE DOES.
#
# Miss distance exists only on competitive swings in 2023-2026 (miss_grade_data.rds).
# The locked bin is 2020-2026, so the miss-distance contrast is the 2023+ subset, not the
# full twenty seasons. The swing-gate sweep uses the 2020-2026 whiff residual so the
# comparison is to the locked number, then repeats 2023+ so miss distance and whiffs
# share a population.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 200)
MDIR <- "data/statcast_model"
LOCK <- readRDS(file.path(MDIR, "locked_spec.rds"))
SPEC <- LOCK$spec
L40 <- LOCK$data; setDT(L40); L40 <- L40[league == "MLB"]

## ---- full MLB pool at every swing count, same other gates --------------------------------------
AS <- readRDS(file.path(MDIR, "active_spin_long.rds")); setDT(AS)
F  <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(F); FC <- F[pitch_type == "CH"]
RW <- readRDS(file.path(MDIR, "mlb_whiff_locaware.rds")); setDT(RW)
S <- FC[, .(axis = mean(axis_diff, na.rm = TRUE), arm = mean(arm_angle, na.rm = TRUE),
            vs = -mean(speed_diff, na.rm = TRUE)), by = .(pitcher, season)]
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, ec = active_spin)], by = c("pitcher","season"))
S <- merge(S, AS[pitch_type == "FF", .(pitcher, season, ef = active_spin)], by = c("pitcher","season"))
S <- merge(S, RW[, .(y_raw = 100*mean(r_all), nsw = .N), by = .(pitcher, season)],
           by = c("pitcher","season"))
S <- S[is.finite(arm) & is.finite(ec) & is.finite(ef)]
S[, id := as.character(pitcher)]
# Usage from the arsenal files, not from the 40-swing lock, or loosening the gate
# cannot add anyone: they were already dropped before this table was built.
AR <- rbindlist(lapply(2020:2026, function(y) {
  x <- fread(sprintf("data/savant/arsenal_%d.csv", y), showProgress = FALSE)
  setnames(x, names(x)[2], "pitcher"); x[, season := y]
  x[, .(id = as.character(pitcher), season,
        ff_use = suppressWarnings(as.numeric(n_ff)),
        si_use = suppressWarnings(as.numeric(n_si)))] }), use.names = TRUE)
AR[is.na(ff_use), ff_use := 0][is.na(si_use), si_use := 0]
S <- merge(S, AR, by = c("id","season"))
S <- S[is.finite(ff_use) & is.finite(si_use) & ff_use >= si_use]
# Arm percentile from the locked 40-swing four-seam-primary pool, not recomputed at each
# gate, so loosening the gate does not quietly change who counts as high-slot.
thr <- L40$arm_thr[1]
S[, `:=`(arm_thr = thr,
         recipe = ec >= SPEC$eff_min & ef >= SPEC$eff_min & axis <= SPEC$axis_max & arm >= thr)]
S[, y := resid(lm(y_raw ~ vs))]

## ---- miss distance, 2023-2026 swings only ------------------------------------------------------
MG <- readRDS(file.path(MDIR, "miss_grade_data.rds")); setDT(MG)
MG <- MG[pitch_type == "CH" & is.finite(miss_distance)]
MD <- MG[, .(miss = mean(miss_distance), nmd = .N, wh = 100*mean(is_whiff)),
         by = .(id = as.character(pitcher), season)]
P <- merge(S, MD, by = c("id","season"), all.x = TRUE)

cat("=== miss distance, locked recipe, 2023-2026 only ===\n")
Z <- P[nsw >= 40 & is.finite(miss)]
cat(sprintf("  seasons with both a 40-swing residual and a miss-distance mean: %d (recipe %d)\n",
            nrow(Z), sum(Z$recipe)))
if (sum(Z$recipe) >= 3) {
  t <- t.test(Z[recipe == TRUE]$miss, Z[recipe == FALSE]$miss)
  cat(sprintf("  miss (inches)  bin %.2f  out %.2f  gap %+.2f  p %.3f\n",
              t$estimate[1], t$estimate[2], diff(rev(t$estimate)), t$p.value))
  tw <- t.test(Z[recipe == TRUE]$wh, Z[recipe == FALSE]$wh)
  cat(sprintf("  raw whiff%%     bin %.1f  out %.1f  gap %+.1f  p %.3f\n",
              tw$estimate[1], tw$estimate[2], diff(rev(tw$estimate)), tw$p.value))
  ty <- t.test(Z[recipe == TRUE]$y, Z[recipe == FALSE]$y)
  cat(sprintf("  whiff residual bin %+.2f  out %+.2f  gap %+.2f  p %.3f\n",
              ty$estimate[1], ty$estimate[2], diff(rev(ty$estimate)), ty$p.value))
  cat(sprintf("  inside recipe, miss vs residual r = %+.3f (n=%d)\n",
              cor(Z[recipe == TRUE]$miss, Z[recipe == TRUE]$y), sum(Z$recipe)))
  cat(sprintf("  outside recipe, miss vs residual r = %+.3f\n",
              cor(Z[recipe == FALSE]$miss, Z[recipe == FALSE]$y)))
}
nm <- unique(FC[, .(id = as.character(pitcher), player_name)])
cat("\n  2023+ recipe seasons with miss distance:\n")
print(merge(Z[recipe == TRUE], nm, by = "id")[order(-miss),
      .(player_name, season, nsw, miss = round(miss,2), wh = round(wh,1), y = round(y,1))],
      row.names = FALSE)

## ---- swing-gate sweep --------------------------------------------------------------------------
show <- function(D, lab) {
  A <- D[recipe == TRUE]; B <- D[recipe == FALSE]
  if (nrow(A) < 3) {
    cat(sprintf("  %-22s  recipe %2d  (too small)\n", lab, nrow(A))); return(invisible(NULL))
  }
  e <- mean(A$y) - mean(B$y)
  se <- sqrt(var(A$y)/nrow(A) + var(B$y)/nrow(B))
  extra <- ""
  if (D[is.finite(miss), .N] > 10 && sum(A$is.finite <- is.finite(A$miss)) >= 3) {
    Am <- A[is.finite(miss)]; Bm <- B[is.finite(miss)]
    extra <- sprintf("   miss %+0.2f (n=%d vs %d)", mean(Am$miss) - mean(Bm$miss),
                     nrow(Am), nrow(Bm))
  }
  cat(sprintf("  %-22s  recipe %2d / %2d arms   whiff res %+6.2f +/- %.2f   p %.3f%s\n",
              lab, nrow(A), uniqueN(A$id), e, se, 2*pnorm(-abs(e/se)), extra))
}

cat("\n=== swing-gate sweep, 2020-2026, recipe otherwise unchanged ===\n")
for (g in c(10L, 20L, 30L, 40L, 60L, 80L)) show(P[nsw >= g], sprintf("%d+ swings", g))

cat("\n=== same sweep, 2023-2026 only (the miss-distance window) ===\n")
for (g in c(10L, 20L, 30L, 40L, 60L)) {
  D <- P[season >= 2023 & nsw >= g]
  A <- D[recipe == TRUE]; B <- D[recipe == FALSE]
  if (nrow(A) < 3) { cat(sprintf("  %d+  recipe %d  too small\n", g, nrow(A))); next }
  ey <- mean(A$y) - mean(B$y)
  sey <- sqrt(var(A$y)/nrow(A) + var(B$y)/nrow(B))
  Am <- A[is.finite(miss)]; Bm <- B[is.finite(miss)]
  em <- if (nrow(Am) >= 3) mean(Am$miss) - mean(Bm$miss) else NA_real_
  cat(sprintf("  %d+  recipe %2d / %2d arms   residual %+6.2f +/- %.2f   miss gap %+0.2f in (n=%d)\n",
              g, nrow(A), uniqueN(A$id), ey, sey, em, nrow(Am)))
}
