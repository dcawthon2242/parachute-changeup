#!/usr/bin/env Rscript

# M1: DOES THE BIN EFFECT SCALE WITH FOUR-SEAM USAGE?
#
# The confound and the mechanism make opposite predictions here. If these arms are simply good, the
# advantage is a property of the pitcher and should not care how often the hitter actually sees the
# fastball the changeup is supposedly imitating. If the changeup works by being mistaken for the
# four-seamer, the advantage needs the four-seamer to be present, and should grow with its usage.
#
# The test is the bin-by-usage interaction, not the usage main effect. Usage on its own is
# confounded six ways - four-seam-heavy pitchers differ from sinker-heavy ones in velocity, slot
# and role - so the identifying comparison is between high- and low-usage arms *within* the bin,
# benchmarked against the same contrast outside it.
#
# MLB usage comes from Savant's arsenal leaderboard, which reports the share of each pitch type
# over a pitcher's whole season. That is the right denominator and it reaches back to 2020, so no
# bin season is lost. D1 usage is computed directly from the TrackMan pitch table.

# bit64 is mandatory here: see the note at the top of spec_lock.R. The NCAA PitcherId is integer64
# and stringifies to garbage without it, which joins cleanly against other garbage and not at all
# against the real digits.
suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 200); set.seed(11)
MDIR <- "data/statcast_model"
L <- readRDS(file.path(MDIR, "locked_spec.rds"))$data; setDT(L)

## ---- MLB four-seam usage from the Savant arsenal leaderboard -----------------------------------
AR <- rbindlist(lapply(2020:2026, function(y) {
  x <- fread(sprintf("data/savant/arsenal_%d.csv", y), showProgress = FALSE)
  setnames(x, names(x)[2], "pitcher")
  x[, season := y]
  x[, .(pitcher = as.integer(pitcher), season,
        ff_use = suppressWarnings(as.numeric(n_ff)),
        ch_use = suppressWarnings(as.numeric(n_ch)),
        si_use = suppressWarnings(as.numeric(n_si)))] }), use.names = TRUE)
AR[is.na(ff_use), ff_use := 0][is.na(ch_use), ch_use := 0][is.na(si_use), si_use := 0]
cat(sprintf("MLB arsenal rows 2020-2026: %s\n", format(nrow(AR), big.mark = ",")))

## ---- D1 four-seam usage from the pitch table ---------------------------------------------------
UD <- file.path(MDIR, "ncaa_usage.rds")
if (!file.exists(UD)) {
  D <- readRDS(file.path(MDIR, "ncaa_d1_seq_pitches.rds")); setDT(D)
  D[, PitcherId := as.character(PitcherId)]
  u <- D[, .(tot = .N, nff = sum(pt == "Four-Seam"), nch = sum(pt == "Changeup"),
             nsi = sum(pt == "Sinker")), by = .(PitcherId, season)]
  u[, `:=`(ff_use = 100*nff/tot, ch_use = 100*nch/tot, si_use = 100*nsi/tot)]
  saveRDS(u, UD); rm(D); gc()
}
UN <- readRDS(UD); setDT(UN)

## ---- attach usage to the locked pool -----------------------------------------------------------
# After the four-seam-primary lock, usage already lives on locked_spec. Re-merging it would
# produce ff_use.x / ff_use.y and then fail. The arsenal files are still read above so this
# script can rebuild usage if it is ever run against an older spec that lacks the columns.
L[, id := as.character(id)]
if (all(c("ff_use","si_use") %in% names(L))) {
  X <- copy(L)[is.finite(ff_use) & is.finite(si_use)]
} else {
  LM <- merge(L[league == "MLB"], AR[, .(id = as.character(pitcher), season, ff_use, ch_use, si_use)],
              by = c("id","season"), all.x = TRUE)
  LD <- merge(L[league == "D1"], UN[, .(id = PitcherId, season, ff_use, ch_use, si_use)],
              by = c("id","season"), all.x = TRUE)
  X <- rbind(LM, LD, fill = TRUE)[is.finite(ff_use)]
}
cat(sprintf("pooled seasons with usage attached: %d of %d (%d in bin)\n",
            nrow(X), nrow(L), sum(X$bin)))
cat("\nfour-seam usage %, in bin versus out, by league:\n")
print(X[, .(seasons = .N, ff_use = round(mean(ff_use),1), median = round(median(ff_use),1)),
        by = .(league, bin)][order(league, bin)], row.names = FALSE)

## ---- the interaction test ----------------------------------------------------------------------
# Usage is centred within league so the bin main effect stays interpretable as the effect at an
# average-usage pitcher, and so the two leagues' different usage distributions do not induce a
# spurious interaction when pooled.
X[, ffc := ff_use - mean(ff_use), by = league]
X[, lg := factor(league)]
X[, w := pmin(nsw, 400)]

fit <- lm(y ~ lg + bin*ffc, data = X, weights = w)
S <- summary(fit)$coefficients
cat("\n=== M1: bin x four-seam usage, pooled, weighted by swings ===\n")
for (r in c("binTRUE","ffc","binTRUE:ffc")) if (r %in% rownames(S))
  cat(sprintf("  %-14s %+7.3f  se %5.3f  t %+5.2f  p %.4f\n", r, S[r,1], S[r,2], S[r,3], S[r,4]))
cat(sprintf("\n  interaction reads as: each +1 pct point of four-seam usage changes the bin\n"))
cat(sprintf("  advantage by %+.3f whiff points (%.2f per +10 points of usage)\n",
            S["binTRUE:ffc",1], 10*S["binTRUE:ffc",1]))

cat("\n=== the same test run inside each league ===\n")
for (lgv in c("MLB","D1")) {
  Z <- X[league == lgv]
  f2 <- lm(y ~ bin*ffc, data = Z, weights = w); s2 <- summary(f2)$coefficients
  cat(sprintf("  %-4s n=%4d (bin %2d)  bin %+6.2f  interaction %+6.3f se %.3f  p %.3f\n",
              lgv, nrow(Z), sum(Z$bin), s2["binTRUE",1], s2["binTRUE:ffc",1],
              s2["binTRUE:ffc",2], s2["binTRUE:ffc",4]))
}

## ---- the same thing without a functional form --------------------------------------------------
# A linear interaction assumes the dose response is a straight line. Splitting the bin at its own
# median usage asks the question without that assumption, and it is the version worth reporting if
# the two disagree.
cat("\n=== bin effect split at the bin's own median four-seam usage ===\n")
for (lgv in c("MLB","D1")) {
  Z <- X[league == lgv]; cut <- median(Z[bin == TRUE]$ff_use)
  for (hi in c(FALSE, TRUE)) {
    A <- Z[bin == TRUE & (ff_use >= cut) == hi]; B <- Z[bin == FALSE]
    e <- weighted.mean(A$y, A$w) - weighted.mean(B$y, B$w)
    se <- sqrt(var(A$y)/nrow(A) + var(B$y)/nrow(B))
    cat(sprintf("  %-4s %-9s usage (cut %4.1f%%)  n=%2d  %+6.2f +/- %.2f\n",
                lgv, if (hi) "high" else "low", cut, nrow(A), e, se))
  }
}
saveRDS(X, file.path(MDIR, "mech_usage.rds"))
cat("\nwrote mech_usage.rds\n")
