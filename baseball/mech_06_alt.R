#!/usr/bin/env Rscript

# ALTERNATIVE MECHANISM SEARCH.
#
# All three pre-specified mechanism tests came back null, so the plan's branch is taken: look for
# some other structure in where the advantage sits. This is exploratory and is labelled as such -
# nothing here was pre-registered, three families are examined at once, and any p-value below should
# be read against that. The purpose is to generate a hypothesis for a future season of data, not to
# rescue the current one.
#
# A1  Batter handedness. A changeup's platoon behaviour is its most basic property. If the bin's
#     advantage is concentrated against opposite-handed hitters that is at least consistent with a
#     pitch-level story; if it is flat across handedness it is more likely a pitcher trait.
# A2  Count state. Deception should matter most where the hitter is protecting - two strikes.
# A3  Release consistency. The most likely confound left standing. Pitchers who match spin axis
#     between two pitches may simply be pitchers who repeat their delivery, and repeatability is a
#     well-known correlate of command and of quality. If controlling for it removes the bin effect,
#     the bin was never about axis.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 200); set.seed(29)
MDIR <- "data/statcast_model"
LOCK <- readRDS(file.path(MDIR, "locked_spec.rds"))
L <- LOCK$data; setDT(L); L[, id := as.character(id)]

## ---- pitch-level table with handedness and count -----------------------------------------------
AL <- file.path(MDIR, "mech_alt_pitches.rds")
if (!file.exists(AL)) {
  X <- readRDS(file.path(MDIR, "parachute_rv.rds"))
  setDT(X); X <- X[pitch_type == "CH", .(game_pk, at_bat_number, pitch_number, same_hand,
                                          balls, strikes)]
  MS <- readRDS(file.path(MDIR, "mlb_whiff_resid_seq.rds")); setDT(MS)
  MS <- merge(MS, X, by = c("game_pk","at_bat_number","pitch_number"))
  MS[, `:=`(league = "MLB", id = as.character(pitcher))]
  DS <- readRDS(file.path(MDIR, "ncaa_whiff_resid_seq.rds")); setDT(DS)
  DS[, `:=`(league = "D1", id = as.character(PitcherId), balls = Balls, strikes = Strikes)]
  # The residual table already carries the batter's side; the pitcher's hand comes from the
  # pitcher-season pair table rather than re-reading the 444 MB pitch file for one column.
  TH <- as.data.table(readRDS(file.path(MDIR, "ncaa_spineff_pairs.rds")))
  DS <- merge(DS, TH[, .(id = as.character(PitcherId), season, throws)],
              by = c("id","season"))
  DS[, same_hand := as.integer((throws == "Left") == (BatterSide == "Left"))]
  saveRDS(rbind(MS[, .(league, id, season, same_hand, balls, strikes, r_all)],
                DS[, .(league, id, season, same_hand, balls, strikes, r_all)]), AL)
}
K <- readRDS(AL); setDT(K)
K <- merge(K, L[, .(league, id, season, bin)], by = c("league","id","season"))
K <- K[!is.na(id) & !is.na(same_hand)]
K[, y := 100*r_all]

# Pitches from one arm are not independent, so the bin mean is formed over arms rather than over
# pitches: each pitcher contributes one number and the spread of those numbers gives the interval.
# Averaging the pitches directly would treat a 300-changeup season as 300 observations and report an
# interval several times too narrow.
split_eff <- function(D, lab) {
  A <- D[bin == TRUE]; B <- D[bin == FALSE]
  if (nrow(A) < 50 || uniqueN(A$id) < 3) return(invisible(NULL))
  am <- A[, .(m = mean(y)), by = id]$m
  e <- mean(am) - mean(B$y); s <- sd(am)/sqrt(length(am))
  cat(sprintf("    %-22s bin %5s pitches / %2d arms   %+6.2f +/- %.2f\n", lab,
              format(nrow(A), big.mark = ","), uniqueN(A$id), e, s))
}

cat("=== A1: does the advantage depend on batter handedness? ===\n")
for (lgv in c("MLB","D1")) { cat(sprintf("  %s\n", lgv))
  split_eff(K[league == lgv & same_hand == 0], "opposite-handed")
  split_eff(K[league == lgv & same_hand == 1], "same-handed") }

cat("\n=== A2: does it depend on the count? ===\n")
K[, cs := fifelse(strikes == 2, "two strikes", fifelse(balls > strikes, "hitter ahead", "even/behind"))]
for (lgv in c("MLB","D1")) { cat(sprintf("  %s\n", lgv))
  for (c in c("two strikes","even/behind","hitter ahead")) split_eff(K[league == lgv & cs == c], c) }

## ---- A3: release consistency --------------------------------------------------------------------
# Two quantities. rel_sd is how tightly a pitcher repeats his own changeup release, and arm_gap is
# how far the changeup's slot sits from the four-seamer's. Both are look-alike cues that have
# nothing to do with spin axis, and both plausibly correlate with it.
RC <- file.path(MDIR, "mech_release.rds")
if (!file.exists(RC)) {
  FC <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(FC); FC <- FC[pitch_type == "CH"]
  M <- FC[, .(rel_sd = sqrt(var(release_pos_x, na.rm = TRUE) + var(release_pos_z, na.rm = TRUE)),
              arm_gap = abs(mean(arm_diff, na.rm = TRUE))),
          by = .(id = as.character(pitcher), season)]
  M[, league := "MLB"]
  D <- readRDS(file.path(MDIR, "ncaa_d1_seq_pitches.rds")); setDT(D)
  D[, id := as.character(PitcherId)]
  ch <- D[pt == "Changeup", .(rel_sd = sqrt(var(RelSide, na.rm=TRUE) + var(RelHeight, na.rm=TRUE)),
                              cx = mean(RelSide, na.rm=TRUE), cz = mean(RelHeight, na.rm=TRUE)),
          by = .(id, season)]
  ff <- D[pt == "Four-Seam", .(fx = mean(RelSide, na.rm=TRUE), fz = mean(RelHeight, na.rm=TRUE)),
          by = .(id, season)]
  d <- merge(ch, ff, by = c("id","season"))
  d[, `:=`(league = "D1", arm_gap = sqrt((cx-fx)^2 + (cz-fz)^2))]
  saveRDS(rbind(M[, .(league, id, season, rel_sd, arm_gap)],
                d[, .(league, id, season, rel_sd, arm_gap)]), RC); rm(D); gc()
}
R <- readRDS(RC); setDT(R)
P <- merge(L, R, by = c("league","id","season"))
P <- P[is.finite(rel_sd) & is.finite(arm_gap)]
cat(sprintf("\n=== A3: release consistency (%d seasons, %d in bin) ===\n", nrow(P), sum(P$bin)))
cat("  is the bin selecting on repeatability?\n")
print(P[, .(seasons = .N, rel_sd = round(mean(rel_sd),3), arm_gap = round(mean(arm_gap),2)),
        by = .(league, bin)][order(league, -bin)], row.names = FALSE)
for (lgv in c("MLB","D1")) { Z <- P[league == lgv]
  cat(sprintf("  %-4s axis gap vs rel_sd r = %+.3f | axis gap vs arm_gap r = %+.3f\n", lgv,
              cor(Z$axis, Z$rel_sd), cor(Z$axis, Z$arm_gap))) }

cat("\n  does the bin effect survive controlling for both?\n")
P[, `:=`(zs = scale(rel_sd)[,1], za = scale(arm_gap)[,1]), by = league]
P[, w := pmin(nsw, 400)]
for (mod in list(c("bin"), c("bin","zs","za"))) {
  f <- lm(reformulate(c("league", mod), "y"), data = P, weights = w); s <- summary(f)$coefficients
  cat(sprintf("    %-24s bin %+6.2f  se %.2f  p %.3f\n",
              paste(mod, collapse = " + "), s["binTRUE",1], s["binTRUE",2], s["binTRUE",4])) }
cat("\n  and the reverse - do the release cues predict anything on their own?\n")
f <- lm(y ~ league + zs + za, data = P, weights = P$w); s <- summary(f)$coefficients
for (r in c("zs","za")) cat(sprintf("    %-4s %+6.2f  se %.2f  p %.3f\n", r, s[r,1], s[r,2], s[r,4]))
