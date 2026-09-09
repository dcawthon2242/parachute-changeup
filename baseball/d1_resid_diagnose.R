#!/usr/bin/env Rscript

# WHY DOES VELOCITY SURVIVE IN THE D1 WHIFF RESIDUAL?
#
# speed_diff is a feature of the D1 whiff model, so the residual should be flat in it. At the
# pitcher-season level it is not: r = +0.080, p = 0.008, a slope of +0.23 whiff points per mph.
# MLB's equivalent model, built the same way, is flat at +0.014. Something is specific to D1.
#
# There are two candidate explanations and they call for different fixes, so the diagnosis has to
# come before the repair.
#
# If the leakage is present AT THE PITCH LEVEL, the model is simply underfitting velocity and the
# answer is a better model - more capacity, or velocity entered in a form the trees can use.
#
# If the pitch level is clean and the leakage only appears AFTER AGGREGATION, the cause is
# different and more interesting: the model is absorbing the within-pitcher variation in separation
# while missing the between-pitcher component. That happens when a pitcher's mean separation
# carries information his individual pitches do not - it proxies for something about the pitcher
# rather than the pitch - and no amount of extra tree capacity on pitch-level features will remove
# it. The fix there is to residualise at the level the analysis actually uses.

suppressPackageStartupMessages(library(data.table)); options(width = 200)
MDIR <- "data/statcast_model"

R <- readRDS(file.path(MDIR, "ncaa_whiff_resid.rds")); setDT(R)
# PitcherId arrives as integer64, which will not join against the double-typed ids in the pair
# tables. Character is the safe common type and the id is only ever used as a key.
R[, PitcherId := as.character(PitcherId)]
cat(sprintf("D1 whiff residual rows: %s from %d pitcher-seasons\n",
            format(nrow(R), big.mark = ","), uniqueN(paste(R$PitcherId, R$season))))

cat("\n=== 1. where does the leakage live? ===\n")
for (v in c("r_tj","r_all")) {
  p <- cor.test(R$speed_diff, R[[v]])
  S <- R[, .(y = mean(get(v)), vs = -mean(speed_diff), n = .N), by = .(PitcherId, season)][n >= 40]
  s <- cor.test(S$vs, S$y)
  cat(sprintf("  %-6s pitch level (n = %s): r = %+.4f (p = %.3f)   |   season level (n = %d): r = %+.4f (p = %.4f)\n",
              v, format(nrow(R), big.mark = ","), p$estimate, p$p.value, nrow(S), s$estimate, s$p.value))
}
cat("  A clean pitch level with a dirty season level is the aggregation case, not underfitting.\n")

## ---- 2. is it between-pitcher or within-pitcher? -----------------------------------------------
# Splitting the season-level correlation into its between- and within-pitcher parts localises it
# exactly. A purely between-pitcher slope means mean separation is standing in for pitcher quality.
S <- R[, .(y = mean(r_all), vs = -mean(speed_diff), n = .N), by = .(PitcherId, season)][n >= 40]
S[, `:=`(vs_bar = mean(vs), y_bar = mean(y)), by = PitcherId]
S[, `:=`(vs_w = vs - vs_bar, y_w = y - y_bar)]
BW <- unique(S[, .(PitcherId, vs_bar, y_bar)])
cb <- cor.test(BW$vs_bar, BW$y_bar)
MW <- S[, .N, by = PitcherId][N >= 2, PitcherId]
cw <- cor.test(S[PitcherId %in% MW, vs_w], S[PitcherId %in% MW, y_w])
cat(sprintf("\n=== 2. decomposition ===\n  between pitchers (n = %d): r = %+.3f (p = %.4f)\n",
            nrow(BW), cb$estimate, cb$p.value))
cat(sprintf("  within pitchers  (n = %d deviations from %d arms): r = %+.3f (p = %.3f)\n",
            S[PitcherId %in% MW, .N], length(MW), cw$estimate, cw$p.value))

## ---- 3. does velocity separation proxy for sample size or stuff? --------------------------------
# If harder separators are simply better pitchers who also throw more, the leakage is confounding
# rather than model error, and controlling it is mandatory regardless of its source.
cat("\n=== 3. what does mean separation travel with? ===\n")
Q <- merge(S, R[, .(spin = mean(speed_diff), rz = mean(RelHeight)), by = .(PitcherId, season)],
           by = c("PitcherId","season"))
cat(sprintf("  separation vs swings: r = %+.3f | vs release height: r = %+.3f\n",
            cor(Q$vs, Q$n), cor(Q$vs, Q$rz)))
cat(sprintf("  raw whiff-above spread: sd = %.2f over %d seasons\n", sd(S$y)*100, nrow(S)))

## ---- 4. the repair, and what it costs -----------------------------------------------------------
# The analysis consumes pitcher-season means, so the honest correction is applied there: regress the
# season-level residual on season-level mean separation and keep what is left. This is a strictly
# more conservative residual than the raw one and cannot manufacture an effect.
cat("\n=== 4. after residualising at the season level ===\n")
S[, y_adj := resid(lm(y ~ vs, S))]
cat(sprintf("  baseline correlation with separation: %+.4f -> %+.4f\n",
            cor(S$vs, S$y), cor(S$vs, S$y_adj)))
saveRDS(S[, .(PitcherId, season, n, y_raw = 100*y, y_adj = 100*y_adj, velo_sep = vs)],
        file.path(MDIR, "ncaa_whiff_season_adj.rds"))
cat(sprintf("  wrote %s\n", file.path(MDIR, "ncaa_whiff_season_adj.rds")))

## ---- 5. the bin, measured both ways ---------------------------------------------------------------
ARM <- as.data.table(readRDS(file.path(MDIR, "ncaa_armangle.rds")))
PR  <- as.data.table(readRDS(file.path(MDIR, "ncaa_spineff_pairs.rds")))
PR  <- merge(PR, ARM[, .(PitcherId, season, arm_hat)], by = c("PitcherId","season"))
setnames(PR, "axis_gap", "axis", skip_absent = TRUE)
PR[, PitcherId := as.character(PitcherId)]
D <- merge(PR, S[, .(PitcherId, season, y_raw = 100*y, y_adj = 100*y_adj)],
           by = c("PitcherId","season"))
thr <- as.numeric(quantile(D$arm_hat, 2/3))
i <- D[, eff_ch >= .85 & eff_ff >= .85 & axis <= 10 & arm_hat >= thr]
cat(sprintf("\n=== 5. the locked cell on D1 (slot threshold %.1f deg) ===\n", thr))
for (v in c("y_raw","y_adj")) { t <- t.test(D[[v]][i], D[[v]][!i])
  cat(sprintf("  %-6s n = %d from %d arms | %+.2f vs %+.2f | diff %+.2f, se %.2f, p = %.3f\n",
              v, sum(i), uniqueN(D$PitcherId[i]), mean(D[[v]][i]), mean(D[[v]][!i]),
              diff(rev(t$estimate)), sd(D[[v]][i])/sqrt(sum(i)), t$p.value)) }
