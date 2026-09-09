#!/usr/bin/env Rscript

# sd_sep HAS SIGNAL BUT DID NOT IMPROVE RMSE. WHICH OF THOSE IS THE REAL ANSWER?
#
# The screen says changeup velocity consistency matters: -2.81 whiff points per mph of season
# standard deviation, p = .005 clustered on the arm. The model says adding it as a feature makes
# honest out-of-sample season RMSE 2 percent WORSE. Both cannot be the practical answer, so this
# gives the feature the two fairest shots it can get.
#
#   1  a LINEAR correction instead of a tree feature. Gradient boosting can spend a season-level
#      aggregate on memorising training arms; a single coefficient cannot. Fit leave-one-arm-out so
#      the correction applied to a pitcher never saw that pitcher.
#   2  the PRIOR season's consistency predicting THIS season's residual. This is the deployable
#      form and it is immune to same-season leakage by construction, since the feature is measured
#      before the outcome exists.
#
# Also checks the obvious confounds: consistency is estimated from a finite sample, so it is
# mechanically noisier for low-usage arms, and a tiring starter is more variable than a one-inning
# reliever for reasons that have nothing to do with pitch design.

suppressPackageStartupMessages({ library(data.table) })
set.seed(47); options(width = 205)
MDIR <- "data/statcast_model"

F <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(F)
F <- F[is.finite(speed_diff)]
PS <- F[, .(np = .N, sd_sep = sd(speed_diff), sep = -mean(speed_diff),
            ch_velo = mean(release_speed)), by = .(pitcher, season)]

L <- readRDS(file.path(MDIR, "velo_sep_resid.rds"))
S <- merge(as.data.table(L$W)[, .(nsw = .N, wh_aware = 100*mean(whiff - p_aware),
                                  wh_blind = 100*mean(whiff - p_blind)), by = .(pitcher, season)],
           as.data.table(L$R)[, .(rv_aware = 100*mean(rv - q_aware)), by = .(pitcher, season)],
           by = c("pitcher","season"))
S <- merge(S, PS, by = c("pitcher","season"))[nsw >= 40 & is.finite(sd_sep)]
S[, id := as.character(pitcher)]
cat(sprintf("%d pitcher-seasons, %d arms\n", nrow(S), uniqueN(S$id)))

## ---- 0. confounds --------------------------------------------------------------------------------
cat("\n=== 0. what is sd_sep correlated with? ===\n")
for (v in c("np","sep","ch_velo","nsw"))
  cat(sprintf("  sd_sep vs %-8s r = %+.3f\n", v, cor(S$sd_sep, S[[v]])))
cat("\n  the screen again, with usage controlled:\n")
crob <- function(f, D, k = "sd_sep") {
  m <- lm(f, D); u <- residuals(m); X <- model.matrix(m); nc <- uniqueN(D$id)
  b <- solve(crossprod(X)); V <- b %*% crossprod(rowsum(X * u, D$id)) %*% b * (nc/(nc-1))
  e <- coef(m)[k]; s <- sqrt(diag(V))[k]
  sprintf("%+.3f (se %.3f) | p %.4f", e, s, 2*pt(-abs(e/s), nc-1))
}
cat(sprintf("    wh_aware ~ sd_sep + sep                : %s\n", crob(wh_aware ~ sd_sep + sep, S)))
cat(sprintf("    wh_aware ~ sd_sep + sep + log(np)      : %s\n", crob(wh_aware ~ sd_sep + sep + log(np), S)))
cat(sprintf("    wh_aware ~ sd_sep + sep + log(np) + velo: %s\n",
            crob(wh_aware ~ sd_sep + sep + log(np) + ch_velo, S)))

## ---- 1. leave-one-arm-out linear correction --------------------------------------------------------
cat("\n=== 1. linear correction, leave-one-arm-out ===\n")
cat("   season RMSE of the aware residual, with and without an sd_sep correction applied\n")
cat("   (a correction that helps should shrink the residual toward zero)\n\n")
arms <- unique(S$id); K <- 10; fa <- sample(rep(1:K, length.out = length(arms)))
S[, fold := fa[match(id, arms)]]
for (v in c("wh_aware","rv_aware")) {
  pred_sd <- rep(NA_real_, nrow(S)); pred_null <- rep(NA_real_, nrow(S))
  for (f in 1:K) {
    tr <- S[fold != f]; te <- which(S$fold == f)
    m1 <- lm(as.formula(paste(v, "~ sd_sep + sep")), tr)
    m0 <- lm(as.formula(paste(v, "~ sep")), tr)
    pred_sd[te]   <- predict(m1, S[te]); pred_null[te] <- predict(m0, S[te])
  }
  y <- S[[v]]; w <- S$nsw
  r_raw <- sqrt(weighted.mean(y^2, w))
  r_nul <- sqrt(weighted.mean((y - pred_null)^2, w))
  r_sd  <- sqrt(weighted.mean((y - pred_sd)^2, w))
  cat(sprintf("  %-9s uncorrected %.4f | separation only %.4f | + sd_sep %.4f  (%+.3f%% vs sep only)\n",
              v, r_raw, r_nul, r_sd, 100*(r_sd - r_nul)/r_nul))
}

## ---- 2. prior season, the deployable form ------------------------------------------------------------
cat("\n=== 2. prior-season consistency predicting this season ===\n")
Y <- merge(S[, .(id, season, wh_aware, rv_aware, nsw, sep, sd_sep)],
           S[, .(id, season = season + 1L, sd_p = sd_sep, sep_p = sep, wh_p = wh_aware)],
           by = c("id","season"))
cat(sprintf("  %d consecutive-season pairs from %d arms\n", nrow(Y), uniqueN(Y$id)))
cat(sprintf("  consistency is itself only moderately stable: sd_sep year to year r = %+.3f\n",
            cor(Y$sd_p, Y$sd_sep)))
for (v in c("wh_aware","rv_aware")) {
  z <- cor.test(Y$sd_p, Y[[v]])
  Y2 <- copy(Y); Y2[, id := id]
  cat(sprintf("  prior sd_sep -> current %-9s r = %+.3f (p %.4f) | slope %s\n", v, z$estimate,
              z$p.value, crob(as.formula(paste(v, "~ sd_p + sep_p")), Y2, "sd_p")))
}

## ---- 3. how big is it, in the units that matter -----------------------------------------------------
cat("\n=== 3. practical size ===\n")
q <- quantile(S$sd_sep, c(.1,.25,.5,.75,.9))
b <- coef(lm(wh_aware ~ sd_sep + sep, S))["sd_sep"]
cat(sprintf("  sd_sep deciles: p10 %.2f  p25 %.2f  p50 %.2f  p75 %.2f  p90 %.2f mph\n",
            q[1], q[2], q[3], q[4], q[5]))
cat(sprintf("  slope %+.2f whiff points per mph, so p10 to p90 (%.2f mph) is worth %.2f whiff points\n",
            b, q[5]-q[1], abs(b)*(q[5]-q[1])))
cat(sprintf("  for comparison, the season-to-season SD of the aware whiff residual is %.2f points,\n",
            sd(S$wh_aware)))
cat(sprintf("  and mean separation across its own p10-p90 range is worth %.2f points.\n",
            abs(coef(lm(wh_blind ~ sep, S))["sep"]) * diff(quantile(S$sep, c(.1,.9)))))
print(S[, .(seasons = .N, sd_sep = round(mean(sd_sep),2), whiff_over_aware = round(mean(wh_aware),2),
            se = round(sd(wh_aware)/sqrt(.N),2)),
        by = .(bin = cut(sd_sep, quantile(sd_sep, 0:5/5), labels = c("most consistent","2","3","4",
               "most variable"), include.lowest = TRUE))][order(bin)], row.names = FALSE)
