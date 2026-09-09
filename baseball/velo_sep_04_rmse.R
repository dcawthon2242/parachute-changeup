#!/usr/bin/env Rscript

# DOES CHANGEUP VELOCITY CONSISTENCY ADD ANYTHING? RMSE, HONESTLY MEASURED.
#
# Baseline RMSE for the two models, then the same models with the pitcher-season standard deviation
# of velocity separation added as a feature.
#
# ONE THING TO KNOW BEFORE READING THE FEATURE. speed_diff is computed against a season-constant
# fastball velocity, so sd(speed_diff) and sd(release_speed) correlate at exactly 1.000. This is
# changeup velocity consistency, not variability in the gap against the fastball he actually threw
# that inning. Median 1.20 mph, IQR 1.04 to 1.40.
#
# AND ONE THING ABOUT HOW IT IS TESTED. The feature is a pitcher-season aggregate broadcast onto
# every pitch in that season. With pitch-level random folds the model can recover pitcher-season
# identity from it and effectively look up that season's own whiff rate, which shows up as an
# improvement that would not survive contact with a new pitcher. This project has been burned by
# exactly that aggregation trap before, so every model is fit twice:
#
#   random   folds assigned per pitch - the leaky version, reported to size the leak
#   grouped  folds assigned per PITCHER, so no arm appears in both train and test
#
# Pitch-level RMSE is also nearly all irreducible noise, so a real but small feature moves it in
# the fourth decimal. Season-level RMSE - aggregate the out-of-fold predictions to pitcher-season
# and compare against the observed rate - is far more sensitive and is reported alongside.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(43); options(width = 205)
MDIR <- "data/statcast_model"

F <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(F)
F <- F[is.finite(speed_diff) & is.finite(release_speed) & is.finite(release_spin_rate) &
       is.finite(ax) & is.finite(az) & is.finite(release_extension) &
       is.finite(release_pos_x) & is.finite(release_pos_z) & is.finite(rv)]
F[, lh := p_throws == "L"]
F[, `:=`(tj_ax = fifelse(lh, -ax, ax), tj_x0 = fifelse(lh, -release_pos_x, release_pos_x),
         tj_ax_diff = fifelse(lh, -ax_diff, ax_diff), sep = -speed_diff)]
F[, `:=`(sd_sep = sd(speed_diff), n_ps = .N), by = .(pitcher, season)]
F <- F[n_ps >= 40 & is.finite(sd_sep)]

BLIND <- c("release_speed","release_spin_rate","release_extension","tj_ax","az","tj_x0","release_pos_z")
AWARE <- c(BLIND, "speed_diff", "tj_ax_diff", "az_diff")
WSD   <- c(AWARE, "sd_sep")

cat(sprintf("%s changeups, %d pitcher-seasons, %d arms\n",
            format(nrow(F), big.mark = ","), uniqueN(F[, .(pitcher, season)]), uniqueN(F$pitcher)))
cat(sprintf("sd of separation: median %.2f mph, IQR %.2f-%.2f\n\n",
            median(F$sd_sep), quantile(F$sd_sep,.25), quantile(F$sd_sep,.75)))

## ---- fit one configuration ------------------------------------------------------------------------
run <- function(D, feats, label, obj, scheme, tag) {
  K <- 4
  if (scheme == "random") fold <- sample(rep(1:K, length.out = nrow(D)))
  else { arms <- unique(D$pitcher); fa <- sample(rep(1:K, length.out = length(arms)))
         fold <- fa[match(D$pitcher, arms)] }
  y <- D[[label]]; p <- rep(NA_real_, nrow(D))
  for (f in 1:K) {
    tri <- which(fold != f); n <- length(tri); vi <- sample(n, floor(.12*n))
    dtr <- lgb.Dataset(as.matrix(D[tri[-vi], ..feats]), label = y[tri[-vi]])
    dva <- lgb.Dataset.create.valid(dtr, as.matrix(D[tri[vi], ..feats]), label = y[tri[vi]])
    m <- lgb.train(params = list(objective = obj,
                   metric = if (obj == "binary") "binary_logloss" else "l2",
                   learning_rate = .06, num_leaves = 31, min_data_in_leaf = 300,
                   feature_fraction = .8, bagging_fraction = .8, bagging_freq = 1),
                   data = dtr, nrounds = 1500, valids = list(v = dva),
                   early_stopping_rounds = 50, verbose = -1)
    p[fold == f] <- predict(m, as.matrix(D[fold == f, ..feats]))
  }
  S <- data.table(pitcher = D$pitcher, season = D$season, y = y, p = p
        )[, .(n = .N, ay = mean(y), ap = mean(p)), by = .(pitcher, season)]
  data.table(model = tag, scheme = scheme, feats = length(feats),
             rmse = sqrt(mean((y - p)^2)), r2 = 1 - var(y - p)/var(y),
             season_rmse = sqrt(weighted.mean((S$ay - S$ap)^2, S$n)))
}

W <- F[is_swing == 1 & is.finite(whiff)]
cat(sprintf("whiff fit on %s swings | run value fit on %s pitches\n\n",
            format(nrow(W), big.mark = ","), format(nrow(F), big.mark = ",")))

JOBS <- list(list("blind", BLIND), list("aware", AWARE), list("aware + sd_sep", WSD))
OUT <- rbindlist(lapply(c("random","grouped"), function(sc)
  rbindlist(lapply(JOBS, function(j) {
    set.seed(43); a <- run(W, j[[2]], "whiff", "binary", sc, j[[1]])[, out := "whiff"]
    set.seed(43); b <- run(F, j[[2]], "rv", "regression", sc, j[[1]])[, out := "rv"]
    rbind(a, b) }))))

report <- function(o, ttl, dig) {
  X <- OUT[out == o]
  cat(sprintf("\n=== %s ===\n", ttl))
  for (sc in c("random","grouped")) {
    Y <- X[scheme == sc]; base <- Y[model == "blind"]
    cat(sprintf("  %s folds\n", sc))
    for (i in 1:nrow(Y)) cat(sprintf("    %-16s %d feats | RMSE %.*f | R2 %+.5f | season RMSE %.5f\n",
        Y$model[i], Y$feats[i], dig, Y$rmse[i], Y$r2[i], Y$season_rmse[i]))
    aw <- Y[model == "aware"]; sd <- Y[model == "aware + sd_sep"]
    cat(sprintf("    -> adding sd_sep to aware: RMSE %+.6f (%+.4f%%), season RMSE %+.6f (%+.3f%%)\n",
        sd$rmse - aw$rmse, 100*(sd$rmse - aw$rmse)/aw$rmse,
        sd$season_rmse - aw$season_rmse, 100*(sd$season_rmse - aw$season_rmse)/aw$season_rmse))
  }
}
report("whiff", "WHIFF (binary; RMSE on 0/1 is the root Brier score)", 6)
report("rv",    "RUN VALUE (per pitch)", 6)

## ---- is there any signal in sd_sep at all? --------------------------------------------------------
# Cheap screen independent of the models: does consistency correlate with the residual the
# separation-aware model leaves behind, at the pitcher-season level, clustered on the arm?
cat("\n=== direct screen: sd_sep against the aware model's residual, pitcher-season ===\n")
L <- readRDS(file.path(MDIR, "velo_sep_resid.rds"))
S <- merge(as.data.table(L$W)[, .(nsw = .N, wh_aware = 100*mean(whiff - p_aware)),
                              by = .(pitcher, season)],
           as.data.table(L$R)[, .(rv_aware = 100*mean(rv - q_aware)), by = .(pitcher, season)],
           by = c("pitcher","season"))
S <- merge(S, unique(F[, .(pitcher, season, sd_sep, sep = mean(sep)), by = .(pitcher, season)][
             , .(pitcher, season, sd_sep, sep)]), by = c("pitcher","season"))[nsw >= 40]
S[, id := as.character(pitcher)]
for (v in c("wh_aware","rv_aware")) {
  f <- as.formula(paste(v, "~ sd_sep + sep"))
  m <- lm(f, S); u <- residuals(m); X <- model.matrix(m)
  b <- solve(crossprod(X)); nc <- uniqueN(S$id)
  V <- b %*% crossprod(rowsum(X * u, S$id)) %*% b * (nc/(nc-1))
  e <- coef(m)["sd_sep"]; s <- sqrt(diag(V))["sd_sep"]
  z <- suppressWarnings(cor.test(S$sd_sep, S[[v]], method = "spearman", exact = FALSE))
  cat(sprintf("  %-9s slope %+.4f per mph of sd (se %.4f) | p %.4f | raw spearman %+.3f (p %.3f)\n",
              v, e, s, 2*pt(-abs(e/s), nc-1), z$estimate, z$p.value))
}
cat(sprintf("  (cor between sd_sep and mean separation = %+.3f, so they are not redundant)\n",
            cor(S$sd_sep, S$sep)))

fwrite(OUT, file.path(MDIR, "article_assets", "ext_velo_sep_rmse.csv"))
cat("\nwrote ext_velo_sep_rmse.csv\n")
