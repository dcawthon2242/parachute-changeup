#!/usr/bin/env Rscript

# DOES LAST SEASON'S DRIFT PREDICT THIS SEASON'S RESULTS?
#
# Everything in the drift finding so far is same-season: the drift measure and the residual are
# computed from the same pitches. That does not make it circular - a velocity dispersion and a
# whiff outcome have no arithmetic link - but it leaves open that a pitcher having a bad year is
# scattered AND ineffective for a common reason, which is a description rather than a prediction.
#
# Lagging the measure closes it. Drift measured in season t-1 cannot be caused by anything that
# happens in season t, so if it still predicts, the relationship runs from consistency to results.
#
# Three specifications, increasingly demanding:
#
#   1  current residual on prior drift. The basic forecast.
#   2  the same, adding the PRIOR residual. Drift has to beat "he was already good last year",
#      which is the obvious rival explanation and the one that would make this uninteresting.
#   3  the same, adding prior separation, since the level and the variability are correlated.
#
# The residuals come from models fit with folds grouped by pitcher, so no arm's residual was ever
# informed by its own pitches. Split by role throughout, because the same-season result was a
# starters-only effect.

suppressPackageStartupMessages({ library(data.table) })
options(width = 205)
MDIR <- "data/statcast_model"
C <- readRDS(file.path(MDIR, "anchor_model_table2.rds")); setDT(C)
S <- readRDS(file.path(MDIR, "velo_sep_within_seasons.rds")); setDT(S)
G <- C[, .(games = uniqueN(game_pk), np = .N), by = .(pitcher, season)]
S <- merge(S, G, by = c("pitcher","season"))
S[, `:=`(id = as.character(pitcher), ch_per_game = np/games)]

Y <- merge(S[, .(id, season, wh_aware, ch_aware, sep, sd_within, sd_between, ch_per_game, nsw, nz)],
           S[, .(id, season = season + 1L, d_p = sd_between, w_p = sd_within, sep_p = sep,
                 wh_p = wh_aware, ch_p = ch_aware, cpg_p = ch_per_game)],
           by = c("id","season"))
cat(sprintf("%d consecutive-season pairs from %d arms\n", nrow(Y), uniqueN(Y$id)))

crob <- function(f, D, k) {
  D <- D[complete.cases(D[, c(all.vars(f), "id"), with = FALSE])]
  m <- lm(f, D); u <- residuals(m); X <- model.matrix(m); nc <- uniqueN(D$id)
  b <- solve(crossprod(X)); V <- b %*% crossprod(rowsum(X*u, D$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[k]); s <- unname(sqrt(diag(V))[k])
  sprintf("%+.3f (se %.3f) | p %.4f | n=%d, %d arms", e, s, 2*pt(-abs(e/s), nc-1), nrow(D), nc)
}

## ---- is drift even a stable trait? -----------------------------------------------------------------
cat("\n=== 0. drift has to persist before it can forecast ===\n")
for (v in list(c("d_p","sd_between","between-outing drift"),
               c("w_p","sd_within","within-outing repeatability"),
               c("sep_p","sep","separation level, for scale"))) {
  z <- cor.test(Y[[v[1]]], Y[[v[2]]])
  cat(sprintf("  %-28s year to year r = %+.3f (p %.4f)\n", v[3], z$estimate, z$p.value))
}

## ---- the forecast ------------------------------------------------------------------------------------
Y[, `:=`(zd = as.numeric(scale(d_p)), zw = as.numeric(scale(w_p)))]
cat("\n=== 1-3. prior-season drift predicting the current season, per 1 SD of drift ===\n")
cat("   negative = arms that were scattered last year underperform their grade this year\n")
for (out in c("wh_aware","ch_aware")) {
  cat(sprintf("\n  ---- %s ----\n", out))
  cat(sprintf("    prior drift only                 %s\n",
              crob(as.formula(paste(out, "~ zd")), Y, "zd")))
  cat(sprintf("    + prior within-outing            %s\n",
              crob(as.formula(paste(out, "~ zd + zw")), Y, "zd")))
  cat(sprintf("    + PRIOR RESIDUAL (the rival)     %s\n",
              crob(as.formula(paste(out, "~ zd + zw +",
                   if (out == "wh_aware") "wh_p" else "ch_p")), Y, "zd")))
  cat(sprintf("    + prior separation               %s\n",
              crob(as.formula(paste(out, "~ zd + zw + sep_p +",
                   if (out == "wh_aware") "wh_p" else "ch_p")), Y, "zd")))
  cat(sprintf("    starters only (>=4 CH/outing)    %s\n",
              crob(as.formula(paste(out, "~ zd + zw")), Y[cpg_p >= 4], "zd")))
  cat(sprintf("    relievers only                   %s\n",
              crob(as.formula(paste(out, "~ zd + zw")), Y[cpg_p < 4], "zd")))
}

## ---- how much of it is just persistence? ---------------------------------------------------------------
cat("\n=== how the rival explanation performs on its own ===\n")
for (out in c("wh_aware","ch_aware")) {
  pv <- if (out == "wh_aware") "wh_p" else "ch_p"
  z <- cor.test(Y[[pv]], Y[[out]])
  cat(sprintf("  prior %s residual -> current: r = %+.3f (p %.5f)\n", out, z$estimate, z$p.value))
}

## ---- the deployable table --------------------------------------------------------------------------
cat("\n=== sorted on LAST season's drift, showing THIS season's results ===\n")
Y[, bin := cut(d_p, quantile(d_p, 0:5/5), labels = c("steadiest","2","3","4","most drift"),
               include.lowest = TRUE)]
print(Y[, .(pairs = .N, prior_drift = round(mean(d_p),2),
            whiff_over = round(mean(wh_aware),2), se_w = round(sd(wh_aware)/sqrt(.N),2),
            chase_over = round(mean(ch_aware),2), se_c = round(sd(ch_aware)/sqrt(.N),2)),
        by = bin][order(bin)], row.names = FALSE)
cat("\n  starters only:\n")
print(Y[cpg_p >= 4][, .(pairs = .N, prior_drift = round(mean(d_p),2),
            whiff_over = round(mean(wh_aware),2), se_w = round(sd(wh_aware)/sqrt(.N),2),
            chase_over = round(mean(ch_aware),2), se_c = round(sd(ch_aware)/sqrt(.N),2)),
        by = .(bin = cut(d_p, quantile(d_p, 0:5/5),
               labels = c("steadiest","2","3","4","most drift"), include.lowest = TRUE))][order(bin)],
      row.names = FALSE)
saveRDS(Y, file.path(MDIR, "velo_sep_prior.rds"))
cat("\nwrote velo_sep_prior.rds\n")
