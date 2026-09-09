#!/usr/bin/env Rscript

# DO TRADITIONAL STUFF MODELS UNDERRATE BIG-SEPARATION CHANGEUPS?
#
# A stuff model grades a pitch on what the pitch is: velocity, spin, movement, extension, release
# point. It does not know what the pitcher's fastball looks like. A changeup thrown 12 mph off a
# 96 mph heater and the same changeup thrown 5 mph off an 89 mph heater get the same grade, and
# only one of them is actually hard to hit. If that is the gap Martinez lives in, then velocity
# separation should predict what a fastball-blind model misses.
#
# The project's own whiff model cannot answer this. It already carries speed_diff as a feature, so
# separation is priced in by construction and its residual is silent on the question. So two models
# are fit here on identical rows and identical folds, differing only in whether they can see the
# fastball:
#
#   BLIND  release_speed, spin rate, extension, ax, az, release x/z   - a stuff model
#   AWARE  the same, plus speed_diff, ax_diff, az_diff                - separation priced in
#
# Neither sees location or count. A stuff model grades the pitch, not the execution, and leaving
# location out of both means it cancels in the comparison rather than contaminating one side.
#
# Two outcomes, because "outperforming the model" can mean either. Whiff is the project's currency;
# run value is what actually shows up in the standings.
#
# The prediction if the hypothesis is right: separation correlates with the BLIND residual and that
# correlation shrinks toward zero under AWARE. If it survives AWARE too, then linear separation is
# not the missing ingredient and something else is going on.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(29); options(width = 205)
MDIR <- "data/statcast_model"; CACHE <- file.path(MDIR, "velo_sep_resid.rds")

F <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(F)
F <- F[is.finite(speed_diff) & is.finite(release_speed) & is.finite(release_spin_rate) &
       is.finite(ax) & is.finite(az) & is.finite(release_extension) &
       is.finite(release_pos_x) & is.finite(release_pos_z) & is.finite(rv)]
F[, lh := p_throws == "L"]
F[, `:=`(tj_ax = fifelse(lh, -ax, ax), tj_x0 = fifelse(lh, -release_pos_x, release_pos_x),
         tj_ax_diff = fifelse(lh, -ax_diff, ax_diff), sep = -speed_diff)]

BLIND <- c("release_speed","release_spin_rate","release_extension","tj_ax","az","tj_x0","release_pos_z")
AWARE <- c(BLIND, "speed_diff", "tj_ax_diff", "az_diff")

cat(sprintf("changeups with a four-seam anchor: %s across %d pitcher-seasons, %d arms\n",
            format(nrow(F), big.mark = ","), uniqueN(F[, .(pitcher, season)]), uniqueN(F$pitcher)))
cat(sprintf("velocity separation: median %.1f mph, IQR %.1f-%.1f, 90th pct %.1f, max %.1f\n\n",
            median(F$sep), quantile(F$sep,.25), quantile(F$sep,.75),
            quantile(F$sep,.90), max(F$sep)))

fit_oof <- function(D, feats, label, obj) {
  K <- 4; fold <- sample(rep(1:K, length.out = nrow(D))); p <- rep(NA_real_, nrow(D))
  y <- D[[label]]
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
  p
}

if (!file.exists(CACHE)) {
  W <- F[is_swing == 1 & is.finite(whiff)]
  cat(sprintf("fitting on %s swings (whiff) and %s pitches (run value)\n",
              format(nrow(W), big.mark = ","), format(nrow(F), big.mark = ",")))
  set.seed(29); W[, p_blind := fit_oof(W, BLIND, "whiff", "binary")]
  set.seed(29); W[, p_aware := fit_oof(W, AWARE, "whiff", "binary")]
  set.seed(31); F[, q_blind := fit_oof(F, BLIND, "rv", "regression")]
  set.seed(31); F[, q_aware := fit_oof(F, AWARE, "rv", "regression")]
  cat(sprintf("  whiff  blind R2 %.4f | aware R2 %.4f  (gain %+.4f)\n",
              1 - var(W$whiff - W$p_blind)/var(W$whiff),
              1 - var(W$whiff - W$p_aware)/var(W$whiff),
              var(W$whiff - W$p_blind)/var(W$whiff) - var(W$whiff - W$p_aware)/var(W$whiff)))
  cat(sprintf("  rv     blind R2 %.5f | aware R2 %.5f\n",
              1 - var(F$rv - F$q_blind)/var(F$rv), 1 - var(F$rv - F$q_aware)/var(F$rv)))
  saveRDS(list(
    W = W[, .(pitcher, player_name, season, sep, whiff, p_blind, p_aware)],
    R = F[, .(pitcher, player_name, season, sep, rv, q_blind, q_aware)]), CACHE)
} else cat("(using cached fits)\n")

L <- readRDS(CACHE); W <- L$W; R <- L$R; setDT(W); setDT(R)

## ---- pitcher-season table -----------------------------------------------------------------------
S <- merge(
  W[, .(nsw = .N, sep = mean(sep), wh = 100*mean(whiff),
        wh_blind = 100*mean(whiff - p_blind), wh_aware = 100*mean(whiff - p_aware)),
    by = .(pitcher, player_name, season)],
  R[, .(np = .N, rv100 = 100*mean(rv),
        rv_blind = 100*mean(rv - q_blind), rv_aware = 100*mean(rv - q_aware)),
    by = .(pitcher, season)], by = c("pitcher","season"))
S <- S[nsw >= 40]
S[, id := as.character(pitcher)]
cat(sprintf("\n%d pitcher-seasons with 40+ changeup swings, %d arms\n", nrow(S), uniqueN(S$id)))

crob <- function(f, D) {
  D <- D[complete.cases(D[, c(all.vars(f), "id"), with = FALSE])]
  m <- lm(f, D); u <- residuals(m); X <- model.matrix(m); cl <- D$id
  nc <- uniqueN(cl); b <- solve(crossprod(X))
  V <- b %*% crossprod(rowsum(X * u, cl)) %*% b * (nc/(nc-1))
  k <- "sep"; e <- coef(m)[k]; s <- sqrt(diag(V))[k]
  sprintf("%+.4f per mph (se %.4f) | p %.4f | %d arms", e, s, 2*pt(-abs(e/s), nc-1), nc)
}

## ---- 1. does separation predict what each model misses? -----------------------------------------
cat("\n=== 1. separation vs residual, one row per pitcher-season, clustered on the arm ===\n")
# rv here is signed from the PITCHER's side: whiffs average +0.112, line drives -0.274. So both
# outcomes point the same way and a positive residual is the pitch beating its grade.
cat("   positive = the model UNDERRATES big-separation changeups, on both outcomes\n\n")
for (v in c("wh_blind","wh_aware","rv_blind","rv_aware")) {
  z <- suppressWarnings(cor.test(S$sep, S[[v]], method = "spearman", exact = FALSE))
  cat(sprintf("  %-9s  spearman %+.3f (p %.4f)   slope %s\n",
              v, z$estimate, z$p.value, crob(as.formula(paste(v, "~ sep")), S)))
}

## ---- 2. is it the separation, or just a slow changeup? -------------------------------------------
# A stuff model already sees absolute changeup velocity, and slow changeups tend to be big-
# separation changeups, so some of the effect above could be the model simply mishandling slow
# pitches. Separation is fb_velo - ch_velo exactly, so all three cannot go in the same regression -
# that is rank deficient and produces meaningless standard errors. Holding the FASTBALL fixed is
# the identified version and it is also the question worth asking: given two pitchers with the same
# heater, does the one who takes more off it beat the grade by more?
cat("\n=== 2. separation with the fastball held fixed (changeup velo is then implied) ===\n")
CV <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(CV)
CV <- CV[is.finite(speed_diff), .(ch_velo = mean(release_speed),
                                  fb_velo = mean(release_speed - speed_diff)),
         by = .(pitcher, season)]
S <- merge(S, CV, by = c("pitcher","season"))
cat(sprintf("  check: cor(sep, fb_velo - ch_velo) = %.4f, so only two of the three are free\n\n",
            cor(S$sep, S$fb_velo - S$ch_velo)))
for (v in c("wh_blind","wh_aware","rv_blind","rv_aware"))
  cat(sprintf("  %-9s  %s\n", v, crob(as.formula(paste(v, "~ sep + fb_velo")), S)))
cat("\n   and the mirror image, holding the CHANGEUP fixed (fastball velo implied):\n")
for (v in c("wh_blind","rv_blind"))
  cat(sprintf("  %-9s  %s\n", v, crob(as.formula(paste(v, "~ sep + ch_velo")), S)))

## ---- 3. deciles ------------------------------------------------------------------------------------
cat("\n=== 3. by separation decile (pitcher-seasons, weighted by swings) ===\n")
S[, dec := cut(sep, quantile(sep, 0:10/10), labels = 1:10, include.lowest = TRUE)]
print(S[, .(seasons = .N, sep = round(mean(sep),1), swings = sum(nsw),
            whiff = round(weighted.mean(wh, nsw),1),
            blind_over = round(weighted.mean(wh_blind, nsw),2),
            aware_over = round(weighted.mean(wh_aware, nsw),2),
            rv_blind = round(weighted.mean(rv_blind, np),2),
            rv_aware = round(weighted.mean(rv_aware, np),2)), by = dec][order(dec)],
      row.names = FALSE)

## ---- 4. the top of the distribution ------------------------------------------------------------
cat("\n=== 4. the biggest separators, 10+ mph, vs everyone else ===\n")
for (v in c("wh_blind","wh_aware","rv_blind","rv_aware")) {
  a <- S[sep >= 10]; b <- S[sep < 10]
  t <- t.test(a[[v]], b[[v]])
  cat(sprintf("  %-9s  10+ mph %+.3f (n=%d)  rest %+.3f (n=%d)  diff %+.3f  p %.4f\n",
              v, mean(a[[v]]), nrow(a), mean(b[[v]]), nrow(b), diff(rev(t$estimate)), t$p.value))
}

## ---- 5. Martinez ---------------------------------------------------------------------------------
cat("\n=== 5. Nick Martinez, every qualifying season ===\n")
print(S[grepl("Martinez, Nick", player_name)][order(season),
        .(season, swings = nsw, sep = round(sep,1), ch = round(ch_velo,1), fb = round(fb_velo,1),
          whiff = round(wh,1), blind_over = round(wh_blind,2), aware_over = round(wh_aware,2),
          rv_blind = round(rv_blind,2), rv_aware = round(rv_aware,2))], row.names = FALSE)

cat("\n=== the 15 biggest separators with 200+ swings ===\n")
print(S[nsw >= 200][order(-sep)][1:15,
        .(player_name = sub(",.*","",player_name), season, swings = nsw, sep = round(sep,1),
          whiff = round(wh,1), blind_over = round(wh_blind,2), aware_over = round(wh_aware,2),
          rv_blind = round(rv_blind,2))], row.names = FALSE)

## ---- 6. how much of the blind gap does separation-awareness actually recover? --------------------
cat("\n=== 6. share of the blind model's miss that pricing in separation recovers ===\n")
for (p in list(c("wh_blind","wh_aware","whiff pts per mph"), c("rv_blind","rv_aware","rv100 per mph"))) {
  b <- coef(lm(as.formula(paste(p[1], "~ sep")), S))["sep"]
  a <- coef(lm(as.formula(paste(p[2], "~ sep")), S))["sep"]
  cat(sprintf("  %-18s blind %+.4f -> aware %+.4f   %.0f%% absorbed, %.0f%% still unexplained\n",
              p[3], b, a, 100*(1 - a/b), 100*a/b))
}

## ---- 7. is this exploitable, or just a description? ----------------------------------------------
# A correction is only worth applying if the thing driving it is stable. Separation is a pitcher
# trait and should carry over year to year; the residual it predicts has to carry over as well, or
# the model is being blamed for noise.
cat("\n=== 7. year-over-year persistence, same arm in consecutive seasons ===\n")
Y <- merge(S[, .(id, season, sep, wh_blind, wh_aware)],
           S[, .(id, season = season + 1L, sep_p = sep, blind_p = wh_blind, aware_p = wh_aware)],
           by = c("id","season"))
cat(sprintf("  %d consecutive-season pairs from %d arms\n", nrow(Y), uniqueN(Y$id)))
for (v in list(c("sep_p","sep","velocity separation itself"),
               c("blind_p","wh_blind","blind-model overperformance"),
               c("aware_p","wh_aware","aware-model overperformance"))) {
  z <- cor.test(Y[[v[1]]], Y[[v[2]]])
  cat(sprintf("  %-30s r = %+.3f (p %.4f)\n", v[3], z$estimate, z$p.value))
}
cat("\n  does last year's separation predict THIS year's blind overperformance?\n")
z <- cor.test(Y$sep_p, Y$wh_blind)
cat(sprintf("  prior-season separation vs current blind residual: r = %+.3f (p %.5f, n = %d)\n",
            z$estimate, z$p.value, nrow(Y)))

saveRDS(S, file.path(MDIR, "velo_sep_seasons.rds"))
cat(sprintf("\nwrote %s\n", file.path(MDIR, "velo_sep_seasons.rds")))
