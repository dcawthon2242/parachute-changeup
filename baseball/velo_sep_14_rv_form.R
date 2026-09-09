#!/usr/bin/env Rscript

# A BETTER REPRESENTATION OF VELOCITY SEPARATION FOR RUN VALUE.
#
# On whiff, adding a linear speed_diff absorbs 77 percent of what the fastball-blind model misses
# and the leftover slope is not significant. On run value it absorbs 43 percent and the leftover is
# +0.040 per mph at p = .061. That gap is the whole motivation here: either run value carries
# something about separation that a raw mph difference does not encode, or the leftover is noise.
#
# LightGBM can already bend speed_diff nonlinearly, so this is not about nonlinearity per se. Three
# things a tree on pitch-level speed_diff genuinely cannot represent cheaply:
#
#   sep_season   the pitcher's characteristic separation. Pitch-level speed_diff carries release
#                noise; the trait a hitter prepares for is the season-level tendency. Run value is
#                so noisy per pitch that the model may never resolve the trait from the pitch.
#   dt_ms        arrival-time difference. Timing is the hitter's actual task and time is not linear
#                in velocity - the same 5 mph gap buys more milliseconds at 85 than at 95. A tree
#                can approximate this only by interacting release_speed with speed_diff.
#   sep_pct      separation as a fraction of the fastball, the scale-free version.
#
# EVALUATION. Run value has a blind R2 near .0003, so pitch-level RMSE differences are tiny and
# seed noise is large. Three criteria are reported and they answer different questions:
#   1. RMSE under pitcher-grouped folds, paired across seeds - does it predict better out of sample
#   2. absorption of the blind separation slope - does it capture the specific thing being missed
#   3. season-level correlation - does it rank pitcher-seasons better, which is where the trait lives
# Criterion 2 is the one that decides whether the representation is right; 1 is whether it matters.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
options(width = 205); MDIR <- "data/statcast_model"
CACHE <- file.path(MDIR, "rv_form_results.rds")

F <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(F)
F <- F[is.finite(speed_diff) & is.finite(release_speed) & is.finite(release_spin_rate) &
       is.finite(ax) & is.finite(az) & is.finite(release_extension) &
       is.finite(release_pos_x) & is.finite(release_pos_z) & is.finite(rv) &
       is.finite(ax_diff) & is.finite(az_diff)]
F[, lh := p_throws == "L"]
F[, `:=`(tj_ax = fifelse(lh, -ax, ax), tj_x0 = fifelse(lh, -release_pos_x, release_pos_x),
         tj_ax_diff = fifelse(lh, -ax_diff, ax_diff), sep = -speed_diff,
         fb_velo = release_speed - speed_diff)]
F[, np_ps := .N, by = .(pitcher, season)]; F <- F[np_ps >= 40]

## ---- candidate encodings ---------------------------------------------------------------------------
F[, sep_season := mean(sep), by = .(pitcher, season)]
# Arrival time to the plate, in milliseconds. Average flight speed runs about 92 percent of release
# speed; the constant scales everything equally so it does not affect the model, only readability.
D_PLATE <- 60.5 - 6.0
F[, dt_ms := 1000 * D_PLATE / 1.467 / 0.92 * (1/release_speed - 1/fb_velo)]
F[, sep_pct := 100 * sep / fb_velo]
# Total stuff separation: how far the changeup sits from the fastball in velocity AND movement.
F[, sepvec := sqrt(scale(sep)^2 + scale(tj_ax_diff)^2 + scale(az_diff)^2)]
F[, dt_season := mean(dt_ms), by = .(pitcher, season)]

cat(sprintf("%s changeups | %d pitcher-seasons | %d arms\n", format(nrow(F), big.mark = ","),
            uniqueN(F[, .(pitcher, season)]), uniqueN(F$pitcher)))
cat("\ncandidate encodings, and how much they differ from raw separation:\n")
print(F[, .(sep = round(mean(sep),2), sep_season = round(mean(sep_season),2),
            dt_ms = round(mean(dt_ms),1), sep_pct = round(mean(sep_pct),1),
            sepvec = round(mean(sepvec),2))], row.names = FALSE)
cat(sprintf("  cor(sep, dt_ms) = %+.3f   cor(sep, sep_pct) = %+.3f   cor(sep, sep_season) = %+.3f\n",
            cor(F$sep, F$dt_ms), cor(F$sep, F$sep_pct), cor(F$sep, F$sep_season)))
cat("  dt_ms is NOT a relabelling of sep: the same mph gap buys different time at different velocity\n")
print(F[, .(pitches = .N, mean_dt = round(mean(dt_ms),1)),
        by = .(sep_bucket = cut(sep, c(-Inf,6,8,10,Inf), labels = c("<6","6-8","8-10","10+")),
               fb = cut(fb_velo, c(-Inf,91,94,Inf), labels = c("soft <91","mid 91-94","hard 94+")))][
        order(sep_bucket, fb)], row.names = FALSE)

BLIND <- c("release_speed","release_spin_rate","release_extension","tj_ax","az","tj_x0","release_pos_z")
AWARE <- c(BLIND, "speed_diff", "tj_ax_diff", "az_diff")
VAR <- list(
  blind        = BLIND,
  aware        = AWARE,
  `+sep_season`= c(AWARE, "sep_season"),
  `+dt_ms`     = c(AWARE, "dt_ms"),
  `+sep_pct`   = c(AWARE, "sep_pct"),
  `+sepvec`    = c(AWARE, "sepvec"),
  everything   = c(AWARE, "sep_season", "dt_ms", "sep_pct", "sepvec", "dt_season"))

## ---- grouped-fold evaluation, several seeds ---------------------------------------------------------
oof <- function(feats, seed) {
  set.seed(seed)
  arms <- unique(F$pitcher); fa <- sample(rep(1:4, length.out = length(arms)))
  fold <- fa[match(F$pitcher, arms)]; p <- rep(NA_real_, nrow(F))
  for (f in 1:4) {
    tri <- which(fold != f); vi <- sample(length(tri), floor(.12*length(tri)))
    dtr <- lgb.Dataset(as.matrix(F[tri[-vi], ..feats]), label = F$rv[tri[-vi]])
    dva <- lgb.Dataset.create.valid(dtr, as.matrix(F[tri[vi], ..feats]), label = F$rv[tri[vi]])
    m <- lgb.train(params = list(objective = "regression", metric = "l2", learning_rate = .06,
                   num_leaves = 31, min_data_in_leaf = 300, feature_fraction = .8,
                   bagging_fraction = .8, bagging_freq = 1), data = dtr, nrounds = 800,
                   valids = list(v = dva), early_stopping_rounds = 40, verbose = -1)
    p[fold == f] <- predict(m, as.matrix(F[fold == f, ..feats]))
  }
  p
}
SEEDS <- c(11, 23, 37)
if (!file.exists(CACHE)) {
  RES <- list()
  for (nm in names(VAR)) for (s in SEEDS) {
    p <- oof(VAR[[nm]], s)
    RES[[paste(nm, s)]] <- list(variant = nm, seed = s, rmse = sqrt(mean((F$rv - p)^2)),
                                res = F$rv - p, pred = p)
    cat(sprintf("  %-12s seed %2d  rmse %.7f\n", nm, s, sqrt(mean((F$rv - p)^2))))
  }
  saveRDS(RES, CACHE)
} else cat("\n(using cached fits)\n")
RES <- readRDS(CACHE)

crob <- function(D, f, k) {
  m <- lm(f, D); u <- residuals(m); X <- model.matrix(m); nc <- uniqueN(D$id)
  b <- solve(crossprod(X)); V <- b %*% crossprod(rowsum(X*u, D$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[k]); s <- unname(sqrt(diag(V))[k])
  c(est = e, se = s, p = 2*pt(-abs(e/s), nc-1))
}
KEY <- F[, .(pitcher, season, sep, rv)]
score <- function(r) {
  S <- cbind(KEY, resid = r$res, pred = r$pred)[
       , .(np = .N, sepm = mean(sep), rvres = 100*mean(resid),
           act = 100*mean(rv), prd = 100*mean(pred)), by = .(pitcher, season)]
  S[, id := as.character(pitcher)]
  c(crob(S, rvres ~ sepm, "sepm"), season_r = cor(S$prd, S$act))
}
T <- rbindlist(lapply(RES, function(r) {
  sc <- score(r)
  data.table(variant = r$variant, seed = r$seed, rmse = r$rmse, slope = sc["est"],
             p = sc["p"], season_r = sc["season_r"]) }))
A <- T[, .(rmse = mean(rmse), slope = mean(slope), p = mean(p), season_r = mean(season_r)),
       by = variant]
base <- A[variant == "aware", rmse]; bslope <- A[variant == "blind", slope]
A[, `:=`(d_rmse = 1e6*(rmse - base), absorbed = paste0(round(100*(1 - slope/bslope)), "%"))]
setorder(A, rmse)
cat("\n=== results, averaged over 3 pitcher-grouped seeds ===\n")
cat("   d_rmse is RMSE minus the aware baseline, in units of 1e-6; negative is better\n")
cat("   absorbed is the share of the BLIND model's separation slope the variant removes\n\n")
print(A[, .(variant, rmse = round(rmse,7), d_rmse = round(d_rmse,2),
            resid_slope = round(slope,4), p = round(p,4), absorbed,
            season_r = round(season_r,4))], row.names = FALSE)

cat("\n=== paired across seeds: is any RMSE gain bigger than seed noise? ===\n")
W <- dcast(T, seed ~ variant, value.var = "rmse")
for (nm in setdiff(names(VAR), c("blind","aware"))) {
  d <- W[["aware"]] - W[[nm]]
  cat(sprintf("  %-12s mean gain %+.2e over 3 seeds  (per-seed: %s)\n",
              nm, mean(d), paste(sprintf("%+.1e", d), collapse = " ")))
}
