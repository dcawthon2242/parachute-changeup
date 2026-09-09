#!/usr/bin/env Rscript

# HOW BIG IS THE MISS, AND CAN A MODEL BE TAUGHT TO STOP MAKING IT?
#
# Three questions, in the order that matters.
#
# 1. SIZE. Every archetype effect so far has been reported as a residual in percentage points,
#    which is not a unit anyone makes decisions in. Here the same gaps are expressed as predicted
#    versus actual whiff rate, as a fraction of the model's own predicted spread, and as runs per
#    season at the pitcher's real changeup usage.
#
# 2. FORECASTABILITY. A model can only absorb a miss that is a stable property of the pitcher. If
#    the residual is one-season noise, no feature set recovers it and the honest answer is that
#    the model is already doing as well as it can. Year-over-year persistence and split-half
#    reliability settle this before any refitting happens.
#
# 3. ABSORPTION. Three nested feature sets are fit with PITCHER-GROUPED folds, so no arm appears
#    in both training and validation. That grouping is the whole point: random pitch-level folds
#    let the model memorise a pitcher and report a gain that evaporates on a new arm.
#
#      STUFF    what the pitch is - the standard stuff-model surface
#      ARSENAL  + the fastball-relative deltas (separation, movement gaps)
#      SHAPE    + spin-axis gap, seam deviation magnitude and sign, arm slot, arm-angle gap
#
#    Judged three ways, because they can disagree: global pitch-level loss, pitcher-season RMSE,
#    and - the one that actually matters here - whether the archetype's residual gap closes.
#
# 4. TEMPORAL HOLDOUT. Trained on 2020-2023, tested on 2024-2026, which is the situation a team
#    is really in. A feature that helps in-sample and dies forward is worthless.

suppressPackageStartupMessages({ library(data.table); library(lightgbm); library(bit64) })
set.seed(414); options(width = 215); MDIR <- "data/statcast_model"
CACHE <- file.path(MDIR, "absorb_preds.rds")

## ---------------------------------------------------------------------------------------------
## Feature construction
P <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(P)
P <- P[is.finite(release_speed) & is.finite(release_spin_rate) & is.finite(release_extension) &
       is.finite(ax) & is.finite(az) & is.finite(release_pos_x) & is.finite(release_pos_z) &
       is.finite(sax) & is.finite(cax) & is.finite(speed_diff) & is.finite(rv)]
P[, lh := p_throws == "L"]
P[, `:=`(tj_ax = fifelse(lh, -ax, ax), tj_x0 = fifelse(lh, -release_pos_x, release_pos_x),
         tj_ax_diff = fifelse(lh, -ax_diff, ax_diff), sep = -speed_diff,
         mx = fifelse(lh, -ax, ax), mz = az + 32.174)]
P[, spin_axis := (atan2(sax, cax) * 180/pi) %% 360]
P[, spin_axis := fifelse(lh, (360 - spin_axis) %% 360, spin_axis)]
wrap <- function(d) ((d + 180) %% 360) - 180
P[, dev_raw := wrap(atan2(mx, mz) * 180/pi + (spin_axis - 180))]
P[, dev := wrap(dev_raw - median(dev_raw, na.rm = TRUE))][, absdev := abs(dev)]
# out of zone, for the chase model
P[, oz := abs(plate_x) > 0.83 | plate_z < z_rel_bot | plate_z > z_rel_top]
P[, id := as.character(pitcher)]

STUFF   <- c("release_speed","release_spin_rate","release_extension","tj_ax","az","tj_x0","release_pos_z")
ARSENAL <- c(STUFF, "speed_diff", "tj_ax_diff", "az_diff")
SHAPE   <- c(ARSENAL, "axis_diff", "dev", "absdev", "arm_angle", "arm_diff")
SETS <- list(STUFF = STUFF, ARSENAL = ARSENAL, SHAPE = SHAPE)

## ---------------------------------------------------------------------------------------------
## Archetype membership, rebuilt on this row set
M <- readRDS(file.path(MDIR, "archetype_roster.rds")); setDT(M)
TAGS <- c("A2_wide","A4_broad","A5_seam","A6_extreme","A7_trad","A8_mismatch")
LAB <- c(A2_wide = "matched axis (optimized)", A4_broad = "axis + separation",
         A5_seam = "seam-shifted", A6_extreme = "extreme seam shift",
         A7_trad = "traditional (low seam dev)", A8_mismatch = "mismatch")
KEY <- M[, c("id","season","nm","nsw", TAGS), with = FALSE]

grp_oof <- function(D, feats, label, obj, K = 5) {
  arms <- unique(D$id); fa <- data.table(id = arms, fold = sample(rep(1:K, length.out = length(arms))))
  D <- merge(D, fa, by = "id", sort = FALSE); y <- D[[label]]; p <- rep(NA_real_, nrow(D))
  for (f in 1:K) {
    tri <- which(D$fold != f)
    va <- sample(tri, floor(.10 * length(tri))); tr <- setdiff(tri, va)
    dtr <- lgb.Dataset(as.matrix(D[tr, ..feats]), label = y[tr])
    dva <- lgb.Dataset.create.valid(dtr, as.matrix(D[va, ..feats]), label = y[va])
    m <- lgb.train(params = list(objective = obj,
                     metric = if (obj == "binary") "binary_logloss" else "l2",
                     learning_rate = .06, num_leaves = 31, min_data_in_leaf = 300,
                     feature_fraction = .8, bagging_fraction = .8, bagging_freq = 1),
                   data = dtr, nrounds = 2000, valids = list(v = dva),
                   early_stopping_rounds = 60, verbose = -1)
    p[D$fold == f] <- predict(m, as.matrix(D[D$fold == f, ..feats]))
  }
  list(p = p, D = D)
}

if (!file.exists(CACHE)) {
  W <- P[is_swing == 1 & is.finite(whiff)]
  C <- P[oz == TRUE & is.finite(is_swing)]
  OUT <- list()
  for (s in names(SETS)) {
    cat(sprintf("fitting %s ...\n", s))
    rw <- grp_oof(W, SETS[[s]], "whiff", "binary")
    rc <- grp_oof(C, SETS[[s]], "is_swing", "binary")
    rr <- grp_oof(P, SETS[[s]], "rv", "regression")
    OUT[[s]] <- list(
      w = data.table(id = rw$D$id, season = rw$D$season, y = rw$D$whiff, p = rw$p),
      c = data.table(id = rc$D$id, season = rc$D$season, y = rc$D$is_swing, p = rc$p),
      r = data.table(id = rr$D$id, season = rr$D$season, y = rr$D$rv, p = rr$p))
  }
  saveRDS(OUT, CACHE)
} else OUT <- readRDS(CACHE)

agg <- function(X) X[, .(n = .N, act = mean(y), pred = mean(p)), by = .(id, season)]

## =============================================================================================
cat("\n=== 1. HOW FAR OFF IS THE STUFF MODEL, IN UNITS PEOPLE USE ===\n\n")
AW <- merge(agg(OUT$STUFF$w), KEY, by = c("id","season"))
AC <- merge(agg(OUT$STUFF$c), KEY, by = c("id","season"))
AR <- merge(agg(OUT$STUFF$r), KEY, by = c("id","season"))
sdpred <- AW[n >= 40, sd(pred)] * 100

SZ <- rbindlist(lapply(TAGS, function(t) {
  w <- AW[get(t) == TRUE & n >= 40]; cc <- AC[get(t) == TRUE & n >= 40]; r <- AR[get(t) == TRUE]
  # runs per season: the season's miss per pitch, multiplied back up by that season's pitch count
  runs <- r$n * (r$act - r$pred)
  data.table(archetype = LAB[[t]], seasons = nrow(w), pitches = round(mean(r$n)),
             pred_wh = 100*weighted.mean(w$pred, w$n), act_wh = 100*weighted.mean(w$act, w$n),
             pred_ch = 100*weighted.mean(cc$pred, cc$n), act_ch = 100*weighted.mean(cc$act, cc$n),
             runs_med = median(runs), runs_mean = mean(runs))
}))
SZ[, `:=`(gap_wh = act_wh - pred_wh, gap_ch = act_ch - pred_ch)]
SZ[, `:=`(rel_wh = 100*gap_wh/pred_wh, sd_wh = gap_wh/sdpred)]
print(SZ[, .(archetype, seasons, pred_wh = round(pred_wh,1), act_wh = round(act_wh,1),
             gap_wh = round(gap_wh,2), `rel%` = round(rel_wh,1), `in SDs` = round(sd_wh,2),
             pred_ch = round(pred_ch,1), act_ch = round(act_ch,1), gap_ch = round(gap_ch,2),
             ch_per_yr = pitches, runs_yr = round(runs_mean,2),
             runs_yr_med = round(runs_med,2))], row.names = FALSE)
cat(sprintf("\n  For scale: the stuff model's predicted whiff rate has an SD of %.1f points across\n", sdpred))
cat("  pitcher-seasons, so the 'in SDs' column says how much of the model's own working range\n  the miss represents. Runs are summed over that season's actual changeups, sign flipped\n  so positive means the pitcher gained runs the model did not credit.\n")

## =============================================================================================
cat("\n=== 2. IS THE MISS A STABLE TRAIT OR ONE-YEAR NOISE? ===\n\n")
RS <- AW[n >= 40][, r := 100*(act - pred)]
setorder(RS, id, season)
RS[, `:=`(r_next = shift(r, -1L), s_next = shift(season, -1L), n_next = shift(n, -1L)), by = id]
PR <- RS[!is.na(r_next) & s_next == season + 1]
memb <- function(t) PR[get(t) == TRUE]
cat(sprintf("  %-30s %5s  %7s  %7s\n", "population", "pairs", "r(t,t+1)", "p"))
rows <- c(list(all = PR), setNames(lapply(TAGS, memb), TAGS))
for (nmx in names(rows)) {
  D <- rows[[nmx]]
  if (nrow(D) < 8) next
  ct <- cor.test(D$r, D$r_next)
  cat(sprintf("  %-30s %5d  %+7.3f  %7.4f\n",
              if (nmx == "all") "every changeup season" else LAB[[nmx]], nrow(D),
              unname(ct$estimate), ct$p.value))
}
# reliability: how much of the observed season-to-season variance is signal
bin_noise <- AW[n >= 40, mean(100^2 * act*(1-act)/n)]
obs_var <- RS[, var(r)]
cat(sprintf("\n  observed variance of the season residual   %6.1f\n", obs_var))
cat(sprintf("  variance expected from binomial noise alone %6.1f  (%.0f%% of observed)\n",
            bin_noise, 100*bin_noise/obs_var))
cat(sprintf("  implied signal share                        %6.1f  (%.0f%%)\n",
            obs_var - bin_noise, 100*(obs_var - bin_noise)/obs_var))

## =============================================================================================
cat("\n=== 3. DOES A RICHER FEATURE SET ABSORB IT? (pitcher-grouped folds) ===\n\n")
ll <- function(y, p) { p <- pmin(pmax(p, 1e-9), 1-1e-9); -mean(y*log(p) + (1-y)*log(1-p)) }
cat(sprintf("  %-9s | %-24s | %-24s | %s\n", "", "whiff", "chase", "run value"))
cat(sprintf("  %-9s | %10s %12s | %10s %12s | %10s %12s\n", "features",
            "logloss", "season RMSE", "logloss", "season RMSE", "RMSE", "season RMSE"))
for (s in names(SETS)) {
  w <- OUT[[s]]$w; cc <- OUT[[s]]$c; r <- OUT[[s]]$r
  aw <- agg(w)[n >= 40]; ac <- agg(cc)[n >= 40]; ar <- agg(r)[n >= 40]
  cat(sprintf("  %-9s | %10.5f %12.3f | %10.5f %12.3f | %10.5f %12.4f\n", s,
              ll(w$y, w$p), sqrt(mean((100*(aw$act - aw$pred))^2)),
              ll(cc$y, cc$p), sqrt(mean((100*(ac$act - ac$pred))^2)),
              sqrt(mean((r$y - r$p)^2)), sqrt(mean((100*(ar$act - ar$pred))^2))))
}

cat("\n  --- and the part that matters: does the archetype gap close? ---\n")
crob <- function(D, f, k) {
  environment(f) <- environment()
  D <- D[complete.cases(D[, c(all.vars(f), "id", "n"), with = FALSE])]
  wt <- as.numeric(D$n); m <- lm(f, D, weights = wt)
  u <- residuals(m)*sqrt(wt); X <- model.matrix(m)*sqrt(wt)
  nc <- uniqueN(D$id); b <- solve(crossprod(X))
  V <- b %*% crossprod(rowsum(X*u, D$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[k]); se <- unname(sqrt(diag(V))[k]); c(e, se, 2*pt(-abs(e/se), nc-1))
}
GAP <- rbindlist(lapply(TAGS, function(t) {
  out <- list(archetype = LAB[[t]])
  for (s in names(SETS)) {
    A <- merge(agg(OUT[[s]]$w), KEY, by = c("id","season"))[n >= 40]
    A[, `:=`(r = 100*(act - pred), gg = get(t))]
    e <- crob(A, r ~ gg, "ggTRUE")
    out[[paste0(s, "_e")]] <- e[1]; out[[paste0(s, "_p")]] <- e[3]
  }
  as.data.table(out)
}))
GAP[, absorbed := 100*(1 - SHAPE_e/STUFF_e)]
print(GAP[, .(archetype, stuff = round(STUFF_e,2), p1 = round(STUFF_p,4),
              arsenal = round(ARSENAL_e,2), p2 = round(ARSENAL_p,4),
              shape = round(SHAPE_e,2), p3 = round(SHAPE_p,4),
              `absorbed%` = round(absorbed,0))], row.names = FALSE)

## =============================================================================================
cat("\n=== 4. TRAIN ON 2020-2023, TEST ON 2024-2026 ===\n\n")
TC <- file.path(MDIR, "absorb_temporal.rds")
if (!file.exists(TC)) {
  W <- P[is_swing == 1 & is.finite(whiff)]
  tr <- W[season <= 2023]; te <- W[season >= 2024]
  TE <- list()
  for (s in names(SETS)) {
    ft <- SETS[[s]]; vi <- sample(nrow(tr), floor(.1*nrow(tr)))
    dtr <- lgb.Dataset(as.matrix(tr[-vi, ..ft]), label = tr$whiff[-vi])
    dva <- lgb.Dataset.create.valid(dtr, as.matrix(tr[vi, ..ft]), label = tr$whiff[vi])
    m <- lgb.train(params = list(objective = "binary", metric = "binary_logloss",
                     learning_rate = .06, num_leaves = 31, min_data_in_leaf = 300,
                     feature_fraction = .8, bagging_fraction = .8, bagging_freq = 1),
                   data = dtr, nrounds = 2000, valids = list(v = dva),
                   early_stopping_rounds = 60, verbose = -1)
    TE[[s]] <- data.table(id = te$id, season = te$season, y = te$whiff,
                          p = predict(m, as.matrix(te[, ..ft])))
  }
  saveRDS(TE, TC)
} else TE <- readRDS(TC)

cat(sprintf("  %-9s %10s %13s %14s\n", "features", "logloss", "season RMSE", "seam-shift gap"))
for (s in names(SETS)) {
  A <- merge(agg(TE[[s]]), KEY, by = c("id","season"))[n >= 40]
  A[, `:=`(r = 100*(act - pred), gg = A5_seam)]
  g <- crob(A, r ~ gg, "ggTRUE")
  cat(sprintf("  %-9s %10.5f %13.3f %8.2f (p=%.3f)\n", s,
              ll(TE[[s]]$y, TE[[s]]$p), sqrt(mean(A[, (100*(act-pred))^2])), g[1], g[3]))
}

cat("\n  --- a post-hoc correction instead of a feature: fit the adjustment on 2020-2023 only ---\n")
A_tr <- merge(agg(OUT$STUFF$w), KEY, by = c("id","season"))[n >= 40 & season <= 2023]
A_tr[, r := 100*(act - pred)]
adj <- sapply(TAGS, function(t) A_tr[get(t) == TRUE, weighted.mean(r, n)])
A_te <- merge(agg(TE$STUFF), KEY, by = c("id","season"))[n >= 40]
A_te[, r := 100*(act - pred)]
base <- sqrt(mean(A_te$r^2))
cat(sprintf("  %-30s %10s %10s %8s\n", "adjustment learned on 2020-23", "value", "test RMSE", "change"))
cat(sprintf("  %-30s %10s %10.3f %8s\n", "no adjustment", "-", base, "-"))
for (t in TAGS) {
  A_te[, r2 := r - fifelse(get(t) == TRUE, adj[[t]], 0)]
  nr <- sqrt(mean(A_te$r2^2))
  cat(sprintf("  %-30s %+10.2f %10.3f %+8.3f\n", LAB[[t]], adj[[t]], nr, nr - base))
}
# all at once, additive, no double counting for overlaps
A_te[, tot := 0]
for (t in TAGS) A_te[get(t) == TRUE, tot := tot + adj[[t]]]
cat(sprintf("  %-30s %10s %10.3f %+8.3f\n", "all archetypes additively", "-",
            sqrt(mean((A_te$r - A_te$tot)^2)), sqrt(mean((A_te$r - A_te$tot)^2)) - base))
