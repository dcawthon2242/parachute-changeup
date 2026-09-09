#!/usr/bin/env Rscript

# WHO THE SEASON-MEAN COMPRESSION ACTUALLY RE-GRADES.
#
# Feeding the model season-mean versions of the arsenal features was the single biggest improvement
# found, 8.1% off pitcher-season RMSE, and it costs nothing to compute. This puts names to it.
#
# Both models are scored with the SAME pitcher-grouped folds, so no arm is graded by a model that
# trained on it. Predictions are then converted to a Stuff+-style index with a SINGLE linear map -
# the baseline model's own mean and standard deviation - so the two indices sit on one scale and a
# difference between them is a real re-grade rather than an artefact of rescaling each separately.
#
# The test of whether a move is any good is not its size. It is whether it points at the residual
# it was supposed to fix, so every mover is shown against what the pitcher actually did.

suppressPackageStartupMessages({ library(data.table); library(lightgbm); library(bit64) })
set.seed(29); options(width = 225); MDIR <- "data/statcast_model"
CACHE <- file.path(MDIR, "smean_oof.rds")

P <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(P)
P <- P[is.finite(release_speed) & is.finite(release_spin_rate) & is.finite(release_extension) &
       is.finite(ax) & is.finite(az) & is.finite(release_pos_x) & is.finite(release_pos_z) &
       is.finite(sax) & is.finite(cax) & is.finite(speed_diff) & is.finite(rv) &
       is.finite(axis_diff) & is.finite(arm_angle) & is.finite(arm_diff)]
P[, lh := p_throws == "L"]
P[, `:=`(tj_ax = fifelse(lh, -ax, ax), tj_x0 = fifelse(lh, -release_pos_x, release_pos_x),
         tj_ax_diff = fifelse(lh, -ax_diff, ax_diff), sep = -speed_diff,
         mx = fifelse(lh, -ax, ax), mz = az + 32.174)]
P[, spin_axis := (atan2(sax, cax) * 180/pi) %% 360]
P[, spin_axis := fifelse(lh, (360 - spin_axis) %% 360, spin_axis)]
wrap <- function(d) ((d + 180) %% 360) - 180
P[, dev_raw := wrap(atan2(mx, mz) * 180/pi + (spin_axis - 180))]
P[, dev := wrap(dev_raw - median(dev_raw, na.rm = TRUE))][, absdev := abs(dev)]
P[, id := as.character(pitcher)]
SM <- P[, .(m_axis = mean(axis_diff), m_dev = mean(dev), m_abs = mean(absdev), m_sep = mean(sep),
            m_arm = mean(arm_angle), m_ad = mean(arm_diff)), by = .(id, season)]
P <- merge(P, SM, by = c("id","season"), sort = FALSE)

SHAPE <- c("release_speed","release_spin_rate","release_extension","tj_ax","az","tj_x0",
           "release_pos_z","speed_diff","tj_ax_diff","az_diff","axis_diff","dev","absdev",
           "arm_angle","arm_diff")
SMEAN <- c(SHAPE, "m_axis","m_dev","m_abs","m_sep","m_arm","m_ad")

grp_oof <- function(D, feats, label, obj, K = 5) {
  arms <- unique(D$id); fa <- data.table(id = arms, fold = sample(rep(1:K, length.out = length(arms))))
  D <- merge(D, fa, by = "id", sort = FALSE); y <- D[[label]]; p <- rep(NA_real_, nrow(D))
  for (f in 1:K) {
    tri <- which(D$fold != f); va <- sample(tri, floor(.10*length(tri))); tr <- setdiff(tri, va)
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
  data.table(id = D$id, season = D$season, y = y, p = p)
}

if (!file.exists(CACHE)) {
  W <- P[is_swing == 1 & is.finite(whiff)]
  cat("fitting SMEAN whiff and run value with pitcher-grouped folds ...\n")
  saveRDS(list(w = grp_oof(W, SMEAN, "whiff", "binary"),
               r = grp_oof(P, SMEAN, "rv", "regression")), CACHE)
}
S <- readRDS(CACHE); B <- readRDS(file.path(MDIR, "absorb_preds.rds"))$STUFF

agg <- function(X) X[, .(n = .N, act = mean(y), pred = mean(p)), by = .(id, season)]
NM <- unique(P[, .(id, season, nm = player_name)])
M <- readRDS(file.path(MDIR, "archetype_roster.rds")); setDT(M)
TAGS <- c("A2_wide","A4_broad","A5_seam","A6_extreme","A7_trad","A8_mismatch")
SHORT <- c(A2_wide = "matched", A4_broad = "broad", A5_seam = "seam", A6_extreme = "extreme",
           A7_trad = "trad", A8_mismatch = "mismatch")

W0 <- agg(B$w)[, .(id, season, nsw = n, act = 100*act, s0 = 100*pred)]
W1 <- agg(S$w)[, .(id, season, s1 = 100*pred)]
R0 <- agg(B$r)[, .(id, season, np = n, rvact = 100*act, r0 = 100*pred)]
R1 <- agg(S$r)[, .(id, season, r1 = 100*pred)]
D <- Reduce(function(a, b) merge(a, b, by = c("id","season")), list(W0, W1, R0, R1))
D <- merge(D, NM, by = c("id","season"))
D <- merge(D, M[, c("id","season", TAGS), with = FALSE], by = c("id","season"), all.x = TRUE)
D <- D[nsw >= 75]

# One linear map, taken from the baseline model, applied to both. 100 = average, 10 = one SD.
mu <- D[, mean(s0)]; sdv <- D[, sd(s0)]
D[, `:=`(x0 = 100 + 10*(s0 - mu)/sdv, x1 = 100 + 10*(s1 - mu)/sdv)]
muR <- D[, mean(r0)]; sdR <- D[, sd(r0)]
D[, `:=`(v0 = 100 + 10*(r0 - muR)/sdR, v1 = 100 + 10*(r1 - muR)/sdR)]
D[, `:=`(mv = x1 - x0, mvR = v1 - v0, e0 = act - s0, e1 = act - s1)]
D[, tag := apply(as.matrix(D[, ..TAGS]), 1, function(r) {
  k <- SHORT[TAGS[which(r == TRUE)]]; if (!length(k)) "" else paste(k, collapse = "/") })]

cat(sprintf("\n%d pitcher-seasons, %d arms, 75-swing floor.\n", nrow(D), uniqueN(D$id)))
cat(sprintf("index scale: 1 SD = 10 points = %.2f whiff points. one whiff point = %.1f index points.\n",
            sdv, 10/sdv))
cat(sprintf("mean absolute move %.1f index points; %.0f%% of seasons move more than 5 points, %.0f%% more than 10.\n",
            D[, mean(abs(mv))], 100*D[, mean(abs(mv) > 5)], 100*D[, mean(abs(mv) > 10)]))

## =============================================================================================
cat("\n=== 1. IS THE RE-GRADE POINTED AT THE RIGHT THING? ===\n\n")
cat(sprintf("  correlation of the move with the baseline residual  %+.3f\n", D[, cor(mv, e0)]))
cat(sprintf("  share of moves that reduce the absolute error       %.0f%%\n", 100*D[, mean(abs(e1) < abs(e0))]))
cat(sprintf("  RMSE  baseline %.2f  ->  compressed %.2f  (%.1f%%)\n",
            D[, sqrt(mean(e0^2))], D[, sqrt(mean(e1^2))],
            100*(D[, sqrt(mean(e1^2))] - D[, sqrt(mean(e0^2))])/D[, sqrt(mean(e0^2))]))
Q <- D[, .(seasons = .N, mean_move = mean(mv), err_before = sqrt(mean(e0^2)), err_after = sqrt(mean(e1^2))),
       by = .(bucket = cut(mv, c(-Inf,-10,-5,5,10,Inf),
                           labels = c("down 10+","down 5-10","little change","up 5-10","up 10+")))][order(bucket)]
Q[, `:=`(improved = round(100*(err_before - err_after)/err_before))]
print(Q[, .(bucket, seasons, mean_move = round(mean_move,1), err_before = round(err_before,2),
            err_after = round(err_after,2), `error cut %` = improved)], row.names = FALSE)

## =============================================================================================
cat("\n=== 2. BIGGEST UPGRADES ===\n\n")
show <- function(X) X[, .(pitcher = nm, season, swings = nsw, stuff_plus = round(x0,1),
                          compressed = round(x1,1), move = round(mv,1),
                          actual_whiff = round(act,1), miss_before = round(e0,1),
                          miss_after = round(e1,1), rv_move = round(mvR,1), archetype = tag)]
print(show(D[order(-mv)][1:20]), row.names = FALSE)
cat("\n=== 3. BIGGEST DOWNGRADES ===\n\n")
print(show(D[order(mv)][1:20]), row.names = FALSE)

## =============================================================================================
cat("\n=== 4. BIGGEST MOVERS AMONG SEASONS WITH REAL SAMPLE (200+ SWINGS) ===\n\n")
BIG <- D[nsw >= 200]
cat(sprintf("  %d seasons clear 200 swings. mean absolute move %.1f index points.\n\n",
            nrow(BIG), BIG[, mean(abs(mv))]))
print(show(rbind(BIG[order(-mv)][1:10], BIG[order(mv)][1:10])), row.names = FALSE)

## =============================================================================================
cat("\n=== 5. CAREER-LEVEL: ARMS THE COMPRESSION SYSTEMATICALLY RE-RATES ===\n\n")
CA <- D[, .(seasons = .N, swings = sum(nsw), stuff_plus = weighted.mean(x0, nsw),
            compressed = weighted.mean(x1, nsw), move = weighted.mean(mv, nsw),
            actual = weighted.mean(act, nsw), miss_before = weighted.mean(e0, nsw),
            miss_after = weighted.mean(e1, nsw), yrs_up = sum(mv > 0)),
        by = .(id, nm)][seasons >= 3]
cat("  ten arms lifted most, of", nrow(CA), "with three or more qualifying seasons:\n")
pr <- function(X) X[, .(pitcher = nm, seasons, swings, stuff_plus = round(stuff_plus,1),
                        compressed = round(compressed,1), move = round(move,1),
                        actual_whiff = round(actual,1), miss_before = round(miss_before,1),
                        miss_after = round(miss_after,1), yrs_up = paste0(yrs_up,"/",seasons))]
print(pr(CA[order(-move)][1:10]), row.names = FALSE)
cat("\n  ten arms cut most:\n")
print(pr(CA[order(move)][1:10]), row.names = FALSE)

## =============================================================================================
cat("\n=== 6. MOVE BY ARCHETYPE ===\n\n")
AR <- rbindlist(lapply(TAGS, function(t) {
  X <- D[get(t) == TRUE]
  if (!nrow(X)) return(NULL)
  data.table(archetype = SHORT[[t]], seasons = nrow(X), move = mean(X$mv),
             miss_before = mean(X$e0), miss_after = mean(X$e1),
             pct_up = 100*mean(X$mv > 0))
}))
AR <- rbind(data.table(archetype = "everyone", seasons = nrow(D), move = mean(D$mv),
                       miss_before = mean(D$e0), miss_after = mean(D$e1),
                       pct_up = 100*mean(D$mv > 0)), AR)
print(AR[, .(archetype, seasons, mean_move = round(move,2), miss_before = round(miss_before,2),
             miss_after = round(miss_after,2), `% moved up` = round(pct_up))], row.names = FALSE)
saveRDS(D, file.path(MDIR, "smean_movers.rds"))
