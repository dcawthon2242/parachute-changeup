#!/usr/bin/env Rscript

# DOES THE SEPARATION BIAS SHOW UP IN RUN VALUE?
#
# The whiff answer was clean: a fastball-blind model underrates big-separation changeups by 0.70
# whiff points per mph, and pricing separation in absorbs 77 percent of it. Run value said yes too,
# +0.071 RV/100 per mph, but pricing separation in absorbed only 43 percent. That asymmetry is the
# interesting part. If the run-value edge were purely the whiff edge cashed in, the same feature
# that explains the whiffs should explain the runs, and it does not.
#
# So this asks four things:
#
#   1  how big is it in runs a team would actually notice
#   2  which channel carries it - takes, whiffs, fouls, or balls in play
#   3  is the in-play piece real contact quality or BABIP luck, using xwOBA-implied RV
#   4  is it a count artifact. Neither model sees the count, and run value depends on the count
#      enormously, so if big-separation arms deploy the changeup in better counts they would look
#      good for a reason that has nothing to do with the pitch
#
# The cached per-pitch predictions carry no pitch key, so the source table is rebuilt under the
# identical filter and the alignment is asserted on four columns before anything is joined.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(37); options(width = 205)
MDIR <- "data/statcast_model"
L <- readRDS(file.path(MDIR, "velo_sep_resid.rds")); R <- as.data.table(L$R)

F <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(F)
F <- F[is.finite(speed_diff) & is.finite(release_speed) & is.finite(release_spin_rate) &
       is.finite(ax) & is.finite(az) & is.finite(release_extension) &
       is.finite(release_pos_x) & is.finite(release_pos_z) & is.finite(rv)]
F[, lh := p_throws == "L"]
F[, `:=`(tj_ax = fifelse(lh, -ax, ax), tj_x0 = fifelse(lh, -release_pos_x, release_pos_x),
         tj_ax_diff = fifelse(lh, -ax_diff, ax_diff), sep = -speed_diff)]
stopifnot(nrow(F) == nrow(R), identical(F$pitcher, R$pitcher), identical(F$season, R$season),
          isTRUE(all.equal(F$rv, R$rv)), isTRUE(all.equal(F$sep, R$sep)))
F[, `:=`(q_blind = R$q_blind, q_aware = R$q_aware)]
cat(sprintf("aligned %s changeups with cached blind/aware run-value grades\n\n",
            format(nrow(F), big.mark = ",")))

F[, ch := factor(fifelse(!is_swing, "take", fifelse(whiff == 1, "whiff",
                 fifelse(is_bip == TRUE, "bip", "foul"))), levels = c("take","whiff","foul","bip"))]
F[, big := sep >= 10]

## ---- 1. how many runs is this worth? -------------------------------------------------------------
S <- F[, .(np = .N, sep = mean(sep), rv100 = 100*mean(rv),
           rv_blind = 100*mean(rv - q_blind), rv_aware = 100*mean(rv - q_aware)),
       by = .(pitcher, player_name, season)]
NS <- R[, .N, by = .(pitcher, season)]  # swing gate needs the swing table, use pitches as proxy
SW <- as.data.table(L$W)[, .(nsw = .N), by = .(pitcher, season)]
S <- merge(S, SW, by = c("pitcher","season"))[nsw >= 40]
S[, id := as.character(pitcher)]
sl <- coef(lm(rv_blind ~ sep, S))["sep"]
cat("=== 1. size of the effect in runs ===\n")
cat(sprintf("  slope %+.4f RV/100 per mph of separation\n", sl))
cat(sprintf("  a pitcher 5 mph above the league mean, throwing the median %d changeups in a season,\n",
            round(median(S$np))))
cat(sprintf("  is misgraded by %.2f runs over the year. At the 95th percentile of separation (%.1f mph,\n",
            5*sl*median(S$np)/100, quantile(S$sep,.95)))
cat(sprintf("  %+.1f above the mean) and a heavy changeup workload of %d, it is %.2f runs.\n",
            quantile(S$sep,.95) - mean(S$sep), round(quantile(S$np,.95)),
            (quantile(S$sep,.95) - mean(S$sep))*sl*quantile(S$np,.95)/100))
cat("  For scale, the same arms are underrated by roughly 2.7 whiff points, which is a real\n")
cat("  pitch-quality difference but a fraction of a win over a season.\n")

## ---- 2. which channel carries it? ----------------------------------------------------------------
# The residual is decomposed the same way the total was in mech_10: each channel's contribution is
# its share of pitches times its mean residual, so the four contributions sum to the total.
cat("\n=== 2. which channel carries the underrating (10+ mph vs the rest) ===\n")
dec <- function(D) {
  n <- nrow(D)
  D[, .(share = .N/n, resid100 = 100*mean(rv - q_blind),
        contrib = (.N/n) * 100*mean(rv - q_blind)), by = ch][order(ch)]
}
A <- dec(F[big == TRUE]); B <- dec(F[big == FALSE])
C <- merge(A, B, by = "ch", suffixes = c("_big","_rest"))
C[, `:=`(d_share = 100*(share_big - share_rest), d_contrib = contrib_big - contrib_rest)]
print(C[, .(ch, share_big = round(100*share_big,1), share_rest = round(100*share_rest,1),
            d_share_pts = round(d_share,2), resid_big = round(resid100_big,2),
            resid_rest = round(resid100_rest,2), d_contrib = round(d_contrib,3))],
      row.names = FALSE)
cat(sprintf("\n  channel contributions sum to %+.3f RV/100, the total blind-model gap\n",
            sum(C$d_contrib)))

## ---- 2b. the take channel is the surprise, so test it rather than eyeball it ---------------------
# A take is a called strike or a ball, and which one it is depends on the zone. z_rel_bot is the
# height above the bottom of the zone and z_rel_top the height above the top, so a pitch is in the
# zone when the first is positive, the second negative, and it is within a plate-half of centre.
F[, in_zone := is.finite(z_rel_bot) & is.finite(z_rel_top) & is.finite(plate_x) &
               z_rel_bot > 0 & z_rel_top < 0 & abs(plate_x) <= 0.83]
cat("\n  is the take advantage a called-strike advantage?\n")
print(F[, .(zone_rate = round(100*mean(in_zone),1),
            take_rate = round(100*mean(ch == "take"),1),
            called_strike_rate = round(100*mean(ch == "take" & in_zone),1),
            of_takes_in_zone = round(100*mean(in_zone[ch == "take"]),1)),
        by = .(group = fifelse(big, "10+ mph", "rest"))], row.names = FALSE)

TK <- F[ch == "take", .(tk = 100*mean(rv - q_blind), n = .N, sep = mean(sep)),
        by = .(pitcher, season)][n >= 40]
TK[, id := as.character(pitcher)]
mt <- lm(tk ~ sep, TK); ut <- residuals(mt); Xt <- model.matrix(mt)
bt <- solve(crossprod(Xt))
Vt <- bt %*% crossprod(rowsum(Xt * ut, TK$id)) %*% bt * (uniqueN(TK$id)/(uniqueN(TK$id)-1))
cat(sprintf("\n  take-channel residual vs separation: %+.4f per mph (se %.4f) | p %.4f | %d arms\n",
            coef(mt)["sep"], sqrt(diag(Vt))["sep"],
            2*pt(-abs(coef(mt)["sep"]/sqrt(diag(Vt))["sep"]), uniqueN(TK$id)-1), uniqueN(TK$id)))

## ---- 3. is the contact real or lucky? --------------------------------------------------------------
cat("\n=== 3. balls in play: actual vs xwOBA-implied ===\n")
BIP <- F[ch == "bip" & is.finite(estimated_woba_using_speedangle)]
fit <- lm(rv ~ estimated_woba_using_speedangle, BIP)
BIP[, xrv := as.numeric(predict(fit, .SD))]
cat(sprintf("  xwOBA -> RV map: slope %.3f, R2 %.3f on %s changeups in play\n",
            coef(fit)[2], summary(fit)$r.squared, format(nrow(BIP), big.mark = ",")))
for (g in c(TRUE, FALSE)) {
  X <- BIP[big == g]
  cat(sprintf("  %-8s n=%6s  xwOBA %.3f  actual RV/100 %+6.2f  expected RV/100 %+6.2f  gap %+.2f\n",
              if (g) "10+ mph" else "rest", format(nrow(X), big.mark = ","),
              mean(X$estimated_woba_using_speedangle), 100*mean(X$rv), 100*mean(X$xrv),
              100*(mean(X$rv) - mean(X$xrv))))
}
cat(sprintf("  xwOBA difference p = %.4f | actual BIP RV difference p = %.4f\n",
            t.test(BIP[big == TRUE]$estimated_woba_using_speedangle,
                   BIP[big == FALSE]$estimated_woba_using_speedangle)$p.value,
            t.test(BIP[big == TRUE]$rv, BIP[big == FALSE]$rv)$p.value))

## ---- 4. is it a count artifact? ---------------------------------------------------------------------
cat("\n=== 4. count deployment ===\n")
cat("  how often the changeup is thrown with two strikes, and in a hitter's count:\n")
print(F[, .(two_strike = round(100*mean(strikes == 2),1),
            hitter_count = round(100*mean(balls > strikes),1),
            rv100 = round(100*mean(rv),2)),
        by = .(group = fifelse(big, "10+ mph", "rest"))], row.names = FALSE)

CFEAT <- c("release_speed","release_spin_rate","release_extension","tj_ax","az","tj_x0",
           "release_pos_z","balls","strikes")
CACHE2 <- file.path(MDIR, "velo_sep_rv_count.rds")
if (!file.exists(CACHE2)) {
  K <- 4; fold <- sample(rep(1:K, length.out = nrow(F))); p <- rep(NA_real_, nrow(F))
  for (f in 1:K) {
    tri <- which(fold != f); n <- length(tri); vi <- sample(n, floor(.12*n))
    dtr <- lgb.Dataset(as.matrix(F[tri[-vi], ..CFEAT]), label = F$rv[tri[-vi]])
    dva <- lgb.Dataset.create.valid(dtr, as.matrix(F[tri[vi], ..CFEAT]), label = F$rv[tri[vi]])
    m <- lgb.train(params = list(objective = "regression", metric = "l2", learning_rate = .06,
                   num_leaves = 31, min_data_in_leaf = 300, feature_fraction = .8,
                   bagging_fraction = .8, bagging_freq = 1), data = dtr, nrounds = 1500,
                   valids = list(v = dva), early_stopping_rounds = 50, verbose = -1)
    p[fold == f] <- predict(m, as.matrix(F[fold == f, ..CFEAT]))
  }
  saveRDS(p, CACHE2)
} else cat("  (using cached count-aware fit)\n")
F[, q_count := readRDS(CACHE2)]
cat(sprintf("  blind R2 %.5f -> blind+count R2 %.5f\n",
            1 - var(F$rv - F$q_blind)/var(F$rv), 1 - var(F$rv - F$q_count)/var(F$rv)))
S2 <- F[, .(rv_count = 100*mean(rv - q_count)), by = .(pitcher, season)]
S <- merge(S, S2, by = c("pitcher","season"))

crob <- function(f, D) {
  D <- D[complete.cases(D[, c(all.vars(f), "id"), with = FALSE])]
  m <- lm(f, D); u <- residuals(m); X <- model.matrix(m); cl <- D$id
  nc <- uniqueN(cl); b <- solve(crossprod(X))
  V <- b %*% crossprod(rowsum(X * u, cl)) %*% b * (nc/(nc-1))
  e <- coef(m)["sep"]; s <- sqrt(diag(V))["sep"]
  sprintf("%+.4f RV/100 per mph (se %.4f) | p %.4f", e, s, 2*pt(-abs(e/s), nc-1))
}
## ---- 4b. pitch quality or usage? the decisive test -----------------------------------------------
# The take channel turned out to be a zone-rate story, and where a pitcher chooses to throw the ball
# is not something a stuff model is supposed to grade. So this gives the model plate location. If
# the separation slope survives, big-separation changeups really do beat their grade as pitches. If
# it collapses, the run-value edge was a usage advantage that separation buys rather than a
# mispriced pitch, and a stuff model is not wrong to miss it.
LFEAT <- c("release_speed","release_spin_rate","release_extension","tj_ax","az","tj_x0",
           "release_pos_z","tj_px","plate_z")
F[, tj_px := fifelse(lh, -plate_x, plate_x)]
CACHE3 <- file.path(MDIR, "velo_sep_rv_loc.rds")
FL <- F[is.finite(tj_px) & is.finite(plate_z)]
if (!file.exists(CACHE3)) {
  K <- 4; fold <- sample(rep(1:K, length.out = nrow(FL))); p <- rep(NA_real_, nrow(FL))
  for (f in 1:K) {
    tri <- which(fold != f); n <- length(tri); vi <- sample(n, floor(.12*n))
    dtr <- lgb.Dataset(as.matrix(FL[tri[-vi], ..LFEAT]), label = FL$rv[tri[-vi]])
    dva <- lgb.Dataset.create.valid(dtr, as.matrix(FL[tri[vi], ..LFEAT]), label = FL$rv[tri[vi]])
    m <- lgb.train(params = list(objective = "regression", metric = "l2", learning_rate = .06,
                   num_leaves = 31, min_data_in_leaf = 300, feature_fraction = .8,
                   bagging_fraction = .8, bagging_freq = 1), data = dtr, nrounds = 1500,
                   valids = list(v = dva), early_stopping_rounds = 50, verbose = -1)
    p[fold == f] <- predict(m, as.matrix(FL[fold == f, ..LFEAT]))
  }
  saveRDS(p, CACHE3)
} else cat("\n  (using cached location-aware fit)\n")
FL[, q_loc := readRDS(CACHE3)]
cat(sprintf("  blind R2 %.5f -> blind+location R2 %.5f\n",
            1 - var(FL$rv - FL$q_blind)/var(FL$rv), 1 - var(FL$rv - FL$q_loc)/var(FL$rv)))
S <- merge(S, FL[, .(rv_loc = 100*mean(rv - q_loc)), by = .(pitcher, season)],
           by = c("pitcher","season"))

cat("\n  separation slope against each run-value residual, clustered on the arm:\n")
for (v in c("rv_blind","rv_count","rv_aware","rv_loc"))
  cat(sprintf("    %-9s %s\n", v, crob(as.formula(paste(v, "~ sep")), S)))
cat("\n  and the whiff comparison, for reference: blind +0.70/mph (p<1e-4), aware +0.16 (p=.16)\n")

## ---- 5. persistence ---------------------------------------------------------------------------------
cat("\n=== 5. does the run-value residual persist? run value is far noisier than whiff ===\n")
Y <- merge(S[, .(id, season, sep, rv_blind, rv_aware)],
           S[, .(id, season = season + 1L, sep_p = sep, rvb_p = rv_blind, rva_p = rv_aware)],
           by = c("id","season"))
for (v in list(c("rvb_p","rv_blind","blind RV residual, year to year"),
               c("rva_p","rv_aware","aware RV residual, year to year"))) {
  z <- cor.test(Y[[v[1]]], Y[[v[2]]])
  cat(sprintf("  %-36s r = %+.3f (p %.4f, n = %d)\n", v[3], z$estimate, z$p.value, nrow(Y)))
}
z <- cor.test(Y$sep_p, Y$rv_blind)
cat(sprintf("  prior-season separation -> current RV residual   r = %+.3f (p %.4f)\n",
            z$estimate, z$p.value))

saveRDS(list(seasons = S, channels = C), file.path(MDIR, "velo_sep_rv.rds"))
cat("\nwrote velo_sep_rv.rds\n")
