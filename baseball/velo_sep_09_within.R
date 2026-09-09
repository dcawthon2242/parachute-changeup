#!/usr/bin/env Rscript

# WITHIN-OUTING REPEATABILITY, AND THE WHOLE THING AGAIN ON CHASE.
#
# TWO JOBS.
#
# 1. SPLIT THE CONSISTENCY MEASURE. sd_game removed day-to-day fastball drift but still pools every
#    changeup in the season, so a pitcher whose changeup alone fades in September still scores as
#    inconsistent. The variance of the gap decomposes cleanly:
#
#      within   pooled standard deviation of the gap INSIDE an outing - can he repeat the pitch
#      between  standard deviation of his game means - does the pitch drift across the season
#
#    Note that inside a single outing the game anchor is a constant, so the within term is exactly
#    within-game changeup velocity spread. That is the point: it is repeatability with the day's
#    conditions differenced out. Games need at least MIN_CH changeups to contribute a variance.
#
# 2. CHASE. Everything so far has been whiff and run value. Chase asks a different question - not
#    whether the hitter was beaten once committed, but whether he was fooled into offering at a
#    ball. A velocity gap is a plausible mechanism for that and there is no reason it has to give
#    the same answer, so it is a genuine second test rather than a re-slice.
#
#    Out of zone is the Statcast strike zone with a baseball's width of margin: more than half a
#    plate off centre horizontally, or above sz_top / below sz_bot.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(59); options(width = 205)
MDIR <- "data/statcast_model"; TBL <- file.path(MDIR, "anchor_model_table2.rds")
MIN_FF <- 5L; MIN_CH <- 5L

if (!file.exists(TBL)) {
  COLS <- c("pitch_type","release_speed","pitcher","player_name","p_throws","release_spin_rate",
            "release_extension","release_pos_x","release_pos_z","ax","az","plate_x","plate_z",
            "sz_top","sz_bot","description","delta_pitcher_run_exp","game_pk","at_bat_number",
            "pitch_number","balls","strikes")
  A <- rbindlist(lapply(2020:2026, function(y) {
    x <- fread(sprintf("data/statcast_%d/statcast_%d_all.csv", y, y), select = COLS,
               showProgress = FALSE)
    x <- x[pitch_type %chin% c("FF","CH") & is.finite(release_speed) & is.finite(ax) & is.finite(az)]
    x[, season := y]; cat(sprintf("  %d: %s rows\n", y, format(nrow(x), big.mark = ","))); x }))

  FS <- A[pitch_type == "FF"]
  GB <- FS[, .(fb_v_game = mean(release_speed), n_ff = .N),
           by = .(pitcher, season, game_pk)][n_ff >= MIN_FF]
  SB <- FS[, .(fb_v_season = mean(release_speed), fb_ax = mean(ax), fb_az = mean(az)),
           by = .(pitcher, season)]

  C <- A[pitch_type == "CH"]
  C <- merge(C, GB, by = c("pitcher","season","game_pk"))
  C <- merge(C, SB, by = c("pitcher","season"))
  C[, `:=`(lh = p_throws == "L",
           sep_season = fb_v_season - release_speed, sep_game = fb_v_game - release_speed)]
  C[, `:=`(tj_ax = fifelse(lh, -ax, ax), tj_x0 = fifelse(lh, -release_pos_x, release_pos_x),
           tj_px = fifelse(lh, -plate_x, plate_x),
           tj_ax_diff = fifelse(lh, -(ax - fb_ax), ax - fb_ax), az_diff = az - fb_az,
           is_swing = as.integer(description %chin%
             c("swinging_strike","swinging_strike_blocked","foul_tip","missed_bunt","foul",
               "foul_bunt","hit_into_play","hit_into_play_score","hit_into_play_no_out")),
           whiff = as.integer(description %chin%
             c("swinging_strike","swinging_strike_blocked","foul_tip","missed_bunt")),
           rv = delta_pitcher_run_exp)]
  C[, ooz := as.integer(abs(plate_x) > 0.83 | plate_z > sz_top | plate_z < sz_bot)]

  # Variance decomposition of the gap. Only outings with enough changeups can carry a within term.
  GV <- C[, .(ng = .N, mg = mean(sep_game), vg = var(sep_game)), by = .(pitcher, season, game_pk)]
  DEC <- GV[ng >= MIN_CH, .(games = .N, np_dec = sum(ng),
                            sd_within = sqrt(sum((ng-1)*vg)/sum(ng-1)),
                            sd_between = sd(mg)), by = .(pitcher, season)]
  C <- merge(C, DEC, by = c("pitcher","season"), all.x = TRUE)
  C[, `:=`(sd_season = sd(sep_season), sd_game = sd(sep_game), n_ps = .N), by = .(pitcher, season)]
  C <- C[n_ps >= 40 & is.finite(release_spin_rate) & is.finite(release_extension) &
         is.finite(tj_ax_diff) & is.finite(rv) & is.finite(plate_z)]
  saveRDS(C[, .(pitcher, player_name, season, game_pk, release_speed, release_spin_rate,
                release_extension, tj_ax, az, tj_x0, release_pos_z, tj_ax_diff, az_diff,
                tj_px, plate_z, sep_season, sep_game, sd_season, sd_game, sd_within, sd_between,
                games, whiff, is_swing, ooz, rv, balls, strikes)], TBL)
} else cat("(using cached table)\n")

C <- readRDS(TBL); setDT(C)
cat(sprintf("\n%s changeups | %d pitcher-seasons | %d arms | out of zone %.1f%%\n",
            format(nrow(C), big.mark = ","), uniqueN(C[, .(pitcher, season)]),
            uniqueN(C$pitcher), 100*mean(C$ooz)))

## ---- 1. the variance decomposition -----------------------------------------------------------------
D <- unique(C[is.finite(sd_within) & is.finite(sd_between),
              .(pitcher, season, sd_season, sd_game, sd_within, sd_between, games)])
cat(sprintf("\n=== 1. where does the variability live? %d pitcher-seasons with 5+ usable outings ===\n",
            nrow(D)))
cat(sprintf("  median WITHIN-outing sd   %.3f mph   (can he repeat it)\n", median(D$sd_within)))
cat(sprintf("  median BETWEEN-outing sd  %.3f mph   (does it drift across the year)\n", median(D$sd_between)))
cat(sprintf("  the within term is %.0f%% of total variance - repeatability dominates drift\n",
            100*median(D$sd_within^2/(D$sd_within^2 + D$sd_between^2))))
cat(sprintf("  cor(within, between) = %+.3f | cor(within, season sd) = %+.3f\n",
            cor(D$sd_within, D$sd_between), cor(D$sd_within, D$sd_season)))

## ---- fit helper --------------------------------------------------------------------------------------
BASE <- c("release_speed","release_spin_rate","release_extension","tj_ax","az","tj_x0",
          "release_pos_z","tj_ax_diff","az_diff")
run <- function(D, feats, label, seed = 59) {
  set.seed(seed)
  arms <- unique(D$pitcher); fa <- sample(rep(1:4, length.out = length(arms)))
  fold <- fa[match(D$pitcher, arms)]; y <- D[[label]]; p <- rep(NA_real_, nrow(D))
  for (f in 1:4) {
    tri <- which(fold != f); n <- length(tri); vi <- sample(n, floor(.12*n))
    dtr <- lgb.Dataset(as.matrix(D[tri[-vi], ..feats]), label = y[tri[-vi]])
    dva <- lgb.Dataset.create.valid(dtr, as.matrix(D[tri[vi], ..feats]), label = y[tri[vi]])
    m <- lgb.train(params = list(objective = "binary", metric = "binary_logloss",
                   learning_rate = .06, num_leaves = 31, min_data_in_leaf = 300,
                   feature_fraction = .8, bagging_fraction = .8, bagging_freq = 1),
                   data = dtr, nrounds = 1500, valids = list(v = dva),
                   early_stopping_rounds = 50, verbose = -1)
    p[fold == f] <- predict(m, as.matrix(D[fold == f, ..feats]))
  }
  list(p = p, rmse = sqrt(mean((y-p)^2)), r2 = 1 - var(y-p)/var(y))
}
crob <- function(f, D, k) {
  D <- D[complete.cases(D[, c(all.vars(f), "id"), with = FALSE])]
  m <- lm(f, D); u <- residuals(m); X <- model.matrix(m); nc <- uniqueN(D$id)
  b <- solve(crossprod(X)); V <- b %*% crossprod(rowsum(X*u, D$id)) %*% b * (nc/(nc-1))
  e <- coef(m)[k]; s <- sqrt(diag(V))[k]
  sprintf("%+.3f (se %.3f) | p %.4f", e, s, 2*pt(-abs(e/s), nc-1))
}

## ---- 2. CHASE: blind vs separation-aware --------------------------------------------------------------
Z <- C[ooz == 1]
cat(sprintf("\n=== 2. CHASE: %s out-of-zone changeups, chase rate %.3f ===\n",
            format(nrow(Z), big.mark = ","), mean(Z$is_swing)))
BL <- setdiff(BASE, c("tj_ax_diff","az_diff"))
rb <- run(Z, BL, "is_swing"); ra <- run(Z, c(BASE, "sep_season"), "is_swing")
cat(sprintf("  blind (7 feats)  RMSE %.6f | R2 %+.5f\n", rb$rmse, rb$r2))
cat(sprintf("  aware (10 feats) RMSE %.6f | R2 %+.5f\n", ra$rmse, ra$r2))
Z[, `:=`(ch_blind = 100*(is_swing - rb$p), ch_aware = 100*(is_swing - ra$p))]

## ---- 3. all three outcomes, one table ------------------------------------------------------------------
W <- C[is_swing == 1]
rwb <- run(W, BL, "whiff"); rwa <- run(W, c(BASE, "sep_season"), "whiff")
W[, `:=`(wh_blind = 100*(whiff - rwb$p), wh_aware = 100*(whiff - rwa$p))]

S <- merge(Z[, .(nz = .N, ch_blind = mean(ch_blind), ch_aware = mean(ch_aware)),
             by = .(pitcher, season)],
           W[, .(nsw = .N, wh_blind = mean(wh_blind), wh_aware = mean(wh_aware)),
             by = .(pitcher, season)], by = c("pitcher","season"))
S <- merge(S, unique(C[, .(pitcher, season, sep = mean(sep_game), sd_season, sd_game,
                           sd_within, sd_between), by = .(pitcher, season)][
             , .(pitcher, season, sep, sd_season, sd_game, sd_within, sd_between)]),
           by = c("pitcher","season"))
S <- S[nsw >= 40 & nz >= 40]; S[, id := as.character(pitcher)]
cat(sprintf("\n%d pitcher-seasons with 40+ swings and 40+ out-of-zone changeups, %d arms\n",
            nrow(S), uniqueN(S$id)))

cat("\n=== does SEPARATION predict what each model misses? ===\n")
cat("   positive = the model underrates big-separation changeups\n")
for (v in c("wh_blind","wh_aware","ch_blind","ch_aware"))
  cat(sprintf("  %-9s %s per mph of separation\n", v, crob(as.formula(paste(v,"~ sep")), S, "sep")))

cat("\n=== does CONSISTENCY predict the residual, and which component of it? ===\n")
for (v in c("wh_aware","ch_aware")) {
  cat(sprintf("  -- %s --\n", v))
  for (k in c("sd_season","sd_game","sd_within","sd_between"))
    cat(sprintf("     %-11s %s\n", k, crob(as.formula(paste(v,"~",k,"+ sep")), S, k)))
  cat(sprintf("     %-11s %s   (both in together)\n", "within|betw",
              crob(as.formula(paste(v,"~ sd_within + sd_between + sep")), S, "sd_within")))
  cat(sprintf("     %-11s %s\n", "  between",
              crob(as.formula(paste(v,"~ sd_within + sd_between + sep")), S, "sd_between")))
}
saveRDS(S, file.path(MDIR, "velo_sep_within_seasons.rds"))
cat("\nwrote velo_sep_within_seasons.rds\n")
