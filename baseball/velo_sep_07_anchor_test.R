#!/usr/bin/env Rscript

# DOES THE OUTING ANCHOR BEAT THE SEASON CONSTANT?
#
# Built self-contained from the raw season CSVs, because parachute_ff carries no game key and so
# cannot be joined to an outing-level anchor at all.
#
# Every model below shares the same seven changeup-only features and the same season-anchored
# movement differences. The ONLY thing that varies is where the velocity separation is measured
# from, so any RMSE difference is attributable to the anchor and nothing else:
#
#   season   ch velo minus the pitcher's season mean four-seam         - what the project used
#   game     ch velo minus his four-seam mean in that outing
#   recent   ch velo minus his last ten four-seamers earlier that outing - causal within the start
#   game+sd  the game anchor plus the season sd of the game-anchored gap
#
# Folds are grouped by PITCHER throughout. The season-anchored result is the one to beat, and it is
# refit here rather than carried over so that both sides see identical rows - the anchor is only
# available for 80 percent of changeups and comparing against a model fit on the other 20 percent
# too would be a sample difference masquerading as a modelling result.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(53); options(width = 205)
MDIR <- "data/statcast_model"; TBL <- file.path(MDIR, "anchor_model_table.rds")
MIN_FF <- 5L

if (!file.exists(TBL)) {
  COLS <- c("pitch_type","release_speed","pitcher","player_name","p_throws","release_spin_rate",
            "release_extension","release_pos_x","release_pos_z","ax","az","plate_x","plate_z",
            "description","delta_pitcher_run_exp","game_pk","at_bat_number","pitch_number",
            "balls","strikes")
  A <- rbindlist(lapply(2020:2026, function(y) {
    x <- fread(sprintf("data/statcast_%d/statcast_%d_all.csv", y, y), select = COLS,
               showProgress = FALSE)
    x <- x[pitch_type %chin% c("FF","CH") & is.finite(release_speed) & is.finite(ax) & is.finite(az)]
    x[, season := y]; cat(sprintf("  %d: %s rows\n", y, format(nrow(x), big.mark = ","))); x }))

  # Four-seam reference at both resolutions.
  FS <- A[pitch_type == "FF"]
  GB <- FS[, .(fb_v_game = mean(release_speed), n_ff = .N), by = .(pitcher, season, game_pk)][n_ff >= MIN_FF]
  SB <- FS[, .(fb_v_season = mean(release_speed), fb_ax = mean(ax), fb_az = mean(az)),
           by = .(pitcher, season)]

  setorder(A, pitcher, game_pk, at_bat_number, pitch_number)
  A[, seq := seq_len(.N), by = .(pitcher, game_pk)]
  FF <- A[pitch_type == "FF", .(pitcher, game_pk, seq, release_speed)]
  setorder(FF, pitcher, game_pk, seq)
  FF[, ci := seq_len(.N), by = .(pitcher, game_pk)]
  FF[, cs := cumsum(release_speed), by = .(pitcher, game_pk)]
  FF[, fb_v_recent := fifelse(ci >= 10L, (cs - shift(cs, 10L, fill = 0))/10, cs/ci),
     by = .(pitcher, game_pk)]

  C <- A[pitch_type == "CH"]
  C[, rid := .I]
  setkey(FF, pitcher, game_pk, seq)
  C <- merge(C, FF[C, on = .(pitcher, game_pk, seq), roll = TRUE, .(rid = i.rid, fb_v_recent)],
             by = "rid", all.x = TRUE)
  C <- merge(C, GB, by = c("pitcher","season","game_pk"), all.x = TRUE)
  C <- merge(C, SB, by = c("pitcher","season"), all.x = TRUE)

  C[, `:=`(whiff = as.integer(description %chin% c("swinging_strike","swinging_strike_blocked",
                                                   "foul_tip","missed_bunt")),
           is_swing = as.integer(description %chin%
             c("swinging_strike","swinging_strike_blocked","foul_tip","missed_bunt","foul",
               "foul_bunt","hit_into_play","hit_into_play_score","hit_into_play_no_out")),
           rv = delta_pitcher_run_exp, lh = p_throws == "L")]
  C[, `:=`(sep_season = fb_v_season - release_speed,
           sep_game   = fb_v_game   - release_speed,
           sep_recent = fb_v_recent - release_speed,
           tj_ax = fifelse(lh, -ax, ax), tj_x0 = fifelse(lh, -release_pos_x, release_pos_x),
           tj_ax_diff = fifelse(lh, -(ax - fb_ax), ax - fb_ax), az_diff = az - fb_az)]
  C <- C[is.finite(sep_season) & is.finite(sep_game) & is.finite(sep_recent) & is.finite(rv) &
         is.finite(release_spin_rate) & is.finite(release_extension) & is.finite(tj_ax_diff)]
  C[, `:=`(sd_season = sd(sep_season), sd_game = sd(sep_game), n_ps = .N), by = .(pitcher, season)]
  C <- C[n_ps >= 40]
  saveRDS(C[, .(pitcher, player_name, season, game_pk, release_speed, release_spin_rate,
                release_extension, tj_ax, az, tj_x0, release_pos_z, tj_ax_diff, az_diff,
                sep_season, sep_game, sep_recent, sd_season, sd_game, whiff, is_swing, rv,
                balls, strikes)], TBL)
} else cat("(using cached anchor model table)\n")

C <- readRDS(TBL); setDT(C)
cat(sprintf("\n%s changeups, %d pitcher-seasons, %d arms, all with an outing anchor\n",
            format(nrow(C), big.mark = ","), uniqueN(C[, .(pitcher, season)]), uniqueN(C$pitcher)))
cat(sprintf("median season sd of the gap: season anchor %.3f, outing anchor %.3f mph\n",
            median(C$sd_season), median(C$sd_game)))
cat(sprintf("cor(sep_season, sep_game) at the pitch level = %.4f\n\n",
            cor(C$sep_season, C$sep_game)))

BASE <- c("release_speed","release_spin_rate","release_extension","tj_ax","az","tj_x0",
          "release_pos_z","tj_ax_diff","az_diff")
SETS <- list(season = c(BASE,"sep_season"), game = c(BASE,"sep_game"),
             recent = c(BASE,"sep_recent"),
             `game + sd` = c(BASE,"sep_game","sd_game"),
             `season + sd` = c(BASE,"sep_season","sd_season"))

run <- function(D, feats, label, obj) {
  K <- 4; arms <- unique(D$pitcher); fa <- sample(rep(1:K, length.out = length(arms)))
  fold <- fa[match(D$pitcher, arms)]; y <- D[[label]]; p <- rep(NA_real_, nrow(D))
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
  list(rmse = sqrt(mean((y - p)^2)), r2 = 1 - var(y - p)/var(y),
       srmse = sqrt(weighted.mean((S$ay - S$ap)^2, S$n)), p = p)
}

W <- C[is_swing == 1]
cat(sprintf("whiff on %s swings | run value on %s pitches | grouped folds by pitcher\n",
            format(nrow(W), big.mark = ","), format(nrow(C), big.mark = ",")))

for (job in list(list("WHIFF", W, "whiff", "binary"), list("RUN VALUE", C, "rv", "regression"))) {
  cat(sprintf("\n=== %s ===\n", job[[1]]))
  base <- NULL
  for (nm in names(SETS)) {
    set.seed(53); r <- run(job[[2]], SETS[[nm]], job[[3]], job[[4]])
    if (is.null(base)) base <- r
    cat(sprintf("  %-12s RMSE %.6f | R2 %+.5f | season RMSE %.5f | vs season anchor: RMSE %+.6f, season RMSE %+.3f%%\n",
                nm, r$rmse, r$r2, r$srmse, r$rmse - base$rmse,
                100*(r$srmse - base$srmse)/base$srmse))
  }
}

## ---- the residual screen, both anchors side by side -----------------------------------------------
cat("\n=== consistency against the residual: does the outing anchor sharpen the signal? ===\n")
set.seed(53); rw <- run(W, SETS$game, "whiff", "binary")
W[, res := 100*(whiff - rw$p)]
S <- W[, .(nsw = .N, res = mean(res), sd_season = sd_season[1], sd_game = sd_game[1],
           sep = mean(sep_game)), by = .(pitcher, season)][nsw >= 40]
S[, id := as.character(pitcher)]
for (v in c("sd_season","sd_game")) {
  f <- as.formula(paste("res ~", v, "+ sep"))
  m <- lm(f, S); u <- residuals(m); X <- model.matrix(m); nc <- uniqueN(S$id)
  b <- solve(crossprod(X)); V <- b %*% crossprod(rowsum(X * u, S$id)) %*% b * (nc/(nc-1))
  e <- coef(m)[v]; s <- sqrt(diag(V))[v]
  cat(sprintf("  %-10s %+.3f whiff pts per mph (se %.3f) | p %.4f | %d arms\n",
              v, e, s, 2*pt(-abs(e/s), nc-1), nc))
}
cat(sprintf("  the two sd measures correlate at %+.3f, so this is close to the same question twice\n",
            cor(S$sd_season, S$sd_game)))
