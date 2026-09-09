#!/usr/bin/env Rscript

# The only surviving result was ground-ball rate: +5.8 pp, p = .005. Before it counts as
# evidence for the parachute IDEA it has to survive the obvious objection, which is that the
# bin is defined half by IVB kill and a changeup that drops more gets grounders for reasons
# that have nothing to do with spin matching the fastball. So: hold drop constant and ask
# whether the axis half of the definition contributes anything on its own.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(1); options(width = 200)
MDIR <- "data/statcast_model"
d <- readRDS(file.path(MDIR, "parachute_rv.rds"))

sw <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip","hit_into_play",
        "foul_bunt","missed_bunt","bunt_foul_tip")
CH <- d[pitch_type == "CH"]
CH[, `:=`(is_swing = description %in% sw, is_bip = !is.na(launch_speed) & bb_type != "")]
B <- CH[is_bip == TRUE & is.finite(axis_diff) & is.finite(az_diff)]
B[, gb := as.integer(bb_type == "ground_ball")]

# Per-pitch ground-ball model on everything EXCEPT axis gap, so the residual is the part of
# grounder tendency that drop, velocity separation, location and handedness cannot explain.
FEAT <- c("release_speed","release_spin_rate","release_extension","release_pos_x",
          "release_pos_z","ax","az","speed_diff","ax_diff","az_diff","plate_x","plate_z",
          "plate_x_in","plate_x_arm","z_rel_bot","z_rel_top","VAA","HAA","HAA_in",
          "stand_R","throws_R","same_hand","balls","strikes")
B <- B[stats::complete.cases(B[, ..FEAT])]
K <- 4; B[, fold := sample(rep(1:K, length.out = .N))]; p <- rep(NA_real_, nrow(B))
for (f in 1:K) {
  tr <- B[fold != f]; n <- nrow(tr); vi <- sample(n, floor(.12*n))
  dtr <- lgb.Dataset(as.matrix(tr[-vi, ..FEAT]), label = tr$gb[-vi])
  dva <- lgb.Dataset.create.valid(dtr, as.matrix(tr[vi, ..FEAT]), label = tr$gb[vi])
  m <- lgb.train(params = list(objective="binary", metric="binary_logloss", learning_rate=.06,
                 num_leaves=31, min_data_in_leaf=300, feature_fraction=.8, bagging_fraction=.8,
                 bagging_freq=1), data = dtr, nrounds = 1500, valids = list(val=dva),
                 early_stopping_rounds = 50, verbose = -1)
  p[B$fold == f] <- predict(m, as.matrix(B[fold == f, ..FEAT]))
}
B[, gb_res := gb - p]
cat(sprintf("GB model: n=%d  base GB%%=%.1f  AUC-ish R2=%.4f\n\n",
            nrow(B), 100*mean(B$gb), 1 - var(B$gb - p)/var(B$gb)))

S <- B[, .(nbip = .N, axis = mean(axis_diff), ivb_kill = -mean(az_diff),
           gbpct = 100*mean(gb), gb_res = 100*mean(gb_res),
           ev = mean(launch_speed, na.rm = TRUE)), by = .(pitcher, player_name, season)]
S <- S[nbip >= 25]
S[, `:=`(hi_kill = ivb_kill >= 13.9, lo_axis = axis <= 15)]

cat("=== 2x2: does the spin-match half do anything once drop is held fixed? ===\n")
q <- S[, .(n = .N, gb = round(mean(gbpct),1), gb_over_model = round(mean(gb_res),2),
           ev = round(mean(ev),1)), by = .(hi_kill, lo_axis)][order(-hi_kill, -lo_axis)]
print(q, row.names = FALSE)

hk <- S[hi_kill == TRUE]
cat(sprintf("\nWithin HIGH-KILL changeups only (n=%d), low-axis vs high-axis:\n", nrow(hk)))
for (v in c("gbpct","gb_res","ev")) {
  t <- t.test(hk[lo_axis == TRUE][[v]], hk[lo_axis == FALSE][[v]])
  cat(sprintf("  %-8s  low-axis %7.2f   high-axis %7.2f   diff %+.2f   p=%.3f\n",
              v, mean(hk[lo_axis == TRUE][[v]]), mean(hk[lo_axis == FALSE][[v]]),
              diff(rev(t$estimate)), t$p.value))
}
a <- suppressWarnings(cor.test(S$axis, S$gb_res, method = "spearman"))
b <- suppressWarnings(cor.test(S$ivb_kill, S$gbpct, method = "spearman"))
cat(sprintf("\nContinuous, all changeups (n=%d):\n  axis gap vs GB%% over model   r=%+.3f  p=%.3g\n  IVB kill vs raw GB%%          r=%+.3f  p=%.3g\n",
            nrow(S), a$estimate, a$p.value, b$estimate, b$p.value))

ano <- summary(aov(gb_res ~ hi_kill * lo_axis, data = S))[[1]]
cat("\nANOVA on GB over model:\n"); print(round(ano[, c("F value","Pr(>F)")], 4))
