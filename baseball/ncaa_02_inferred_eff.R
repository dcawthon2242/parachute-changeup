#!/usr/bin/env Rscript

# CAN SPIN EFFICIENCY BE INFERRED FROM WHAT TRACKMAN ACTUALLY REPORTS?
#
# TrackMan has no efficiency column. The standard workaround - what Rapsodo and most college
# labs do - is to exploit the fact that gyro spin produces no Magnus force. Total spin is
# measured; the movement a pitch WOULD have if all of that spin were transverse can be
# modelled; and the ratio of observed movement to that ceiling estimates the transverse
# fraction. Note this is a magnitude argument, not a direction one: gyro spin shrinks the break
# without rotating it, which is also why the axis-versus-break test in the previous script
# speaks to seam-shifted wake rather than to efficiency.
#
# Whether that inference is good enough is an empirical question, and MLB can answer it because
# it has both sides: Statcast reports the same physical inputs TrackMan does AND publishes
# measured active spin from 3D seam tracking. So the estimator is built from TrackMan-available
# inputs only, then scored against the measured truth. Whatever accuracy it achieves on MLB is
# the ceiling for what it can achieve on NCAA data.
#
# An earlier pass in this project found inferred efficiency unreliable and switched to the
# measured leaderboard. This quantifies how unreliable, which is what decides the NCAA port.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(5); options(width = 200)
MDIR <- "data/statcast_model"

# parachute_ff.rds holds changeups only; the gate needs the four-seamer too, so the wider file
# is used here.
CH <- readRDS(file.path(MDIR, "parachute_rv.rds"))
AS <- readRDS(file.path(MDIR, "active_spin_long.rds"))
CH <- CH[pitch_type %in% c("CH","FF")]

# Magnus acceleration magnitude, gravity removed. This is the "observed break" side.
CH[, `:=`(mag_x = ax, mag_z = az + 32.174)]
CH[, magnus := sqrt(mag_x^2 + mag_z^2)]

# The physical ceiling. Magnus force scales with the transverse spin rate and with velocity, so
# a spin rate that buys a lot of break at 95 mph buys less at 80. Rather than impose a lift
# coefficient, the ratio is formed and the mapping to efficiency is learned.
CH[, ratio := magnus / (release_spin_rate * release_speed / 1000)]

P <- CH[is.finite(ratio) & is.finite(release_spin_rate) & release_spin_rate > 500,
        .(n = .N, ratio = mean(ratio), magnus = mean(magnus), spin = mean(release_spin_rate),
          velo = mean(release_speed), ext = mean(release_extension, na.rm = TRUE)),
        by = .(pitcher, season, pitch_type)]
P <- merge(P, AS[, .(pitcher, season, pitch_type, measured = active_spin)],
           by = c("pitcher","season","pitch_type"))[n >= 50 & is.finite(measured)]
cat(sprintf("pitcher-season-pitch-types with both inferred inputs and measured active spin: %d\n",
            nrow(P)))
cat(sprintf("pitch types: %s\n\n", paste(sort(unique(P$pitch_type)), collapse = ", ")))

cat("=== simple ratio estimator vs measured active spin ===\n")
for (pt in c("CH","FF")) { D <- P[pitch_type == pt]
  if (nrow(D) < 30) next
  ct <- cor.test(D$ratio, D$measured)
  cat(sprintf("  %s  n=%4d  r = %+.3f  [%+.3f, %+.3f]  R2 = %.3f\n", pt, nrow(D), ct$estimate,
              ct$conf.int[1], ct$conf.int[2], ct$estimate^2)) }
ct <- cor.test(P$ratio, P$measured)
cat(sprintf("  all n=%4d  r = %+.3f  R2 = %.3f\n", nrow(P), ct$estimate, ct$estimate^2))

# A learned version, given every input TrackMan exposes, is the best case for the inference.
cat("\n=== best case: gradient boosting on all TrackMan-available inputs, out of fold ===\n")
F <- c("ratio","magnus","spin","velo","ext")
K <- 5; fold <- sample(rep(1:K, length.out = nrow(P))); P[, pred := NA_real_]
for (f in 1:K) {
  tr <- P[fold != f]
  m <- lgb.train(params = list(objective = "regression", learning_rate = .05, num_leaves = 15,
                 min_data_in_leaf = 40, feature_fraction = .9), verbose = -1, nrounds = 400,
                 data = lgb.Dataset(as.matrix(tr[, ..F]), label = tr$measured))
  P$pred[fold == f] <- predict(m, as.matrix(P[fold == f, ..F]))
}
cat(sprintf("  overall  r = %+.3f  R2 = %.3f  RMSE = %.3f (measured active spin sd = %.3f)\n",
            cor(P$pred, P$measured), cor(P$pred, P$measured)^2,
            sqrt(mean((P$pred - P$measured)^2)), sd(P$measured)))
for (pt in c("CH","FF")) { D <- P[pitch_type == pt]
  cat(sprintf("  %s       r = %+.3f  R2 = %.3f  RMSE = %.3f\n", pt, cor(D$pred, D$measured),
              cor(D$pred, D$measured)^2, sqrt(mean((D$pred - D$measured)^2)))) }

# The question is not correlation in the abstract but whether the inference can reproduce the
# gate that actually matters: does a pitcher clear .85 active spin on both pitches?
cat("\n=== can the inference reproduce the .85 gate the Core bin needs? ===\n")
W <- dcast(P, pitcher + season ~ pitch_type, value.var = c("measured","pred"))
W <- W[is.finite(measured_CH) & is.finite(measured_FF) & is.finite(pred_CH) & is.finite(pred_FF)]
truth <- W$measured_CH >= .85 & W$measured_FF >= .85
guess <- W$pred_CH >= .85 & W$pred_FF >= .85
cat(sprintf("  %d pitcher-seasons with both pitches\n", nrow(W)))
cat(sprintf("  truly above the gate: %d   inferred above: %d   agree: %d\n",
            sum(truth), sum(guess), sum(truth == guess)))
cat(sprintf("  precision %.1f%% (of those the inference admits, this share really belong)\n",
            100*sum(truth & guess)/max(sum(guess),1)))
cat(sprintf("  recall    %.1f%% (of those that really belong, this share are found)\n",
            100*sum(truth & guess)/max(sum(truth),1)))
cat(sprintf("  agreement %.1f%%, versus %.1f%% from always guessing the majority class\n",
            100*mean(truth == guess), 100*max(mean(truth), 1-mean(truth))))

cat("\n=== and the efficiency GAP between the two pitches, which the Core bin uses ===\n")
ct <- cor.test(W$pred_CH - W$pred_FF, W$measured_CH - W$measured_FF)
cat(sprintf("  inferred gap vs measured gap: r = %+.3f  R2 = %.3f\n", ct$estimate, ct$estimate^2))
tg <- abs(W$measured_CH - W$measured_FF) <= .10; gg <- abs(W$pred_CH - W$pred_FF) <= .10
cat(sprintf("  |gap| <= .10 test: precision %.1f%%, recall %.1f%%, agreement %.1f%%\n",
            100*sum(tg & gg)/max(sum(gg),1), 100*sum(tg & gg)/max(sum(tg),1), 100*mean(tg == gg)))
