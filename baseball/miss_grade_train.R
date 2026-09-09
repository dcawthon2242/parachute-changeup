#!/usr/bin/env Rscript

# Train + compare the three architectures for the miss-distance pitch grade model:
#   A) single model, all pitch types (TJStuff+ default: shape only, no type label)
#   B) one model per pitch type
#   C) three grouped models: fastball / breaking / offspeed
# Selection metric: 2026 holdout RMSE (+ R2, Spearman). Writes the winner choice.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(1)
MDIR <- file.path("data","statcast_model")
d <- readRDS(file.path(MDIR, "miss_grade_data.rds"))

FEAT <- c("release_speed","release_spin_rate","release_extension",
          "release_pos_x","release_pos_z","sax","cax","ax","az",
          "speed_diff","ax_diff","az_diff")

train <- d[set=="train"]; hold <- d[set=="holdout"]
cat(sprintf("train rows=%d  holdout rows=%d\n", nrow(train), nrow(hold)))

train_lgb <- function(dtr, feats, nrounds=2000){
  dtr <- dtr[stats::complete.cases(dtr[, ..feats]) & is.finite(miss_distance)]
  n <- nrow(dtr); vi <- sample(n, floor(0.15*n))
  Xtr <- as.matrix(dtr[-vi, ..feats]); Xv <- as.matrix(dtr[vi, ..feats])
  dtrain <- lgb.Dataset(Xtr, label=dtr$miss_distance[-vi])
  dval   <- lgb.Dataset.create.valid(dtrain, Xv, label=dtr$miss_distance[vi])
  lgb.train(params=list(objective="regression", metric="rmse", learning_rate=0.05,
    num_leaves=31, min_data_in_leaf=200, feature_fraction=0.8,
    bagging_fraction=0.8, bagging_freq=1),
    data=dtrain, nrounds=nrounds, valids=list(val=dval),
    early_stopping_rounds=60, verbose=-1)
}
pred_lgb <- function(m, dd, feats) predict(m, as.matrix(dd[, ..feats]))
met <- function(pred, act){
  ok <- is.finite(pred)&is.finite(act); pred<-pred[ok]; act<-act[ok]
  c(rmse=sqrt(mean((pred-act)^2)),
    r2=1-sum((act-pred)^2)/sum((act-mean(act))^2),
    spearman=suppressWarnings(cor(pred,act,method="spearman")), n=length(pred))
}

# ---- A: all types ----
mA <- train_lgb(train, FEAT)
hold[, predA := pred_lgb(mA, hold, FEAT)]

# ---- C: grouped fastball/breaking/offspeed ----
hold[, predC := NA_real_]
for (gg in c("fastball","breaking","offspeed")) {
  mg <- train_lgb(train[grp==gg], FEAT)
  hold[grp==gg, predC := pred_lgb(mg, hold[grp==gg], FEAT)]
}
hold[is.na(predC), predC := pred_lgb(mA, hold[is.na(predC)], FEAT)]  # 'other' fallback

# ---- B: per pitch type (fallback to A for sparse types) ----
hold[, predB := NA_real_]
types <- train[, .N, by=pitch_type][N>=1000]$pitch_type
for (pt in types) {
  mt <- train_lgb(train[pitch_type==pt], FEAT)
  hold[pitch_type==pt, predB := pred_lgb(mt, hold[pitch_type==pt], FEAT)]
}
hold[is.na(predB), predB := pred_lgb(mA, hold[is.na(predB)], FEAT)]

# ---- evaluate ----
cat("\n=== 2026 HOLDOUT metrics (target = miss_distance | swing) ===\n")
res <- rbind(
  data.table(arch="A: all-types",       t(met(hold$predA, hold$miss_distance))),
  data.table(arch="B: per-pitch-type",  t(met(hold$predB, hold$miss_distance))),
  data.table(arch="C: FB/BRK/OFF group",t(met(hold$predC, hold$miss_distance))))
print(res[, .(arch, rmse=round(rmse,4), r2=round(r2,4), spearman=round(spearman,4), n)])

cat("\n=== Holdout RMSE by pitch type (each arch) ===\n")
bt <- hold[, .(n=.N,
  A=sqrt(mean((predA-miss_distance)^2)),
  B=sqrt(mean((predB-miss_distance)^2)),
  C=sqrt(mean((predC-miss_distance)^2))), by=pitch_type][order(-n)]
print(bt[, .(pitch_type, n, A=round(A,3), B=round(B,3), C=round(C,3))])

winner <- res$arch[which.min(res$rmse)]
cat(sprintf("\nWINNER (lowest holdout RMSE): %s\n", winner))
saveRDS(list(winner=winner, res=res, by_type=bt, FEAT=FEAT),
        file.path(MDIR, "arch_compare.rds"))
fwrite(hold[, .(season, pitcher, player_name, pitch_type, grp, miss_distance, predA, predB, predC)],
       file.path(MDIR, "holdout_preds.csv"))
cat(sprintf("Wrote %s and holdout_preds.csv\n", file.path(MDIR,"arch_compare.rds")))
