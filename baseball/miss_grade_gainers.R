#!/usr/bin/env Rscript

# Approach A (all-types), PUT-AWAY counts (0-2/1-2/2-2): which pitchers x pitch
# types gain the most predicted miss distance from adding the new features
# (path_ratio + axis_diff + spin_eff_diff) on top of the TJStuff+ shape set?
# gain = mean(pred_augmented - pred_base) per pitcher x pitch type on 2026 holdout.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(1)
MDIR <- file.path("data","statcast_model")
d  <- readRDS(file.path(MDIR, "miss_grade_data.rds"))
FEAT  <- readRDS(file.path(MDIR, "arch_compare.rds"))$FEAT
EXTRA <- c("path_ratio","axis_diff","spin_eff_diff")
BREAKING <- c("SL","ST","CU","KC","SV","CS"); OFFSPEED <- c("CH","FS","FO")

d <- d[strikes==2 & balls <= 2]              # put-away counts
train <- d[set=="train"]; hold <- d[set=="holdout"]

train_lgb <- function(dtr, feats, nrounds=2000){
  dtr <- dtr[stats::complete.cases(dtr[, ..FEAT]) & is.finite(miss_distance)]  # same rows, base + aug
  n <- nrow(dtr); vi <- sample(n, floor(0.15*n))
  dtrain <- lgb.Dataset(as.matrix(dtr[-vi, ..feats]), label=dtr$miss_distance[-vi])
  dval   <- lgb.Dataset.create.valid(dtrain, as.matrix(dtr[vi, ..feats]), label=dtr$miss_distance[vi])
  lgb.train(params=list(objective="regression", metric="rmse", learning_rate=0.05,
    num_leaves=31, min_data_in_leaf=150, feature_fraction=0.8, bagging_fraction=0.8, bagging_freq=1),
    data=dtrain, nrounds=nrounds, valids=list(val=dval), early_stopping_rounds=60, verbose=-1)
}

mB <- train_lgb(train, FEAT)
mA <- train_lgb(train, c(FEAT, EXTRA))
hold[, pred_base := predict(mB, as.matrix(.SD[, FEAT, with=FALSE]))]
hold[, pred_aug  := predict(mA, as.matrix(.SD[, c(FEAT,EXTRA), with=FALSE]))]
hold[, gain := pred_aug - pred_base]

agg <- hold[, .(n=.N, act_miss=round(mean(miss_distance),2),
  base=round(mean(pred_base),2), aug=round(mean(pred_aug),2), gain=round(mean(gain),3),
  path_ratio=round(mean(path_ratio,na.rm=TRUE),2), axis_diff=round(mean(axis_diff,na.rm=TRUE),1),
  eff_diff_k=round(mean(spin_eff_diff,na.rm=TRUE)*1000,2)),
  by=.(player_name, pitch_type, grp)][n>=25]

cat("=== TOP 10 BREAKING BALLS that gain most from the new features (approach A, put-away) ===\n")
cat("    (gain = predicted miss WITH new features minus base; path_ratio lower = tighter tunnel)\n")
print(head(agg[grp=="breaking"][order(-gain),
  .(player_name, pitch_type, n, act_miss, base, aug, gain, path_ratio)], 10))

cat("\n=== TOP 10 OFFSPEED that gain most from the new features (approach A, put-away) ===\n")
cat("    (axis_diff = spin-direction gap vs FB in deg; eff_diff_k = spin-eff gap vs FB, x1000)\n")
print(head(agg[grp=="offspeed"][order(-gain),
  .(player_name, pitch_type, n, act_miss, base, aug, gain, axis_diff, eff_diff_k)], 10))

fwrite(agg[order(-gain)], file.path(MDIR, "putaway_gainers_A.csv"))
cat(sprintf("\nWrote %s\n", file.path(MDIR,"putaway_gainers_A.csv")))
