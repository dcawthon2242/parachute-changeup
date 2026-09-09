#!/usr/bin/env Rscript

# NON-2-strike counts (strikes < 2): confirm the optimal-RMSE approach, then list
# the biggest gainers from the new features (path_ratio + spin similarity).
# Optimal RMSE config tested = augmented all-types (approach A + new features).

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(1)
MDIR <- file.path("data","statcast_model")
d  <- readRDS(file.path(MDIR, "miss_grade_data.rds"))
FEAT  <- readRDS(file.path(MDIR, "arch_compare.rds"))$FEAT
EXTRA <- c("path_ratio","axis_diff","spin_eff_diff")

d <- d[strikes < 2]                          # non-2-strike counts
train <- d[set=="train"]; hold <- d[set=="holdout"]
cat(sprintf("NON-2-STRIKE counts: train=%d holdout=%d (holdout whiff%%=%.1f, mean miss=%.2f)\n",
    nrow(train), nrow(hold), 100*mean(hold$is_whiff), mean(hold$miss_distance)))

train_lgb <- function(dtr, feats, nrounds=2000){
  dtr <- dtr[stats::complete.cases(dtr[, ..FEAT]) & is.finite(miss_distance)]
  n <- nrow(dtr); vi <- sample(n, floor(0.15*n))
  dtrain <- lgb.Dataset(as.matrix(dtr[-vi, ..feats]), label=dtr$miss_distance[-vi])
  dval   <- lgb.Dataset.create.valid(dtrain, as.matrix(dtr[vi, ..feats]), label=dtr$miss_distance[vi])
  lgb.train(params=list(objective="regression", metric="rmse", learning_rate=0.05,
    num_leaves=31, min_data_in_leaf=150, feature_fraction=0.8, bagging_fraction=0.8, bagging_freq=1),
    data=dtrain, nrounds=nrounds, valids=list(val=dval), early_stopping_rounds=60, verbose=-1)
}
rmse <- function(p,a){ok<-is.finite(p)&is.finite(a); sqrt(mean((p[ok]-a[ok])^2))}

# ---- architecture RMSE check on this subset (base) + augmented A ----
mB <- train_lgb(train, FEAT)
mA <- train_lgb(train, c(FEAT, EXTRA))
hold[, pred_base := predict(mB, as.matrix(.SD[, FEAT, with=FALSE]))]
hold[, pred_aug  := predict(mA, as.matrix(.SD[, c(FEAT,EXTRA), with=FALSE]))]
predC <- rep(NA_real_, nrow(hold))
for (gg in c("fastball","breaking","offspeed")) { mg <- train_lgb(train[grp==gg], FEAT)
  idx <- which(hold$grp==gg); predC[idx] <- predict(mg, as.matrix(hold[idx, FEAT, with=FALSE])) }
predC[is.na(predC)] <- predict(mB, as.matrix(hold[is.na(predC), FEAT, with=FALSE]))
cat(sprintf("\nHoldout RMSE:  A base=%.4f | C grouped base=%.4f | A augmented=%.4f\n",
    rmse(hold$pred_base,hold$miss_distance), rmse(predC,hold$miss_distance), rmse(hold$pred_aug,hold$miss_distance)))

# ---- gainers from the new features (approach A augmented) ----
hold[, gain := pred_aug - pred_base]
agg <- hold[, .(n=.N, act_miss=round(mean(miss_distance),2),
  base=round(mean(pred_base),2), aug=round(mean(pred_aug),2), gain=round(mean(gain),3),
  path_ratio=round(mean(path_ratio,na.rm=TRUE),2), axis_diff=round(mean(axis_diff,na.rm=TRUE),1)),
  by=.(player_name, pitch_type, grp)][n>=40]

cat("\n=== TOP 10 BREAKING BALLS gaining most from new features (non-2-strike, approach A) ===\n")
print(head(agg[grp=="breaking"][order(-gain), .(player_name,pitch_type,n,act_miss,base,aug,gain,path_ratio)],10))
cat("\n=== TOP 10 OFFSPEED gaining most from new features (non-2-strike, approach A) ===\n")
print(head(agg[grp=="offspeed"][order(-gain), .(player_name,pitch_type,n,act_miss,base,aug,gain,axis_diff)],10))
fwrite(agg[order(-gain)], file.path(MDIR, "non2k_gainers_A.csv"))
cat(sprintf("\nWrote %s\n", file.path(MDIR,"non2k_gainers_A.csv")))
