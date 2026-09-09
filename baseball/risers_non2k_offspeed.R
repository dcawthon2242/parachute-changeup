#!/usr/bin/env Rscript

# NON-2-strike offspeed risers, gated Approach A (all-types).
# Spin-similarity features (axis_diff + measured active-spin gap) gated to offspeed only,
# tunneling path_ratio gated to breaking only. gain = pred WITH new features - WITHOUT.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(1)
MDIR <- file.path("data","statcast_model")
FEAT <- readRDS(file.path(MDIR, "arch_compare.rds"))$FEAT
as_long <- readRDS(file.path(MDIR, "active_spin_long.rds"))
d <- readRDS(file.path(MDIR, "miss_grade_data_activespin.rds"))

fb <- as_long[pitch_type %in% c("FF","SI","FC")][, pr:=match(pitch_type,c("FF","SI","FC"))][
  order(pitcher,season,pr)][, .SD[1], by=.(pitcher,season)][, .(pitcher,season,fb_active=active_spin)]
d <- merge(d, fb, by=c("pitcher","season"), all.x=TRUE)
d[, meas_eff_diff := active_spin - fb_active]
d[, path_ratio_g := ifelse(grp=="breaking", path_ratio, NA_real_)]
d[, axis_diff_g  := ifelse(grp=="offspeed", axis_diff,   NA_real_)]
d[, eff_diff_g   := ifelse(grp=="offspeed", meas_eff_diff, NA_real_)]
A_X <- c("path_ratio_g","axis_diff_g","eff_diff_g")

train_lgb <- function(dtr, feats, nr=1500){
  dtr <- dtr[is.finite(miss_distance)]; n <- nrow(dtr); vi <- sample(n, floor(0.15*n))
  dtrain <- lgb.Dataset(as.matrix(dtr[-vi, ..feats]), label=dtr$miss_distance[-vi])
  dval   <- lgb.Dataset.create.valid(dtrain, as.matrix(dtr[vi, ..feats]), label=dtr$miss_distance[vi])
  lgb.train(params=list(objective="regression",metric="rmse",learning_rate=0.05,num_leaves=31,
    min_data_in_leaf=150,feature_fraction=0.8,bagging_fraction=0.8,bagging_freq=1),
    data=dtrain, nrounds=nr, valids=list(val=dval), early_stopping_rounds=50, verbose=-1)
}
rmse <- function(p,a){ok<-is.finite(p)&is.finite(a); sqrt(mean((p[ok]-a[ok])^2))}

dd <- d[strikes < 2]
tr <- dd[set=="train"]; ho <- dd[set=="holdout" & is.finite(miss_distance)]
cat(sprintf("NON-2K: train=%d holdout=%d (mean miss=%.2f)\n",
    nrow(tr[is.finite(miss_distance)]), nrow(ho), mean(ho$miss_distance)))

mB <- train_lgb(tr, FEAT); mA <- train_lgb(tr, c(FEAT, A_X))
ho[, pred_base := predict(mB, as.matrix(.SD[, FEAT, with=FALSE]))]
ho[, pred_aug  := predict(mA, as.matrix(.SD[, c(FEAT,A_X), with=FALSE]))]
cat(sprintf("holdout RMSE  base=%.4f  aug=%.4f  delta=%+.4f\n",
    rmse(ho$pred_base,ho$miss_distance), rmse(ho$pred_aug,ho$miss_distance),
    rmse(ho$pred_aug,ho$miss_distance)-rmse(ho$pred_base,ho$miss_distance)))
ho[, gain := pred_aug - pred_base]

agg <- ho[grp=="offspeed", .(n=.N, act=round(mean(miss_distance),2),
  base=round(mean(pred_base),2), aug=round(mean(pred_aug),2), gain=round(mean(gain),3),
  axis_diff=round(mean(axis_diff,na.rm=TRUE),1), as_gap=round(mean(meas_eff_diff,na.rm=TRUE),3)),
  by=.(player_name, pitch_type)][n>=40]
cat("\n=== TOP 10 OFFSPEED RISERS (non-2K, gated approach A) ===\n")
print(head(agg[order(-gain), .(player_name,pitch_type,n,act,base,aug,gain,axis_diff,as_gap)],10))
fwrite(agg[order(-gain)], file.path(MDIR, "risers_non2K_offspeed.csv"))
cat("\nsaved risers_non2K_offspeed.csv\n")
