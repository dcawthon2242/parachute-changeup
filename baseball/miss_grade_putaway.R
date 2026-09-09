#!/usr/bin/env Rscript

# Restrict the miss-distance model to PUT-AWAY counts (0-2, 1-2, 2-2; two strikes,
# full count excluded) and re-run all three architectures, then re-test the two
# insight features across each architecture on this subset.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(1)
MDIR <- file.path("data","statcast_model")
d  <- readRDS(file.path(MDIR, "miss_grade_data.rds"))
FEAT <- readRDS(file.path(MDIR, "arch_compare.rds"))$FEAT

# put-away counts only
d <- d[strikes==2 & balls <= 2]
train <- d[set=="train"]; hold <- d[set=="holdout"]
cat(sprintf("PUT-AWAY counts (0-2/1-2/2-2): train=%d  holdout=%d  (holdout whiff%%=%.1f, mean miss=%.2f)\n",
    nrow(train), nrow(hold), 100*mean(hold$is_whiff), mean(hold$miss_distance)))

train_lgb <- function(dtr, feats, req=feats, nrounds=2000){
  dtr <- dtr[stats::complete.cases(dtr[, ..req]) & is.finite(miss_distance)]
  n <- nrow(dtr); vi <- sample(n, floor(0.15*n))
  dtrain <- lgb.Dataset(as.matrix(dtr[-vi, ..feats]), label=dtr$miss_distance[-vi])
  dval   <- lgb.Dataset.create.valid(dtrain, as.matrix(dtr[vi, ..feats]), label=dtr$miss_distance[vi])
  lgb.train(params=list(objective="regression", metric="rmse", learning_rate=0.05,
    num_leaves=31, min_data_in_leaf=150, feature_fraction=0.8, bagging_fraction=0.8, bagging_freq=1),
    data=dtrain, nrounds=nrounds, valids=list(val=dval), early_stopping_rounds=60, verbose=-1)
}
rmse <- function(p,a){ok<-is.finite(p)&is.finite(a); sqrt(mean((p[ok]-a[ok])^2))}
r2   <- function(p,a){ok<-is.finite(p)&is.finite(a); p<-p[ok];a<-a[ok];1-sum((a-p)^2)/sum((a-mean(a))^2)}
sp   <- function(p,a){ok<-is.finite(p)&is.finite(a); suppressWarnings(cor(p[ok],a[ok],method="spearman"))}

# ================= architecture comparison (base TJStuff+) =================
mA <- train_lgb(train, FEAT)
hold[, predA := predict(mA, as.matrix(.SD[, FEAT, with=FALSE]))]
hold[, predC := NA_real_]
for (gg in c("fastball","breaking","offspeed")) {
  mg <- train_lgb(train[grp==gg], FEAT)
  hold[grp==gg, predC := predict(mg, as.matrix(.SD[, FEAT, with=FALSE]))]
}
hold[is.na(predC), predC := predict(mA, as.matrix(.SD[, FEAT, with=FALSE]))]
hold[, predB := NA_real_]
types <- train[, .N, by=pitch_type][N>=500]$pitch_type
for (pt in types) { mt <- train_lgb(train[pitch_type==pt], FEAT)
  hold[pitch_type==pt, predB := predict(mt, as.matrix(.SD[, FEAT, with=FALSE]))] }
hold[is.na(predB), predB := predict(mA, as.matrix(.SD[, FEAT, with=FALSE]))]

cat("\n=== PUT-AWAY holdout: architecture comparison (base TJStuff+) ===\n")
res <- data.table(
  arch=c("A: all-types","B: per-pitch-type","C: FB/BRK/OFF group"),
  rmse=c(rmse(hold$predA,hold$miss_distance), rmse(hold$predB,hold$miss_distance), rmse(hold$predC,hold$miss_distance)),
  r2  =c(r2(hold$predA,hold$miss_distance),   r2(hold$predB,hold$miss_distance),   r2(hold$predC,hold$miss_distance)),
  spearman=c(sp(hold$predA,hold$miss_distance), sp(hold$predB,hold$miss_distance), sp(hold$predC,hold$miss_distance)))
print(res[, .(arch, rmse=round(rmse,4), r2=round(r2,4), spearman=round(spearman,4))])
cat(sprintf("WINNER: %s\n", res$arch[which.min(res$rmse)]))

# ================= insight features across the 3 approaches =================
evalfit <- function(lab, dtr, dte, extra){
  base <- FEAT; aug <- c(FEAT, extra)
  mb <- train_lgb(dtr, base, FEAT); ma <- train_lgb(dtr, aug, FEAT)
  pb <- predict(mb, as.matrix(dte[, ..base])); pa <- predict(ma, as.matrix(dte[, ..aug]))
  imp <- lgb.importance(ma); gr <- match(extra, imp$Feature)
  cat(sprintf("  %-26s n_te=%6d | RMSE %.4f->%.4f (%+.4f) | R2 %+.4f (%+.4f) | gain %s\n",
    lab, nrow(dte), rmse(pb,dte$miss_distance), rmse(pa,dte$miss_distance),
    rmse(pa,dte$miss_distance)-rmse(pb,dte$miss_distance), r2(pa,dte$miss_distance),
    r2(pa,dte$miss_distance)-r2(pb,dte$miss_distance),
    paste(sprintf("%s=#%s(%.1f%%)",extra,gr,100*imp$Gain[gr]),collapse=" ")))
}
extraA <- c("path_ratio","axis_diff","spin_eff_diff")
cat("\n=== A) all-types + path_ratio + spin-sim ===\n")
evalfit("all-types overall", train, hold, extraA)
cat("\n=== C) grouped: breaking+path_ratio, offspeed+spin-sim ===\n")
evalfit("breaking + path_ratio", train[grp=="breaking"], hold[grp=="breaking"], "path_ratio")
evalfit("offspeed + spin-sim",  train[grp=="offspeed"], hold[grp=="offspeed"], c("axis_diff","spin_eff_diff"))
cat("\n=== B) per-pitch-type ===\n")
for (pt in c("SL","ST","CU","KC")) if (train[pitch_type==pt,.N]>=500)
  evalfit(sprintf("%s + path_ratio",pt), train[pitch_type==pt], hold[pitch_type==pt], "path_ratio")
for (pt in c("CH","FS")) if (train[pitch_type==pt,.N]>=500)
  evalfit(sprintf("%s + spin-sim",pt), train[pitch_type==pt], hold[pitch_type==pt], c("axis_diff","spin_eff_diff"))
cat("\nDone.\n")
