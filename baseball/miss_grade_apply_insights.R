#!/usr/bin/env Rscript

# Apply the two hypothesis features to the OTHER two architectures (beyond the
# grouped winner) and see if they help there:
#   A) all-types single model  -> add path_ratio + axis_diff + spin_eff_diff
#   B) per-pitch-type models    -> path_ratio on breaking types, spin similarity on offspeed types
#
# Marginal value = 2026 holdout RMSE/R2 delta (+ permutation importance, gain rank).
# NOTE: path_ratio is NA on fastballs, so augmented models keep the SAME rows as
# their base (completeness enforced on the base TJStuff+ feats only) and let
# LightGBM handle NA in the added features natively -- a fair marginal test.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(1)
MDIR <- file.path("data","statcast_model")
d  <- readRDS(file.path(MDIR, "miss_grade_data.rds"))
FEAT <- readRDS(file.path(MDIR, "arch_compare.rds"))$FEAT
train <- d[set=="train"]; hold <- d[set=="holdout"]

BREAKING <- c("SL","ST","CU","KC","SV","CS"); OFFSPEED <- c("CH","FS","FO")

# train on rows complete in `req`; model matrix uses `feats` (may contain NA cols -> LGBM handles)
train_lgb <- function(dtr, feats, req, nrounds=2000){
  dtr <- dtr[stats::complete.cases(dtr[, ..req]) & is.finite(miss_distance)]
  n <- nrow(dtr); vi <- sample(n, floor(0.15*n))
  dtrain <- lgb.Dataset(as.matrix(dtr[-vi, ..feats]), label=dtr$miss_distance[-vi])
  dval   <- lgb.Dataset.create.valid(dtrain, as.matrix(dtr[vi, ..feats]), label=dtr$miss_distance[vi])
  lgb.train(params=list(objective="regression", metric="rmse", learning_rate=0.05,
    num_leaves=31, min_data_in_leaf=200, feature_fraction=0.8, bagging_fraction=0.8, bagging_freq=1),
    data=dtrain, nrounds=nrounds, valids=list(val=dval), early_stopping_rounds=60, verbose=-1)
}
rmse <- function(p,a){ok<-is.finite(p)&is.finite(a); sqrt(mean((p[ok]-a[ok])^2))}
r2   <- function(p,a){ok<-is.finite(p)&is.finite(a); p<-p[ok];a<-a[ok];1-sum((a-p)^2)/sum((a-mean(a))^2)}

evalfit <- function(lab, dtr, dte, extra, req=FEAT){
  base <- FEAT; aug <- c(FEAT, extra)
  mb <- train_lgb(dtr, base, req); ma <- train_lgb(dtr, aug, req)
  pb <- predict(mb, as.matrix(dte[, ..base])); pa <- predict(ma, as.matrix(dte[, ..aug]))
  hp <- copy(dte); for(e in extra) hp[[e]] <- sample(hp[[e]])
  pp <- predict(ma, as.matrix(hp[, ..aug]))
  imp <- lgb.importance(ma); gr <- match(extra, imp$Feature)
  cat(sprintf("  %-28s n_te=%6d | RMSE %.4f->%.4f (%+.4f) | R2 %+.4f->%+.4f (%+.4f) | permRise %+.4f\n",
    lab, nrow(dte), rmse(pb,dte$miss_distance), rmse(pa,dte$miss_distance),
    rmse(pa,dte$miss_distance)-rmse(pb,dte$miss_distance),
    r2(pb,dte$miss_distance), r2(pa,dte$miss_distance), r2(pa,dte$miss_distance)-r2(pb,dte$miss_distance),
    rmse(pp,dte$miss_distance)-rmse(pa,dte$miss_distance)))
  cat(sprintf("      added-feat gain rank: %s\n",
    paste(sprintf("%s=#%s(%.1f%%)", extra, gr, 100*imp$Gain[gr]), collapse="  ")))
}

# ================= ARCHITECTURE A: all-types single model =================
cat("=== A) ALL-TYPES single model + path_ratio + spin similarity ===\n")
extraA <- c("path_ratio","axis_diff","spin_eff_diff")
evalfit("all-types (overall holdout)", train, hold, extraA)
cat("  same model, evaluated on subsets:\n")
mb <- train_lgb(train, FEAT, FEAT); ma <- train_lgb(train, c(FEAT,extraA), FEAT)
for (gg in c("fastball","breaking","offspeed")) {
  he <- hold[grp==gg]
  pb <- predict(mb, as.matrix(he[, ..FEAT])); pa <- predict(ma, as.matrix(he[, c(FEAT,extraA), with=FALSE]))
  cat(sprintf("    %-9s n=%6d | RMSE %.4f->%.4f (%+.4f)\n", gg, nrow(he),
    rmse(pb,he$miss_distance), rmse(pa,he$miss_distance), rmse(pa,he$miss_distance)-rmse(pb,he$miss_distance)))
}

# ================= ARCHITECTURE B: per-pitch-type models =================
cat("\n=== B) PER-PITCH-TYPE models + insight for that pitch's family ===\n")
cat("  breaking types + path_ratio:\n")
for (pt in c("SL","ST","CU","KC","SV")) {
  if (train[pitch_type==pt, .N] < 1000) next
  evalfit(sprintf("  %s + path_ratio", pt), train[pitch_type==pt], hold[pitch_type==pt], "path_ratio")
}
cat("  offspeed types + spin similarity:\n")
for (pt in c("CH","FS")) {
  if (train[pitch_type==pt, .N] < 1000) next
  evalfit(sprintf("  %s + spin-sim", pt), train[pitch_type==pt], hold[pitch_type==pt], c("axis_diff","spin_eff_diff"))
}
cat("\nDone.\n")
