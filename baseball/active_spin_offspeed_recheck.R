#!/usr/bin/env Rscript

# Rebuild the offspeed "FB spin-similarity" feature with MEASURED active spin and
# test whether it tracks miss distance better than my inferred spin_eff_diff.
#  measured_eff_diff = active_spin(offspeed) - active_spin(primary FB)   [per pitcher-season]
# Compare: (a) correlation with miss_distance, (b) marginal RMSE in the offspeed model.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(1)
MDIR <- file.path("data","statcast_model")
as_long <- readRDS(file.path(MDIR, "active_spin_long.rds"))
d <- readRDS(file.path(MDIR, "miss_grade_data_activespin.rds"))   # already has per-pitch active_spin

# ---- primary fastball measured active spin per pitcher-season (FF>SI>FC) ----
fb <- as_long[pitch_type %in% c("FF","SI","FC")]
fb[, pr := match(pitch_type, c("FF","SI","FC"))]
fb <- fb[order(pitcher, season, pr)][, .SD[1], by=.(pitcher, season)][, .(pitcher, season, fb_active=active_spin)]
d <- merge(d, fb, by=c("pitcher","season"), all.x=TRUE)
d[, meas_eff_diff := active_spin - fb_active]     # measured similarity vs FB

off <- d[grp=="offspeed" & !is.na(miss_distance)]
cat(sprintf("offspeed swings: %d | with meas_eff_diff: %d | with inferred spin_eff_diff: %d\n",
    nrow(off), sum(!is.na(off$meas_eff_diff)), sum(!is.na(off$spin_eff_diff))))

# ---- (a) correlation with miss distance ----
cc <- off[!is.na(meas_eff_diff) & !is.na(spin_eff_diff)]
cat(sprintf("\n=== corr with miss_distance (offspeed, n=%d) ===\n", nrow(cc)))
cat(sprintf("  inferred spin_eff_diff : Pearson %.4f  Spearman %.4f\n",
    cor(cc$spin_eff_diff, cc$miss_distance), cor(cc$spin_eff_diff, cc$miss_distance, method="spearman")))
cat(sprintf("  MEASURED meas_eff_diff : Pearson %.4f  Spearman %.4f\n",
    cor(cc$meas_eff_diff, cc$miss_distance), cor(cc$meas_eff_diff, cc$miss_distance, method="spearman")))
cat(sprintf("  measured axis_diff     : Pearson %.4f  Spearman %.4f\n",
    cor(cc$axis_diff, cc$miss_distance), cor(cc$axis_diff, cc$miss_distance, method="spearman")))

# pitcher x pitch aggregate (deception is a pitch-level trait)
ag <- cc[, .(miss=mean(miss_distance), meas=mean(meas_eff_diff), inf=mean(spin_eff_diff),
             ax=mean(axis_diff), n=.N), by=.(pitcher, pitch_type)][n>=40]
cat(sprintf("\naggregate (pitcher x pt, n>=40: %d)\n", nrow(ag)))
cat(sprintf("  inferred vs miss  r=%.3f | MEASURED vs miss  r=%.3f | axis_diff vs miss r=%.3f\n",
    cor(ag$inf, ag$miss), cor(ag$meas, ag$miss), cor(ag$ax, ag$miss)))

# ---- (b) marginal RMSE in the offspeed model ----
FEAT <- readRDS(file.path(MDIR, "arch_compare.rds"))$FEAT
train_lgb <- function(dtr, feats, nr=2000){
  dtr <- dtr[stats::complete.cases(dtr[, ..feats]) & is.finite(miss_distance)]
  n<-nrow(dtr); vi<-sample(n, floor(0.15*n))
  dtrain<-lgb.Dataset(as.matrix(dtr[-vi,..feats]),label=dtr$miss_distance[-vi])
  dval<-lgb.Dataset.create.valid(dtrain,as.matrix(dtr[vi,..feats]),label=dtr$miss_distance[vi])
  lgb.train(params=list(objective="regression",metric="rmse",learning_rate=0.05,num_leaves=31,
    min_data_in_leaf=150,feature_fraction=0.8,bagging_fraction=0.8,bagging_freq=1),
    data=dtrain,nrounds=nr,valids=list(val=dval),early_stopping_rounds=60,verbose=-1)
}
rmse<-function(p,a){ok<-is.finite(p)&is.finite(a);sqrt(mean((p[ok]-a[ok])^2))}

# same rows for a fair base/inferred/measured comparison (need both diffs present)
d[, feat_ok := stats::complete.cases(.SD), .SDcols = FEAT]
dd <- d[grp=="offspeed" & is.finite(miss_distance) & feat_ok &
        !is.na(spin_eff_diff) & !is.na(meas_eff_diff) & !is.na(axis_diff)]
tr <- dd[set=="train"]; ho <- dd[set=="holdout"]
cat(sprintf("\n=== offspeed model, common-row holdout n=%d ===\n", nrow(ho)))
mk <- function(extra){ m<-train_lgb(tr, c(FEAT, extra)); rmse(predict(m, as.matrix(ho[, c(FEAT,extra), with=FALSE])), ho$miss_distance) }
b   <- mk(character(0))
inf <- mk(c("axis_diff","spin_eff_diff"))
mea <- mk(c("axis_diff","meas_eff_diff"))
both<- mk(c("axis_diff","spin_eff_diff","meas_eff_diff"))
cat(sprintf("  base (shape only)                : RMSE %.4f\n", b))
cat(sprintf("  + axis_diff + INFERRED eff_diff  : RMSE %.4f  (delta %+.4f)\n", inf, inf-b))
cat(sprintf("  + axis_diff + MEASURED eff_diff  : RMSE %.4f  (delta %+.4f)\n", mea, mea-b))
cat(sprintf("  + axis_diff + BOTH eff diffs     : RMSE %.4f  (delta %+.4f)\n", both, both-b))
