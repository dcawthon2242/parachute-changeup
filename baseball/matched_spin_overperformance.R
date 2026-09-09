#!/usr/bin/env Rscript

# Prove: offspeed pitches that SPIN LIKE THE FASTBALL overperform their shape
# expectation. Overperformance = actual miss - expected miss from a shape-only
# (TJStuff+) model trained on 2023-25 and applied to 2026. Then relate that
# residual to spin SIMILARITY (small axis gap AND small active-spin gap vs FB).

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(1)
MDIR <- file.path("data","statcast_model")
FEAT <- readRDS(file.path(MDIR, "arch_compare.rds"))$FEAT
as_long <- readRDS(file.path(MDIR, "active_spin_long.rds"))
d <- readRDS(file.path(MDIR, "miss_grade_data_activespin.rds"))
fb <- as_long[pitch_type %in% c("FF","SI","FC")][, pr:=match(pitch_type,c("FF","SI","FC"))][
  order(pitcher,season,pr)][, .SD[1], by=.(pitcher,season)][, .(pitcher,season,fb_active=active_spin)]
d <- merge(d, fb, by=c("pitcher","season"), all.x=TRUE)
d[, as_gap := active_spin - fb_active]          # signed active-spin gap vs FB
d[, abs_as := abs(as_gap)]                        # active-spin DISSIMILARITY
# spin SIMILARITY score, driven by MEASURED ACTIVE SPIN (what Cease & Vesia share),
# with a softer axis term. 1 = spins like the FB, 0 = very different.
d[, spin_sim := exp(-(abs_as/0.10)^2) * exp(-(axis_diff/45)^2)]

# ---- shape-only expectation: train 2023-25, predict 2026 ----
tr <- d[set=="train" & is.finite(miss_distance)]
n <- nrow(tr); vi <- sample(n, floor(0.15*n))
dtrain <- lgb.Dataset(as.matrix(tr[-vi, ..FEAT]), label=tr$miss_distance[-vi])
dval   <- lgb.Dataset.create.valid(dtrain, as.matrix(tr[vi, ..FEAT]), label=tr$miss_distance[vi])
m <- lgb.train(params=list(objective="regression",metric="rmse",learning_rate=0.05,num_leaves=31,
  min_data_in_leaf=150,feature_fraction=0.8,bagging_fraction=0.8,bagging_freq=1),
  data=dtrain, nrounds=2000, valids=list(val=dval), early_stopping_rounds=60, verbose=-1)

ho <- d[set=="holdout" & is.finite(miss_distance) & grp=="offspeed"]
ho[, exp_miss := predict(m, as.matrix(.SD[, FEAT, with=FALSE]))]
ho[, resid := miss_distance - exp_miss]           # + = overperforms shape

ag <- ho[, .(n=.N, miss=round(mean(miss_distance),2), exp=round(mean(exp_miss),2),
  resid=round(mean(resid),3), axis_diff=round(mean(axis_diff,na.rm=TRUE),1),
  as_gap=round(mean(as_gap,na.rm=TRUE),3), spin_sim=round(mean(spin_sim,na.rm=TRUE),3)),
  by=.(player_name, pitch_type)][n>=40]

CO <- function(a,b,...) cor(a,b,use="complete.obs",...)
cat("=== does SPIN SIMILARITY predict overperformance? (2026 offspeed, pitcher x pt, n>=40) ===\n")
cat(sprintf("  n groups=%d\n", nrow(ag)))
cat(sprintf("  resid vs spin_sim (match->overperform) : r=%+.3f  (Spearman %+.3f)\n",
    CO(ag$resid, ag$spin_sim), CO(ag$resid, ag$spin_sim, method="spearman")))
cat(sprintf("  resid vs axis_diff (bigger axis gap)   : r=%+.3f\n", CO(ag$resid, ag$axis_diff)))
cat(sprintf("  resid vs |active-spin gap| (dissimilar): r=%+.3f\n", CO(ag$resid, abs(ag$as_gap))))
cat(sprintf("  resid vs signed active-spin gap        : r=%+.3f\n", CO(ag$resid, ag$as_gap)))

cat("\n--- matched-active-spin cohort (|as_gap|<=0.05) vs the rest ---\n")
print(ag[!is.na(as_gap), .(n_pitchers=.N, mean_resid=round(mean(resid),3), mean_miss=round(mean(miss),2),
    mean_exp=round(mean(exp),2)), by=.(matched = abs(as_gap)<=0.05)][order(-matched)])
cat("\n--- by active-spin similarity tier ---\n")
print(ag[!is.na(as_gap), .(n=.N, mean_resid=round(mean(resid),3)),
    by=.(tier=cut(abs(as_gap), c(-.01,.03,.07,.15,1), labels=c("<=.03 (near-identical)",".03-.07",".07-.15",">.15")))][order(tier)])

cat("\n=== CHECK: Cease & Vesia ===\n")
print(ho[grepl("Cease|Vesia", player_name), .(player_name, pitch_type,
  miss=round(mean(miss_distance),2), exp=round(mean(exp_miss),2), resid=round(mean(resid),3),
  axis_diff=round(mean(axis_diff),1), as_gap=round(mean(as_gap),3), spin_sim=round(mean(spin_sim),3)),
  by=.(player_name,pitch_type)])

cat("\n=== TOP 15 MATCHED-SPIN OVERPERFORMERS (spin_sim>=0.4), by residual over shape ===\n")
print(head(ag[spin_sim>=0.4][order(-resid),
  .(player_name,pitch_type,n,miss,exp,resid,axis_diff,as_gap,spin_sim)], 15))
fwrite(ag[order(-resid)], file.path(MDIR, "matched_spin_overperformance.csv"))
cat("\nsaved matched_spin_overperformance.csv\n")
