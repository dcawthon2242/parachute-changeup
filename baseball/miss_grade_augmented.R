#!/usr/bin/env Rscript

# Step 4 + 5: test whether the two hypothesis features add signal beyond the
# TJStuff+ shape set, then produce 2026 pitch grades.
#
#   Breaking balls: + path_ratio (mean path_to_location_ratio off the fastball)
#   Offspeed:       + axis_diff + spin_eff_diff (spin similarity to the fastball)
#
# Marginal value is measured at the group level (where the features live) via
# 2026 holdout RMSE/R2 delta + permutation importance. The final grade model uses
# each added feature only for the group where it helped on holdout.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(1)
MDIR <- file.path("data","statcast_model")
d  <- readRDS(file.path(MDIR, "miss_grade_data.rds"))
ac <- readRDS(file.path(MDIR, "arch_compare.rds"))
FEAT <- ac$FEAT
cat(sprintf("Architecture winner from step 3: %s\n", ac$winner))

train <- d[set=="train"]; hold <- d[set=="holdout"]

train_lgb <- function(dtr, feats, nrounds=2000){
  dtr <- dtr[stats::complete.cases(dtr[, ..feats]) & is.finite(miss_distance)]
  n <- nrow(dtr); vi <- sample(n, floor(0.15*n))
  dtrain <- lgb.Dataset(as.matrix(dtr[-vi, ..feats]), label=dtr$miss_distance[-vi])
  dval   <- lgb.Dataset.create.valid(dtrain, as.matrix(dtr[vi, ..feats]), label=dtr$miss_distance[vi])
  lgb.train(params=list(objective="regression", metric="rmse", learning_rate=0.05,
    num_leaves=31, min_data_in_leaf=200, feature_fraction=0.8, bagging_fraction=0.8, bagging_freq=1),
    data=dtrain, nrounds=nrounds, valids=list(val=dval), early_stopping_rounds=60, verbose=-1)
}
rmse <- function(p,a){ok<-is.finite(p)&is.finite(a); sqrt(mean((p[ok]-a[ok])^2))}
r2   <- function(p,a){ok<-is.finite(p)&is.finite(a); p<-p[ok];a<-a[ok];1-sum((a-p)^2)/sum((a-mean(a))^2)}

test_add <- function(gg, extra){
  base <- FEAT; aug <- c(FEAT, extra)
  hb <- hold[grp==gg]
  cov <- mean(stats::complete.cases(hb[, ..extra]))
  mb <- train_lgb(train[grp==gg], base)
  ma <- train_lgb(train[grp==gg], aug)
  pb <- predict(mb, as.matrix(hb[, ..base]))
  pa <- predict(ma, as.matrix(hb[, ..aug]))
  # permutation importance of the added feature(s) on holdout
  hp <- copy(hb); for (e in extra) hp[[e]] <- sample(hp[[e]])
  pperm <- predict(ma, as.matrix(hp[, ..aug]))
  imp <- lgb.importance(ma)
  cat(sprintf("\n--- %s : + %s  (holdout n=%d, feature coverage=%.0f%%) ---\n",
              gg, paste(extra,collapse="+"), nrow(hb), 100*cov))
  cat(sprintf("  RMSE  base=%.4f  aug=%.4f  delta=%+.4f\n", rmse(pb,hb$miss_distance),
              rmse(pa,hb$miss_distance), rmse(pa,hb$miss_distance)-rmse(pb,hb$miss_distance)))
  cat(sprintf("  R2    base=%.4f  aug=%.4f  delta=%+.4f\n", r2(pb,hb$miss_distance),
              r2(pa,hb$miss_distance), r2(pa,hb$miss_distance)-r2(pb,hb$miss_distance)))
  cat(sprintf("  permutation-shuffle RMSE rise from added feature: %+.4f\n",
              rmse(pperm,hb$miss_distance)-rmse(pa,hb$miss_distance)))
  cat("  added-feature gain rank (of ", nrow(imp), " feats): ",
      paste(sprintf("%s=#%d(%.1f%%)", extra,
        match(extra, imp$Feature), 100*imp$Gain[match(extra, imp$Feature)]), collapse="  "), "\n", sep="")
  list(help = rmse(pa,hb$miss_distance) < rmse(pb,hb$miss_distance), aug=aug)
}

cat("\n=== STEP 4: marginal value of the hypothesis features (2026 holdout) ===\n")
brk <- test_add("breaking", "path_ratio")
off <- test_add("offspeed", c("axis_diff","spin_eff_diff"))

# ---- STEP 5: final grade model (per-group, add features where they helped) ----
cat("\n=== STEP 5: final grades ===\n")
feat_fb  <- FEAT
feat_brk <- if (brk$help) brk$aug else FEAT
feat_off <- if (off$help) off$aug else FEAT
cat(sprintf("  fastball feats: base | breaking: %s | offspeed: %s\n",
  if(brk$help) "base + path_ratio" else "base",
  if(off$help) "base + spin-similarity" else "base"))

d[, pred_miss := NA_real_]
for (cfg in list(list(g="fastball",f=feat_fb), list(g="breaking",f=feat_brk), list(g="offspeed",f=feat_off))) {
  m <- train_lgb(train[grp==cfg$g], cfg$f)
  d[grp==cfg$g, pred_miss := predict(m, as.matrix(.SD[, cfg$f, with=FALSE]))]
}
d[is.na(pred_miss) & grp=="other", pred_miss := mean(train$miss_distance)]

# scale per-pitch predictions to 100 +/- 10 over 2026 (higher pred miss = better)
hp <- d[set=="holdout" & is.finite(pred_miss)]
mu <- mean(hp$pred_miss); sdv <- sd(hp$pred_miss)
hp[, grade100 := 100 + 10*(pred_miss-mu)/sdv]

grades <- hp[, .(n=.N, act_miss=round(mean(miss_distance),2),
  pred_miss=round(mean(pred_miss),2), grade=round(mean(grade100),1)),
  by=.(pitcher, player_name, pitch_type, grp)][n>=25]
# 20-80 scouting grade within pitch type
grades[, g2080 := {z <- (grade-mean(grade))/sd(grade); pmax(20, pmin(80, round(50+10*z)))}, by=pitch_type]
setorder(grades, -grade)
fwrite(grades, file.path(MDIR, "miss_grade_2026.csv"))

cat(sprintf("\nWrote %s (%d pitcher-pitchtype grades, min 25 swings)\n",
            file.path(MDIR,"miss_grade_2026.csv"), nrow(grades)))
cat("\nTop 15 miss-distance pitch grades (2026):\n")
print(head(grades[, .(player_name, pitch_type, n, act_miss, pred_miss, grade, g2080)], 15))
cat("\nTop breaking balls:\n")
print(head(grades[grp=="breaking", .(player_name, pitch_type, n, act_miss, grade, g2080)], 8))
cat("\nTop offspeed:\n")
print(head(grades[grp=="offspeed", .(player_name, pitch_type, n, act_miss, grade, g2080)], 8))

cat("\n=== VERDICT ===\n")
cat(sprintf("  Architecture: %s\n", ac$winner))
cat(sprintf("  path_to_location_ratio (breaking): %s\n",
    if(brk$help) "IMPROVES holdout -> kept" else "no holdout improvement -> dropped"))
cat(sprintf("  spin similarity (offspeed): %s\n",
    if(off$help) "IMPROVES holdout -> kept" else "no holdout improvement -> dropped"))
