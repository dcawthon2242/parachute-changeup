#!/usr/bin/env Rscript

# FULL GRID:
#   Approaches  A = all-types, B = per-pitch-type, C = grouped (FB/BRK/OFF)
#   Config      base (TJStuff+ shape only)  vs  aug (+ tunneling path_ratio for
#               breaking, + spin similarity [axis_diff + MEASURED active-spin gap] for offspeed)
#   Count set   ALL competitive swings  vs  2K (0-2/1-2/2-2 put-away)
# Then: for the winning approach x set, top-10 breaking & offspeed RISERS
#   (gain = predicted miss WITH new features - WITHOUT).

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(1)
MDIR <- file.path("data","statcast_model")
FEAT <- readRDS(file.path(MDIR, "arch_compare.rds"))$FEAT
as_long <- readRDS(file.path(MDIR, "active_spin_long.rds"))
d <- readRDS(file.path(MDIR, "miss_grade_data_activespin.rds"))   # has per-pitch active_spin

# measured spin-similarity vs primary FB (FF>SI>FC), season-level
fb <- as_long[pitch_type %in% c("FF","SI","FC")]
fb[, pr := match(pitch_type, c("FF","SI","FC"))]
fb <- fb[order(pitcher,season,pr)][, .SD[1], by=.(pitcher,season)][, .(pitcher,season,fb_active=active_spin)]
d <- merge(d, fb, by=c("pitcher","season"), all.x=TRUE)
d[, meas_eff_diff := active_spin - fb_active]

BREAKING <- c("SL","ST","CU","KC","SV","CS"); OFFSPEED <- c("CH","FS","FO")
# GATED features for the all-types model: tunneling only lives on breaking balls,
# spin-similarity only on offspeed. NA elsewhere (LightGBM treats NA as missing).
d[, path_ratio_g := ifelse(grp=="breaking", path_ratio, NA_real_)]
d[, axis_diff_g  := ifelse(grp=="offspeed", axis_diff,   NA_real_)]
d[, eff_diff_g   := ifelse(grp=="offspeed", meas_eff_diff, NA_real_)]
BRK_X <- "path_ratio"; OFF_X <- c("axis_diff","meas_eff_diff"); A_X <- c("path_ratio_g","axis_diff_g","eff_diff_g")

train_lgb <- function(dtr, feats, nr=1500){
  dtr <- dtr[is.finite(miss_distance)]; n <- nrow(dtr)
  vi <- sample(n, max(1, floor(0.15*n)))
  dtrain <- lgb.Dataset(as.matrix(dtr[-vi, ..feats]), label=dtr$miss_distance[-vi])
  dval   <- lgb.Dataset.create.valid(dtrain, as.matrix(dtr[vi, ..feats]), label=dtr$miss_distance[vi])
  lgb.train(params=list(objective="regression",metric="rmse",learning_rate=0.05,num_leaves=31,
    min_data_in_leaf=150,feature_fraction=0.8,bagging_fraction=0.8,bagging_freq=1),
    data=dtrain, nrounds=nr, valids=list(val=dval), early_stopping_rounds=50, verbose=-1)
}
rmse <- function(p,a){ok<-is.finite(p)&is.finite(a); sqrt(mean((p[ok]-a[ok])^2))}
r2   <- function(p,a){ok<-is.finite(p)&is.finite(a); 1-sum((a[ok]-p[ok])^2)/sum((a[ok]-mean(a[ok]))^2)}

# ---- approach predictors: return base & aug holdout prediction vectors ----
pred_A <- function(tr, ho){
  mb <- train_lgb(tr, FEAT); ma <- train_lgb(tr, c(FEAT, A_X))
  list(base=predict(mb, as.matrix(ho[, FEAT, with=FALSE])),
       aug =predict(ma, as.matrix(ho[, c(FEAT,A_X), with=FALSE])))
}
pred_C <- function(tr, ho){
  pb <- pa <- rep(NA_real_, nrow(ho))
  for (g in c("fastball","breaking","offspeed")){
    trg <- tr[grp==g]; idx <- which(ho$grp==g); if(!length(idx)) next
    ax <- if(g=="breaking") BRK_X else if(g=="offspeed") OFF_X else character(0)
    mb <- train_lgb(trg, FEAT); pb[idx] <- predict(mb, as.matrix(ho[idx, FEAT, with=FALSE]))
    fa <- c(FEAT, ax); ma <- train_lgb(trg, fa); pa[idx] <- predict(ma, as.matrix(ho[idx, fa, with=FALSE]))
  }
  list(base=pb, aug=pa)
}
pred_B <- function(tr, ho){
  types <- names(which(table(tr$pitch_type) >= 1000))
  pb <- pa <- rep(NA_real_, nrow(ho))
  for (pt in types){
    trp <- tr[pitch_type==pt]; idx <- which(ho$pitch_type==pt); if(!length(idx)) next
    ax <- if(pt %in% BREAKING) BRK_X else if(pt %in% OFFSPEED) OFF_X else character(0)
    mb <- train_lgb(trp, FEAT); pb[idx] <- predict(mb, as.matrix(ho[idx, FEAT, with=FALSE]))
    fa <- c(FEAT, ax); ma <- train_lgb(trp, fa); pa[idx] <- predict(ma, as.matrix(ho[idx, fa, with=FALSE]))
  }
  left <- tr[!pitch_type %in% types]; lidx <- which(!ho$pitch_type %in% types)
  if(length(lidx) && nrow(left)>50){ m <- train_lgb(left, FEAT)
    pb[lidx] <- pa[lidx] <- predict(m, as.matrix(ho[lidx, FEAT, with=FALSE])) }
  list(base=pb, aug=pa)
}

run_set <- function(dd, label, min_n){
  tr <- dd[set=="train"]; ho <- dd[set=="holdout" & is.finite(miss_distance)]
  cat(sprintf("\n########## COUNT SET: %s  (train=%d holdout=%d, holdout mean miss=%.2f) ##########\n",
      label, nrow(tr[is.finite(miss_distance)]), nrow(ho), mean(ho$miss_distance)))
  res <- list(A=pred_A(tr,ho), B=pred_B(tr,ho), C=pred_C(tr,ho))
  tab <- rbindlist(lapply(names(res), function(k){
    ok <- is.finite(res[[k]]$base) & is.finite(res[[k]]$aug)
    data.table(approach=k,
      base_rmse=rmse(res[[k]]$base[ok], ho$miss_distance[ok]),
      aug_rmse =rmse(res[[k]]$aug[ok],  ho$miss_distance[ok]),
      base_r2  =r2(res[[k]]$base[ok],   ho$miss_distance[ok]),
      aug_r2   =r2(res[[k]]$aug[ok],    ho$miss_distance[ok]))
  }))
  tab[, delta_rmse := aug_rmse - base_rmse]
  tab[, `:=`(base_rmse=round(base_rmse,4), aug_rmse=round(aug_rmse,4),
             delta_rmse=round(delta_rmse,4), base_r2=round(base_r2,4), aug_r2=round(aug_r2,4))]
  lab <- c(A="A all-types", B="B per-pitch-type", C="C grouped FB/BRK/OFF")
  tab[, approach := lab[approach]]
  cat("\n--- RMSE grid (with vs without new features) ---\n"); print(tab)

  # winner = lowest aug_rmse
  win <- names(res)[which.min(sapply(res, function(x){ ok<-is.finite(x$base)&is.finite(x$aug); rmse(x$aug[ok], ho$miss_distance[ok]) }))]
  cat(sprintf("\n>>> winning combination for %s: approach %s (augmented)\n", label, win))
  ho[, `:=`(pb=res[[win]]$base, pa=res[[win]]$aug)]
  ho[, gain := pa - pb]
  agg <- ho[is.finite(gain), .(n=.N, act=round(mean(miss_distance),2),
      base=round(mean(pb),2), aug=round(mean(pa),2), gain=round(mean(gain),3),
      path_ratio=round(mean(path_ratio,na.rm=TRUE),2), axis_diff=round(mean(axis_diff,na.rm=TRUE),1),
      as_gap=round(mean(meas_eff_diff,na.rm=TRUE),3)),
      by=.(player_name, pitch_type, grp)][n>=min_n]
  cat(sprintf("\n=== TOP 10 BREAKING RISERS (%s, approach %s) ===\n", label, win))
  print(head(agg[grp=="breaking"][order(-gain), .(player_name,pitch_type,n,act,base,aug,gain,path_ratio)],10))
  cat(sprintf("\n=== TOP 10 OFFSPEED RISERS (%s, approach %s) ===\n", label, win))
  print(head(agg[grp=="offspeed"][order(-gain), .(player_name,pitch_type,n,act,base,aug,gain,axis_diff,as_gap)],10))
  fwrite(agg[order(-gain)], file.path(MDIR, sprintf("risers_%s.csv", label)))
  tab[, set := label]; tab
}

grid <- rbindlist(list(
  run_set(d[is.finite(miss_distance)], "ALL", 40),
  run_set(d[strikes==2 & balls<=2 & is.finite(miss_distance)], "2K", 25)
))
cat("\n\n================ RMSE GRID SUMMARY ================\n"); print(grid[, .(set,approach,base_rmse,aug_rmse,delta_rmse,base_r2,aug_r2)])
fwrite(grid, file.path(MDIR, "grid_rmse_summary.csv"))
cat("\nsaved grid_rmse_summary.csv, risers_ALL.csv, risers_2K.csv\n")
