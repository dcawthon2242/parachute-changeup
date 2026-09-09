#!/usr/bin/env Rscript

# Add LOCATION to the breaking-ball miss-distance model and see what happens to the
# tunneling metric (path_ratio). Four configs on identical rows:
#   base = TJStuff+ shape | +tunnel | +location | +both
# Report holdout RMSE, the marginal value of path_ratio WITHOUT vs WITH location,
# and permutation importance. Target = miss_distance (E[miss|swing]).

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(1)
MDIR <- file.path("data","statcast_model")
FEAT <- readRDS(file.path(MDIR, "arch_compare.rds"))$FEAT
d <- readRDS(file.path(MDIR, "miss_grade_data_activespin.rds"))

kin <- rbindlist(lapply(2023:2026, function(y) fread(sprintf("data/statcast_%d/statcast_%d_all.csv",y,y),
  select=c("game_pk","at_bat_number","pitch_number","plate_x","plate_z","sz_bot","vx0","vy0","vz0","ax","ay","az"))))
kin <- unique(kin, by=c("game_pk","at_bat_number","pitch_number"))
yf<-17/12; kin[,vyf:=-sqrt(pmax(vy0^2-2*ay*(50-yf),0))]; kin[,tf:=(vyf-vy0)/ay]
kin[,vzf:=vz0+az*tf][,vxf:=vx0+ax*tf]
kin[, VAA := atan2(vzf,abs(vyf))*180/pi][, HAA := atan2(vxf,abs(vyf))*180/pi][, z_rel_bot := plate_z - sz_bot]
d <- merge(d, kin[, .(game_pk,at_bat_number,pitch_number,plate_x,plate_z,z_rel_bot,VAA,HAA)],
           by=c("game_pk","at_bat_number","pitch_number"), all.x=TRUE)

LOC <- c("plate_x","plate_z","z_rel_bot","VAA","HAA"); TUN <- "path_ratio"
brk <- d[grp=="breaking"]
need <- c(FEAT, LOC, TUN)
brk <- brk[stats::complete.cases(brk[, ..need]) & is.finite(miss_distance)]

train_lgb <- function(dtr, feats, nr=2000){ n<-nrow(dtr); vi<-sample(n,floor(0.15*n))
  dtrain<-lgb.Dataset(as.matrix(dtr[-vi,..feats]),label=dtr$miss_distance[-vi])
  dval<-lgb.Dataset.create.valid(dtrain,as.matrix(dtr[vi,..feats]),label=dtr$miss_distance[vi])
  lgb.train(params=list(objective="regression",metric="rmse",learning_rate=0.05,num_leaves=31,
    min_data_in_leaf=150,feature_fraction=0.8,bagging_fraction=0.8,bagging_freq=1),
    data=dtrain,nrounds=nr,valids=list(val=dval),early_stopping_rounds=60,verbose=-1) }
rmse<-function(p,a) sqrt(mean((p-a)^2))
perm_imp <- function(m, ho, feats, target){ base<-rmse(predict(m,as.matrix(ho[,..feats])), ho[[target]])
  sapply(feats, function(f){ h<-copy(ho); h[[f]]<-sample(h[[f]]); rmse(predict(m,as.matrix(h[,..feats])), h[[target]])-base }) }

run <- function(sub, label){
  tr<-sub[set=="train"]; ho<-sub[set=="holdout"]
  cat(sprintf("\n########## BREAKING: %s (train=%d holdout=%d) ##########\n", label, nrow(tr), nrow(ho)))
  cfg <- list(base=FEAT, `+tunnel`=c(FEAT,TUN), `+location`=c(FEAT,LOC), `+both`=c(FEAT,LOC,TUN))
  R <- sapply(cfg, function(fe){ m<-train_lgb(tr,fe); rmse(predict(m,as.matrix(ho[,..fe])), ho$miss_distance) })
  tab <- data.table(config=names(R), rmse=round(R,4))
  print(tab)
  cat(sprintf("  tunneling marginal WITHOUT location: %+.4f  (base -> +tunnel)\n", R["+tunnel"]-R["base"]))
  cat(sprintf("  tunneling marginal WITH location   : %+.4f  (+location -> +both)\n", R["+both"]-R["+location"]))
  cat(sprintf("  location marginal (base -> +location): %+.4f\n", R["+location"]-R["base"]))
  # permutation importance in the full (+both) model
  mfull <- train_lgb(tr, c(FEAT,LOC,TUN))
  imp <- perm_imp(mfull, ho, c(FEAT,LOC,TUN), "miss_distance")
  cat("  permutation importance (RMSE rise when shuffled), full model, top of interest:\n")
  pr <- sort(imp[c(TUN,LOC)], decreasing=TRUE)
  for(nm in names(pr)) cat(sprintf("      %-12s %+.4f\n", nm, pr[nm]))
  cat(sprintf("  corr(path_ratio, plate_z)=%.3f  corr(path_ratio, z_rel_bot)=%.3f  corr(path_ratio, VAA)=%.3f\n",
      cor(ho$path_ratio,ho$plate_z), cor(ho$path_ratio,ho$z_rel_bot), cor(ho$path_ratio,ho$VAA)))
}

run(brk, "ALL counts")
run(brk[strikes==2 & balls<=2], "2-strike (put-away)")
