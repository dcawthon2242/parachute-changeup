#!/usr/bin/env Rscript

# THE proper test of "matched spin overperforms": build the whiff/miss expectation
# from MOVEMENT + velo + release ONLY (no spin rate, no spin axis). If pitches whose
# spin LOOKS like the fastball beat that movement-only expectation, the spin *look*
# adds deception beyond the physical movement -> proves the hypothesis.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(1)
MDIR <- file.path("data","statcast_model")
FEAT <- readRDS(file.path(MDIR, "arch_compare.rds"))$FEAT
FEAT_MOVE <- setdiff(FEAT, c("release_spin_rate","sax","cax"))   # drop all spin-appearance
cat("movement-only features: ", paste(FEAT_MOVE, collapse=", "), "\n")

as_long <- readRDS(file.path(MDIR, "active_spin_long.rds"))
d <- readRDS(file.path(MDIR, "miss_grade_data_activespin.rds"))
fb <- as_long[pitch_type %in% c("FF","SI","FC")][, pr:=match(pitch_type,c("FF","SI","FC"))][
  order(pitcher,season,pr)][, .SD[1], by=.(pitcher,season)][, .(pitcher,season,fb_active=active_spin)]
d <- merge(d, fb, by=c("pitcher","season"), all.x=TRUE)
d[, as_gap := active_spin - fb_active]; d[, abs_as := abs(as_gap)]
d[, spin_sim := exp(-(abs_as/0.10)^2) * exp(-(axis_diff/45)^2)]

tr <- d[set=="train" & is.finite(miss_distance)]
mk <- function(obj,y,feats){ n<-nrow(tr); vi<-sample(n,floor(0.15*n))
  dtrain<-lgb.Dataset(as.matrix(tr[-vi,..feats]),label=tr[[y]][-vi])
  dval<-lgb.Dataset.create.valid(dtrain,as.matrix(tr[vi,..feats]),label=tr[[y]][vi])
  lgb.train(params=list(objective=obj,metric=if(obj=="binary")"binary_logloss" else "rmse",
    learning_rate=0.05,num_leaves=31,min_data_in_leaf=150,feature_fraction=0.8,
    bagging_fraction=0.8,bagging_freq=1),data=dtrain,nrounds=2000,valids=list(val=dval),
    early_stopping_rounds=60,verbose=-1) }

mw_move <- mk("binary","is_whiff", FEAT_MOVE)
mw_full <- mk("binary","is_whiff", FEAT)
ho <- d[set=="holdout" & grp=="offspeed" & !is.na(is_whiff)]
ho[, exp_move := predict(mw_move, as.matrix(.SD[, FEAT_MOVE, with=FALSE]))]
ho[, exp_full := predict(mw_full, as.matrix(.SD[, FEAT, with=FALSE]))]
ho[, wres_move := is_whiff - exp_move]   # whiffs above MOVEMENT-only expectation
ho[, wres_full := is_whiff - exp_full]   # whiffs above FULL-shape expectation (incl spin)

ag <- ho[, .(n=.N, whiff=round(mean(is_whiff),3),
  wres_move=round(mean(wres_move),3), wres_full=round(mean(wres_full),3),
  as_gap=round(mean(as_gap,na.rm=TRUE),3), axis=round(mean(axis_diff,na.rm=TRUE),1),
  spin_sim=round(mean(spin_sim,na.rm=TRUE),3)), by=.(player_name,pitch_type)][n>=40]
CO <- function(a,b) cor(a,b,use="complete.obs")

cat("\n=== does SPIN SIMILARITY predict whiffs ABOVE a movement-only expectation? ===\n")
cat(sprintf("  (if hypothesis true: spin_sim should predict wres_move but NOT wres_full)\n"))
cat(sprintf("  wres_MOVE vs spin_sim          : r=%+.3f\n", CO(ag$wres_move, ag$spin_sim)))
cat(sprintf("  wres_MOVE vs |active-spin gap| : r=%+.3f\n", CO(ag$wres_move, abs(ag$as_gap))))
cat(sprintf("  wres_FULL vs spin_sim          : r=%+.3f   (spin already priced in)\n", CO(ag$wres_full, ag$spin_sim)))
cat("\n  whiffs-above-MOVEMENT by active-spin similarity tier:\n")
print(ag[!is.na(as_gap), .(n=.N, wres_move=round(mean(wres_move),3), wres_full=round(mean(wres_full),3)),
  by=.(tier=cut(abs(as_gap),c(-.01,.03,.07,.15,1),labels=c("<=.03 match","0.03-.07","0.07-.15",">.15 different")))][order(tier)])

cat("\n=== archetypes: whiffs above MOVEMENT-only vs FULL-shape ===\n")
print(ag[grepl("Cease|Ragans|Vesia|Chivilli",player_name)][order(-wres_move),
  .(player_name,pitch_type,n,whiff,wres_move,wres_full,as_gap,axis,spin_sim)])
fwrite(ag[order(-wres_move)], file.path(MDIR,"movement_only_spin_test.csv"))
cat("\nsaved movement_only_spin_test.csv\n")
