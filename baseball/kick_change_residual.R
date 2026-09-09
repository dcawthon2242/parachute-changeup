#!/usr/bin/env Rscript

# Mirror test: do KICK CHANGES (very low spin rate + spin very DIFFERENT from the FB)
# UNDERPERFORM their expectation -- i.e., are they easy for hitters to pick up?
# Residual = actual whiff - expected whiff, over (a) full-shape and (b) movement-only
# expectations, 2026 offspeed. Look at the low-spin / high-dissimilarity end.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(1)
MDIR <- file.path("data","statcast_model")
FEAT <- readRDS(file.path(MDIR, "arch_compare.rds"))$FEAT
FEAT_MOVE <- setdiff(FEAT, c("release_spin_rate","sax","cax"))
as_long <- readRDS(file.path(MDIR, "active_spin_long.rds"))
d <- readRDS(file.path(MDIR, "miss_grade_data_activespin.rds"))

# FB anchors: active spin + raw spin rate (FF>SI>FC)
fb <- as_long[pitch_type %in% c("FF","SI","FC")][, pr:=match(pitch_type,c("FF","SI","FC"))][
  order(pitcher,season,pr)][, .SD[1], by=.(pitcher,season)][, .(pitcher,season,fb_active=active_spin)]
d <- merge(d, fb, by=c("pitcher","season"), all.x=TRUE)
fbsp <- d[pitch_type %in% c("FF","SI","FC")][, pr:=match(pitch_type,c("FF","SI","FC"))][
  order(pitcher,season,pr), .(fb_spin=mean(release_spin_rate,na.rm=TRUE)), by=.(pitcher,season)]
d <- merge(d, fbsp, by=c("pitcher","season"), all.x=TRUE)
d[, as_gap := active_spin - fb_active]; d[, abs_as := abs(as_gap)]
d[, spin_drop := fb_spin - release_spin_rate]      # + = changeup spins slower than FB (kick = large)

# IVB for 2026 (to confirm kick = low vert)
raw26 <- fread("data/statcast_2026/statcast_2026_all.csv",
               select=c("game_pk","at_bat_number","pitch_number","pfx_z"))
raw26[, ivb := pfx_z*12]
d <- merge(d, raw26[, .(game_pk,at_bat_number,pitch_number,ivb)],
           by=c("game_pk","at_bat_number","pitch_number"), all.x=TRUE)

tr <- d[set=="train" & is.finite(miss_distance)]
mk <- function(feats){ n<-nrow(tr); vi<-sample(n,floor(0.15*n))
  dtrain<-lgb.Dataset(as.matrix(tr[-vi,..feats]),label=tr$is_whiff[-vi])
  dval<-lgb.Dataset.create.valid(dtrain,as.matrix(tr[vi,..feats]),label=tr$is_whiff[vi])
  lgb.train(params=list(objective="binary",metric="binary_logloss",learning_rate=0.05,num_leaves=31,
    min_data_in_leaf=150,feature_fraction=0.8,bagging_fraction=0.8,bagging_freq=1),
    data=dtrain,nrounds=2000,valids=list(val=dval),early_stopping_rounds=60,verbose=-1) }
mw_full <- mk(FEAT); mw_move <- mk(FEAT_MOVE)

ho <- d[set=="holdout" & pitch_type %in% c("CH","FS") & !is.na(is_whiff)]
ho[, exp_full := predict(mw_full, as.matrix(.SD[, FEAT, with=FALSE]))]
ho[, exp_move := predict(mw_move, as.matrix(.SD[, FEAT_MOVE, with=FALSE]))]
ho[, wres_full := is_whiff - exp_full][, wres_move := is_whiff - exp_move]

CO <- function(a,b) cor(a,b,use="complete.obs")
cat("=== do LOWER-SPIN / more-different changeups underperform? (2026 CH+FS, pitch-level) ===\n")
cat("    (hypothesis: easy to pick up -> NEGATIVE residual)\n")
cat(sprintf("  wres_full vs release_spin_rate (abs) : r=%+.3f  (want +: low spin -> low resid)\n", CO(ho$wres_full, ho$release_spin_rate)))
cat(sprintf("  wres_full vs spin_drop vs FB         : r=%+.3f  (want -: big drop -> low resid)\n", CO(ho$wres_full, ho$spin_drop)))
cat(sprintf("  wres_full vs |active-spin gap|       : r=%+.3f  (want -)\n", CO(ho$wres_full, ho$abs_as)))
cat(sprintf("  wres_full vs axis gap vs FB          : r=%+.3f  (want -)\n", CO(ho$wres_full, ho$axis_diff)))

cat("\n--- whiff residual by RAW spin-rate bin (CH+FS) ---\n")
ho[, spin_bin := cut(release_spin_rate, c(0,1200,1500,1800,2100,2400,9000),
   labels=c("<1200 (kick)","1200-1500","1500-1800","1800-2100","2100-2400",">2400"))]
print(ho[!is.na(spin_bin), .(pitches=.N, whiff=round(mean(is_whiff),3), exp=round(mean(exp_full),3),
  wres_full=round(mean(wres_full),3), wres_move=round(mean(wres_move),3),
  mean_ivb=round(mean(ivb,na.rm=TRUE),1), mean_as_gap=round(mean(as_gap,na.rm=TRUE),3)), by=spin_bin][order(spin_bin)])

# explicit kick-change cohort: low spin AND low IVB (the classic definition)
ho[, kick := release_spin_rate<=1500 & ivb<=4]
cat("\n--- KICK-CHANGE cohort (spin<=1500 & IVB<=4) vs other changeups ---\n")
print(ho[!is.na(kick), .(pitches=.N, whiff=round(mean(is_whiff),3), exp=round(mean(exp_full),3),
  wres_full=round(mean(wres_full),3), wres_move=round(mean(wres_move),3)), by=kick][order(kick)])

ag <- ho[, .(n=.N, whiff=round(mean(is_whiff),3), exp=round(mean(exp_full),3),
  wres=round(mean(wres_full),3), spin=round(mean(release_spin_rate),0),
  spin_drop=round(mean(spin_drop),0), ivb=round(mean(ivb,na.rm=TRUE),1),
  as_gap=round(mean(as_gap,na.rm=TRUE),3)), by=.(player_name,pitch_type)][n>=40]
cat("\n=== known kick-change guys (Munoz / Birdsong / others low-spin) ===\n")
print(ag[grepl("Munoz|Muñoz|Birdsong",player_name)])
cat("\n=== lowest-spin changeups in 2026 (n>=40): do they underperform? ===\n")
print(head(ag[order(spin), .(player_name,pitch_type,n,spin,spin_drop,ivb,as_gap,whiff,exp,wres)],15))
fwrite(ag[order(spin)], file.path(MDIR,"kick_change_residual.csv"))
cat("\nsaved kick_change_residual.csv\n")
