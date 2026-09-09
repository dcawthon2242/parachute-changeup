#!/usr/bin/env Rscript

# Kick-change underperformance test across 2023-2026 with OUT-OF-FOLD expectations
# (honest residuals for every season). Does very-low-spin / spin-different offspeed
# whiff LESS than expected? Pull named kick-change guys (Munoz, Birdsong).

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(1)
MDIR <- file.path("data","statcast_model")
FEAT <- readRDS(file.path(MDIR, "arch_compare.rds"))$FEAT
FEAT_MOVE <- setdiff(FEAT, c("release_spin_rate","sax","cax"))
as_long <- readRDS(file.path(MDIR, "active_spin_long.rds"))
d <- readRDS(file.path(MDIR, "miss_grade_data_activespin.rds"))

fb <- as_long[pitch_type %in% c("FF","SI","FC")][, pr:=match(pitch_type,c("FF","SI","FC"))][
  order(pitcher,season,pr)][, .SD[1], by=.(pitcher,season)][, .(pitcher,season,fb_active=active_spin)]
d <- merge(d, fb, by=c("pitcher","season"), all.x=TRUE)
fbsp <- d[pitch_type %in% c("FF","SI","FC")][, pr:=match(pitch_type,c("FF","SI","FC"))][
  order(pitcher,season,pr), .(fb_spin=mean(release_spin_rate,na.rm=TRUE)), by=.(pitcher,season)]
d <- merge(d, fbsp, by=c("pitcher","season"), all.x=TRUE)
d[, as_gap := active_spin - fb_active]; d[, abs_as := abs(as_gap)]
d[, spin_drop := fb_spin - release_spin_rate]

# IVB (pfx_z) for every season
ivb_all <- rbindlist(lapply(2023:2026, function(y){
  f <- sprintf("data/statcast_%d/statcast_%d_all.csv", y, y)
  x <- fread(f, select=c("game_pk","at_bat_number","pitch_number","pfx_z"))
  x[, ivb := pfx_z*12][, .(game_pk,at_bat_number,pitch_number,ivb)] }))
d <- merge(d, unique(ivb_all, by=c("game_pk","at_bat_number","pitch_number")),
           by=c("game_pk","at_bat_number","pitch_number"), all.x=TRUE)

# ---- out-of-fold whiff expectations (train all pitch types, predict offspeed) ----
D <- d[!is.na(is_whiff) & stats::complete.cases(d[, ..FEAT])]
K <- 4; D[, fold := sample(rep(1:K, length.out=.N))]
lgb_fit <- function(dtr, feats){ n<-nrow(dtr); vi<-sample(n,floor(0.12*n))
  dtrain<-lgb.Dataset(as.matrix(dtr[-vi,..feats]),label=dtr$is_whiff[-vi])
  dval<-lgb.Dataset.create.valid(dtrain,as.matrix(dtr[vi,..feats]),label=dtr$is_whiff[vi])
  lgb.train(params=list(objective="binary",metric="binary_logloss",learning_rate=0.06,num_leaves=31,
    min_data_in_leaf=200,feature_fraction=0.8,bagging_fraction=0.8,bagging_freq=1),
    data=dtrain,nrounds=1200,valids=list(val=dval),early_stopping_rounds=50,verbose=-1) }
D[, `:=`(exp_full=NA_real_, exp_move=NA_real_)]
for (f in 1:K){ cat("fold",f,"..\n")
  mf <- lgb_fit(D[fold!=f], FEAT);      D[fold==f, exp_full := predict(mf, as.matrix(.SD[, FEAT, with=FALSE]))]
  mm <- lgb_fit(D[fold!=f], FEAT_MOVE); D[fold==f, exp_move := predict(mm, as.matrix(.SD[, FEAT_MOVE, with=FALSE]))]
}
D[, wres_full := is_whiff - exp_full][, wres_move := is_whiff - exp_move]
off <- D[pitch_type %in% c("CH","FS")]
CO <- function(a,b) cor(a,b,use="complete.obs")

cat(sprintf("\n=== 2023-2026 CH+FS, OOF residuals (n=%d) ===\n", nrow(off)))
cat(sprintf("  wres_full vs spin rate (abs)  : r=%+.3f (want +)\n", CO(off$wres_full, off$release_spin_rate)))
cat(sprintf("  wres_full vs spin_drop vs FB  : r=%+.3f (want -)\n", CO(off$wres_full, off$spin_drop)))
cat(sprintf("  wres_full vs |active-spin gap|: r=%+.3f (want -)\n", CO(off$wres_full, off$abs_as)))
cat(sprintf("  wres_full vs axis gap         : r=%+.3f (want -)\n", CO(off$wres_full, off$axis_diff)))

off[, spin_bin := cut(release_spin_rate, c(0,1200,1500,1800,2100,2400,9000),
   labels=c("<1200 (kick)","1200-1500","1500-1800","1800-2100","2100-2400",">2400"))]
cat("\n--- whiff residual by RAW spin-rate bin (CH+FS, all yrs) ---\n")
print(off[!is.na(spin_bin), .(pitches=.N, whiff=round(mean(is_whiff),3), exp=round(mean(exp_full),3),
  wres_full=round(mean(wres_full),3), wres_move=round(mean(wres_move),3),
  ivb=round(mean(ivb,na.rm=TRUE),1), as_gap=round(mean(as_gap,na.rm=TRUE),3)), by=spin_bin][order(spin_bin)])

off[, kick := release_spin_rate<=1500 & ivb<=4]
cat("\n--- KICK cohort (spin<=1500 & IVB<=4) vs other, by pitch type ---\n")
print(off[!is.na(kick), .(pitches=.N, whiff=round(mean(is_whiff),3), exp=round(mean(exp_full),3),
  wres_full=round(mean(wres_full),3), wres_move=round(mean(wres_move),3)), by=.(pitch_type,kick)][order(pitch_type,kick)])

cat("\n=== NAMED kick-change guys (all their offspeed, pooled 2023-2026, n>=25) ===\n")
nk <- D[grepl("Mu.oz, Andr|Birdsong|Bird song", player_name) &
        pitch_type %in% c("CH","FS","FO","SL")]
print(nk[, .(n=.N, spin=round(mean(release_spin_rate),0), spin_drop=round(mean(spin_drop),0),
  ivb=round(mean(ivb,na.rm=TRUE),1), as_gap=round(mean(as_gap,na.rm=TRUE),3),
  whiff=round(mean(is_whiff),3), exp=round(mean(exp_full),3),
  wres_full=round(mean(wres_full),3), wres_move=round(mean(wres_move),3)),
  by=.(player_name,pitch_type,season)][n>=25][order(player_name,pitch_type,season)])

ag <- off[, .(n=.N, spin=round(mean(release_spin_rate),0), spin_drop=round(mean(spin_drop),0),
  ivb=round(mean(ivb,na.rm=TRUE),1), as_gap=round(mean(as_gap,na.rm=TRUE),3),
  whiff=round(mean(is_whiff),3), exp=round(mean(exp_full),3), wres=round(mean(wres_full),3)),
  by=.(player_name,pitch_type,season)][n>=40]
cat("\n=== lowest-spin CHANGEUPS (CH only, all yrs, n>=40): do they underperform? ===\n")
print(head(ag[pitch_type=="CH"][order(spin), .(player_name,season,n,spin,spin_drop,ivb,as_gap,whiff,exp,wres)],20))
fwrite(ag[order(spin)], file.path(MDIR,"kick_change_allyears.csv"))
cat("\nsaved kick_change_allyears.csv\n")
