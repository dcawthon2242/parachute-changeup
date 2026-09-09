#!/usr/bin/env Rscript

# Who overperforms AFTER controlling for shape AND location/approach angle?
# Build OOF whiff expectation on FEAT + plate location + VAA/HAA. Residual = whiffs
# beyond what shape+location predict. Rank remaining overperformers; test whether
# this location-adjusted residual finally links to SPIN SIMILARITY.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(1)
MDIR <- file.path("data","statcast_model")
FEAT <- readRDS(file.path(MDIR, "arch_compare.rds"))$FEAT
as_long <- readRDS(file.path(MDIR, "active_spin_long.rds"))
d <- readRDS(file.path(MDIR, "miss_grade_data_activespin.rds"))
fb <- as_long[pitch_type %in% c("FF","SI","FC")][, pr:=match(pitch_type,c("FF","SI","FC"))][
  order(pitcher,season,pr)][, .SD[1], by=.(pitcher,season)][, .(pitcher,season,fb_active=active_spin)]
d <- merge(d, fb, by=c("pitcher","season"), all.x=TRUE)
d[, as_gap := active_spin - fb_active]
d[, spin_sim := exp(-(abs(as_gap)/0.10)^2)*exp(-(axis_diff/45)^2)]

kin <- rbindlist(lapply(2023:2026, function(y) fread(sprintf("data/statcast_%d/statcast_%d_all.csv",y,y),
  select=c("game_pk","at_bat_number","pitch_number","plate_x","plate_z","sz_bot","vy0","vz0","ay","az"))))
kin <- unique(kin, by=c("game_pk","at_bat_number","pitch_number"))
yf<-17/12; kin[,vyf:=-sqrt(pmax(vy0^2-2*ay*(50-yf),0))]; kin[,tf:=(vyf-vy0)/ay]; kin[,vzf:=vz0+az*tf]
kin[, VAA := atan2(vzf,abs(vyf))*180/pi]; kin[, HAA := atan2((vx0<-NA),1)]  # HAA set below
kin2 <- rbindlist(lapply(2023:2026, function(y) fread(sprintf("data/statcast_%d/statcast_%d_all.csv",y,y),
  select=c("game_pk","at_bat_number","pitch_number","vx0","ax"))))
kin2 <- unique(kin2, by=c("game_pk","at_bat_number","pitch_number"))
kin <- merge(kin, kin2, by=c("game_pk","at_bat_number","pitch_number"))
kin[, vxf := vx0 + ax*tf][, HAA := atan2(vxf, abs(vyf))*180/pi]
kin[, z_rel_bot := plate_z - sz_bot]
d <- merge(d, kin[, .(game_pk,at_bat_number,pitch_number,plate_x,plate_z,z_rel_bot,VAA,HAA)],
           by=c("game_pk","at_bat_number","pitch_number"), all.x=TRUE)

LOC <- c("plate_x","plate_z","z_rel_bot","VAA","HAA")
FULL <- c(FEAT, LOC)
D <- d[!is.na(is_whiff) & stats::complete.cases(d[, ..FULL])]
K<-4; D[, fold := sample(rep(1:K, length.out=.N))]
fit <- function(dtr, feats){ n<-nrow(dtr); vi<-sample(n,floor(0.12*n))
  dtrain<-lgb.Dataset(as.matrix(dtr[-vi,..feats]),label=dtr$is_whiff[-vi])
  dval<-lgb.Dataset.create.valid(dtrain,as.matrix(dtr[vi,..feats]),label=dtr$is_whiff[vi])
  lgb.train(params=list(objective="binary",metric="binary_logloss",learning_rate=0.06,num_leaves=31,
    min_data_in_leaf=200,feature_fraction=0.8,bagging_fraction=0.8,bagging_freq=1),
    data=dtrain,nrounds=1200,valids=list(val=dval),early_stopping_rounds=50,verbose=-1) }
D[, `:=`(exp_shape=NA_real_, exp_full=NA_real_)]
for (f in 1:K){ cat("fold",f,"..\n")
  D[fold==f, exp_shape := predict(fit(D[fold!=f], FEAT), as.matrix(.SD[, FEAT, with=FALSE]))]
  D[fold==f, exp_full  := predict(fit(D[fold!=f], FULL), as.matrix(.SD[, FULL, with=FALSE]))] }
D[, res_shape := is_whiff - exp_shape][, res_loc := is_whiff - exp_full]

off <- D[pitch_type %in% c("CH","FS")]
CO <- function(a,b) cor(a,b,use="complete.obs")
cat("\n=== does location-adjusted overperformance link to spin similarity? (offspeed) ===\n")
ag <- off[, .(n=.N, whiff=round(mean(is_whiff),3),
   res_shape=mean(res_shape), res_loc=mean(res_loc),
   spin_sim=mean(spin_sim,na.rm=TRUE), as_gap=mean(as_gap,na.rm=TRUE), axis=mean(axis_diff,na.rm=TRUE)),
   by=.(player_name,pitch_type)][n>=40]
cat(sprintf("  aggregate (pooled seasons, n>=40: %d pitchers)\n", nrow(ag)))
cat(sprintf("  res_SHAPE-only vs spin_sim   : r=%+.3f\n", CO(ag$res_shape, ag$spin_sim)))
cat(sprintf("  res_LOC-adjusted vs spin_sim : r=%+.3f\n", CO(ag$res_loc,  ag$spin_sim)))
cat(sprintf("  res_LOC-adjusted vs |as_gap| : r=%+.3f\n", CO(ag$res_loc,  abs(ag$as_gap))))
cat(sprintf("  res_LOC-adjusted vs axis gap : r=%+.3f\n", CO(ag$res_loc,  ag$axis)))
cat(sprintf("  var explained by location: shape-resid SD=%.3f -> loc-resid SD=%.3f (%.0f%% shrink)\n",
    sd(ag$res_shape), sd(ag$res_loc), 100*(1-sd(ag$res_loc)/sd(ag$res_shape))))

cat("\n=== TOP 20 REMAINING overperformers after shape+location (pooled, n>=40) ===\n")
print(head(ag[order(-res_loc), .(player_name,pitch_type,n, whiff,
   res_shape=round(res_shape,3), res_loc=round(res_loc,3),
   spin_sim=round(spin_sim,3), as_gap=round(as_gap,3))], 20))

cat("\n=== archetypes ===\n")
print(ag[grepl("Cease|Vesia|Chivilli|Ragans|Mu.oz|Birdsong",player_name),
   .(player_name,pitch_type,n,whiff, res_shape=round(res_shape,3), res_loc=round(res_loc,3),
     spin_sim=round(spin_sim,3), as_gap=round(as_gap,3))][order(-res_loc)])

cat("\n--- is the remaining top cohort more spin-similar than average? ---\n")
top <- ag[order(-res_loc)][1:30]
cat(sprintf("  top-30 remaining: mean spin_sim=%.3f, mean|as_gap|=%.3f\n", mean(top$spin_sim,na.rm=TRUE), mean(abs(top$as_gap),na.rm=TRUE)))
cat(sprintf("  all offspeed    : mean spin_sim=%.3f, mean|as_gap|=%.3f\n", mean(ag$spin_sim,na.rm=TRUE), mean(abs(ag$as_gap),na.rm=TRUE)))
fwrite(ag[order(-res_loc)], file.path(MDIR,"remaining_residual_offspeed.csv"))
cat("\nsaved remaining_residual_offspeed.csv\n")
