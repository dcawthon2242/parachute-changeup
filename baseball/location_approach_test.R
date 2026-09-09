#!/usr/bin/env Rscript

# Does LOCATION/COMMAND or APPROACH ANGLE explain offspeed whiff overperformance
# (the residual spin & tunneling couldn't)? Build OOF whiff residual over full shape,
# then merge plate location (below-zone / how low), VAA/HAA (approach angle), and
# release slot from raw kinematics (2023-2026). Correlate + multivariate + archetypes.

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

# ---- raw kinematics + location for VAA/HAA + plate loc, all seasons ----
kin <- rbindlist(lapply(2023:2026, function(y){
  fread(sprintf("data/statcast_%d/statcast_%d_all.csv", y, y),
    select=c("game_pk","at_bat_number","pitch_number","plate_x","plate_z","sz_top","sz_bot",
             "vx0","vy0","vz0","ax","ay","az")) }))
kin <- unique(kin, by=c("game_pk","at_bat_number","pitch_number"))
yf <- 17/12
kin[, vyf := -sqrt(pmax(vy0^2 - 2*ay*(50-yf), 0))]
kin[, tf := (vyf - vy0)/ay]
kin[, vzf := vz0 + az*tf][, vxf := vx0 + ax*tf]
kin[, VAA := atan2(vzf, vyf)*180/pi]         # negative = steeper downward
kin[, HAA := atan2(vxf, vyf)*180/pi]
kin[, below_zone := as.integer(plate_z < sz_bot)]
kin[, z_rel_bot := plate_z - sz_bot]          # how far above (or below, -) zone bottom
d <- merge(d, kin[, .(game_pk,at_bat_number,pitch_number,plate_x,plate_z,sz_bot,
                      VAA,HAA,below_zone,z_rel_bot)],
           by=c("game_pk","at_bat_number","pitch_number"), all.x=TRUE)

# ---- OOF whiff expectation over full shape (cache per-pitch residual) ----
D <- d[!is.na(is_whiff) & stats::complete.cases(d[, ..FEAT])]
K <- 4; D[, fold := sample(rep(1:K, length.out=.N))]
lgb_fit <- function(dtr){ n<-nrow(dtr); vi<-sample(n,floor(0.12*n))
  dtrain<-lgb.Dataset(as.matrix(dtr[-vi,..FEAT]),label=dtr$is_whiff[-vi])
  dval<-lgb.Dataset.create.valid(dtrain,as.matrix(dtr[vi,..FEAT]),label=dtr$is_whiff[vi])
  lgb.train(params=list(objective="binary",metric="binary_logloss",learning_rate=0.06,num_leaves=31,
    min_data_in_leaf=200,feature_fraction=0.8,bagging_fraction=0.8,bagging_freq=1),
    data=dtrain,nrounds=1200,valids=list(val=dval),early_stopping_rounds=50,verbose=-1) }
D[, exp_full := NA_real_]
for (f in 1:K){ cat("fold",f,"..\n"); mf<-lgb_fit(D[fold!=f])
  D[fold==f, exp_full := predict(mf, as.matrix(.SD[, FEAT, with=FALSE]))] }
D[, wres := is_whiff - exp_full]
saveRDS(D[, .(game_pk,at_bat_number,pitch_number,season,pitcher,player_name,pitch_type,grp,
   is_whiff,exp_full,wres,plate_x,plate_z,z_rel_bot,below_zone,VAA,HAA,as_gap,spin_sim,
   release_pos_x,release_pos_z)], file.path(MDIR,"oof_whiff_resid.rds"))

off <- D[pitch_type %in% c("CH","FS")]
CO <- function(a,b) cor(a,b,use="complete.obs")
cat(sprintf("\n=== offspeed CH+FS (n=%d): whiff overperformance vs location/approach ===\n", nrow(off)))
cat(sprintf("  pitch-level  wres vs below_zone : r=%+.3f\n", CO(off$wres, off$below_zone)))
cat(sprintf("  pitch-level  wres vs z_rel_bot  : r=%+.3f  (lower in zone -> ?)\n", CO(off$wres, off$z_rel_bot)))
cat(sprintf("  pitch-level  wres vs VAA        : r=%+.3f  (steeper=more negative)\n", CO(off$wres, off$VAA)))
cat(sprintf("  pitch-level  wres vs |HAA|      : r=%+.3f\n", CO(off$wres, abs(off$HAA))))

ag <- off[, .(n=.N, whiff=round(mean(is_whiff),3), exp=round(mean(exp_full),3), wres=mean(wres),
  below=mean(below_zone,na.rm=TRUE), z_rel=mean(z_rel_bot,na.rm=TRUE), VAA=mean(VAA,na.rm=TRUE),
  cmd_z=sd(plate_z,na.rm=TRUE), cmd_x=sd(plate_x,na.rm=TRUE),
  as_gap=mean(as_gap,na.rm=TRUE), spin_sim=mean(spin_sim,na.rm=TRUE)),
  by=.(player_name,pitch_type,season)][n>=40]
cat(sprintf("\n=== pitcher x pt x season aggregate (n>=40: %d) ===\n", nrow(ag)))
cat(sprintf("  wres vs below-zone rate     : r=%+.3f\n", CO(ag$wres, ag$below)))
cat(sprintf("  wres vs mean z above bottom : r=%+.3f\n", CO(ag$wres, ag$z_rel)))
cat(sprintf("  wres vs VAA (steeper)       : r=%+.3f\n", CO(ag$wres, ag$VAA)))
cat(sprintf("  wres vs command (SD plate_z): r=%+.3f  (tighter=lower SD)\n", CO(ag$wres, ag$cmd_z)))
cat(sprintf("  wres vs command (SD plate_x): r=%+.3f\n", CO(ag$wres, ag$cmd_x)))
cat(sprintf("  wres vs spin_sim (recap)    : r=%+.3f\n", CO(ag$wres, ag$spin_sim)))

# multivariate: variance in overperformance explained by location+approach
mv <- ag[stats::complete.cases(ag[, .(wres,below,z_rel,VAA,cmd_z,cmd_x)])]
fit <- lm(wres ~ below + z_rel + VAA + cmd_z + cmd_x, data=mv)
cat(sprintf("\n  multivariate R^2 (location+approach explain overperformance) = %.3f\n", summary(fit)$r.squared))

cat("\n=== archetypes: location & approach ===\n")
print(ag[grepl("Cease|Vesia|Chivilli|Mu.oz|Birdsong",player_name)][order(-wres),
  .(player_name,pitch_type,season,n, wres=round(wres,3), below=round(below,2),
    z_rel=round(z_rel,2), VAA=round(VAA,1), cmd_z=round(cmd_z,2))])
fwrite(ag[order(-wres)], file.path(MDIR,"location_approach_offspeed.csv"))
cat("\nsaved oof_whiff_resid.rds + location_approach_offspeed.csv\n")
