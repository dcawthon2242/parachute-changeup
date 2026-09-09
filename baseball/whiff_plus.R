#!/usr/bin/env Rscript
# tscore+ : predict a pitcher's timing-disruption score (tscore, runs saved /100 swings from
# contact depth) from pitch characteristics. Four feature approaches, 2026 holdout.
suppressPackageStartupMessages({ library(data.table); library(splines); library(lightgbm) }); set.seed(41); options(width=230)
cols <- c("game_year","game_type","player_name","pitcher","batter","stand","p_throws","pitch_type","description","balls","strikes",
          "plate_x","plate_z","sz_top","sz_bot","release_speed","pfx_x","pfx_z","release_pos_x","release_pos_z","release_extension","arm_angle","spin_axis","release_spin_rate",
          "delta_run_exp","intercept_ball_minus_batter_pos_y_inches","vy0","vz0","ay","az")
dt <- rbindlist(lapply(2024:2026, function(yr) fread(sprintf("data/statcast_%d/statcast_%d_all.csv",yr,yr), showProgress=FALSE, select=cols)))
setnames(dt, "intercept_ball_minus_batter_pos_y_inches", "depth")
dt <- dt[game_type=="R" & pitch_type!="" & balls<=3 & strikes<=2 & is.finite(delta_run_exp) & is.finite(release_speed) & is.finite(plate_x) & is.finite(plate_z) & is.finite(pfx_x) & is.finite(pfx_z)]
FB<-c("FF","SI","FC"); BR<-c("SL","ST","CU","KC","SV","CS"); OS<-c("CH","FS","FO")
dt[, pgrp := fifelse(pitch_type %in% FB,"FB", fifelse(pitch_type %in% BR,"BR", fifelse(pitch_type %in% OS,"OS",NA_character_)))]; dt <- dt[!is.na(pgrp)]
dt[, pg2 := fifelse(pgrp=="FB","FB","OFF")]
dt[, hb_arm := fifelse(p_throws=="R", -pfx_x*12, pfx_x*12)]; dt[, ivb := pfx_z*12]
dt[, rel_x_arm := fifelse(p_throws=="R", -release_pos_x, release_pos_x)]
dt[, px_bat := fifelse(stand=="R",-plate_x,plate_x)]; dt[, pz_rel := (plate_z-sz_bot)/pmax(sz_top-sz_bot,0.1)]
dt[, tt := (-vy0 - sqrt(pmax(vy0^2 - 2*ay*(50-17/12),0)))/ay]; dt[, vaa := atan2(vz0+az*tt, abs(vy0+ay*tt))*180/pi]
# fastball reference per pitcher-season (FF > SI > FC by usage)
fbref <- dt[pitch_type %in% FB, .(n=.N, fb_velo=mean(release_speed), fb_ivb=mean(ivb), fb_hb=mean(hb_arm), fb_axis=atan2(mean(sin(spin_axis*pi/180),na.rm=TRUE),mean(cos(spin_axis*pi/180),na.rm=TRUE))*180/pi,
                                  fb_spin=mean(release_spin_rate,na.rm=TRUE), fb_relz=mean(release_pos_z,na.rm=TRUE), fb_vaa=mean(vaa,na.rm=TRUE)), by=.(pitcher,game_year,pitch_type)][order(pitcher,game_year,-n)][, .SD[1], by=.(pitcher,game_year)]
dt <- merge(dt, fbref[, -c("n","pitch_type")], by=c("pitcher","game_year"))
dt[, sep := fb_velo-release_speed]; dt[, ivb_gap := fb_ivb-ivb]; dt[, hb_gap := fb_hb-hb_arm]
dt[, axis_diff := abs(((spin_axis-fb_axis+180) %% 360)-180)]; dt[, spin_gap := fb_spin-release_spin_rate]; dt[, vaa_gap := fb_vaa-vaa]
usage <- dt[, .(n_pt=.N), by=.(pitcher,game_year,pitch_type)]; usage[, use := n_pt/sum(n_pt), by=.(pitcher,game_year)]; usage[, n_pt := NULL]
dt <- merge(dt, usage, by=c("pitcher","game_year","pitch_type"))
fbsh <- dt[, .(fbsh=mean(pg2=="FB")), by=.(pitcher,game_year)]; dt <- merge(dt, fbsh, by=c("pitcher","game_year"))
# ---- whiff+ : shape-only probability of a whiff given a swing (model-B features), OOF by season
SWING <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip","hit_into_play")
sw <- dt[description %in% SWING]
sw[, y := as.numeric(description %in% c("swinging_strike","swinging_strike_blocked","foul_tip"))]
sw[, hand := as.integer(p_throws=="R")]; sw[, pgi := match(pgrp,c("FB","BR","OS"))]; sw[, pti := match(pitch_type, c("FF","SI","FC","SL","ST","CU","KC","SV","CH","FS","FO","CS"))]
FEAT <- c("release_speed","ivb","hb_arm","rel_x_arm","release_pos_z","release_extension","arm_angle","release_spin_rate","spin_axis","vaa","hand","pgi","pti","sep","ivb_gap","hb_gap","axis_diff","spin_gap","vaa_gap","fb_velo","fb_ivb","fb_hb","fb_relz","use","fbsh")
fitw <- function(tr, te) { vi <- sample(nrow(tr), floor(0.1*nrow(tr)))
  dtr <- lgb.Dataset(as.matrix(tr[-vi, ..FEAT]), label=tr$y[-vi]); dva <- lgb.Dataset.create.valid(dtr, as.matrix(tr[vi, ..FEAT]), label=tr$y[vi])
  m <- lgb.train(params=list(objective="binary", metric="binary_logloss", learning_rate=0.05, num_leaves=63, min_data_in_leaf=400, feature_fraction=0.8, bagging_fraction=0.8, bagging_freq=1, lambda_l2=5), data=dtr, nrounds=3000, valids=list(v=dva), early_stopping_rounds=100, verbose=-1)
  list(m=m, pred=predict(m, as.matrix(te[, ..FEAT]))) }
sw[, pred := NA_real_]; for (yr in 2024:2026) { g <- fitw(sw[game_year!=yr], sw[game_year==yr]); sw[game_year==yr, pred := g$pred]; cat("oof", yr, "\n"); if (yr==2026) print(lgb.importance(g$m)[1:12, .(Feature, Gain=round(Gain,3))]) }
saveRDS(sw[, .(pitcher, game_year, pitch_type, pgrp, pg2, y, pred)], "data/swing_timing/whiff_plus_oof_swings.rds")
P <- sw[, .(swings=.N, whiff=mean(y), whiff_plus_raw=mean(pred)), by=.(pitcher,game_year)][swings>=300]
cat(sprintf("pitcher-seasons %d | OOF cor(whiff+, whiff) = %.3f\n", nrow(P), cor(P$whiff_plus_raw, P$whiff)))
a <- P; b <- copy(P)[, game_year := game_year-1L]; q <- merge(a,b,by=c("pitcher","game_year"),suffixes=c("","_n"))
cat(sprintf("yoy: whiff %.3f | whiff+ %.3f | whiff+ -> next whiff %.3f | whiff -> next whiff %.3f (n=%d)\n", cor(q$whiff,q$whiff_n), cor(q$whiff_plus_raw,q$whiff_plus_raw_n), cor(q$whiff_plus_raw,q$whiff_n), cor(q$whiff,q$whiff_n), nrow(q)))
fwrite(P, "data/swing_timing/whiff_plus_pitchers.csv"); cat("DONE\n")
