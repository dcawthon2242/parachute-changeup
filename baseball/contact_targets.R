#!/usr/bin/env Rscript
# tscore+ : predict a pitcher's timing-disruption score (tscore, runs saved /100 swings from
# contact depth) from pitch characteristics. Four feature approaches, 2026 holdout.
suppressPackageStartupMessages({ library(data.table); library(splines); library(lightgbm) }); set.seed(41); options(width=230)
cols <- c("game_year","game_type","player_name","pitcher","batter","stand","p_throws","pitch_type","description","balls","strikes",
          "plate_x","plate_z","sz_top","sz_bot","release_speed","pfx_x","pfx_z","release_pos_x","release_pos_z","release_extension","arm_angle","spin_axis","release_spin_rate",
          "delta_run_exp","intercept_ball_minus_batter_pos_y_inches","vy0","vz0","ay","az","launch_speed","launch_angle","launch_speed_angle","estimated_woba_using_speedangle","bb_type","events")
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
# ---- alternative contact targets, all modeled from the same model-B shape features, OOF by season
dt[, hand := as.integer(p_throws=="R")]; dt[, pgi := match(pgrp,c("FB","BR","OS"))]; dt[, pti := match(pitch_type, c("FF","SI","FC","SL","ST","CU","KC","SV","CH","FS","FO","CS"))]
FEAT <- c("release_speed","ivb","hb_arm","rel_x_arm","release_pos_z","release_extension","arm_angle","release_spin_rate","spin_axis","vaa","hand","pgi","pti","sep","ivb_gap","hb_gap","axis_diff","spin_gap","vaa_gap","fb_velo","fb_ivb","fb_hb","fb_relz","use","fbsh")
extra <- rbindlist(lapply(2024:2026, function(yr) fread(sprintf("data/statcast_%d/statcast_%d_all.csv",yr,yr), select=c("game_pk","at_bat_number","pitch_number","launch_speed","launch_angle","launch_speed_angle","estimated_woba_using_speedangle","bb_type","events","woba_value"), showProgress=FALSE)))
SWING <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip","hit_into_play")
dt[, swung := description %in% SWING]; dt[, bip := description=="hit_into_play" & is.finite(launch_speed) & is.finite(launch_angle)]
B <- dt[bip==TRUE]
B[, `:=`(t_xw = estimated_woba_using_speedangle, t_ev = launch_speed, t_hard = as.numeric(launch_speed>=95), t_barrel = as.numeric(launch_speed_angle==6),
         t_la = launch_angle, t_gb = as.numeric(launch_angle < 10), t_hr = as.numeric(events=="home_run"), t_air_ev = fifelse(launch_angle>=10, launch_speed, NA_real_),
         t_under = as.numeric(launch_speed_angle==3), t_topped = as.numeric(launch_speed_angle %in% 1:2), t_rv_con = -delta_run_exp)]
# per-pitch run value (all pitches) and per-swing contact-or-whiff run value as broader targets
dt[, t_rv_pitch := -delta_run_exp]; S <- dt[swung==TRUE]; S[, t_rv_swing := -delta_run_exp]
TG <- list(xw="t_xw", ev="t_ev", hard="t_hard", barrel="t_barrel", la="t_la", gb="t_gb", hr="t_hr", air_ev="t_air_ev", under="t_under", topped="t_topped", rv_con="t_rv_con")
fit1 <- function(d, y, binary) { d <- d[is.finite(get(y))]; d[, pred := NA_real_]
  for (yr in 2024:2026) { tr <- d[game_year!=yr]; te <- d[game_year==yr]; vi <- sample(nrow(tr), floor(0.1*nrow(tr)))
    dtr <- lgb.Dataset(as.matrix(tr[-vi, ..FEAT]), label=tr[[y]][-vi]); dva <- lgb.Dataset.create.valid(dtr, as.matrix(tr[vi, ..FEAT]), label=tr[[y]][vi])
    m <- lgb.train(params=list(objective=if (binary) "binary" else "regression", metric=if (binary) "binary_logloss" else "rmse", learning_rate=0.05, num_leaves=63, min_data_in_leaf=400, feature_fraction=0.8, bagging_fraction=0.8, bagging_freq=1, lambda_l2=5), data=dtr, nrounds=3000, valids=list(v=dva), early_stopping_rounds=100, verbose=-1)
    d[game_year==yr, pred := predict(m, as.matrix(te[, ..FEAT]))] }
  d[, .(n=.N, actual=mean(get(y)), modeled=mean(pred)), by=.(pitcher,game_year)] }
out <- list()
for (nm in names(TG)) { y <- TG[[nm]]; bin <- nm %in% c("hard","barrel","gb","hr","under","topped"); r <- fit1(B, y, bin); setnames(r, c("n","actual","modeled"), c(paste0("n_",nm), paste0("act_",nm), paste0("mod_",nm))); out[[nm]] <- r; cat("done", nm, "\n") }
r <- fit1(dt, "t_rv_pitch", FALSE); setnames(r, c("n","actual","modeled"), c("n_rvp","act_rvp","mod_rvp")); out[["rvp"]] <- r; cat("done rv_pitch\n")
r <- fit1(S, "t_rv_swing", FALSE); setnames(r, c("n","actual","modeled"), c("n_rvs","act_rvs","mod_rvs")); out[["rvs"]] <- r; cat("done rv_swing\n")
X <- Reduce(function(a,b) merge(a,b,by=c("pitcher","game_year"),all=TRUE), out)
fwrite(X, "data/swing_timing/contact_targets_pitchers.csv"); cat("DONE\n")
