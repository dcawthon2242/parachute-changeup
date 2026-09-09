#!/usr/bin/env Rscript
# bscore+ : shape-only model of the 3-axis barrel-positioning score (bscore), model-B feature set, OOF by season.
suppressPackageStartupMessages({ library(data.table); library(splines); library(lightgbm) }); set.seed(41); options(width=240)
cols <- c("game_year","game_type","player_name","pitcher","batter","stand","p_throws","pitch_type","description","bb_type","balls","strikes","plate_x","plate_z","sz_top","sz_bot",
          "release_speed","pfx_x","pfx_z","release_pos_x","release_pos_z","release_extension","arm_angle","spin_axis","release_spin_rate","launch_speed","launch_speed_angle","estimated_woba_using_speedangle","delta_run_exp","vy0","vz0","ay","az",
          "intercept_ball_minus_batter_pos_x_inches","intercept_ball_minus_batter_pos_y_inches","attack_angle")
dt <- rbindlist(lapply(2024:2026, function(yr) fread(sprintf("data/statcast_%d/statcast_%d_all.csv",yr,yr), showProgress=FALSE, select=cols)))
setnames(dt, c("intercept_ball_minus_batter_pos_x_inches","intercept_ball_minus_batter_pos_y_inches"), c("ix","depth"))
dt <- dt[game_type=="R" & pitch_type!="" & balls<=3 & strikes<=2 & is.finite(delta_run_exp) & is.finite(release_speed) & is.finite(plate_x) & is.finite(plate_z) & is.finite(pfx_x) & is.finite(pfx_z)]
FB<-c("FF","SI","FC"); BR<-c("SL","ST","CU","KC","SV","CS"); OS<-c("CH","FS","FO")
dt[, pgrp := fifelse(pitch_type %in% FB,"FB", fifelse(pitch_type %in% BR,"BR", fifelse(pitch_type %in% OS,"OS",NA_character_)))]; dt <- dt[!is.na(pgrp)]; dt[, pg2 := fifelse(pgrp=="FB","FB","OFF")]
dt[, hb_arm := fifelse(p_throws=="R", -pfx_x*12, pfx_x*12)]; dt[, ivb := pfx_z*12]; dt[, rel_x_arm := fifelse(p_throws=="R", -release_pos_x, release_pos_x)]
dt[, px_bat := fifelse(stand=="R",-plate_x,plate_x)]; dt[, pz_rel := (plate_z-sz_bot)/pmax(sz_top-sz_bot,0.1)]
dt[, tt := (-vy0 - sqrt(pmax(vy0^2 - 2*ay*(50-17/12),0)))/ay]; dt[, vaa := atan2(vz0+az*tt, abs(vy0+ay*tt))*180/pi]
fbref <- dt[pitch_type %in% FB, .(n=.N, fb_velo=mean(release_speed), fb_ivb=mean(ivb), fb_hb=mean(hb_arm), fb_axis=atan2(mean(sin(spin_axis*pi/180),na.rm=TRUE),mean(cos(spin_axis*pi/180),na.rm=TRUE))*180/pi,
                                  fb_spin=mean(release_spin_rate,na.rm=TRUE), fb_relz=mean(release_pos_z,na.rm=TRUE), fb_vaa=mean(vaa,na.rm=TRUE)), by=.(pitcher,game_year,pitch_type)][order(pitcher,game_year,-n)][, .SD[1], by=.(pitcher,game_year)]
dt <- merge(dt, fbref[, -c("n","pitch_type")], by=c("pitcher","game_year"))
dt[, sep := fb_velo-release_speed]; dt[, ivb_gap := fb_ivb-ivb]; dt[, hb_gap := fb_hb-hb_arm]; dt[, axis_diff := abs(((spin_axis-fb_axis+180) %% 360)-180)]; dt[, spin_gap := fb_spin-release_spin_rate]; dt[, vaa_gap := fb_vaa-vaa]
usage <- dt[, .(n_pt=.N), by=.(pitcher,game_year,pitch_type)]; usage[, use := n_pt/sum(n_pt), by=.(pitcher,game_year)]; usage[, n_pt := NULL]; dt <- merge(dt, usage, by=c("pitcher","game_year","pitch_type"))
dt <- merge(dt, dt[, .(fbsh=mean(pg2=="FB")), by=.(pitcher,game_year)], by=c("pitcher","game_year"))
# ---- bscore target
SWING <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip","hit_into_play")
sw <- dt[description %in% SWING & is.finite(depth) & is.finite(ix) & is.finite(attack_angle) & is.finite(vaa)]
sw[, ix_bat := fifelse(stand=="R", ix, -ix)]; sw[, zraw := attack_angle + vaa]
ra <- function(v) residuals(lm(as.formula(paste(v, "~ ns(px_bat,5)*pgrp + ns(pz_rel,5) + ns(release_speed,4) + stand")), data=sw))
sw[, r_t := ra("depth")]; sw[, r_x := ra("ix_bat")]; sw[, r_z := ra("zraw")]; sw[, nb := .N, by=batter]; sw <- sw[nb>=200]
sw[, tdev := r_t-mean(r_t), by=batter]; sw[, xdev := r_x-mean(r_x), by=batter]; sw[, zdev := r_z-mean(r_z), by=batter]
sw[, cnti := balls*3+strikes]; sw[, pgi := match(pgrp,c("FB","BR","OS"))]
Xs <- as.matrix(sw[, .(tdev, xdev, zdev, pgi, cnti)])
ms <- lgb.train(params=list(objective="regression", learning_rate=0.05, num_leaves=15, min_data_in_leaf=2000, lambda_l2=10, verbose=-1), data=lgb.Dataset(Xs, label=sw$delta_run_exp), nrounds=400)
pj <- predict(ms, Xs); pjb <- pj - ave(pj, sw$batter); sw[, y := -100*(pjb - ave(pjb, sw$pg2))]      # per-swing bscore, runs/100, pitcher-positive
sw[, whiff := description %in% c("swinging_strike","swinging_strike_blocked","foul_tip")]; sw[, isbip := description=="hit_into_play" & bb_type!="" & is.finite(launch_speed) & is.finite(launch_speed_angle)]
sw[, hand := as.integer(p_throws=="R")]; sw[, pti := match(pitch_type, c("FF","SI","FC","SL","ST","CU","KC","SV","CH","FS","FO","CS"))]
FEAT <- c("release_speed","ivb","hb_arm","rel_x_arm","release_pos_z","release_extension","arm_angle","release_spin_rate","spin_axis","vaa","hand","pgi","pti",
          "sep","ivb_gap","hb_gap","axis_diff","spin_gap","vaa_gap","fb_velo","fb_ivb","fb_hb","fb_relz","use","fbsh")
fitgbm <- function(tr, te) { vi <- sample(nrow(tr), floor(0.1*nrow(tr)))
  dtr <- lgb.Dataset(as.matrix(tr[-vi, ..FEAT]), label=tr$y[-vi]); dva <- lgb.Dataset.create.valid(dtr, as.matrix(tr[vi, ..FEAT]), label=tr$y[vi])
  m <- lgb.train(params=list(objective="regression", metric="rmse", learning_rate=0.05, num_leaves=63, min_data_in_leaf=400, feature_fraction=0.8, bagging_fraction=0.8, bagging_freq=1, lambda_l2=5), data=dtr, nrounds=3000, valids=list(v=dva), early_stopping_rounds=100, verbose=-1)
  list(m=m, pred=predict(m, as.matrix(te[, ..FEAT]))) }
sw[, pred := NA_real_]; imps <- list(); mods <- list()
for (yr in 2024:2026) { g <- fitgbm(sw[game_year!=yr], sw[game_year==yr]); sw[game_year==yr, pred := g$pred]; imps[[as.character(yr)]] <- lgb.importance(g$m); mods[[as.character(yr)]] <- g$m; cat("oof", yr, "\n") }
saveRDS(sw[, .(pitcher, game_year, pitch_type, pgrp, pg2, y, pred, tdev, xdev, zdev, whiff, isbip, launch_speed, estimated_woba_using_speedangle)], "data/swing_timing/bscore_plus_oof_swings.rds")
IMP <- rbindlist(imps)[, .(Gain=mean(Gain)), by=Feature][order(-Gain)]
cat("\n=== bscore+ feature importance (mean gain share across the 3 OOF fits) ===\n"); print(IMP[, .(Feature, Gain=round(Gain,3))])
# ---- pitcher-season
P <- sw[, .(swings=.N, bscore=mean(y), bplus=mean(pred), whiff=100*mean(whiff), hard=100*sum(isbip & launch_speed>=95)/pmax(sum(isbip),1), ev=mean(launch_speed[isbip]), xw=mean(estimated_woba_using_speedangle[isbip],na.rm=TRUE)), by=.(pitcher,game_year)][swings>=300]
P[, resid := bscore-bplus]; NM <- unique(dt[, .(pitcher,player_name)], by="pitcher"); P <- merge(P, NM, by="pitcher")
cat(sprintf("\npitcher-seasons %d | OOF cor(bscore+, bscore) = %.3f | resid SD %.3f vs bscore SD %.3f\n", nrow(P), cor(P$bplus,P$bscore), sd(P$resid), sd(P$bscore)))
a <- P[, .(pitcher, game_year, bscore, bplus, resid, hard, ev, xw)]; b <- copy(a)[, game_year := game_year-1L]; q <- merge(a,b,by=c("pitcher","game_year"), suffixes=c("","_n"))
cat(sprintf("yoy: bscore %.3f | bscore+ %.3f | residual %.3f (n=%d)\n", cor(q$bscore,q$bscore_n), cor(q$bplus,q$bplus_n), cor(q$resid,q$resid_n), nrow(q)))
cat(sprintf("predict next-year: bscore -> bscore %.3f | bscore+ -> bscore %.3f | both R %.3f\n", cor(q$bscore,q$bscore_n), cor(q$bplus,q$bscore_n), sqrt(summary(lm(bscore_n ~ bscore + bplus, data=q))$r.squared)))
for (v in c("hard","ev","xw")) cat(sprintf("  next %-4s: own R2 %.3f | own+bscore %.3f | own+bscore+ %.3f | bscore+ alone %.3f\n", v, summary(lm(as.formula(paste0(v,"_n ~ ",v)),q))$r.squared, summary(lm(as.formula(paste0(v,"_n ~ ",v,"+bscore")),q))$r.squared, summary(lm(as.formula(paste0(v,"_n ~ ",v,"+bplus")),q))$r.squared, summary(lm(as.formula(paste0(v,"_n ~ bplus")),q))$r.squared))
# ---- partial dependence (2026 model, 40k sample), p5->p95 swing in bscore+ points (10 pts = 1 SD of pitcher bscore)
sdB <- sd(P[game_year<=2025]$bscore); muB <- weighted.mean(P[game_year<=2025]$bscore, P[game_year<=2025]$swings); sc <- 10/sdB
m <- mods[["2026"]]; S <- sw[game_year==2026][sample(.N, 40000)]
pd <- function(v, grp) { s <- S[pgrp %in% grp]; qs <- quantile(s[[v]], c(.05,.25,.5,.75,.95), na.rm=TRUE); base <- as.matrix(s[, ..FEAT]); p <- sapply(qs, function(x) { X <- base; X[, v] <- x; mean(predict(m, X)) }); list(q=qs, p=(p-mean(p))*sc) }
cat("\n=== partial dependence, bscore+ points vs mean at p5/p25/p50/p75/p95 ===\n")
for (grp in list(c("FB"), c("BR"), c("OS"))) { cat("--", grp, "--\n")
  vs <- if (grp=="FB") c("release_speed","ivb","hb_arm","vaa","release_pos_z","rel_x_arm","arm_angle","release_extension","release_spin_rate","fbsh","use") else c("release_speed","sep","ivb","hb_arm","vaa","vaa_gap","axis_diff","spin_gap","release_spin_rate","fb_velo","fb_ivb","use","fbsh")
  for (v in vs) { r <- pd(v, grp); cat(sprintf("  %-18s | %s | p5->p95 %+.1f\n", v, paste(sprintf("%s:%+.1f", formatC(r$q,digits=3,format="fg"), r$p), collapse="  "), r$p[5]-r$p[1])) } }
# ---- residual correlates using the existing pitcher-season feature bank
FBK <- fread("data/swing_timing/tscore_plus_residual_features.csv"); drop <- c("swings","tscore","plus","fb_ts","off_ts","fb_plus","off_plus","resid","fb_resid","off_resid","pitches","g","tdev_sd","tdev_fb","tdev_off","whiff","hard","ev","xw")
FBK <- FBK[, setdiff(names(FBK), drop), with=FALSE]; R <- merge(P, FBK, by=c("pitcher","game_year"))
SK <- setdiff(names(FBK), c("pitcher","game_year"))
RC <- rbindlist(lapply(SK, function(v) { x <- R[[v]]; ok <- is.finite(x); data.table(feature=v, r_resid=round(cor(x[ok],R$resid[ok]),3), r_bscore=round(cor(x[ok],R$bscore[ok]),3), r_bplus=round(cor(x[ok],R$bplus[ok]),3)) }))
cat("\n=== residual correlates (actual bscore - OOF bscore+), |r| >= 0.10 ===\n"); print(RC[abs(r_resid)>=0.10][order(-abs(r_resid))])
blk <- list(count=c("strike1","behind","ahead","ball"), swing_dec=c("swing","chase","zswing"), location=c("fb_locsd","off_locsd","dloc_prev","off_px_in","zone","meat","edge","fb_high","off_low"), sequencing=c("dvelo_prev","same_as_prev","off_after_fb","path_ratio"), opp=c("opp_bat_speed","opp_swing_len","opp_attack","opp_whiff_tend"), release=c("fb_ext","rel_x_spread","rel_z_spread","arm_spread","within_velo_sd","within_relx_sd"))
for (bn in names(blk)) { v <- blk[[bn]]; Z <- R[complete.cases(R[, ..v])]; cat(sprintf("  block %-10s R2 = %.3f\n", bn, summary(lm(as.formula(paste("resid ~", paste(v,collapse="+"))), data=Z, weights=Z$swings))$r.squared)) }
TOP <- RC[order(-abs(r_resid))][1:12]$feature; Z <- R[complete.cases(R[, ..TOP])]; mm <- lm(as.formula(paste("resid ~", paste(TOP,collapse="+"))), data=Z, weights=Z$swings)
cat(sprintf("  top-12 features: in-sample R2 %.3f; leave-one-season-out r: %s\n", summary(mm)$r.squared, paste(round(sapply(2024:2026, function(yr) { f <- lm(as.formula(paste("resid ~", paste(TOP,collapse="+"))), data=Z[game_year!=yr], weights=Z[game_year!=yr]$swings); cor(predict(f, Z[game_year==yr]), Z[game_year==yr]$resid) }),3), collapse=" ")))
# ---- per-pitch grades: bscore+ and tscore+ (from saved OOF), 2026, scaled 100/10 on pitcher-level SDs
TO <- readRDS("data/swing_timing/tscore_plus_oof_swings.rds"); TP <- TO[, .(swings=.N, tscore=mean(y), tplus=mean(pred)), by=.(pitcher,game_year)][swings>=300]
sdT <- sd(TP[game_year<=2025]$tscore); muT <- weighted.mean(TP[game_year<=2025]$tscore, TP[game_year<=2025]$swings)
gB <- function(x) round(100 + 10*(x-muB)/sdB); gT <- function(x) round(100 + 10*(x-muT)/sdT)
PP <- merge(sw[game_year==2026, .(sw_b=.N, bscore=mean(y), bplus=mean(pred)), by=.(pitcher,pitch_type)], TO[game_year==2026, .(sw_t=.N, tscore=mean(y), tplus=mean(pred)), by=.(pitcher,pitch_type)], by=c("pitcher","pitch_type"), all.x=TRUE)
PP <- merge(PP, usage[game_year==2026, .(pitcher, pitch_type, use)], by=c("pitcher","pitch_type"))
SH <- dt[game_year==2026, .(velo=mean(release_speed), ivb=mean(ivb), hb=mean(hb_arm), sep=mean(sep)), by=.(pitcher,pitch_type)]; PP <- merge(PP, SH, by=c("pitcher","pitch_type"))
L <- P[game_year==2026 & swings>=400]; L[, `:=`(bscore_plus=gB(bplus), bscore_act=gB(bscore))]; L <- merge(L, TP[game_year==2026, .(pitcher, tscore_plus=gT(tplus), tscore_act=gT(tscore))], by="pitcher", all.x=TRUE)
cat(sprintf("\nscale: bscore+ 100 = %.2f runs/100 swings, 10 pts = %.2f | tscore+ 10 pts = %.2f\n", muB, sdB, sdT))
showp <- function(ids, lab) { cat("\n########", lab, "########\n")
  for (id in ids) { r <- L[pitcher==id]; cat(sprintf("\n%s  | swings %d | bscore+ %d (actual %d) | tscore+ %d (actual %d) | whiff %.1f hard %.1f xw %.3f\n", r$player_name, r$swings, r$bscore_plus, r$bscore_act, r$tscore_plus, r$tscore_act, r$whiff, r$hard, r$xw))
    x <- PP[pitcher==id & use>=0.04][order(-use)]; print(x[, .(pitch=pitch_type, use=paste0(round(100*use),"%"), velo=round(velo,1), ivb=round(ivb), hb=round(hb), sep=round(sep,1), sw=sw_b, `bscore+`=gB(bplus), bscore=gB(bscore), `tscore+`=gT(tplus), tscore=gT(tscore))], row.names=FALSE) } }
showp(L[order(-bplus)][1:10]$pitcher, "10 BEST by bscore+ (2026, >=400 swings)"); showp(L[order(bplus)][1:10]$pitcher, "10 WORST by bscore+")
fwrite(L, "data/swing_timing/bscore_plus_2026.csv"); fwrite(PP, "data/swing_timing/bscore_plus_2026_per_pitch.csv"); cat("\nDONE\n")
