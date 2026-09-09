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
# ---- tscore target per swing
SWING <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip","hit_into_play")
sw <- dt[description %in% SWING & is.finite(depth)]
fitd <- lm(depth ~ ns(px_bat,5)*pgrp + ns(pz_rel,5) + ns(release_speed,4) + stand, data=sw)
sw[, r1 := residuals(fitd)]; sw[, nb := .N, by=batter]; sw <- sw[nb>=200]; sw[, tdev := r1-mean(r1), by=batter]; sw[, cnt := paste0(balls,"-",strikes)]
surf <- lm(delta_run_exp ~ ns(tdev,6)*pgrp + factor(cnt), data=sw); sw[, trv := predict(surf, newdata=sw)]
sw[, trv_b := trv-mean(trv), by=batter]; sw[, trv_ex := trv_b-mean(trv_b), by=pg2]
sw[, y := -100*trv_ex]                                      # per-swing timing run value, pitcher-positive, runs/100
sw[, hand := as.integer(p_throws=="R")]; sw[, same := as.integer(p_throws==stand)]
sw[, pgi := match(pgrp, c("FB","BR","OS"))]; sw[, pti := match(pitch_type, c("FF","SI","FC","SL","ST","CU","KC","SV","CH","FS","FO","CS"))]
sw[, cnti := balls*3+strikes]
actual <- sw[, .(swings=.N, tscore=mean(y)), by=.(pitcher,game_year)]
NM <- unique(dt[, .(pitcher,player_name)], by="pitcher")
# ---- feature sets
FA <- c("release_speed","ivb","hb_arm","rel_x_arm","release_pos_z","release_extension","arm_angle","release_spin_rate","spin_axis","vaa","hand","pgi","pti")          # A: pitch physics only
FB_ <- c(FA, "sep","ivb_gap","hb_gap","axis_diff","spin_gap","vaa_gap","fb_velo","fb_ivb","fb_hb","fb_relz","use","fbsh")                                          # B: + arsenal context
FC_ <- c(FB_, "px_bat","pz_rel","cnti","same")                                                                                                                        # C: + location & count
fitgbm <- function(feats, tr, te) {
  vi <- sample(nrow(tr), floor(0.1*nrow(tr)))
  dtr <- lgb.Dataset(as.matrix(tr[-vi, ..feats]), label=tr$y[-vi]); dva <- lgb.Dataset.create.valid(dtr, as.matrix(tr[vi, ..feats]), label=tr$y[vi])
  m <- lgb.train(params=list(objective="regression", metric="rmse", learning_rate=0.05, num_leaves=63, min_data_in_leaf=400, feature_fraction=0.8, bagging_fraction=0.8, bagging_freq=1, lambda_l2=5),
                 data=dtr, nrounds=3000, valids=list(v=dva), early_stopping_rounds=100, verbose=-1)
  list(m=m, pred=predict(m, as.matrix(te[, ..feats])), imp=lgb.importance(m))
}
evalp <- function(te, pred, label) {
  te <- copy(te)[, p := pred]
  P <- te[, .(swings=.N, tscore=mean(y), plus_raw=mean(p)), by=.(pitcher,game_year)][swings>=300]
  P <- merge(P, actual[game_year==2025, .(pitcher, tscore_prev=tscore, swings_prev=swings)], by="pitcher", all.x=TRUE)
  r <- cor(P$plus_raw, P$tscore); rmse_pitch <- sqrt(mean((te$y-te$p)^2)); rmse_null <- sqrt(mean((te$y-mean(te$y))^2))
  q <- P[swings_prev>=300 & is.finite(tscore_prev)]
  data.table(approach=label, pitch_rmse_lift=round(100*(1-rmse_pitch/rmse_null),2), pitchers=nrow(P), r_same_year=round(r,3), R2_same_year=round(r^2,3),
             r_2025plus_vs_2026actual=round(cor(q$plus_raw, q$tscore),3), r_2025actual_vs_2026actual=round(cor(q$tscore_prev, q$tscore),3), n_pairs=nrow(q))
}
tr <- sw[game_year<=2025]; te <- sw[game_year==2026]
cat(sprintf("train swings %s (2024-25) | test swings %s (2026)\n", format(nrow(tr),big.mark=","), format(nrow(te),big.mark=",")))
res <- list(); preds <- list(); imps <- list()
for (nm in c("A","B","C")) { feats <- get(paste0("F", if (nm=="A") "A" else paste0(nm,"_"))); g <- fitgbm(feats, tr, te); preds[[nm]] <- g$pred; imps[[nm]] <- g$imp
  # for the predictive test we need 2025 tscore+ from a model NOT trained on 2025: train on 2024+2026, predict 2025
  g2 <- fitgbm(feats, sw[game_year!=2025], sw[game_year==2025]); p25 <- copy(sw[game_year==2025])[, p := g2$pred][, .(swings=.N, plus25=mean(p)), by=pitcher][swings>=300]
  te2 <- copy(te)[, p := g$pred]; P26 <- te2[, .(swings=.N, tscore=mean(y), plus_raw=mean(p)), by=pitcher][swings>=300]
  q <- merge(merge(P26, p25[, .(pitcher, plus25)], by="pitcher"), actual[game_year==2025 & swings>=300, .(pitcher, ts25=tscore)], by="pitcher")
  rmse_pitch <- sqrt(mean((te$y-g$pred)^2)); rmse_null <- sqrt(mean((te$y-mean(te$y))^2))
  res[[nm]] <- data.table(approach=nm, pitch_rmse_lift_pct=round(100*(1-rmse_pitch/rmse_null),2), pitchers_2026=nrow(P26), r_same_year=round(cor(P26$plus_raw,P26$tscore),3),
                          n_pairs=nrow(q), r_plus25_vs_ts26=round(cor(q$plus25,q$tscore),3), r_ts25_vs_ts26=round(cor(q$ts25,q$tscore),3),
                          r_both=round(sqrt(summary(lm(tscore ~ plus25 + ts25, data=q))$r.squared),3))
}
# ---- D: pitcher-season arsenal profile, linear, train 2024-25 predict 2026
prof <- sw[, .(swings=.N, tscore=mean(y), fb_velo=fb_velo[1], fb_ivb=fb_ivb[1], fb_hb=fb_hb[1], fb_relz=fb_relz[1], fbsh=fbsh[1], arm=mean(arm_angle,na.rm=TRUE), ext=mean(release_extension,na.rm=TRUE),
                sep=mean(sep[pg2=="OFF"]), off_ivb=mean(ivb[pg2=="OFF"]), off_hb=mean(hb_arm[pg2=="OFF"]), off_axis=mean(axis_diff[pg2=="OFF"],na.rm=TRUE), off_spingap=mean(spin_gap[pg2=="OFF"],na.rm=TRUE),
                br_share=mean(pgrp=="BR"), os_share=mean(pgrp=="OS"), si_share=mean(pitch_type=="SI"), fc_share=mean(pitch_type=="FC"), fb_high=mean(pz_rel[pg2=="FB"]>0.67), off_low=mean(pz_rel[pg2=="OFF"]<0.33),
                velo_sd=sd(release_speed)), by=.(pitcher,game_year)][swings>=300]
prof <- prof[complete.cases(prof)]
DV <- c("fb_velo","fb_ivb","fb_hb","fb_relz","fbsh","arm","ext","sep","off_ivb","off_hb","off_axis","off_spingap","br_share","os_share","si_share","fc_share","fb_high","off_low","velo_sd")
mD <- lm(as.formula(paste("tscore ~", paste(DV,collapse="+"))), data=prof[game_year<=2025], weights=prof[game_year<=2025]$swings)
pD <- prof[game_year==2026]; pD[, plus_raw := predict(mD, newdata=pD)]
mD25 <- lm(as.formula(paste("tscore ~", paste(DV,collapse="+"))), data=prof[game_year!=2025], weights=prof[game_year!=2025]$swings)
p25D <- prof[game_year==2025][, .(pitcher, plus25=predict(mD25, newdata=prof[game_year==2025]), ts25=tscore)]
qD <- merge(pD, p25D, by="pitcher")
res[["D"]] <- data.table(approach="D", pitch_rmse_lift_pct=NA_real_, pitchers_2026=nrow(pD), r_same_year=round(cor(pD$plus_raw,pD$tscore),3), n_pairs=nrow(qD),
                         r_plus25_vs_ts26=round(cor(qD$plus25,qD$tscore),3), r_ts25_vs_ts26=round(cor(qD$ts25,qD$tscore),3), r_both=round(sqrt(summary(lm(tscore ~ plus25 + ts25, data=qD))$r.squared),3))
cat("\n=== Holdout results (2026). A = pitch physics; B = + fastball-relative arsenal context; C = + location/count; D = pitcher-level arsenal profile (linear) ===\n")
print(rbindlist(res))
cat("\n  r_same_year: tscore+ vs actual tscore, same season (descriptive fit)\n  r_plus25_vs_ts26: 2025 tscore+ (model never saw 2025) predicting 2026 actual; compare to r_ts25_vs_ts26 (raw 2025 tscore predicting 2026)\n  r_both: 2025 tscore+ and 2025 tscore together\n")
cat("\n=== Top features, approach B and C (gain share) ===\n"); print(imps$B[1:12, .(Feature, Gain=round(Gain,3))]); print(imps$C[1:12, .(Feature, Gain=round(Gain,3))])
# ---- tscore+ scale: 100 = league average, 10 points = 1 SD of actual pitcher tscore (>=300 swings, 2024-25)
sdA <- sd(actual[game_year<=2025 & swings>=300]$tscore); muA <- weighted.mean(actual[game_year<=2025 & swings>=300]$tscore, actual[game_year<=2025 & swings>=300]$swings)
te2 <- copy(te)[, `:=`(pA=preds$A, pB=preds$B, pC=preds$C)]
LB <- te2[, .(swings=.N, tscore=mean(y), A=mean(pA), B=mean(pB), C=mean(pC)), by=pitcher][swings>=300]
LB <- merge(LB, pD[, .(pitcher, D=plus_raw)], by="pitcher", all.x=TRUE); LB <- merge(LB, NM, by="pitcher")
for (v in c("tscore","A","B","C","D")) LB[, (paste0(v,"_plus")) := round(100 + 10*(get(v)-muA)/sdA)]
cat(sprintf("\nscale: 100 = %.2f runs/100 swings, 10 pts = %.2f runs/100 swings\n", muA, sdA))
cat("\n=== 2026 tscore+ leaderboard by approach B (>=300 swings): top 15 ===\n")
print(LB[order(-B)][1:15, .(player_name, swings, actual=tscore_plus, A=A_plus, B=B_plus, C=C_plus, D=D_plus)])
cat("\n--- bottom 10 ---\n"); print(LB[order(B)][1:10, .(player_name, swings, actual=tscore_plus, A=A_plus, B=B_plus, C=C_plus, D=D_plus)])
cat("\n--- biggest gaps: actual far above model B (doing something the shape doesn't explain) ---\n"); LB[, gap := tscore_plus-B_plus]; print(LB[swings>=600][order(-gap)][1:8, .(player_name, swings, actual=tscore_plus, B=B_plus, gap)])
cat("--- actual far below model B ---\n"); print(LB[swings>=600][order(gap)][1:8, .(player_name, swings, actual=tscore_plus, B=B_plus, gap)])
fwrite(LB, "data/swing_timing/tscore_plus_2026.csv")
