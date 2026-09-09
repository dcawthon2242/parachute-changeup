#!/usr/bin/env Rscript
# bscore: extend tscore (depth axis only) to all three barrel-positioning axes.
#   y-axis  depth      : intercept_ball_minus_batter_pos_y (late/early)      -> tdev
#   x-axis  horizontal : intercept_ball_minus_batter_pos_x (tied-up/flail)   -> xdev
#   z-axis  vertical   : no published per-pitch field; proxy = swing-plane vs pitch-plane mismatch,
#                        attack_angle + VAA (0 = bat travelling along the pitch line)          -> zdev
# Each axis is residualized on location/pitch group/velo/stand and batter-centered, exactly like tdev.
# Scores: per-axis (univariate RV surface) and joint (lightgbm on the three axes + pitch group + count).
suppressPackageStartupMessages({ library(data.table); library(splines); library(lightgbm) }); set.seed(41); options(width=230)
cols <- c("game_year","game_type","player_name","pitcher","batter","stand","p_throws","pitch_type","description","bb_type","balls","strikes","plate_x","plate_z","sz_top","sz_bot",
          "release_speed","launch_speed","launch_angle","launch_speed_angle","estimated_woba_using_speedangle","delta_run_exp","vy0","vz0","ay","az",
          "intercept_ball_minus_batter_pos_x_inches","intercept_ball_minus_batter_pos_y_inches","attack_angle","bat_speed","swing_length","swing_path_tilt")
dt <- rbindlist(lapply(2024:2026, function(yr) fread(sprintf("data/statcast_%d/statcast_%d_all.csv",yr,yr), showProgress=FALSE, select=cols)))
setnames(dt, c("intercept_ball_minus_batter_pos_x_inches","intercept_ball_minus_batter_pos_y_inches"), c("ix","depth"))
dt <- dt[game_type=="R" & pitch_type!="" & balls<=3 & strikes<=2 & is.finite(delta_run_exp) & is.finite(release_speed) & is.finite(plate_x) & is.finite(plate_z)]
FB<-c("FF","SI","FC"); BR<-c("SL","ST","CU","KC","SV","CS"); OS<-c("CH","FS","FO")
dt[, pgrp := fifelse(pitch_type %in% FB,"FB", fifelse(pitch_type %in% BR,"BR", fifelse(pitch_type %in% OS,"OS",NA_character_)))]; dt <- dt[!is.na(pgrp)]; dt[, pg2 := fifelse(pgrp=="FB","FB","OFF")]
dt[, px_bat := fifelse(stand=="R",-plate_x,plate_x)]; dt[, pz_rel := (plate_z-sz_bot)/pmax(sz_top-sz_bot,0.1)]
dt[, tt := (-vy0 - sqrt(pmax(vy0^2 - 2*ay*(50-17/12),0)))/ay]; dt[, vaa := atan2(vz0+az*tt, abs(vy0+ay*tt))*180/pi]
SWING <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip","hit_into_play")
sw <- dt[description %in% SWING & is.finite(depth) & is.finite(ix) & is.finite(attack_angle) & is.finite(vaa)]
sw[, ix_bat := fifelse(stand=="R", ix, -ix)]     # horizontal intercept toward the batter's pull side? keep sign: positive = ball farther from batter (toward plate outside)
sw[, zraw := attack_angle + vaa]                 # + = bat steeper than pitch line (under the ball), - = flatter (over the ball)
cat(sprintf("swings with all three axes: %s (2024 %.0f%%, 2025 %.0f%%, 2026 %.0f%% of swings)\n", format(nrow(sw),big.mark=","),
  100*nrow(sw[game_year==2024])/nrow(dt[game_year==2024 & description %in% SWING]), 100*nrow(sw[game_year==2025])/nrow(dt[game_year==2025 & description %in% SWING]), 100*nrow(sw[game_year==2026])/nrow(dt[game_year==2026 & description %in% SWING])))
# ---- residualize each axis on location/pitch group/velo/stand, then batter-center
resid_axis <- function(v) { f <- lm(as.formula(paste(v, "~ ns(px_bat,5)*pgrp + ns(pz_rel,5) + ns(release_speed,4) + stand")), data=sw); residuals(f) }
sw[, r_t := resid_axis("depth")]; sw[, r_x := resid_axis("ix_bat")]; sw[, r_z := resid_axis("zraw")]
sw[, nb := .N, by=batter]; sw <- sw[nb>=200]
sw[, tdev := r_t-mean(r_t), by=batter]; sw[, xdev := r_x-mean(r_x), by=batter]; sw[, zdev := r_z-mean(r_z), by=batter]
sw[, cnt := paste0(balls,"-",strikes)]; sw[, bat_rv := delta_run_exp]; sw[, pit_rv := -delta_run_exp]
sw[, whiff := description %in% c("swinging_strike","swinging_strike_blocked","foul_tip")]
sw[, isbip := description=="hit_into_play" & bb_type!="" & is.finite(launch_speed) & is.finite(launch_speed_angle)]
# ---- validate the z proxy: does it track over/under contact?
cat("\n=== z-axis proxy validation (balls in play): launch angle and contact type by zdev quintile ===\n")
b <- sw[(isbip)]; b[, zq := cut(zdev, quantile(zdev, seq(0,1,.2)), include.lowest=TRUE, labels=c("Q1 flat/over","Q2","Q3","Q4","Q5 steep/under"))]
print(b[, .(n=.N, launch_angle=round(mean(launch_angle),1), gb=round(mean(bb_type=="ground_ball"),3), popup=round(mean(bb_type=="popup"),3), topped=round(mean(launch_speed_angle %in% 1:2),3), under=round(mean(launch_speed_angle==3),3), barrel=round(mean(launch_speed_angle==6),3), ev=round(mean(launch_speed),1), xw=round(mean(estimated_woba_using_speedangle,na.rm=TRUE),3)), by=zq][order(zq)])
cat("cor(zdev, launch_angle) on BIP =", round(cor(b$zdev, b$launch_angle),3), " | cor(xdev, launch_speed) =", round(cor(b$xdev,b$launch_speed),3), " | cor(tdev, launch_speed) =", round(cor(b$tdev,b$launch_speed),3), "\n")
cat("\n=== x-axis: what horizontal intercept does (BIP) ===\n"); b[, xq := cut(xdev, quantile(xdev, seq(0,1,.2)), include.lowest=TRUE, labels=c("Q1","Q2","Q3","Q4","Q5"))]
print(b[, .(n=.N, ev=round(mean(launch_speed),1), barrel=round(mean(launch_speed_angle==6),3), weak=round(mean(launch_speed_angle==1),3), xw=round(mean(estimated_woba_using_speedangle,na.rm=TRUE),3), pull=round(mean(fifelse(stand=="R", launch_angle*0 + 1, 1))*0,1)), by=xq][order(xq)][, -"pull"])
cat("whiff rate by |xdev| quintile: "); sw[, axq := cut(abs(xdev), quantile(abs(xdev), seq(0,1,.2)), include.lowest=TRUE, labels=1:5)]; print(sw[, .(whiff=round(mean(whiff),3)), by=axq][order(axq)])
# ---- univariate surfaces (same construction as tscore) -> per-axis scores
axis_score <- function(dev, nm) {
  f <- lm(as.formula(paste("bat_rv ~ ns(", dev, ",6)*pgrp + factor(cnt)")), data=sw); p <- predict(f, newdata=sw)
  pb <- p - ave(p, sw$batter); pe <- pb - ave(pb, sw$pg2); sw[, (nm) := -100*pe]
  cat(sprintf("  %s surface R2 = %.4f\n", nm, summary(f)$r.squared)) }
cat("\n=== run-value surfaces ===\n"); axis_score("tdev","ts"); axis_score("xdev","xs"); axis_score("zdev","zs")
# ---- joint surface: lightgbm on the three axes x pitch group x count (smooth-ish: few leaves, big min_data)
X <- as.matrix(sw[, .(tdev, xdev, zdev, pgi=match(pgrp,c("FB","BR","OS")), cnti=balls*3+strikes)])
m <- lgb.train(params=list(objective="regression", learning_rate=0.05, num_leaves=15, min_data_in_leaf=2000, feature_fraction=1, lambda_l2=10, verbose=-1), data=lgb.Dataset(X, label=sw$bat_rv), nrounds=400)
pj <- predict(m, X); cat(sprintf("  joint 3-axis surface R2 = %.4f  (depth-only surface above for comparison)\n", 1-sum((sw$bat_rv-pj)^2)/sum((sw$bat_rv-mean(sw$bat_rv))^2)))
pjb <- pj - ave(pj, sw$batter); pje <- pjb - ave(pjb, sw$pg2); sw[, bs := -100*pje]
cat("  importance:"); print(lgb.importance(m)[, .(Feature, Gain=round(Gain,3))])
# ---- pitcher-season scores
P <- sw[, .(swings=.N, tscore=mean(ts), xscore=mean(xs), zscore=mean(zs), bscore=mean(bs), tdev=mean(tdev), xdev=mean(xdev), zdev=mean(zdev),
            whiff=100*mean(whiff), hard=100*sum(isbip & launch_speed>=95)/pmax(sum(isbip),1), ev=mean(launch_speed[isbip]), xw=mean(estimated_woba_using_speedangle[isbip],na.rm=TRUE),
            barrel=100*sum(isbip & launch_speed_angle==6)/pmax(sum(isbip),1), rv100=100*mean(pit_rv), n_bip=sum(isbip)), by=.(pitcher,game_year)][swings>=300]
NM <- unique(dt[, .(pitcher,player_name)], by="pitcher"); P <- merge(P, NM, by="pitcher")
cat(sprintf("\npitcher-seasons >=300 swings: %d\n", nrow(P)))
cat("\n=== score SDs and inter-correlations ===\n"); print(round(cor(P[, .(tscore,xscore,zscore,bscore)]),2)); print(P[, .(sd_t=round(sd(tscore),3), sd_x=round(sd(xscore),3), sd_z=round(sd(zscore),3), sd_b=round(sd(bscore),3))])
cat("\n=== same-season correlation with outcomes ===\n")
print(rbindlist(lapply(c("tscore","xscore","zscore","bscore"), function(s) data.table(score=s, whiff=round(cor(P[[s]],P$whiff),3), hard=round(cor(P[[s]],P$hard),3), ev=round(cor(P[[s]],P$ev),3), barrel=round(cor(P[[s]],P$barrel),3), xw=round(cor(P[[s]],P$xw,use="complete.obs"),3), rv100=round(cor(P[[s]],P$rv100),3)))))
cat("\n=== year-over-year reliability (>=300 swings both) ===\n")
a <- copy(P); b2 <- copy(P)[, game_year := game_year-1L]; q <- merge(a, b2, by=c("pitcher","game_year"), suffixes=c("","_n"))
print(rbindlist(lapply(c("tscore","xscore","zscore","bscore","tdev","xdev","zdev","whiff","hard","ev","xw"), function(s) data.table(metric=s, r_yoy=round(cor(q[[s]],q[[paste0(s,"_n")]],use="complete.obs"),3), n=nrow(q)))))
cat("\n=== next-season prediction of contact quality: R2 ===\n")
for (y in c("hard","ev","xw")) { yn <- paste0(y,"_n"); f <- function(rhs) summary(lm(as.formula(paste(yn,"~",rhs)), data=q))$r.squared
  cat(sprintf("  %-4s next | own %.3f | own+tscore %.3f | own+bscore %.3f | own+t+x+z %.3f | tscore alone %.3f | bscore alone %.3f | xscore alone %.3f | zscore alone %.3f\n", y, f(y), f(paste(y,"+tscore")), f(paste(y,"+bscore")), f(paste(y,"+tscore+xscore+zscore")), f("tscore"), f("bscore"), f("xscore"), f("zscore"))) }
cat("\n=== 2026 leaders, bscore (>=600 swings) ===\n")
print(P[game_year==2026 & swings>=600][order(-bscore)][1:15, .(player_name, swings, bscore=round(bscore,2), tscore=round(tscore,2), xscore=round(xscore,2), zscore=round(zscore,2), whiff=round(whiff,1), hard=round(hard,1), xw=round(xw,3))])
cat("\n--- biggest disagreements: high bscore, low tscore (barrel-missers who are not timing-disruptors) ---\n"); P[, gap := scale(bscore)-scale(tscore)]
print(P[game_year==2026 & swings>=600][order(-gap)][1:8, .(player_name, bscore=round(bscore,2), tscore=round(tscore,2), xscore=round(xscore,2), zscore=round(zscore,2), hard=round(hard,1), xw=round(xw,3))])
fwrite(P, "data/swing_timing/barrel_scores_2024_2026.csv")
