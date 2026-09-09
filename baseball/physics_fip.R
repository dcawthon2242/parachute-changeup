#!/usr/bin/env Rscript
# pFIP: a physics-modeled FIP from three shape/command legs -- whiff+ (K), OpenCommand miss (BB), bscore+ (contact/HR).
suppressPackageStartupMessages(library(data.table)); options(width=220)
# ---- outcomes per pitcher-season from Statcast events
cols <- c("game_year","game_type","pitcher","player_name","events","description","bat_score","post_bat_score","woba_value","woba_denom","launch_speed","estimated_woba_using_speedangle","bb_type")
dt <- rbindlist(lapply(2024:2026, function(yr) fread(sprintf("data/statcast_%d/statcast_%d_all.csv",yr,yr), select=cols, showProgress=FALSE)))[game_type=="R"]
dt[, runs := post_bat_score - bat_score]
K <- c("strikeout","strikeout_double_play"); BB <- c("walk"); OUT_EV <- c("field_out","strikeout","strikeout_double_play","force_out","grounded_into_double_play","double_play","sac_fly","sac_bunt","fielders_choice_out","triple_play","sac_fly_double_play","sac_bunt_double_play","other_out","caught_stealing_2b","caught_stealing_3b","caught_stealing_home","pickoff_1b","pickoff_2b","pickoff_3b","pickoff_caught_stealing_2b","pickoff_caught_stealing_3b","pickoff_caught_stealing_home","fielders_choice")
dt[, outs := fifelse(events %in% c("grounded_into_double_play","double_play","strikeout_double_play","sac_fly_double_play","sac_bunt_double_play"), 2L, fifelse(events=="triple_play", 3L, fifelse(events %in% OUT_EV, 1L, 0L)))]
O <- dt[, .(pa=sum(events!=""), k=sum(events %in% K), bb=sum(events %in% BB), hbp=sum(events=="hit_by_pitch"), hr=sum(events=="home_run"), ibb=0, outs=sum(outs), runs=sum(runs, na.rm=TRUE),
            bip=sum(description=="hit_into_play" & is.finite(launch_speed)), hard=sum(description=="hit_into_play" & launch_speed>=95, na.rm=TRUE), xw=mean(estimated_woba_using_speedangle[description=="hit_into_play"], na.rm=TRUE),
            fb=sum(bb_type %in% c("fly_ball","popup")), woba=sum(woba_value,na.rm=TRUE)/sum(woba_denom,na.rm=TRUE)), by=.(pitcher,game_year)]
O[, ip := outs/3]; O[, `:=`(k_pct=k/pa, bb_pct=(bb+hbp)/pa, hr_pct=hr/pa, hr_fb=hr/pmax(fb,1), hard_pct=hard/pmax(bip,1), ra9=9*runs/ip)]
# league FIP constant per season so FIP is on the RA9 scale here (we lack ERA; use RA9)
O[, fip_raw := (13*hr + 3*(bb+hbp) - 2*k)/ip]; O[, cFIP := weighted.mean(ra9, ip) - weighted.mean(fip_raw, ip), by=game_year]; O[, fip := fip_raw + cFIP]
O[, xfip_raw := (13*(fb*weighted.mean(hr_fb, fb)) + 3*(bb+hbp) - 2*k)/ip, by=game_year]; O[, xfip := xfip_raw + cFIP]
NM <- unique(dt[, .(pitcher, player_name)], by="pitcher")
# ---- legs
W <- fread("data/swing_timing/whiff_plus_pitchers.csv")[, .(pitcher, game_year, swings, whiff, whiff_plus_raw)]
BS <- readRDS("data/swing_timing/bscore_plus_oof_swings.rds")[, .(bswings=.N, bscore=mean(y), bplus_raw=mean(pred)), by=.(pitcher,game_year)]
TS <- readRDS("data/swing_timing/tscore_plus_oof_swings.rds")[, .(tscore=mean(y), tplus_raw=mean(pred)), by=.(pitcher,game_year)]
CM <- fread("data/swing_timing/command_pitchers_2024_2026.csv")[, .(pitcher, game_year, n_cmd, miss_mean, miss_med, target_zone, big_miss)]
D <- Reduce(function(a,b) merge(a,b,by=c("pitcher","game_year")), list(O, W, BS, TS, CM)); D <- merge(D, NM, by="pitcher")
D <- D[pa>=200 & swings>=300 & n_cmd>=300 & is.finite(ra9)]
cat(sprintf("pitcher-seasons with all legs: %d (2024 %d, 2025 %d, 2026 %d)\n", nrow(D), sum(D$game_year==2024), sum(D$game_year==2025), sum(D$game_year==2026)))
# plus scales: 100 = league avg, 10 = 1 SD (2024-25 pitcher-seasons); command flipped so higher = better
ref <- D[game_year<=2025]; pl <- function(x, r, flip=FALSE) { z <- (x-mean(r))/sd(r); if (flip) z <- -z; 100+10*z }
D[, `:=`(whiff_plus=pl(whiff_plus_raw, ref$whiff_plus_raw), bscore_plus=pl(bplus_raw, ref$bplus_raw), tscore_plus=pl(tplus_raw, ref$tplus_raw), command_plus=pl(miss_mean, ref$miss_mean, flip=TRUE))]
cat("\n=== leg reliability (yoy) and leg -> component validity ===\n")
a <- copy(D); b <- copy(D)[, game_year := game_year-1L]; q <- merge(a, b, by=c("pitcher","game_year"), suffixes=c("","_n")); cat("pairs:", nrow(q), "\n")
for (v in c("whiff_plus","command_plus","bscore_plus","tscore_plus","k_pct","bb_pct","hr_pct","hard_pct","fip","xfip","ra9")) cat(sprintf("  yoy %-13s %.3f\n", v, cor(q[[v]], q[[paste0(v,"_n")]])))
cat("\nsame-season correlations (legs x components):\n")
M <- D[, .(whiff_plus, command_plus, bscore_plus, tscore_plus, k_pct, bb_pct, hr_pct, hard_pct, xw, fip, ra9)]; print(round(cor(M, use="complete.obs")[1:4, 5:11], 3))
cat("\nnext-season: each leg predicting next-year component (r):\n")
for (y in c("k_pct","bb_pct","hr_pct","hard_pct","fip","ra9")) cat(sprintf("  next %-8s | whiff+ %+.3f  command+ %+.3f  bscore+ %+.3f | own %+.3f\n", y, cor(q$whiff_plus, q[[paste0(y,"_n")]]), cor(q$command_plus, q[[paste0(y,"_n")]]), cor(q$bscore_plus, q[[paste0(y,"_n")]]), cor(q[[y]], q[[paste0(y,"_n")]])))
# ---- (1) literal pFIP: model each FIP component from the legs, plug into the FIP formula
tr <- D[game_year<=2025]
mK <- lm(k_pct ~ whiff_plus + command_plus + bscore_plus, data=tr, weights=tr$pa); mB <- lm(bb_pct ~ whiff_plus + command_plus + bscore_plus, data=tr, weights=tr$pa); mH <- lm(hr_pct ~ whiff_plus + command_plus + bscore_plus, data=tr, weights=tr$pa)
cat("\n=== component models (2024-25), standardized effect per +10 leg points ===\n")
for (nm in c("mK","mB","mH")) { m <- get(nm); co <- summary(m)$coefficients; cat(sprintf("  %s R2 %.3f | whiff+ %+.4f (t %+.1f) command+ %+.4f (t %+.1f) bscore+ %+.4f (t %+.1f)\n", c(mK="K%  ",mB="BB% ",mH="HR% ")[nm], summary(m)$r.squared, 10*co["whiff_plus",1], co["whiff_plus",3], 10*co["command_plus",1], co["command_plus",3], 10*co["bscore_plus",1], co["bscore_plus",3])) }
D[, `:=`(k_hat=predict(mK, D), bb_hat=predict(mB, D), hr_hat=predict(mH, D))]
pa_per_ip <- weighted.mean(tr$pa/tr$ip, tr$ip)
D[, pfip_lit := (13*hr_hat + 3*bb_hat - 2*k_hat)*pa_per_ip + weighted.mean(tr$cFIP, tr$ip)]
# ---- (2) regression pFIP: FIP ~ legs directly (and RA9 ~ legs)
mF <- lm(fip ~ whiff_plus + command_plus + bscore_plus, data=tr, weights=tr$ip); mR <- lm(ra9 ~ whiff_plus + command_plus + bscore_plus, data=tr, weights=tr$ip)
co <- summary(mF)$coefficients; cat(sprintf("\n=== regression pFIP (FIP ~ legs, 2024-25): R2 %.3f | per +10: whiff+ %+.3f (t %+.1f) command+ %+.3f (t %+.1f) bscore+ %+.3f (t %+.1f)\n", summary(mF)$r.squared, 10*co[2,1], co[2,3], 10*co[3,1], co[3,3], 10*co[4,1], co[4,3]))
co <- summary(mR)$coefficients; cat(sprintf("    RA9 ~ legs: R2 %.3f | per +10: whiff+ %+.3f (t %+.1f) command+ %+.3f (t %+.1f) bscore+ %+.3f (t %+.1f)\n", summary(mR)$r.squared, 10*co[2,1], co[2,3], 10*co[3,1], co[3,3], 10*co[4,1], co[4,3]))
D[, pfip_reg := predict(mF, D)]; D[, pra9 := predict(mR, D)]
# with tscore+ instead of / in addition to bscore+
mF2 <- lm(fip ~ whiff_plus + command_plus + tscore_plus, data=tr, weights=tr$ip); mF3 <- lm(fip ~ whiff_plus + command_plus + bscore_plus + tscore_plus, data=tr, weights=tr$ip)
cat(sprintf("    FIP ~ whiff+ command+ tscore+: R2 %.3f | + both: %.3f\n", summary(mF2)$r.squared, summary(mF3)$r.squared))
# ---- validation: 2026 holdout same-season, and next-season prediction
h <- D[game_year==2026]
cat(sprintf("\n=== 2026 holdout (n=%d), same-season r with actual FIP: pFIP literal %.3f | pFIP regression %.3f | xFIP %.3f\n", nrow(h), cor(h$pfip_lit,h$fip), cor(h$pfip_reg,h$fip), cor(h$xfip,h$fip)))
cat(sprintf("    same-season r with actual RA9: pFIP literal %.3f | pFIP reg %.3f | pRA9 %.3f | FIP %.3f | xFIP %.3f\n", cor(h$pfip_lit,h$ra9), cor(h$pfip_reg,h$ra9), cor(h$pra9,h$ra9), cor(h$fip,h$ra9), cor(h$xfip,h$ra9)))
q <- merge(copy(D), copy(D)[, game_year := game_year-1L], by=c("pitcher","game_year"), suffixes=c("","_n"))
cat(sprintf("\n=== next-season prediction (n=%d pairs), r with NEXT-year FIP / RA9 ===\n", nrow(q)))
for (v in c("fip","xfip","ra9","pfip_lit","pfip_reg","pra9")) cat(sprintf("  %-9s -> next FIP %+.3f | next RA9 %+.3f\n", v, cor(q[[v]], q$fip_n), cor(q[[v]], q$ra9_n)))
cat(sprintf("  FIP + pFIP_reg -> next FIP R %.3f (FIP alone %.3f) | -> next RA9 R %.3f (FIP alone %.3f, xFIP alone %.3f)\n", sqrt(summary(lm(fip_n ~ fip + pfip_reg, q))$r.squared), cor(q$fip,q$fip_n), sqrt(summary(lm(ra9_n ~ fip + pfip_reg, q))$r.squared), cor(q$fip,q$ra9_n), cor(q$xfip,q$ra9_n)))
cat(sprintf("  yoy: pFIP_reg %.3f | pFIP_lit %.3f | FIP %.3f | xFIP %.3f | RA9 %.3f\n", cor(q$pfip_reg,q$pfip_reg_n), cor(q$pfip_lit,q$pfip_lit_n), cor(q$fip,q$fip_n), cor(q$xfip,q$xfip_n), cor(q$ra9,q$ra9_n)))
# residual of pFIP: what is it correlated with
D[, res := fip - pfip_reg]; cat(sprintf("\npFIP residual (FIP - pFIP): SD %.2f vs FIP SD %.2f | yoy r %.3f\n", sd(D$res), sd(D$fip), cor(q$fip-q$pfip_reg, q$fip_n-q$pfip_reg_n)))
for (v in c("k_pct","bb_pct","hr_pct","hr_fb","hard_pct","xw","target_zone","big_miss","tscore_plus")) cat(sprintf("  cor(residual, %-11s) = %+.3f\n", v, cor(D$res, D[[v]], use="complete.obs")))
# ---- 2026 leaderboard
L <- h[ip>=50][order(pfip_reg)]
cat("\n=== 2026 pFIP (regression) best 12, >=50 IP ===\n"); print(L[1:12, .(player_name, ip=round(ip), pFIP=round(pfip_reg,2), FIP=round(fip,2), xFIP=round(xfip,2), RA9=round(ra9,2), `whiff+`=round(whiff_plus), `command+`=round(command_plus), `bscore+`=round(bscore_plus), K=round(100*k_pct,1), BB=round(100*bb_pct,1), HR=round(100*hr_pct,1))])
cat("\n--- worst 12 ---\n"); print(L[order(-pfip_reg)][1:12, .(player_name, ip=round(ip), pFIP=round(pfip_reg,2), FIP=round(fip,2), xFIP=round(xfip,2), RA9=round(ra9,2), `whiff+`=round(whiff_plus), `command+`=round(command_plus), `bscore+`=round(bscore_plus), K=round(100*k_pct,1), BB=round(100*bb_pct,1), HR=round(100*hr_pct,1))])
cat("\n--- biggest FIP over-performers vs pFIP (FIP much better than shape) ---\n"); print(L[order(fip-pfip_reg)][1:8, .(player_name, ip=round(ip), pFIP=round(pfip_reg,2), FIP=round(fip,2), gap=round(fip-pfip_reg,2), `whiff+`=round(whiff_plus), `command+`=round(command_plus), `bscore+`=round(bscore_plus))])
cat("--- biggest under-performers ---\n"); print(L[order(-(fip-pfip_reg))][1:8, .(player_name, ip=round(ip), pFIP=round(pfip_reg,2), FIP=round(fip,2), gap=round(fip-pfip_reg,2), `whiff+`=round(whiff_plus), `command+`=round(command_plus), `bscore+`=round(bscore_plus))])
fwrite(D, "data/swing_timing/physics_fip_2024_2026.csv"); cat("\nDONE\n")
