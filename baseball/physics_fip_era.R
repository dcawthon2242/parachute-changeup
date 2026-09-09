suppressPackageStartupMessages({ library(data.table); library(jsonlite) }); options(width=200)
E <- rbindlist(lapply(2024:2026, function(y) { s <- fromJSON(sprintf("data/mlb_stats/pitching_%d.json", y))$stats$splits[[1]]
  data.table(pitcher=s$player$id, game_year=y, era=as.numeric(s$stat$era), ip_off=as.numeric(s$stat$outs)/3, er=s$stat$earnedRuns, r_off=s$stat$runs, k_off=s$stat$strikeOuts, bb_off=s$stat$baseOnBalls, hbp_off=s$stat$hitByPitch, hr_off=s$stat$homeRuns, bf=s$stat$battersFaced, gs=s$stat$gamesStarted, g=s$stat$gamesPlayed) }))
E[, ra9_off := 9*r_off/ip_off]
# official-scale FIP and xFIP (ERA scale) per season
E[, fip_raw := (13*hr_off + 3*(bb_off+hbp_off) - 2*k_off)/ip_off]; E[, c_fip := weighted.mean(era, ip_off) - weighted.mean(fip_raw, ip_off), by=game_year]; E[, fip_off := fip_raw + c_fip]
D <- fread("data/swing_timing/physics_fip_2024_2026.csv"); D <- merge(D, E, by=c("pitcher","game_year"))
# xFIP needs FB count: from statcast side (fb) with league HR/FB
D[, xfip_off := (13*fb*weighted.mean(hr_off/pmax(fb,1), fb) + 3*(bb_off+hbp_off) - 2*k_off)/ip_off + c_fip, by=game_year]
cat(sprintf("pitcher-seasons %d | cor(my RA9 from Statcast, official RA9) = %.3f | cor(my FIP, official FIP) = %.3f | cor(ERA, RA9) = %.3f\n", nrow(D), cor(D$ra9, D$ra9_off), cor(D$fip, D$fip_off), cor(D$era, D$ra9_off)))
# refit pFIP on the ERA scale (legs -> ERA) so units are comparable; also legs -> FIP_off
tr <- D[game_year<=2025]
mE <- lm(era ~ whiff_plus + command_plus + bscore_plus, data=tr, weights=tr$ip_off); mF <- lm(fip_off ~ whiff_plus + command_plus + bscore_plus, data=tr, weights=tr$ip_off)
co <- summary(mE)$coefficients; cat(sprintf("\nlegs -> ERA (2024-25): R2 %.3f | per +10: whiff+ %+.3f (t %+.1f) command+ %+.3f (t %+.1f) bscore+ %+.3f (t %+.1f)\n", summary(mE)$r.squared, 10*co[2,1],co[2,3],10*co[3,1],co[3,3],10*co[4,1],co[4,3]))
D[, pERA := predict(mE, D)]; D[, pFIP := predict(mF, D)]
a <- copy(D); b <- copy(D)[, game_year := game_year-1L]; q <- merge(a, b, by=c("pitcher","game_year"), suffixes=c("","_n"))
for (minip in c(30, 60, 100)) { qq <- q[ip_off>=minip & ip_off_n>=minip]
  cat(sprintf("\n=== next-season ERA, both seasons >= %d IP (n=%d) ===\n", minip, nrow(qq)))
  R <- rbindlist(lapply(c("era","ra9_off","fip_off","xfip_off","pFIP","pERA","pfip_reg"), function(v) data.table(predictor=v, r_next_ERA=round(cor(qq[[v]], qq$era_n),3), r_next_RA9=round(cor(qq[[v]], qq$ra9_off_n),3), r_next_FIP=round(cor(qq[[v]], qq$fip_off_n),3), yoy_self=round(cor(qq[[v]], qq[[paste0(v,"_n")]]),3))))
  print(R)
  f <- function(rhs) sqrt(summary(lm(as.formula(paste("era_n ~", rhs)), data=qq, weights=qq$ip_off_n))$r.squared)
  cat(sprintf("  combos -> next ERA (R): FIP+xFIP %.3f | FIP+pFIP %.3f | xFIP+pFIP %.3f | ERA+FIP+xFIP %.3f | ERA+FIP+xFIP+pFIP %.3f | legs raw (whiff+,command+,bscore+) %.3f\n", f("fip_off+xfip_off"), f("fip_off+pFIP"), f("xfip_off+pFIP"), f("era+fip_off+xfip_off"), f("era+fip_off+xfip_off+pFIP"), f("whiff_plus+command_plus+bscore_plus")))
  m <- lm(era_n ~ era + fip_off + xfip_off + pFIP, data=qq, weights=qq$ip_off_n); co <- summary(m)$coefficients; cat("  t-stats in the full model: "); cat(sprintf("%s %+.1f  ", rownames(co)[-1], co[-1,3])); cat("\n")
  # RMSE of next-year ERA
  for (v in c("era","fip_off","xfip_off","pFIP","pERA")) { z <- lm(as.formula(paste("era_n ~", v)), data=qq, weights=qq$ip_off_n); cat(sprintf("  RMSE next ERA via %-8s %.3f\n", v, sqrt(weighted.mean(resid(z)^2, qq$ip_off_n)))) }
}
# starters only
qq <- q[gs>=15 & gs_n>=15]; cat(sprintf("\n=== starters (>=15 GS both seasons, n=%d) -> next ERA ===\n", nrow(qq)))
for (v in c("era","fip_off","xfip_off","pFIP","pERA")) cat(sprintf("  %-8s r = %+.3f\n", v, cor(qq[[v]], qq$era_n)))
fwrite(D, "data/swing_timing/physics_fip_2024_2026_era.csv")
