suppressPackageStartupMessages(library(data.table)); options(width=230)
X <- fread("data/swing_timing/contact_targets_pitchers.csv"); D <- fread("data/swing_timing/physics_fip_2024_2026_era.csv")
M <- merge(D, X, by=c("pitcher","game_year")); M <- M[n_xw>=100]
cat(sprintf("pitcher-seasons with >=100 BIP and all legs: %d\n", nrow(M)))
sc <- c(bscore="bscore_plus", tscore="tscore_plus", xw="mod_xw", ev="mod_ev", hard="mod_hard", barrel="mod_barrel", la="mod_la", gb="mod_gb", hr="mod_hr", air_ev="mod_air_ev", under="mod_under", topped="mod_topped", rv_con="mod_rv_con", rv_pitch="mod_rvp", rv_swing="mod_rvs")
act <- c(xw="act_xw", ev="act_ev", hard="act_hard", barrel="act_barrel", la="act_la", gb="act_gb", hr="act_hr", air_ev="act_air_ev", under="act_under", topped="act_topped", rv_con="act_rv_con", rv_pitch="act_rvp", rv_swing="act_rvs")
a <- copy(M); b <- copy(M)[, game_year := game_year-1L]; q <- merge(a, b, by=c("pitcher","game_year"), suffixes=c("","_n")); q60 <- q[ip_off>=60 & ip_off_n>=60]
cat(sprintf("year pairs: %d (>=60 IP both: %d)\n", nrow(q), nrow(q60)))
sgn <- function(v) if (v %in% c("mod_xw","mod_ev","mod_hard","mod_barrel","mod_hr","mod_air_ev","mod_under","mod_la","act_xw","act_ev","act_hard","act_barrel","act_hr","act_air_ev","act_under","act_la")) -1 else 1  # orient so higher = better for pitcher
R <- rbindlist(lapply(names(sc), function(nm) { v <- sc[[nm]]; s <- sgn(v); x <- s*M[[v]]
  own <- if (nm %in% names(act)) cor(M[[v]], M[[act[[nm]]]], use="complete.obs") else NA_real_
  data.table(target=nm, learnable_r=round(own,3), yoy_modeled=round(cor(q[[v]], q[[paste0(v,"_n")]], use="complete.obs"),3),
             yoy_actual=round(if (nm %in% names(act)) cor(q[[act[[nm]]]], q[[paste0(act[[nm]],"_n")]], use="complete.obs") else NA_real_,3),
             ss_hard=round(cor(x, M$hard_pct, use="complete.obs"),3), ss_xw=round(cor(x, M$xw, use="complete.obs"),3), ss_hr=round(cor(x, M$hr_pct, use="complete.obs"),3), ss_era=round(cor(x, M$era, use="complete.obs"),3),
             nx_hard=round(cor(s*q[[v]], q$hard_pct_n, use="complete.obs"),3), nx_xw=round(cor(s*q[[v]], q$xw_n, use="complete.obs"),3), nx_hr=round(cor(s*q[[v]], q$hr_pct_n, use="complete.obs"),3), nx_era60=round(cor(s*q60[[v]], q60$era_n, use="complete.obs"),3), nx_fip60=round(cor(s*q60[[v]], q60$fip_off_n, use="complete.obs"),3)) }))
cat("\n=== shape-modeled contact targets (all oriented so higher = better for the pitcher). learnable_r = OOF cor(modeled, actual) at pitcher level; ss = same-season r; nx = next-season r ===\n"); print(R)
cat("\n=== incremental value as the contact leg of pFIP: next-season ERA (>=60 IP both) ~ whiff+ + command+ + leg ===\n")
base <- lm(era_n ~ whiff_plus + command_plus, data=q60, weights=q60$ip_off_n); cat(sprintf("  whiff+ + command+ only: R %.3f\n", sqrt(summary(base)$r.squared)))
INC <- rbindlist(lapply(names(sc), function(nm) { v <- sc[[nm]]; qq <- q60[is.finite(get(v))]; m <- lm(as.formula(paste("era_n ~ whiff_plus + command_plus +", v)), data=qq, weights=qq$ip_off_n); mf <- lm(as.formula(paste("fip_off_n ~ whiff_plus + command_plus +", v)), data=qq, weights=qq$ip_off_n); mh <- lm(as.formula(paste("hr_pct_n ~ whiff_plus + command_plus +", v)), data=qq, weights=qq$ip_off_n)
  data.table(leg=nm, R_nextERA=round(sqrt(summary(m)$r.squared),3), t_leg=round(summary(m)$coefficients[4,3],1), R_nextFIP=round(sqrt(summary(mf)$r.squared),3), t_leg_fip=round(summary(mf)$coefficients[4,3],1), R_nextHR=round(sqrt(summary(mh)$r.squared),3), t_leg_hr=round(summary(mh)$coefficients[4,3],1)) }))
print(INC[order(-R_nextERA)])
cat("\n=== do two contact legs beat one? next ERA (>=60 IP) ~ whiff+ + command+ + bscore+ + X ===\n")
for (v in c("mod_xw","mod_hr","mod_gb","mod_air_ev","mod_barrel","mod_rvp")) { qq <- q60[is.finite(get(v))]; m <- lm(as.formula(paste("era_n ~ whiff_plus + command_plus + bscore_plus +", v)), data=qq, weights=qq$ip_off_n); co <- summary(m)$coefficients; cat(sprintf("  + %-10s R %.3f | t bscore+ %+.1f | t %s %+.1f\n", v, sqrt(summary(m)$r.squared), co["bscore_plus",3], v, co[v,3])) }
cat("\ncorrelations among modeled contact scores:\n"); print(round(cor(M[, .(bscore_plus, mod_xw, mod_ev, mod_hard, mod_barrel, mod_gb, mod_hr, mod_air_ev, mod_rvp)], use="complete.obs"),2))
