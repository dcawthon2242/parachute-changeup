#!/usr/bin/env Rscript
# High velo-gap: chase and chase-channel xRV, by secondary family.

suppressPackageStartupMessages({ library(data.table); library(splines) })
options(width = 215)

P <- fread("data/statcast_2026/velo_gap_tunnel_pairs.csv")
cols <- c("game_year","game_type","pitcher","pitch_type","description",
          "balls","strikes","plate_x","plate_z","sz_top","sz_bot",
          "estimated_woba_using_speedangle","delta_run_exp")
dt <- rbindlist(lapply(c(2025, 2026), function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress=FALSE, select=cols)))
setnames(dt, "estimated_woba_using_speedangle", "xw")
dt <- dt[game_type=="R" & pitch_type!="" & balls<=3 & strikes<=2 &
         is.finite(delta_run_exp) & is.finite(plate_x)]
dt[, in_zone := abs(plate_x)<=0.83 & plate_z>=sz_bot & plate_z<=sz_top]
dt[, bat_rv := delta_run_exp]
WH <- c("swinging_strike","swinging_strike_blocked","foul_tip")
SW <- c(WH,"foul","hit_into_play")
dt[, swung := description %in% SW]
mb <- lm(bat_rv ~ ns(xw,5)+factor(paste0(balls,"-",strikes)),
         data=dt[description=="hit_into_play" & is.finite(xw)])
dt[, xrv_p := bat_rv]
dt[description=="hit_into_play" & is.finite(xw),
   xrv_p := predict(mb, newdata=dt[description=="hit_into_play" & is.finite(xw)])]
dt <- dt[!(description=="hit_into_play" & !is.finite(xw))]

sec <- merge(dt, P[, .(pitcher, pitch_type, velo_gap, fb_velo, n, xrv)],
             by=c("pitcher","pitch_type"))
sec[, fam := fifelse(pitch_type=="SL","SL",
              fifelse(pitch_type=="ST","ST",
              fifelse(pitch_type %in% c("CU","KC","SV"),"CU",
              fifelse(pitch_type %in% c("CH","FS","FO"),"CH","OT"))))]
sec <- sec[fam!="OT"]

C <- sec[, .(
  np=.N,
  chase=100*mean(swung[!in_zone]),
  swing_z=100*mean(swung[in_zone]),
  zone=100*mean(in_zone),
  ball=100*mean(description=="ball"),
  xrv=-100*mean(xrv_p),
  ch_chase=-100*mean(fifelse(!in_zone & swung, xrv_p, 0)),
  ch_ooz_take=-100*mean(fifelse(!in_zone & !swung, xrv_p, 0))
), by=.(pitcher, pitch_type, fam, velo_gap, fb_velo)]

C[, g3 := cut(velo_gap, quantile(velo_gap, 0:3/3), include.lowest=TRUE,
              labels=c("small","mid","big")), by=fam]

cat("############ Chase and chase-channel xRV by family × gap tercile ############\n\n")
print(C[, .(
  pairs=.N, gap=round(mean(velo_gap),1), fb=round(mean(fb_velo),1),
  chase=round(mean(chase),1), swing_z=round(mean(swing_z),1),
  zone=round(mean(zone),1), ball=round(mean(ball),1),
  xrv=round(mean(xrv),3),
  ch_chase=round(mean(ch_chase),3),
  ch_ooz_take=round(mean(ch_ooz_take),3)
), by=.(fam, g3)][order(fam, g3)], row.names=FALSE)

cat("\n  Big minus small, by family:\n")
for (f in c("SL","ST","CU","CH")) {
  a <- C[fam==f & g3=="big"]; b <- C[fam==f & g3=="small"]
  if (nrow(a)<8 || nrow(b)<8) next
  cat(sprintf("\n  %s  n=%d/%d\n", f, nrow(a), nrow(b)))
  for (v in c("chase","swing_z","ball","xrv","ch_chase","ch_ooz_take")) {
    t <- t.test(a[[v]], b[[v]])
    cat(sprintf("     %-12s  %+7.3f vs %+7.3f   diff %+7.3f   p = %.3f\n",
                v, mean(a[[v]]), mean(b[[v]]),
                mean(a[[v]])-mean(b[[v]]), t$p.value))
  }
}

cat("\n  Within family, chase ~ gap + FB velo:\n")
for (f in c("SL","ST","CU","CH")) {
  d <- C[fam==f]
  m <- lm(chase ~ velo_gap + fb_velo, data=d)
  cf <- summary(m)$coefficients
  cat(sprintf("     %-3s  β(gap) = %+6.3f  p = %.3f   β(FB) = %+6.3f  p = %.3f   n=%d\n",
              f, cf["velo_gap","Estimate"], cf["velo_gap","Pr(>|t|)"],
              cf["fb_velo","Estimate"], cf["fb_velo","Pr(>|t|)"], nrow(d)))
}
cat("\n  Within family, xRV ~ gap + FB velo:\n")
for (f in c("SL","ST","CU","CH")) {
  d <- C[fam==f]
  m <- lm(xrv ~ velo_gap + fb_velo, data=d)
  cf <- summary(m)$coefficients
  cat(sprintf("     %-3s  β(gap) = %+6.3f  p = %.3f   β(FB) = %+6.3f  p = %.3f\n",
              f, cf["velo_gap","Estimate"], cf["velo_gap","Pr(>|t|)"],
              cf["fb_velo","Estimate"], cf["fb_velo","Pr(>|t|)"]))
}
