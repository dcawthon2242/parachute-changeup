#!/usr/bin/env Rscript

# Does "harder secondary is not better once FB is held" apply beyond sliders?
# Pitcher-seasons 2020–2026 (ST from 2023). Primary FB 96–100 band + full-sample
# continuous (FB as covariate). Type-specific hard/soft velo bands.

suppressPackageStartupMessages({ library(data.table); library(splines) })
set.seed(7)
options(width = 215)

yrs <- 2020:2026
cols <- c("game_year","game_type","pitcher","player_name","p_throws",
          "pitch_type","description","balls","strikes",
          "release_speed","pfx_x","pfx_z","ax","az",
          "plate_x","plate_z","sz_top","sz_bot",
          "launch_speed","launch_speed_angle",
          "estimated_woba_using_speedangle","delta_run_exp")
dt <- rbindlist(lapply(yrs, function(yr) {
  f <- file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr))
  if (!file.exists(f)) return(NULL)
  cat(sprintf("  reading %d\n", yr))
  fread(f, showProgress=FALSE, select=cols)
}), fill=TRUE)
setnames(dt, "estimated_woba_using_speedangle", "xw")
dt <- dt[game_type=="R" & pitch_type!="" & balls<=3 & strikes<=2 & is.finite(delta_run_exp)]
dt[, bat_rv := delta_run_exp]
dt[, in_zone := is.finite(plate_x) & abs(plate_x)<=0.83 &
                is.finite(plate_z) & plate_z>=sz_bot & plate_z<=sz_top]
WH <- c("swinging_strike","swinging_strike_blocked","foul_tip")
SW <- c(WH,"foul","hit_into_play")
dt[, swung := description %in% SW]
dt[, cnt := paste0(balls,"-",strikes)]
dt[, ax_r := fifelse(p_throws=="L", -ax, ax)]
mb <- lm(bat_rv ~ ns(xw,5)+factor(cnt),
         data=dt[description=="hit_into_play" & is.finite(xw)])
dt[, xrv_p := bat_rv]
dt[description=="hit_into_play" & is.finite(xw),
   xrv_p := predict(mb, newdata=dt[description=="hit_into_play" & is.finite(xw)])]
dt <- dt[!(description=="hit_into_play" & !is.finite(xw))]

FB <- dt[pitch_type %in% c("FF","SI","FC") & is.finite(release_speed)]
PRI <- FB[, .N, by=.(game_year, pitcher, pitch_type)][order(-N), .SD[1],
         by=.(game_year, pitcher)]
FBp <- FB[PRI, on=.(game_year, pitcher, pitch_type),
          .(fb_pri=mean(release_speed), fb_pri_n=.N), by=.(game_year, pitcher)]

# map to families
dt[, fam := fcase(
  pitch_type=="SL", "SL",
  pitch_type=="ST", "ST",
  pitch_type %in% c("CU","KC"), "CU",
  pitch_type %in% c("CH","FS"), "CH",
  default=NA_character_)]
sec <- dt[!is.na(fam)]

PP <- sec[, .(
  n=.N,
  name=player_name[1],
  ptype=pitch_type[1],
  velo=mean(release_speed),
  ax=mean(ax_r, na.rm=TRUE),
  az=mean(az, na.rm=TRUE),
  xrv=-100*mean(xrv_p),
  chase=100*mean(swung[!in_zone]),
  whiff=100*mean(description %in% WH),
  hard=100*mean(launch_speed[description=="hit_into_play"]>=95, na.rm=TRUE),
  xwcon=mean(xw[description=="hit_into_play"], na.rm=TRUE)
), by=.(game_year, pitcher, fam)]
PP <- merge(PP, FBp, by=c("game_year","pitcher"))
PP[, gap := fb_pri - velo]
PP <- PP[n>=80 & fb_pri_n>=80 & is.finite(gap) & is.finite(ax) & is.finite(az)]

cat(sprintf("\nPitcher-seasons by family:\n"))
print(PP[, .(seasons=.N, pitchers=uniqueN(pitcher),
             velo=round(mean(velo),1), p10=round(quantile(velo,.1),1),
             p50=round(quantile(velo,.5),1), p90=round(quantile(velo,.9),1),
             fb=round(mean(fb_pri),1), gap=round(mean(gap),1)),
         by=fam][order(fam)], row.names=FALSE)

# type-specific hard/soft: ~p75 band vs ~p25 band, 2 mph wide around those
# so we are not forcing 90 vs 86 onto curves
bands <- list(
  SL=list(hard=c(88.5, 91.5), soft=c(84.5, 87.5)),  # 90 vs 86
  ST=list(hard=c(83.5, 86.5), soft=c(79.5, 82.5)),  # 85 vs 81
  CU=list(hard=c(81.5, 84.5), soft=c(76.5, 80.0)),  # 83 vs 78
  CH=list(hard=c(87.5, 90.5), soft=c(82.5, 85.5))   # 89 vs 84
)

tt <- function(a, b, lab) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a)<8 || length(b)<8) {
    cat(sprintf("     %-8s  n too small (%d/%d)\n", lab, length(a), length(b)))
    return(invisible(NULL))
  }
  t <- t.test(a, b)
  cat(sprintf("     %-8s  hard %+7.3f vs soft %+7.3f   diff %+7.3f   p = %.3f   n=%d/%d\n",
              lab, mean(a), mean(b), mean(a)-mean(b), t$p.value, length(a), length(b)))
}

cat("\n############ 1. Hard-FB 96–100: type-specific hard vs soft secondary ############\n")
for (f in c("SL","ST","CU","CH")) {
  d <- PP[fam==f & fb_pri>=96 & fb_pri<=100]
  if (f=="ST") d <- d[game_year>=2023]
  bh <- bands[[f]]$hard; bs <- bands[[f]]$soft
  H <- d[velo>=bh[1] & velo<bh[2]]
  S <- d[velo>=bs[1] & velo<bs[2]]
  cat(sprintf("\n  --- %s  FB 96–100  hard [%.1f,%.1f) n=%d  soft [%.1f,%.1f) n=%d ---\n",
              f, bh[1], bh[2], nrow(H), bs[1], bs[2], nrow(S)))
  if (!nrow(H) || !nrow(S)) next
  cat(sprintf("     mean velo/FB/gap  hard %.1f / %.1f / %.1f    soft %.1f / %.1f / %.1f\n",
              mean(H$velo), mean(H$fb_pri), mean(H$gap),
              mean(S$velo), mean(S$fb_pri), mean(S$gap)))
  for (v in c("xrv","chase","whiff","hard","xwcon")) tt(H[[v]], S[[v]], v)
}

cat("\n############ 2. Continuous: xRV ~ sec velo + FB + ax + az + year ############\n")
cat("  Hard-FB 96–100 slice:\n")
for (f in c("SL","ST","CU","CH")) {
  d <- PP[fam==f & fb_pri>=96 & fb_pri<=100]
  if (f=="ST") d <- d[game_year>=2023]
  if (nrow(d)<40) { cat(sprintf("     %-3s  n=%d skip\n", f, nrow(d))); next }
  m1 <- lm(xrv ~ velo + fb_pri + ax + az + factor(game_year), data=d)
  m2 <- lm(xrv ~ gap + fb_pri + ax + az + factor(game_year), data=d)
  c1 <- summary(m1)$coefficients; c2 <- summary(m2)$coefficients
  cat(sprintf("     %-3s  n=%3d  β(velo) = %+6.3f p=%.3f    β(gap) = %+6.3f p=%.3f    β(FB) = %+6.3f p=%.3f\n",
              f, nrow(d),
              c1["velo","Estimate"], c1["velo","Pr(>|t|)"],
              c2["gap","Estimate"], c2["gap","Pr(>|t|)"],
              c1["fb_pri","Estimate"], c1["fb_pri","Pr(>|t|)"]))
}

cat("\n  Full sample (any FB velo) — more power, FB still in the model:\n")
for (f in c("SL","ST","CU","CH")) {
  d <- PP[fam==f]
  if (f=="ST") d <- d[game_year>=2023]
  m1 <- lm(xrv ~ velo + fb_pri + ax + az + factor(game_year), data=d)
  m2 <- lm(xrv ~ gap + fb_pri + ax + az + factor(game_year), data=d)
  c1 <- summary(m1)$coefficients; c2 <- summary(m2)$coefficients
  cat(sprintf("     %-3s  n=%4d  β(velo) = %+6.3f p=%.3f    β(gap) = %+6.3f p=%.3f    β(FB) = %+6.3f p=%.3f   cor(velo,FB)=%+.2f\n",
              f, nrow(d),
              c1["velo","Estimate"], c1["velo","Pr(>|t|)"],
              c2["gap","Estimate"], c2["gap","Pr(>|t|)"],
              c1["fb_pri","Estimate"], c1["fb_pri","Pr(>|t|)"],
              cor(d$velo, d$fb_pri)))
}

cat("\n############ 3. Contact specifically: hard-hit / xwOBAcon ~ velo | FB+ax+az+yr ############\n")
for (y in c("hard","xwcon")) {
  cat(sprintf("\n  %s, full sample:\n", y))
  for (f in c("SL","ST","CU","CH")) {
    d <- PP[fam==f]
    if (f=="ST") d <- d[game_year>=2023]
    m <- lm(as.formula(paste(y, "~ velo + fb_pri + ax + az + factor(game_year)")), data=d)
    cf <- summary(m)$coefficients["velo",]
    cat(sprintf("     %-3s  β(velo) = %+7.4f   p = %.3f\n", f, cf[1], cf[4]))
  }
}

cat("\n############ 4. Unmatched velo (no FB hold) — the stuff-model view ############\n")
for (f in c("SL","ST","CU","CH")) {
  d <- PP[fam==f]
  if (f=="ST") d <- d[game_year>=2023]
  m <- lm(xrv ~ velo, data=d)
  cf <- summary(m)$coefficients["velo",]
  cat(sprintf("     %-3s  unmatched β(velo) = %+6.3f  p=%.3f   r(velo,xRV)=%+.3f   r(FB,xRV)=%+.3f\n",
              f, cf[1], cf[4], cor(d$velo, d$xrv), cor(d$fb_pri, d$xrv)))
}
