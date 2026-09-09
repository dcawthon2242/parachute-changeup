#!/usr/bin/env Rscript

# Section 3 at scale: 98 FB + 90 SL vs 98 FB + 86 SL, 2020–2026.
# Unit is pitcher-season so velo changes are not smeared.

suppressPackageStartupMessages({ library(data.table); library(splines) })
set.seed(7)
options(width = 215)

yrs <- 2020:2026
cols <- c("game_year","game_type","pitcher","player_name","p_throws",
          "pitch_type","description","balls","strikes",
          "release_speed","pfx_x","pfx_z",
          "plate_x","plate_z","sz_top","sz_bot",
          "launch_speed","launch_speed_angle",
          "estimated_woba_using_speedangle","delta_run_exp")

dt <- rbindlist(lapply(yrs, function(yr) {
  f <- file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr))
  if (!file.exists(f)) { cat(sprintf("  missing %s\n", f)); return(NULL) }
  cat(sprintf("  reading %d\n", yr))
  fread(f, showProgress=FALSE, select=cols)
}), fill=TRUE)
setnames(dt, "estimated_woba_using_speedangle", "xw")
dt <- dt[game_type=="R" & pitch_type!="" & balls<=3 & strikes<=2 &
         is.finite(delta_run_exp)]
cat(sprintf("Pitches: %s  years: %s\n\n",
            format(nrow(dt), big.mark=","),
            paste(sort(unique(dt$game_year)), collapse=",")))

dt[, bat_rv := delta_run_exp]
dt[, in_zone := is.finite(plate_x) & abs(plate_x)<=0.83 &
                is.finite(plate_z) & plate_z>=sz_bot & plate_z<=sz_top]
WH <- c("swinging_strike","swinging_strike_blocked","foul_tip")
SW <- c(WH,"foul","hit_into_play")
dt[, swung := description %in% SW]
dt[, cnt := paste0(balls,"-",strikes)]
mb <- lm(bat_rv ~ ns(xw,5)+factor(cnt),
         data=dt[description=="hit_into_play" & is.finite(xw)])
dt[, xrv_p := bat_rv]
dt[description=="hit_into_play" & is.finite(xw),
   xrv_p := predict(mb, newdata=dt[description=="hit_into_play" & is.finite(xw)])]
dt <- dt[!(description=="hit_into_play" & !is.finite(xw))]
dt[, pfx_glove := fifelse(p_throws=="R", -pfx_x, pfx_x)]

# primary FB by pitcher-season
FB <- dt[pitch_type %in% c("FF","SI","FC") & is.finite(release_speed)]
PRI <- FB[, .N, by=.(game_year, pitcher, pitch_type)][order(-N), .SD[1],
         by=.(game_year, pitcher)]
FBp <- FB[PRI, on=.(game_year, pitcher, pitch_type),
          .(fb_pri=mean(release_speed), fb_pri_n=.N), by=.(game_year, pitcher)]

sl <- dt[pitch_type=="SL"]
SL <- sl[, .(
  n=.N,
  name=player_name[1],
  sl_velo=mean(release_speed),
  hmov=mean(pfx_glove, na.rm=TRUE),
  vmov=mean(pfx_z, na.rm=TRUE),
  xrv=-100*mean(xrv_p),
  rv=-100*mean(bat_rv),
  zone=100*mean(in_zone, na.rm=TRUE),
  swing=100*mean(swung),
  swing_z=100*mean(swung[in_zone]),
  chase=100*mean(swung[!in_zone]),
  whiff=100*mean(description %in% WH),
  wswing=100*mean((description %in% WH)[swung]),
  cs=100*mean(description=="called_strike"),
  ball=100*mean(description=="ball"),
  hard=100*mean(launch_speed[description=="hit_into_play"]>=95, na.rm=TRUE),
  brl=100*mean(launch_speed_angle[description=="hit_into_play"]==6, na.rm=TRUE),
  xwcon=mean(xw[description=="hit_into_play"], na.rm=TRUE),
  ch_whiff=-100*mean(fifelse(description %in% WH, xrv_p, 0)),
  ch_bip=-100*mean(fifelse(description=="hit_into_play", xrv_p, 0)),
  ch_foul=-100*mean(fifelse(description=="foul", xrv_p, 0)),
  ch_cs=-100*mean(fifelse(description=="called_strike", xrv_p, 0)),
  ch_ball=-100*mean(fifelse(description=="ball", xrv_p, 0)),
  ch_chase=-100*mean(fifelse(!in_zone & swung, xrv_p, 0)),
  ch_ooz_take=-100*mean(fifelse(!in_zone & !swung, xrv_p, 0)),
  ch_iz_swing=-100*mean(fifelse(in_zone & swung, xrv_p, 0)),
  ch_iz_take=-100*mean(fifelse(in_zone & !swung, xrv_p, 0))
), by=.(game_year, pitcher)]
SL <- merge(SL, FBp, by=c("game_year","pitcher"))
SL[, gap := fb_pri - sl_velo]
SL <- SL[n>=80 & fb_pri_n>=80 & is.finite(gap) & is.finite(hmov)]
cat(sprintf("Pitcher-seasons (80+ SL, 80+ primary FB): %d   unique pitchers: %d\n",
            nrow(SL), uniqueN(SL$pitcher)))
print(SL[, .(seasons=.N, pitchers=uniqueN(pitcher)), by=game_year][order(game_year)],
      row.names=FALSE)

tt <- function(a, b, lab) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a)<8 || length(b)<8) { cat(sprintf("     %-12s n too small\n", lab)); return(invisible(NULL)) }
  t <- t.test(a, b)
  cat(sprintf("     %-12s  90 %+7.3f vs 86 %+7.3f   diff %+7.3f   p = %.3f   n=%d/%d\n",
              lab, mean(a), mean(b), mean(a)-mean(b), t$p.value, length(a), length(b)))
}

HARD <- SL[fb_pri>=96 & fb_pri<=100]
A <- HARD[sl_velo>=89 & sl_velo<91.5]
B <- HARD[sl_velo>=85 & sl_velo<87.5]

cat(sprintf("\n############ Hold FB 96–100: ~90 SL vs ~86 SL, 2020–2026 ############\n\n"))
cat(sprintf("  Hard-FB pitcher-seasons: %d  (unique pitchers %d)\n",
            nrow(HARD), uniqueN(HARD$pitcher)))
cat(sprintf("  ~90 SL: %d seasons / %d pitchers\n", nrow(A), uniqueN(A$pitcher)))
cat(sprintf("  ~86 SL: %d seasons / %d pitchers\n", nrow(B), uniqueN(B$pitcher)))
cat(sprintf("  Mean FB: 90-group %.2f    86-group %.2f\n", mean(A$fb_pri), mean(B$fb_pri)))
cat(sprintf("  Mean SL: 90-group %.2f    86-group %.2f\n", mean(A$sl_velo), mean(B$sl_velo)))
cat(sprintf("  Mean gap: 90-group %.2f   86-group %.2f\n\n", mean(A$gap), mean(B$gap)))

cat("  By year:\n")
print(rbind(
  A[, .(grp="90", seasons=.N), by=game_year],
  B[, .(grp="86", seasons=.N), by=game_year]
)[, dcast(.SD, game_year ~ grp, value.var="seasons", fill=0)], row.names=FALSE)

cat("\n  Outcomes (90 minus 86; positive = 90 is better for the pitcher):\n")
for (v in c("xrv","rv","whiff","wswing","chase","swing_z","zone","cs","ball",
            "hard","brl","xwcon","hmov","vmov","fb_pri","gap")) tt(A[[v]], B[[v]], v)

cat("\n  xRV channels:\n")
for (v in c("ch_whiff","ch_bip","ch_foul","ch_cs","ch_ball",
            "ch_chase","ch_ooz_take","ch_iz_swing","ch_iz_take")) tt(A[[v]], B[[v]], v)

# weighted by pitch count
cat(sprintf("\n  Pitch-weighted xRV: 90 %+0.3f vs 86 %+0.3f\n",
            weighted.mean(A$xrv, A$n), weighted.mean(B$xrv, B$n)))

cat("\n############ Wider bands, FB 96–100: 88–92 vs 84–88 ############\n\n")
Aw <- HARD[sl_velo>=88 & sl_velo<92]
Bw <- HARD[sl_velo>=84 & sl_velo<88]
cat(sprintf("  n = %d vs %d  (pitchers %d vs %d)\n",
            nrow(Aw), nrow(Bw), uniqueN(Aw$pitcher), uniqueN(Bw$pitcher)))
for (v in c("xrv","whiff","chase","hard","xwcon","hmov","vmov","fb_pri",
            "ch_chase","ch_ooz_take","ch_iz_swing","ch_whiff","ch_bip"))
  tt(Aw[[v]], Bw[[v]], v)

cat("\n############ Movement-matched 90 vs 86, FB 96–100 ############\n\n")
POOL <- rbind(copy(A)[, tag:="90"], copy(B)[, tag:="86"])
POOL[, hz := as.numeric(scale(hmov))]
POOL[, vz := as.numeric(scale(vmov))]
# match within year so era/ST classification does not leak
nn <- rbindlist(lapply(unique(A$game_year), function(yr) {
  a <- POOL[tag=="90" & game_year==yr]
  cand <- POOL[tag=="86" & game_year==yr]
  if (nrow(a)==0 || nrow(cand)==0) return(NULL)
  rbindlist(lapply(seq_len(nrow(a)), function(i) {
    d <- (cand$hz - a$hz[i])^2 + (cand$vz - a$vz[i])^2
    j <- which.min(d)
    data.table(yr=yr, p90=a$pitcher[i], p86=cand$pitcher[j], dist=sqrt(d[j]),
               h90=a$hmov[i], h86=cand$hmov[j], v90=a$vmov[i], v86=cand$vmov[j])
  }))
}))
setorder(nn, dist)
nn <- nn[!duplicated(paste(yr, p86))]
cat(sprintf("  Within-year matched pairs: %d   mean dist %.2f SD\n", nrow(nn), mean(nn$dist)))
cat(sprintf("  Mean hmov 90/86: %.2f / %.2f    vmov: %.2f / %.2f\n\n",
            mean(nn$h90), mean(nn$h86), mean(nn$v90), mean(nn$v86)))
M90 <- merge(nn, A, by.x=c("yr","p90"), by.y=c("game_year","pitcher"))
M86 <- merge(nn, B, by.x=c("yr","p86"), by.y=c("game_year","pitcher"))
# paired tests
pt <- function(x, y, lab) {
  ok <- is.finite(x) & is.finite(y)
  t <- t.test(x[ok], y[ok], paired=TRUE)
  cat(sprintf("     %-12s  90 %+7.3f vs 86 %+7.3f   diff %+7.3f   p = %.3f   pairs=%d\n",
              lab, mean(x[ok]), mean(y[ok]), mean(x[ok]-y[ok]), t$p.value, sum(ok)))
}
cat("  Paired, same-year movement match (90 minus 86):\n")
for (v in c("xrv","whiff","wswing","chase","swing_z","zone","cs","ball",
            "hard","xwcon","sl_velo","gap","fb_pri","hmov","vmov",
            "ch_chase","ch_ooz_take","ch_iz_swing","ch_whiff","ch_bip","ch_cs","ch_ball"))
  pt(M90[[v]], M86[[v]], v)

cat("\n############ Continuous, FB 96–100, 2020–2026 ############\n\n")
cat("  xRV ~ sl_velo + fb_pri + hmov + vmov + year\n")
print(round(summary(lm(xrv ~ sl_velo + fb_pri + hmov + vmov + factor(game_year),
                       data=HARD))$coefficients[1:5,], 4))
cat("\n  xRV ~ gap + fb_pri + hmov + vmov + year\n")
print(round(summary(lm(xrv ~ gap + fb_pri + hmov + vmov + factor(game_year),
                       data=HARD))$coefficients[1:5,], 4))
cat("\n  Full SL sample (any FB velo), same model:\n")
print(round(summary(lm(xrv ~ sl_velo + fb_pri + hmov + vmov + factor(game_year),
                       data=SL))$coefficients[1:5,], 4))

cat("\n  Chase / whiff / xwOBAcon ~ sl_velo | FB + movement + year, hard-FB only:\n")
for (y in c("chase","whiff","xwcon","hard")) {
  m <- lm(as.formula(paste(y, "~ sl_velo + fb_pri + hmov + vmov + factor(game_year)")), data=HARD)
  cf <- summary(m)$coefficients["sl_velo",]
  cat(sprintf("     %-8s  β(SL velo) = %+7.4f   p = %.3f\n", y, cf[1], cf[4]))
}

cat("\n############ By era (ST classification starts 2023) ############\n\n")
for (lab in c("2020-2022","2023-2026")) {
  yrs2 <- if (lab=="2020-2022") 2020:2022 else 2023:2026
  H <- HARD[game_year %in% yrs2]
  a <- H[sl_velo>=89 & sl_velo<91.5]; b <- H[sl_velo>=85 & sl_velo<87.5]
  cat(sprintf("  %s  90 n=%d  86 n=%d\n", lab, nrow(a), nrow(b)))
  if (nrow(a)>=8 && nrow(b)>=8) {
    for (v in c("xrv","chase","whiff","hard","xwcon","hmov","vmov","fb_pri"))
      tt(a[[v]], b[[v]], v)
    cat(sprintf("     sl_velo β | FB+mov+yr  p from model:\n"))
    m <- lm(xrv ~ sl_velo + fb_pri + hmov + vmov + factor(game_year), data=H)
    cf <- summary(m)$coefficients["sl_velo",]
    cat(sprintf("                 β = %+0.4f  p = %.3f  n=%d\n\n", cf[1], cf[4], nobs(m)))
  }
}

# clustered SE by pitcher (repeat seasons)
cat("############ Clustered by pitcher (repeat seasons) ############\n\n")
if (requireNamespace("sandwich", quietly=TRUE) && requireNamespace("lmtest", quietly=TRUE)) {
  m <- lm(xrv ~ sl_velo + fb_pri + hmov + vmov + factor(game_year), data=HARD)
  vc <- sandwich::vcovCL(m, cluster=HARD$pitcher)
  print(round(lmtest::coeftest(m, vcov=vc)[1:5,], 4))
} else {
  cat("  sandwich/lmtest not installed — using pitcher-mean collapse instead.\n")
  PM <- HARD[, lapply(.SD, mean),
             .SDcols=c("xrv","sl_velo","fb_pri","hmov","vmov","chase","whiff"),
             by=pitcher]
  cat(sprintf("  Unique hard-FB pitchers: %d\n", nrow(PM)))
  print(round(summary(lm(xrv ~ sl_velo + fb_pri + hmov + vmov, data=PM))$coefficients, 4))
}

fwrite(HARD, "data/statcast_2026/hard_slider_vs_gap_2020_2026.csv")
cat("\nWrote hard-FB pitcher-seasons.\n")
