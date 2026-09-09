#!/usr/bin/env Rscript

# 90 vs 86 sliders when acceleration matches: same ax and az.
# Same force, different velo → the 86 gets more observed break (more flight time).
# 2020–2026 pitcher-seasons, primary FB 96–100.

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
dt <- dt[game_type=="R" & pitch_type!="" & balls<=3 & strikes<=2 &
         is.finite(delta_run_exp)]
dt[, bat_rv := delta_run_exp]
dt[, in_zone := is.finite(plate_x) & abs(plate_x)<=0.83 &
                is.finite(plate_z) & plate_z>=sz_bot & plate_z<=sz_top]
WH <- c("swinging_strike","swinging_strike_blocked","foul_tip")
SW <- c(WH,"foul","hit_into_play")
dt[, swung := description %in% SW]
dt[, cnt := paste0(balls,"-",strikes)]
# RHP-equivalent ax (LHP flipped), same as tjStuff+
dt[, ax_r := fifelse(p_throws=="L", -ax, ax)]
dt[, pfx_glove := fifelse(p_throws=="R", -pfx_x, pfx_x)]

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
          .(fb_pri=mean(release_speed), fb_pri_n=.N,
            fb_ax=mean(ax_r, na.rm=TRUE), fb_az=mean(az, na.rm=TRUE)),
          by=.(game_year, pitcher)]

SL <- dt[pitch_type=="SL", .(
  n=.N,
  name=player_name[1],
  sl_velo=mean(release_speed),
  ax=mean(ax_r, na.rm=TRUE),
  az=mean(az, na.rm=TRUE),
  hmov=mean(pfx_glove, na.rm=TRUE),
  vmov=mean(pfx_z, na.rm=TRUE),
  xrv=-100*mean(xrv_p),
  zone=100*mean(in_zone, na.rm=TRUE),
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
  ch_cs=-100*mean(fifelse(description=="called_strike", xrv_p, 0)),
  ch_ball=-100*mean(fifelse(description=="ball", xrv_p, 0)),
  ch_chase=-100*mean(fifelse(!in_zone & swung, xrv_p, 0)),
  ch_ooz_take=-100*mean(fifelse(!in_zone & !swung, xrv_p, 0)),
  ch_iz_swing=-100*mean(fifelse(in_zone & swung, xrv_p, 0))
), by=.(game_year, pitcher)]
SL <- merge(SL, FBp, by=c("game_year","pitcher"))
SL[, gap := fb_pri - sl_velo]
SL[, ax_diff := ax - fb_ax]
SL[, az_diff := az - fb_az]
SL <- SL[n>=80 & fb_pri_n>=80 & is.finite(gap) & is.finite(ax) & is.finite(az)]

HARD <- SL[fb_pri>=96 & fb_pri<=100]
A <- HARD[sl_velo>=89 & sl_velo<91.5]
B <- HARD[sl_velo>=85 & sl_velo<87.5]
cat(sprintf("Hard-FB pool: %d seasons.  ~90: %d    ~86: %d\n",
            nrow(HARD), nrow(A), nrow(B)))
cat(sprintf("  Unmatched ax  90 %+5.2f vs 86 %+5.2f\n", mean(A$ax), mean(B$ax)))
cat(sprintf("  Unmatched az  90 %+5.2f vs 86 %+5.2f\n", mean(A$az), mean(B$az)))
cat(sprintf("  Unmatched pfx 90 h %+5.2f v %+5.2f   86 h %+5.2f v %+5.2f\n\n",
            mean(A$hmov), mean(A$vmov), mean(B$hmov), mean(B$vmov)))

pt <- function(x, y, lab) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok)<8) { cat(sprintf("     %-12s n too small\n", lab)); return() }
  t <- t.test(x[ok], y[ok], paired=TRUE)
  cat(sprintf("     %-12s  90 %+7.3f vs 86 %+7.3f   diff %+7.3f   p = %.3f   n=%d\n",
              lab, mean(x[ok]), mean(y[ok]), mean(x[ok]-y[ok]), t$p.value, sum(ok)))
}
tt <- function(a, b, lab) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  t <- t.test(a, b)
  cat(sprintf("     %-12s  90 %+7.3f vs 86 %+7.3f   diff %+7.3f   p = %.3f   n=%d/%d\n",
              lab, mean(a), mean(b), mean(a)-mean(b), t$p.value, length(a), length(b)))
}

match_pairs <- function(A, B, v1, v2, caliper=Inf, label="") {
  POOL <- rbind(copy(A)[, tag:="90"], copy(B)[, tag:="86"])
  POOL[, z1 := as.numeric(scale(get(v1)))]
  POOL[, z2 := as.numeric(scale(get(v2)))]
  nn <- rbindlist(lapply(unique(A$game_year), function(yr) {
    a <- POOL[tag=="90" & game_year==yr]
    cand <- POOL[tag=="86" & game_year==yr]
    if (!nrow(a) || !nrow(cand)) return(NULL)
    rbindlist(lapply(seq_len(nrow(a)), function(i) {
      d <- sqrt((cand$z1 - a$z1[i])^2 + (cand$z2 - a$z2[i])^2)
      j <- which.min(d)
      data.table(yr=yr, p90=a$pitcher[i], p86=cand$pitcher[j], dist=d[j])
    }))
  }))
  setorder(nn, dist)
  nn <- nn[!duplicated(paste(yr, p86))]
  nn <- nn[dist <= caliper]
  cat(sprintf("\n############ %s ############\n\n", label))
  cat(sprintf("  Pairs: %d   mean dist %.2f SD   caliper %s\n",
              nrow(nn), mean(nn$dist),
              if (is.finite(caliper)) sprintf("%.2f SD", caliper) else "none"))
  M90 <- merge(nn, A, by.x=c("yr","p90"), by.y=c("game_year","pitcher"))
  M86 <- merge(nn, B, by.x=c("yr","p86"), by.y=c("game_year","pitcher"))
  cat(sprintf("  ax   90 %+5.2f vs 86 %+5.2f\n", mean(M90$ax), mean(M86$ax)))
  cat(sprintf("  az   90 %+5.2f vs 86 %+5.2f\n", mean(M90$az), mean(M86$az)))
  cat(sprintf("  pfxh 90 %+5.2f vs 86 %+5.2f\n", mean(M90$hmov), mean(M86$hmov)))
  cat(sprintf("  pfxv 90 %+5.2f vs 86 %+5.2f\n", mean(M90$vmov), mean(M86$vmov)))
  cat(sprintf("  FB   90 %.2f vs 86 %.2f\n\n", mean(M90$fb_pri), mean(M86$fb_pri)))
  cat("  Paired (90 minus 86):\n")
  for (v in c("xrv","whiff","wswing","chase","swing_z","cs","ball","hard","xwcon",
              "sl_velo","gap","fb_pri","ax","az","hmov","vmov",
              "ch_chase","ch_ooz_take","ch_iz_swing","ch_whiff","ch_bip"))
    pt(M90[[v]], M86[[v]], v)
  invisible(list(nn=nn, M90=M90, M86=M86))
}

# 1. Same ax and az (the question)
m1 <- match_pairs(A, B, "ax", "az", Inf, "1. Same-year match on slider ax and az")
m1c <- match_pairs(A, B, "ax", "az", 0.50, "1b. Same ax/az, caliper 0.50 SD")

# 2. Same observed break, for contrast
m2 <- match_pairs(A, B, "hmov", "vmov", Inf, "2. Contrast: match on pfx (observed break)")

# 3. Continuous
cat("\n############ 3. Continuous: SL velo | FB + ax + az + year ############\n\n")
cat("  Hard-FB 96–100:\n")
print(round(summary(lm(xrv ~ sl_velo + fb_pri + ax + az + factor(game_year),
                       data=HARD))$coefficients[1:5,], 4))
cat("\n  All SL pitcher-seasons:\n")
print(round(summary(lm(xrv ~ sl_velo + fb_pri + ax + az + factor(game_year),
                       data=SL))$coefficients[1:5,], 4))
cat("\n  Hard-FB, also control pfx:\n")
print(round(summary(lm(xrv ~ sl_velo + fb_pri + ax + az + hmov + vmov + factor(game_year),
                       data=HARD))$coefficients[1:7,], 4))

cat("\n  Other outcomes ~ sl_velo | FB + ax + az + year, hard-FB:\n")
for (y in c("chase","whiff","xwcon","hard","hmov","vmov")) {
  m <- lm(as.formula(paste(y, "~ sl_velo + fb_pri + ax + az + factor(game_year)")), data=HARD)
  cf <- summary(m)$coefficients["sl_velo",]
  cat(sprintf("     %-8s  β(SL velo) = %+8.4f   p = %.3f\n", y, cf[1], cf[4]))
}

# physics check: among hard-FB, pfx ~ ax/az * flight time proxy
cat("\n############ 4. Physics check: same a, slower pitch = more pfx ############\n\n")
cat(sprintf("  cor(az, vmov) = %+.3f   cor(ax, hmov) = %+.3f\n",
            cor(HARD$az, HARD$vmov), cor(HARD$ax, HARD$hmov)))
cat("  vmov ~ az + sl_velo  (slower should add break at fixed az):\n")
print(round(summary(lm(vmov ~ az + sl_velo, data=HARD))$coefficients, 4))
cat("  hmov ~ ax + sl_velo:\n")
print(round(summary(lm(hmov ~ ax + sl_velo, data=HARD))$coefficients, 4))
