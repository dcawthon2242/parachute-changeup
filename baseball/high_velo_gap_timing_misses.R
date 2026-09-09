#!/usr/bin/env Rscript

# High-velo, high-gap secondaries that upset timing but still lose runs.
#
# Cell: fastball velo and velo gap both in the top third of the pair panel.
# "Timing benefit" = pitch-type-and-usage-adjusted tscore above the cell median.
# "Performs poorly"  = pitch-type-and-usage-adjusted xRV in the bottom third of
# the cell. Residualizing first so we are not just listing curveballs.

suppressPackageStartupMessages({ library(data.table); library(splines) })
set.seed(7)
options(width = 220)

P <- fread("data/statcast_2026/velo_gap_tunnel_pairs.csv")
cat(sprintf("Pair panel: %d pairs, %d pitchers\n", nrow(P), uniqueN(P$pitcher)))

# =============================================================================
# 1. Define the cell and the leftover
# =============================================================================
vg_cut <- quantile(P$velo_gap, 2/3)
fv_cut <- quantile(P$fb_velo,  2/3)
cat(sprintf("  High gap  = velo_gap >= %.2f mph\n", vg_cut))
cat(sprintf("  High velo = fb_velo  >= %.2f mph\n\n", fv_cut))

P[, hi := velo_gap >= vg_cut & fb_velo >= fv_cut]
H <- P[hi == TRUE]
cat(sprintf("############ 1. High-velo + high-gap cell: %d pairs ############\n\n", nrow(H)))

# leftover = performance not explained by timing (and pitch type / usage)
H[, xrv_hat := predict(lm(xrv_a ~ tscore_a, data=H))]
H[, leftover := xrv_a - xrv_hat]

cat(sprintf("  Inside the cell, cor(tscore_a, xrv_a) = %+.3f\n", cor(H$tscore_a, H$xrv_a)))
cat(sprintf("  Residual SD of xRV after timing = %.3f  (raw xRV SD = %.3f)\n\n",
            sd(H$leftover), sd(H$xrv_a)))

H[, ts_hi := tscore_a >= median(tscore_a)]
H[, xv_lo := xrv_a    <= quantile(xrv_a, 1/3)]
H[, quad := fifelse(ts_hi & xv_lo, "timing+ / xRV-",
            fifelse(ts_hi & !xv_lo, "timing+ / xRV+",
            fifelse(!ts_hi & xv_lo, "timing- / xRV-", "timing- / xRV+")))]

Q <- H[, .(pairs=.N,
           fb_velo=mean(fb_velo), velo_gap=mean(velo_gap),
           tscore=mean(tscore), tscore_a=mean(tscore_a),
           xrv=mean(xrv), xrv_a=mean(xrv_a), leftover=mean(leftover),
           tun=mean(tunnel_dist), com=mean(commit_dist),
           plat=mean(plate_dist), whiff=mean(whiff)), by=quad]
setorder(Q, -leftover)
print(Q[, lapply(.SD, function(x) if (is.numeric(x)) round(x,3) else x)], row.names=FALSE)

MISS <- H[quad == "timing+ / xRV-"]
cat(sprintf("\n  Target group: %d pairs that upset timing and still lose runs.\n", nrow(MISS)))

# =============================================================================
# 2. Enrich with contact quality and location from Statcast
# =============================================================================
cat("\n############ 2. Why do they lose? Contact and location on the secondary ############\n")
cols <- c("game_year","game_type","pitcher","batter","stand","pitch_type","description",
          "balls","strikes","plate_x","plate_z","sz_top","sz_bot","release_speed",
          "launch_speed","launch_angle","launch_speed_angle","estimated_woba_using_speedangle",
          "woba_value","delta_run_exp","hit_distance_sc")
dt <- rbindlist(lapply(c(2025, 2026), function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress=FALSE, select=cols)))
dt <- dt[game_type=="R" & pitch_type != "" & balls<=3 & strikes<=2]
dt[, `:=`(px_bat = fifelse(stand=="R", -plate_x, plate_x),
          pz_rel = (plate_z - sz_bot)/pmax(sz_top - sz_bot, 0.1),
          cnt    = paste0(balls,"-",strikes))]
WH <- c("swinging_strike","swinging_strike_blocked","foul_tip")
SW <- c(WH, "foul","hit_into_play")

key <- unique(H[, .(pitcher, pitch_type)])
d <- merge(dt, key, by=c("pitcher","pitch_type"))
d[, bip := description=="hit_into_play"]

CTX <- d[, .(
  zone     = 100*mean(abs(plate_x) <= 0.83 & plate_z >= sz_bot & plate_z <= sz_top, na.rm=TRUE),
  heart    = 100*mean(abs(px_bat) <= 0.55 & pz_rel >= 0.25 & pz_rel <= 0.75, na.rm=TRUE),
  chase    = 100*mean(!(abs(plate_x) <= 0.83 & plate_z >= sz_bot & plate_z <= sz_top), na.rm=TRUE),
  swing    = 100*mean(description %in% SW),
  inzone_whiff = 100*mean(description %in% WH &
                   abs(plate_x)<=0.83 & plate_z>=sz_bot & plate_z<=sz_top),
  foul     = 100*mean(description=="foul"),
  in_play  = 100*mean(bip),
  xw       = mean(estimated_woba_using_speedangle[bip], na.rm=TRUE),
  woba     = mean(woba_value, na.rm=TRUE),
  hard     = 100*mean(launch_speed[bip] >= 95, na.rm=TRUE),
  barrel   = 100*mean(launch_speed_angle[bip]==6, na.rm=TRUE),
  gb       = 100*mean(launch_angle[bip] < 10, na.rm=TRUE),
  air      = 100*mean(launch_angle[bip] >= 25, na.rm=TRUE),
  pop      = 100*mean(launch_angle[bip] >= 50, na.rm=TRUE),
  ev       = mean(launch_speed[bip], na.rm=TRUE)
), by=.(pitcher, pitch_type)]

H <- merge(H, CTX, by=c("pitcher","pitch_type"), all.x=TRUE)
MISS <- H[quad == "timing+ / xRV-"]
OK   <- H[quad == "timing+ / xRV+"]
BADT <- H[quad == "timing- / xRV-"]

cmp <- function(cols) {
  rbindlist(lapply(cols, function(v)
    data.table(metric=v,
      miss=mean(MISS[[v]], na.rm=TRUE),
      timing_ok=mean(OK[[v]], na.rm=TRUE),
      no_timing=mean(BADT[[v]], na.rm=TRUE),
      cell=mean(H[[v]], na.rm=TRUE))))
}
cat("\n  Timing-benefit pairs that lose (n=", nrow(MISS),
    ") vs timing-benefit pairs that win (n=", nrow(OK), ")\n\n")
print(cmp(c("fb_velo","velo_gap","tscore","tscore_a","xrv","xrv_a","leftover",
            "tunnel_dist","commit_dist","plate_dist","break_ratio","dt_plate",
            "zone","heart","chase","swing","whiff","in_play","foul",
            "xw","hard","barrel","gb","air","pop","ev"))[
  , lapply(.SD, function(x) if (is.numeric(x)) round(x,3) else x)], row.names=FALSE)

# =============================================================================
# 3. What actually separates leftover xRV inside the cell?
# =============================================================================
cat("\n############ 3. What predicts leftover xRV once timing is held? ############\n\n")
Z <- copy(H)
vars <- c("tunnel_dist","commit_dist","plate_dist","zone","heart","chase",
          "hard","barrel","xw","ev","in_play","foul","gb","pop")
for (v in vars) set(Z, j=v, value=as.numeric(scale(Z[[v]])))
m <- lm(leftover ~ tunnel_dist + commit_dist + plate_dist + zone + heart +
          hard + barrel + in_play + foul + gb + pop + factor(pitch_type), Z)
s <- summary(m)
cat(sprintf("  R2 = %.3f  (leftover xRV ~ contact/location/tunnel + pitch type)\n\n", s$r.squared))
co <- as.data.table(s$coefficients, keep.rownames="term")
setnames(co, c("term","beta","se","t","p"))
print(co[!grepl("factor|Intercept", term)][order(-abs(beta))][
  , .(term, beta=round(beta,3), t=round(t,2), p=round(p,4),
      sig=fifelse(p<0.01,"**",fifelse(p<0.05,"*",fifelse(p<0.1,".",""))))],
  row.names=FALSE)

cat("\n  Univariate leftover correlations (timing already removed):\n")
U <- rbindlist(lapply(c("hard","barrel","xw","ev","in_play","zone","heart","chase",
                        "tunnel_dist","commit_dist","plate_dist","whiff","foul","gb","pop"),
  function(v) data.table(metric=v, r=cor(H$leftover, H[[v]], use="complete.obs"))))
setorder(U, r)
print(U[, r := round(r,3)], row.names=FALSE)

# =============================================================================
# 4. Leaderboard of the misses
# =============================================================================
cat("\n############ 4. The misses ############\n\n")
SHOW <- c("player_name","fb_type","pitch_type","use","pitches","fb_velo","velo_gap",
          "tscore","tscore_a","xrv","xrv_a","leftover",
          "tunnel_dist","commit_dist","zone","heart","whiff","in_play",
          "xw","hard","barrel","ev")
cat("  Worst leftover (timing predicted they would be fine; they were not):\n")
print(MISS[order(leftover)][, ..SHOW][
  , lapply(.SD, function(x) if (is.numeric(x)) round(x,2) else x)], row.names=FALSE)

cat("\n  Contrast: same cell, timing benefit, BEST leftover:\n")
print(OK[order(-leftover)][1:min(12,.N), ..SHOW][
  , lapply(.SD, function(x) if (is.numeric(x)) round(x,2) else x)], row.names=FALSE)

# pitcher-level: anyone with 2+ of these misses
cat("\n  Pitchers with two or more such pairs:\n")
print(MISS[, .(pairs=.N,
               types=paste(sort(unique(pitch_type)), collapse=","),
               leftover=round(mean(leftover),2),
               tscore=round(mean(tscore),2),
               xrv=round(mean(xrv),2)),
           by=.(pitcher, player_name)][pairs>=2][order(-pairs, leftover)],
      row.names=FALSE)

fwrite(MISS, "data/statcast_2026/high_velo_gap_timing_misses.csv")
cat("\nWrote data/statcast_2026/high_velo_gap_timing_misses.csv\n")
