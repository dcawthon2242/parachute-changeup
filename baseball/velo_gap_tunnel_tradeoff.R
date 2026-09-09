#!/usr/bin/env Rscript

# The velocity-gap tradeoff, measured directly.
#
# Claim under test: a large fastball-to-secondary velocity gap buys worse hitter
# timing but costs tunneling, and the tunnel loss degrades performance.
#
# To test it we need a real tunnel metric, not a proxy. Statcast gives the full
# 9-parameter trajectory (position, velocity, acceleration at y = 50 ft), so we
# can place every pitch anywhere in flight and measure how far apart a pitcher's
# fastball and secondary actually are at the hitter's commit point.
#
#   tunnel_dist  separation of the two mean trajectories 23.8 ft from the plate,
#                the conventional commit point (inches)
#   plate_dist   separation at the front of the plate (inches)
#   break_ratio  plate_dist / tunnel_dist - how much they diverge after commit
#   commit_dist  separation at the instant the FASTBALL would be 175 ms from
#                arriving, which is the version that actually absorbs velo gap
#                (a slower pitch must be released on a higher line to reach the
#                zone, so at the same instant it sits above the fastball)
#   commit_dz    the signed vertical component of commit_dist (+ = secondary above)
#   dt_plate     extra flight time of the secondary (ms)
#   seq_tunnel   INDEPENDENT CHECK, 2026 only: the sequence-level integrated
#                separation from data/statcast_2026/tunnel_metric_2026.csv, which
#                is measured on actual consecutive fastball -> secondary pitches
#                to the same hitter rather than on mean trajectories
#
# Performance is measured with an expected run value (xRV) rather than realized
# run value, because realized RV is ~88% noise at pair sample sizes.
#
# Revision notes (v2, after the first run):
#   * one row per pitcher x secondary. v1 keyed pairs on the fastball type as
#     well, so 171 secondaries whose pitcher uses SI vs one side and FF vs the
#     other appeared twice with identical outcomes.
#   * fastball velocity is controlled in every stage of the mediation, not just
#     the last one.
#   * tunnel_dist and commit_dist are never in the same regression (r = 0.76,
#     VIF 11 in v1's horse race).
#   * added the sequence-level tunnel check, a 2025 -> 2026 out-of-sample test,
#     and a spline location of where the timing benefit flattens.

suppressPackageStartupMessages({ library(data.table); library(splines) })
set.seed(7)
options(width = 215)

cols <- c("game_year","game_type","player_name","pitcher","batter","stand","p_throws",
          "pitch_type","description","events","balls","strikes","plate_x","plate_z",
          "sz_top","sz_bot","release_speed","release_pos_x","release_pos_z",
          "release_extension","pfx_x","pfx_z","estimated_woba_using_speedangle",
          "delta_run_exp","vx0","vy0","vz0","ax","ay","az",
          "intercept_ball_minus_batter_pos_y_inches")

dt <- rbindlist(lapply(c(2025, 2026), function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress=FALSE, select=cols)))
setnames(dt, c("intercept_ball_minus_batter_pos_y_inches",
               "estimated_woba_using_speedangle"), c("depth","xw"))
dt <- dt[game_type=="R" & pitch_type != "" & balls<=3 & strikes<=2 &
         is.finite(delta_run_exp) & is.finite(plate_x) & is.finite(plate_z) &
         is.finite(vy0) & is.finite(ay) & vy0 < 0]
FB<-c("FF","SI","FC"); BR<-c("SL","ST","CU","KC","SV","CS"); OS<-c("CH","FS","FO")
dt[, pgrp := fifelse(pitch_type %in% FB,"FB", fifelse(pitch_type %in% BR,"BR",
             fifelse(pitch_type %in% OS,"OS",NA_character_)))]
dt <- dt[!is.na(pgrp)]
dt[, pg2 := fifelse(pgrp=="FB","FB","OFF")]
dt[, `:=`(px_bat = fifelse(stand=="R", -plate_x, plate_x),
          pz_rel = (plate_z - sz_bot)/pmax(sz_top - sz_bot, 0.1),
          cnt    = paste0(balls,"-",strikes),
          bat_rv = delta_run_exp)]
cat(sprintf("Pitches: %s (2025-2026 pooled)\n\n", format(nrow(dt), big.mark=",")))

# =============================================================================
# 0. Trajectory reconstruction
# =============================================================================
YP <- 17/12   # front of home plate
tof <- function(y_target, vy0, ay)
  (-vy0 - sqrt(pmax(vy0^2 - 2*ay*(50 - y_target), 0)))/ay
dt[, t_plate := tof(YP, vy0, ay)]
# anchor x0/z0 at y = 50 by integrating back from the known plate crossing
dt[, x0 := plate_x - vx0*t_plate - 0.5*ax*t_plate^2]
dt[, z0 := plate_z - vz0*t_plate - 0.5*az*t_plate^2]
posx <- function(t, x0, vx0, ax) x0 + vx0*t + 0.5*ax*t^2
posz <- function(t, z0, vz0, az) z0 + vz0*t + 0.5*az*t^2

dt[, t_tun := tof(23.8, vy0, ay)]
dt[, `:=`(tun_x = posx(t_tun, x0, vx0, ax), tun_z = posz(t_tun, z0, vz0, az))]
dt <- dt[is.finite(t_plate) & is.finite(t_tun) & t_plate > 0.2 & t_plate < 0.8]
cat(sprintf("  Flight time to plate: mean %.3f s, 10th pct %.3f, 90th pct %.3f\n",
    mean(dt$t_plate), quantile(dt$t_plate,.1), quantile(dt$t_plate,.9)))
cat(sprintf("  Commit point 23.8 ft out is reached at %.3f s, leaving %.0f ms of flight\n",
    mean(dt$t_tun), 1000*mean(dt$t_plate - dt$t_tun)))
# sanity: run the trajectory BACK from y = 50 to the release plane and compare
# with the recorded release point
dt[, y_rel := 60.5 - fifelse(is.finite(release_extension), release_extension, 6.5)]
dt[, t_rel := tof(y_rel, vy0, ay)]   # negative: release plane is behind y = 50
dt[, `:=`(x_rel_chk = posx(t_rel, x0, vx0, ax), z_rel_chk = posz(t_rel, z0, vz0, az))]
cat(sprintf("  Reconstruction check vs recorded release point: mean |dx| = %.2f in, mean |dz| = %.2f in\n\n",
    12*mean(abs(dt$x_rel_chk - dt$release_pos_x), na.rm=TRUE),
    12*mean(abs(dt$z_rel_chk - dt$release_pos_z), na.rm=TRUE)))

# =============================================================================
# 1. Expected run value, so performance is measurable
# =============================================================================
cat("############ 1. Building an expected run value ############\n")
dt[, bip := description == "hit_into_play"]
# For a non-batted-ball, delta_run_exp is essentially determined by the count and
# the outcome. Only balls in play carry real luck, so those get replaced by their
# xwOBA-implied value.
mb <- lm(bat_rv ~ ns(xw, 5) + factor(cnt), data=dt[bip == TRUE & is.finite(xw)])
dt[, xrv := bat_rv]
ib <- dt[, bip == TRUE & is.finite(xw)]
dt[ib, xrv := predict(mb, newdata=dt[ib])]
dt <- dt[!(bip == TRUE & !is.finite(xw))]
cat(sprintf("  Balls in play repriced by xwOBA: %s pitches, model R2 = %.3f\n",
    format(sum(ib), big.mark=","), summary(mb)$r.squared))

# =============================================================================
# 2. Timing score machinery (identical to timing_directional_score.R)
# =============================================================================
SWING <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip","hit_into_play")
WHIFF <- c("swinging_strike","swinging_strike_blocked","foul_tip")
sw <- dt[description %in% SWING & is.finite(depth)]
fitd <- lm(depth ~ ns(px_bat,5)*pgrp + ns(pz_rel,5) + ns(release_speed,4) + stand, data=sw)
sw[, r1 := residuals(fitd)]
sw[, nb := .N, by=batter]; sw <- sw[nb >= 200]
sw[, tdev := r1 - mean(r1), by=batter]
surf <- lm(bat_rv ~ ns(tdev,6)*pgrp + factor(cnt), data=sw)
sw[, trv := predict(surf, newdata=sw)]
sw[, trv_b := trv - mean(trv), by=batter]
sw[, trv_ex := trv_b - mean(trv_b), by=pg2]

# =============================================================================
# 3. Pair construction with real tunnel geometry
# =============================================================================
# Trajectories are averaged within batter handedness, since a pitcher aims a
# different spot at a lefty than a righty and mixing the two fakes divergence.
build_pairs <- function(dt, sw, min_fb=60, min_sec=40, min_sw=80, min_p=120) {
  TRJ <- dt[, .(n=.N, velo=mean(release_speed),
                rx=mean(release_pos_x), rz=mean(release_pos_z),
                tx=mean(tun_x), tz=mean(tun_z),
                px=mean(plate_x), pz=mean(plate_z),
                tp=mean(t_plate), tt=mean(t_tun),
                vy=mean(vy0), a_y=mean(ay), a_x=mean(ax), a_z=mean(az),
                vx=mean(vx0), vz=mean(vz0), X0=mean(x0), Z0=mean(z0)),
            by=.(pitcher, pitch_type, pg2, stand)]
  TRJ[, use := 100*n/sum(n), by=.(pitcher, stand)]

  PRIM <- TRJ[pg2=="FB" & n >= min_fb][order(pitcher, stand, -n)][, .SD[1], by=.(pitcher, stand)]
  fbc <- c("pitch_type","velo","rx","rz","tx","tz","px","pz","tp","tt",
           "vy","a_y","a_x","a_z","vx","vz","X0","Z0","n","use")
  setnames(PRIM, fbc, paste0("fb_", fbc))
  PRIM <- PRIM[, c("pitcher","stand", paste0("fb_", fbc)), with=FALSE]

  SEC <- TRJ[pg2=="OFF" & use >= 5 & n >= min_sec]
  PS <- merge(SEC, PRIM, by=c("pitcher","stand"))

  PS[, rel_dist   := 12*sqrt((rx-fb_rx)^2 + (rz-fb_rz)^2)]
  PS[, tunnel_dist:= 12*sqrt((tx-fb_tx)^2 + (tz-fb_tz)^2)]
  PS[, plate_dist := 12*sqrt((px-fb_px)^2 + (pz-fb_pz)^2)]
  PS[, break_ratio:= plate_dist/pmax(tunnel_dist, 0.5)]
  PS[, velo_gap   := fb_velo - velo]
  PS[, dt_plate   := 1000*(tp - fb_tp)]
  # Time-based commit point: where is each pitch when the fastball is 175 ms out?
  PS[, t_c := fb_tp - 0.175]
  PS[, `:=`(cx_f = fb_X0 + fb_vx*t_c + 0.5*fb_a_x*t_c^2,
            cz_f = fb_Z0 + fb_vz*t_c + 0.5*fb_a_z*t_c^2,
            cx_s = X0 + vx*t_c + 0.5*a_x*t_c^2,
            cz_s = Z0 + vz*t_c + 0.5*a_z*t_c^2)]
  PS[, commit_dist := 12*sqrt((cx_s-cx_f)^2 + (cz_s-cz_f)^2)]
  PS[, commit_dz   := 12*(cz_s - cz_f)]

  # ONE row per pitcher x secondary (v2 fix). fb_type = the fastball it was
  # paired with most often.
  P <- PS[, .(n=sum(n), use=weighted.mean(use,n),
              fb_type=fb_pitch_type[which.max(n)],
              velo_gap=weighted.mean(velo_gap,n),
              tunnel_dist=weighted.mean(tunnel_dist,n),
              plate_dist=weighted.mean(plate_dist,n),
              commit_dist=weighted.mean(commit_dist,n),
              commit_dz=weighted.mean(commit_dz,n),
              rel_dist=weighted.mean(rel_dist,n),
              break_ratio=weighted.mean(break_ratio,n),
              dt_plate=weighted.mean(dt_plate,n),
              fb_velo=weighted.mean(fb_velo,n)),
          by=.(pitcher, pitch_type)]

  OUT <- merge(
    sw[, .(swings=.N, tscore=-100*mean(trv_ex), tdev=mean(tdev),
           whiff=100*mean(description %in% WHIFF)), by=.(pitcher, pitch_type)],
    dt[, .(pitches=.N, xrv=-100*mean(xrv), rv=-100*mean(bat_rv)), by=.(pitcher, pitch_type)],
    by=c("pitcher","pitch_type"))
  P <- merge(P, OUT, by=c("pitcher","pitch_type"))
  P <- P[swings >= min_sw & pitches >= min_p & is.finite(tunnel_dist)]
  NAME <- unique(dt[, .(pitcher, player_name)], by="pitcher")
  merge(P, NAME, by="pitcher")
}
P <- build_pairs(dt, sw)
# pitch-type + usage adjusted versions of everything (residualized), used throughout
adj <- function(P, v) residuals(lm(as.formula(paste(v, "~ factor(pitch_type) + use")), P))
for (v in c("xrv","whiff","tscore","tunnel_dist","commit_dist","commit_dz","plate_dist","velo_gap"))
  set(P, j=paste0(v,"_a"), value=adj(P, v))

cat(sprintf("\n############ 2. Pair panel: %d pairs (one per pitcher x secondary), %d pitchers ############\n\n",
            nrow(P), uniqueN(P$pitcher)))
print(rbindlist(lapply(
  c("velo_gap","rel_dist","tunnel_dist","commit_dist","commit_dz","plate_dist","break_ratio","dt_plate"),
  function(v) data.table(metric=v, mean=mean(P[[v]]), sd=sd(P[[v]]),
    p10=quantile(P[[v]],.1), p90=quantile(P[[v]],.9))))[
  , lapply(.SD, function(x) if (is.numeric(x)) round(x,2) else x)], row.names=FALSE)

cat("\n  Secondary pitch type by velocity-gap quintile (why raw quintiles are confounded):\n")
P[, vgQ := cut(velo_gap, quantile(velo_gap, 0:5/5), include.lowest=TRUE,
               labels=c("Q1 smallest","Q2","Q3","Q4","Q5 largest"))]
print(dcast(P[, .N, by=.(vgQ, pitch_type)], vgQ ~ pitch_type, value.var="N", fill=0), row.names=FALSE)

# is xRV actually measurable? split each pitcher-pitch's own pitches in half
cat("\n  Reliability check (split each pitcher-pitch's own pitches in half):\n")
key <- unique(P[, .(pitcher, pitch_type)])
h <- merge(dt, key, by=c("pitcher","pitch_type"))
h[, half := sample(c(1,2), .N, replace=TRUE), by=.(pitcher, pitch_type)]
relf <- function(d, expr, minn, lab) {
  e <- d[, .(v=eval(expr), k=.N), by=.(pitcher, pitch_type, half)]
  w <- dcast(e, pitcher + pitch_type ~ half, value.var=c("v","k"))[k_1>=minn & k_2>=minn]
  r <- cor(w$v_1, w$v_2)
  data.table(outcome=lab, pairs=nrow(w), split_half=round(r,3),
             reliability=round(2*r/(1+r),3))
}
# tunnel geometry reliability: rebuild the per-(pitcher,type,stand) trajectory
# means on each random half and compare the resulting tunnel_dist
fbk <- unique(P[, .(pitcher, pitch_type=fb_type)])
hh2 <- merge(dt, unique(rbind(key, fbk)), by=c("pitcher","pitch_type"))
hh2[, half := sample(c(1,2), .N, replace=TRUE), by=.(pitcher, pitch_type)]
tun_half <- function(hf) {
  T2 <- hh2[half==hf, .(tx=mean(tun_x), tz=mean(tun_z), n=.N), by=.(pitcher, pitch_type, pg2, stand)]
  fb <- T2[pg2=="FB", .(pitcher, stand, fx=tx, fz=tz)]
  m <- merge(T2[pg2=="OFF"], fb, by=c("pitcher","stand"))
  m[, td := 12*sqrt((tx-fx)^2 + (tz-fz)^2)]
  m[, .(td=weighted.mean(td, n), n=sum(n)), by=.(pitcher, pitch_type)]
}
tw <- merge(tun_half(1), tun_half(2), by=c("pitcher","pitch_type"))[n.x>=40 & n.y>=40]
rt <- cor(tw$td.x, tw$td.y)
print(rbindlist(list(
  relf(h, quote(-100*mean(xrv)), 50, "Expected run value (xRV)"),
  relf(h, quote(-100*mean(bat_rv)), 50, "Realized run value"),
  relf(h[description %in% SWING], quote(100*mean(description %in% WHIFF)), 40, "Whiff rate"),
  data.table(outcome="tunnel_dist (mean-trajectory)", pairs=nrow(tw),
             split_half=round(rt,3), reliability=round(2*rt/(1+rt),3))
)), row.names=FALSE)

# =============================================================================
# 4. Does a big velocity gap actually cost tunneling?
# =============================================================================
cat("\n############ 3. Does a big velocity gap degrade tunneling? ############\n\n")
cat("  (partial = net of secondary pitch type and usage)\n")
for (v in c("rel_dist","tunnel_dist","commit_dist","commit_dz","plate_dist","break_ratio","dt_plate"))
  cat(sprintf("     cor(velo gap, %-12s) = %+.3f     partial = %+.3f\n", v,
      cor(P$velo_gap, P[[v]]),
      cor(P$velo_gap_a, adj(P, v))))
cat(sprintf("\n  In inches: +1 mph of velocity gap = %+.2f in of separation at the 23.8 ft point,\n",
    coef(lm(tunnel_dist_a ~ velo_gap_a, P))[2]))
cat(sprintf("             %+.2f in at the 175 ms time-based commit point (of which %+.2f in is vertical).\n",
    coef(lm(commit_dist_a ~ velo_gap_a, P))[2], coef(lm(commit_dz_a ~ velo_gap_a, P))[2]))

cat("\n  Tunnel geometry by velocity-gap quintile (raw, confounded by pitch type):\n")
print(P[, .(pairs=.N, velo_gap=round(mean(velo_gap),1),
            rel_dist=round(mean(rel_dist),2), tunnel=round(mean(tunnel_dist),2),
            commit=round(mean(commit_dist),2), plate=round(mean(plate_dist),2),
            brk_ratio=round(mean(break_ratio),2), dt_ms=round(mean(dt_plate),1),
            tdev=round(mean(tdev),2), tscore=round(mean(tscore),3),
            whiff=round(mean(whiff),1), xrv=round(mean(xrv),3)),
        by=vgQ][order(vgQ)], row.names=FALSE)

# =============================================================================
# 5. Is performance actually degraded? Look for the turn.
# =============================================================================
cat("\n############ 4. Is there a point where more gap starts to hurt? ############\n\n")
cat("  Deciles of velocity gap (pitch-type and usage adjusted outcomes):\n")
P[, vgD := cut(velo_gap, quantile(velo_gap, 0:10/10), include.lowest=TRUE, labels=1:10)]
print(P[, .(pairs=.N, velo_gap=round(mean(velo_gap),1),
            tunnel_adj=round(mean(tunnel_dist_a),2), commit_adj=round(mean(commit_dist_a),2),
            tscore_adj=round(mean(tscore_a),3),
            whiff_adj=round(mean(whiff_a),2), xrv_adj=round(mean(xrv_a),3)),
        by=vgD][order(vgD)], row.names=FALSE)

cat("\n  Testing for curvature (does the benefit turn over?):\n")
for (y in c("tscore_a","whiff_a","xrv_a")) {
  m1 <- lm(as.formula(paste(y,"~ velo_gap + fb_velo")), P)
  m2 <- lm(as.formula(paste(y,"~ poly(velo_gap,2) + fb_velo")), P)
  m3 <- lm(as.formula(paste(y,"~ ns(velo_gap,4) + fb_velo")), P)
  cat(sprintf("     %-9s linear R2 = %.4f | quadratic p = %.3f | spline vs linear p = %.3f\n",
      y, summary(m1)$r.squared, anova(m1,m2)$`Pr(>F)`[2], anova(m1,m3)$`Pr(>F)`[2]))
}
cat("\n  Spline-fitted adjusted outcome along the velocity gap (fb_velo at its mean):\n")
grid <- data.table(velo_gap=seq(4, 17, by=1), fb_velo=mean(P$fb_velo))
for (y in c("tscore_a","whiff_a","xrv_a")) {
  m3 <- lm(as.formula(paste(y,"~ ns(velo_gap,4) + fb_velo")), P)
  set(grid, j=y, value=predict(m3, newdata=grid))
}
grid[, n_pairs := sapply(velo_gap, function(g) sum(abs(P$velo_gap - g) <= 0.5))]
print(grid[, .(velo_gap, n_pairs, tscore_adj=round(tscore_a,3), whiff_adj=round(whiff_a,2), xrv_adj=round(xrv_a,3))], row.names=FALSE)
# bootstrap the argmax of the timing and whiff curves
set.seed(11); B <- 400
am <- rbindlist(lapply(1:B, function(b) {
  Pb <- P[sample(.N, .N, replace=TRUE)]
  g <- data.table(velo_gap=seq(4, 17, by=0.25), fb_velo=mean(P$fb_velo))
  data.table(
    ts = g$velo_gap[which.max(predict(lm(tscore_a ~ ns(velo_gap,4) + fb_velo, Pb), g))],
    wh = g$velo_gap[which.max(predict(lm(whiff_a  ~ ns(velo_gap,4) + fb_velo, Pb), g))],
    xr = g$velo_gap[which.max(predict(lm(xrv_a    ~ ns(velo_gap,4) + fb_velo, Pb), g))])
}))
cat(sprintf("\n  Velocity gap that maximizes the adjusted timing score: %.1f mph (bootstrap 80%% CI %.1f-%.1f)\n",
    grid$velo_gap[which.max(grid$tscore_a)], quantile(am$ts,.1), quantile(am$ts,.9)))
cat(sprintf("  Velocity gap that maximizes the adjusted whiff rate:    %.1f mph (bootstrap 80%% CI %.1f-%.1f)\n",
    grid$velo_gap[which.max(grid$whiff_a)], quantile(am$wh,.1), quantile(am$wh,.9)))
cat(sprintf("  Velocity gap that maximizes adjusted xRV:               %.1f mph (bootstrap 80%% CI %.1f-%.1f)\n",
    grid$velo_gap[which.max(grid$xrv_a)], quantile(am$xr,.1), quantile(am$xr,.9)))
cat(sprintf("  Share of pairs with a velocity gap above 13 mph: %.1f%%\n", 100*mean(P$velo_gap > 13)))

# =============================================================================
# 6. Mediation: does the tunnel loss eat the timing gain?
# =============================================================================
cat("\n############ 5. Decomposition: velo gap -> tunnel -> performance ############\n\n")
cat("  All variables standardized; pitch type, usage and fastball velocity controlled at every stage.\n\n")
Z <- copy(P)
for (v in c("velo_gap","tunnel_dist","commit_dist","commit_dz","rel_dist","plate_dist","use","fb_velo",
            "tscore_a","whiff_a","xrv_a"))
  set(Z, j=v, value=as.numeric(scale(Z[[v]])))

med <- function(y, lab, mvar="tunnel_dist") {
  ctrl <- "+ factor(pitch_type) + use + fb_velo"
  fa <- lm(as.formula(paste(mvar,"~ velo_gap", ctrl)), Z)
  a  <- coef(fa)["velo_gap"]
  fc <- lm(as.formula(paste(y,"~ velo_gap", ctrl)), Z)
  ct <- coef(fc)["velo_gap"]
  fb <- lm(as.formula(paste(y,"~ velo_gap +",mvar, ctrl)), Z)
  cd <- coef(fb)["velo_gap"]; b <- coef(fb)[mvar]
  pb <- summary(fb)$coefficients[mvar,4]
  cat(sprintf("  -- %s, mediator = %s --\n", lab, mvar))
  cat(sprintf("     velo gap -> mediator            a = %+.3f\n", a))
  cat(sprintf("     mediator -> outcome (net)       b = %+.3f (p = %.4f)\n", b, pb))
  cat(sprintf("     total effect of velo gap        c = %+.3f (p = %.4f)\n", ct,
      summary(fc)$coefficients["velo_gap",4]))
  cat(sprintf("     direct effect, mediator held    c'= %+.3f\n", cd))
  cat(sprintf("     indirect through mediator     a*b = %+.3f  (%.0f%% of total)\n\n",
      a*b, 100*(a*b)/ct))
  invisible(data.table(outcome=lab, mediator=mvar, a=a, b=b, p_b=pb, total=ct, direct=cd,
                       indirect=a*b, share=(a*b)/ct))
}
MED <- rbindlist(list(
  med("tscore_a", "Timing score"),
  med("whiff_a",  "Whiff rate"),
  med("xrv_a",    "Expected run value"),
  med("tscore_a", "Timing score",       "commit_dist"),
  med("whiff_a",  "Whiff rate",         "commit_dist"),
  med("xrv_a",    "Expected run value", "commit_dist")))

cat("  The same decomposition in natural units (per +1 mph of velocity gap):\n\n")
nat <- function(y, ylab, unit, mvar) {
  ctrl <- "+ factor(pitch_type) + use + fb_velo"
  net    <- coef(lm(as.formula(paste(y, "~ velo_gap", ctrl)), P))["velo_gap"]
  direct <- coef(lm(as.formula(paste(y, "~ velo_gap +", mvar, ctrl)), P))["velo_gap"]
  cat(sprintf("     %-19s via %-11s if tunnel were free %+.3f %s | given back through tunnel %+.3f | observed net %+.3f\n",
      ylab, mvar, direct, unit, net - direct, net))
}
for (mv in c("tunnel_dist","commit_dist")) {
  nat("tscore_a", "Timing score", "runs/100 sw", mv)
  nat("whiff_a",  "Whiff rate",   "pts        ", mv)
  nat("xrv_a",    "Expected RV",  "runs/100 p ", mv)
}
cat("\n  Is it the vertical 'pop' or horizontal separation at the commit instant that costs whiffs?\n")
P[, commit_dx := sqrt(pmax(commit_dist^2 - commit_dz^2, 0))]
P[, commit_dx_a := adj(P, "commit_dx")]
Zc <- copy(P); for (v in c("velo_gap","commit_dz","commit_dx","use","fb_velo","whiff_a","tscore_a","xrv_a"))
  set(Zc, j=v, value=as.numeric(scale(Zc[[v]])))
Zc[, commit_dz_abs := as.numeric(scale(abs(P$commit_dz)))]
for (y in c("tscore_a","whiff_a","xrv_a")) {
  m <- lm(as.formula(paste(y, "~ velo_gap + commit_dz + commit_dx + factor(pitch_type) + use + fb_velo")), Zc)
  co <- summary(m)$coefficients
  cat(sprintf("     %-9s velo_gap %+.3f | vertical (signed, + = secondary above) %+.3f (p=%.4f) | horizontal %+.3f (p=%.4f)\n",
      y, co["velo_gap",1], co["commit_dz",1], co["commit_dz",4], co["commit_dx",1], co["commit_dx",4]))
}

cat("\n  Which separation predicts performance? (commit_dist and tunnel_dist tested\n")
cat("  separately - they correlate 0.76 and blow up jointly)\n\n")
for (mv in c("tunnel_dist","commit_dist")) {
  hr <- lm(as.formula(paste("xrv_a ~ velo_gap +", mv, "+ rel_dist + plate_dist + fb_velo + use + factor(pitch_type)")), Z)
  co <- as.data.table(summary(hr)$coefficients, keep.rownames="term")
  setnames(co, c("term","beta","se","t","p"))
  cat(sprintf("  model with %s (R2 = %.4f):\n", mv, summary(hr)$r.squared))
  print(co[!grepl("factor|Intercept", term)][order(-abs(beta))][
    , .(term, beta=round(beta,4), t=round(t,2), p=round(p,4),
        sig=fifelse(p<0.01,"**",fifelse(p<0.05,"*",fifelse(p<0.1,".",""))))],
    row.names=FALSE)
  cat("\n")
}

# =============================================================================
# 7. Independent check: the sequence-level tunnel (2026 only)
# =============================================================================
cat("############ 6. Independent check with the sequence-level tunnel metric (2026) ############\n\n")
seqf <- "data/statcast_2026/tunnel_metric_2026.csv"
if (file.exists(seqf)) {
  S <- fread(seqf)
  S <- S[p_pitch_type %in% FB & !(pitch_type %in% FB) & is.finite(tunnel)]
  SQ <- S[, .(seq_n=.N, seq_tunnel=mean(tunnel)), by=.(pitcher, pitch_type)][seq_n >= 40]
  P26 <- build_pairs(dt[game_year==2026], sw[game_year==2026], min_sw=60, min_p=90)
  P26 <- merge(P26, SQ, by=c("pitcher","pitch_type"))
  for (v in c("xrv","whiff","tscore","tunnel_dist","commit_dist","velo_gap","seq_tunnel"))
    set(P26, j=paste0(v,"_a"), value=adj(P26, v))
  cat(sprintf("  %d 2026 pairs with >= 40 actual fastball -> secondary sequences\n", nrow(P26)))
  cat(sprintf("  cor(sequence tunnel, mean-trajectory tunnel_dist) = %+.3f  |  with commit_dist = %+.3f  (adjusted)\n",
      cor(P26$seq_tunnel_a, P26$tunnel_dist_a), cor(P26$seq_tunnel_a, P26$commit_dist_a)))
  cat(sprintf("  cor(velo gap, sequence tunnel)                    = %+.3f  (adjusted)\n",
      cor(P26$velo_gap_a, P26$seq_tunnel_a)))
  cat("\n  adjusted correlations with outcomes:\n")
  for (x in c("velo_gap_a","seq_tunnel_a","tunnel_dist_a","commit_dist_a"))
    cat(sprintf("     %-14s tscore %+.3f | whiff %+.3f | xrv %+.3f\n", x,
        cor(P26[[x]], P26$tscore_a), cor(P26[[x]], P26$whiff_a), cor(P26[[x]], P26$xrv_a)))
  # The sequence metric integrates the 3D gap over the fastball's reaction window,
  # depth included (WY = 0.4). A slower pitch lags in depth throughout, so the
  # metric grows with velocity gap by construction. How much of it is just velo gap?
  r2 <- summary(lm(seq_tunnel_a ~ velo_gap_a, P26))$r.squared
  cat(sprintf("\n  R2 of sequence tunnel on velo gap alone (both adjusted): %.3f\n", r2))
  P26[, seq_resid := residuals(lm(seq_tunnel_a ~ velo_gap_a))]
  cat("  What is left of the sequence tunnel after removing velo gap, vs outcomes:\n")
  cat(sprintf("     seq_tunnel | velo_gap   tscore %+.3f | whiff %+.3f | xrv %+.3f\n",
      cor(P26$seq_resid, P26$tscore_a), cor(P26$seq_resid, P26$whiff_a), cor(P26$seq_resid, P26$xrv_a)))
  cat("  => the sequence-level tunnel cannot adjudicate the tradeoff; it is mostly a\n")
  cat("     re-expression of the velocity gap. The x/z-only geometry above is the fair test.\n")
} else cat("  (sequence tunnel file not found - skipped)\n")

# =============================================================================
# 8. Out of sample: 2025 pair geometry -> 2026 performance
# =============================================================================
cat("\n############ 7. Out of sample: does 2025 velo gap / tunnel forecast 2026 performance? ############\n\n")
P25  <- build_pairs(dt[game_year==2025], sw[game_year==2025], min_sw=60, min_p=90)
P26b <- build_pairs(dt[game_year==2026], sw[game_year==2026], min_sw=60, min_p=90)
for (v in c("xrv","whiff","tscore","tunnel_dist","commit_dist","velo_gap","rv")) {
  set(P25,  j=paste0(v,"_a"), value=adj(P25, v))
  set(P26b, j=paste0(v,"_a"), value=adj(P26b, v))
}
YY <- merge(P25[, .(pitcher, pitch_type, use, fb_velo, velo_gap_a, tunnel_dist_a, commit_dist_a,
                    tscore25=tscore_a, whiff25=whiff_a, xrv25=xrv_a, rv25=rv_a)],
            P26b[, .(pitcher, pitch_type, tscore26=tscore_a, whiff26=whiff_a, xrv26=xrv_a, rv26=rv_a)],
            by=c("pitcher","pitch_type"))
cat(sprintf("  %d pairs qualify in both seasons\n\n", nrow(YY)))
cat("  year-over-year stability of the pair outcomes (adjusted):\n")
for (v in c("tscore","whiff","xrv","rv"))
  cat(sprintf("     %-7s r(2025, 2026) = %+.3f\n", v, cor(YY[[paste0(v,"25")]], YY[[paste0(v,"26")]])))
cat("\n  2025 geometry -> 2026 outcome (adjusted, standardized betas):\n")
for (y in c("tscore26","whiff26","xrv26","rv26")) {
  m <- lm(as.formula(paste("scale(",y,") ~ scale(velo_gap_a) + scale(tunnel_dist_a) + scale(fb_velo)")), YY)
  co <- summary(m)$coefficients
  cat(sprintf("     %-8s velo_gap %+.3f (p=%.3f) | tunnel_dist %+.3f (p=%.3f) | fb_velo %+.3f | R2 %.3f\n", y,
      co[2,1], co[2,4], co[3,1], co[3,4], co[4,1], summary(m)$r.squared))
  m2 <- lm(as.formula(paste("scale(",y,") ~ scale(velo_gap_a) + scale(commit_dist_a) + scale(fb_velo)")), YY)
  co2 <- summary(m2)$coefficients
  cat(sprintf("     %-8s velo_gap %+.3f (p=%.3f) | commit_dist %+.3f (p=%.3f) | fb_velo %+.3f | R2 %.3f\n", "",
      co2[2,1], co2[2,4], co2[3,1], co2[3,4], co2[4,1], summary(m2)$r.squared))
}

# =============================================================================
# 9. Who beats the tradeoff?
# =============================================================================
cat("\n############ 8. Pairs with a big gap that tunnel anyway ############\n\n")
big <- P$velo_gap >= quantile(P$velo_gap, 2/3)
tight <- P$tunnel_dist_a <= quantile(P$tunnel_dist_a, 1/3)
cat(sprintf("  %d of %d big-gap pairs also tunnel better than expected for their pitch type.\n",
            sum(big & tight), sum(big)))
cmp <- P[, .(pairs=.N, velo_gap=round(mean(velo_gap),1),
             tunnel=round(mean(tunnel_dist),2), tunnel_adj=round(mean(tunnel_dist_a),2),
             tscore_adj=round(mean(tscore_a),3),
             whiff_adj=round(mean(whiff_a),2), xrv_adj=round(mean(xrv_a),3)),
         by=.(grp = fifelse(big & tight, "big gap + tight tunnel",
              fifelse(big, "big gap only",
              fifelse(tight, "tight tunnel only", "neither"))))]
print(cmp[order(-xrv_adj)], row.names=FALSE)

SHOW <- c("player_name","fb_type","pitch_type","use","pitches","velo_gap",
          "tunnel_dist","commit_dist","plate_dist","break_ratio","dt_plate",
          "tdev","tscore","whiff","xrv")
cat("\n  Best 15 of the big-gap + tight-tunnel group by expected run value:\n")
print(P[big & tight][order(-xrv)][1:15, ..SHOW][
  , lapply(.SD, function(x) if (is.numeric(x)) round(x,2) else x)], row.names=FALSE)
cat("\n  Best 15 big-gap pairs that do NOT tunnel (big gap only), by expected run value:\n")
print(P[big & !tight][order(-xrv)][1:15, ..SHOW][
  , lapply(.SD, function(x) if (is.numeric(x)) round(x,2) else x)], row.names=FALSE)

fwrite(P, "data/statcast_2026/velo_gap_tunnel_pairs.csv")
fwrite(MED, "data/statcast_2026/velo_gap_tunnel_mediation.csv")
fwrite(grid, "data/statcast_2026/velo_gap_tunnel_curve.csv")
cat("\nWrote data/statcast_2026/velo_gap_tunnel_pairs.csv, _mediation.csv, _curve.csv\n")
