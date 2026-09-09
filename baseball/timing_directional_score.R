#!/usr/bin/env Rscript

# A pitch-type-aware timing disruption metric.
#
# Prior work established that the damaging direction on the contact-depth axis DEPENDS on
# pitch type: on a fastball, run value rises monotonically as the hitter gets out front, so
# the pitcher wants DEEP contact; on breaking balls and changeups the curve inverts, so the
# pitcher wants the hitter EARLY. A magnitude-only disruption metric averages those two
# opposite goals together and throws the signal away.
#
# The metric built here scores every swing by the run value the league actually earns at that
# (pitch group, contact depth, count) position, then credits the pitcher with the difference
# from what the HITTER he faced normally produces. Units are runs saved per 100 swings.
#
#   tscore = -100 * mean( trv_swing - batter's own mean trv )
#
# The run-value surface is a smooth spline in depth interacted with pitch group rather than
# raw cell means, so thin tails do not drive the ranking. Positive tscore = the pitcher puts
# hitters in worse places on the axis than those hitters normally reach.
#
# The metric is then decomposed into its fastball and offspeed halves, tested for
# repeatability, and raced against the old direction-agnostic disruption metric.

suppressPackageStartupMessages({ library(data.table); library(splines) })
set.seed(41)

cols <- c("game_year","game_type","player_name","pitcher","batter","stand","pitch_type",
          "description","bb_type","events","balls","strikes","plate_x","plate_z",
          "sz_top","sz_bot","release_speed","launch_speed","launch_angle",
          "launch_speed_angle","estimated_woba_using_speedangle","woba_value","woba_denom",
          "delta_run_exp","intercept_ball_minus_batter_pos_y_inches")

dt <- rbindlist(lapply(c(2025, 2026), function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress = FALSE, select = cols)))
setnames(dt, "intercept_ball_minus_batter_pos_y_inches", "depth")
dt <- dt[game_type == "R" & pitch_type != "" & balls <= 3 & strikes <= 2 &
         is.finite(delta_run_exp)]

SWING <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip","hit_into_play")
sw <- dt[description %in% SWING & is.finite(depth) & is.finite(plate_x) &
         is.finite(plate_z) & is.finite(release_speed)]
FB<-c("FF","SI","FC"); BR<-c("SL","ST","CU","KC","SV","CS"); OS<-c("CH","FS","FO")
sw[, pgrp := fifelse(pitch_type %in% FB,"FB", fifelse(pitch_type %in% BR,"BR",
             fifelse(pitch_type %in% OS,"OS",NA_character_)))]
sw <- sw[!is.na(pgrp)]
sw[, pg2 := fifelse(pgrp=="FB","FB","OFF")]
sw[, px_bat := fifelse(stand=="R", -plate_x, plate_x)]
sw[, pz_rel := (plate_z - sz_bot)/pmax(sz_top - sz_bot, 0.1)]
fitd <- lm(depth ~ ns(px_bat,5)*pgrp + ns(pz_rel,5) + ns(release_speed,4) + stand, data=sw)
sw[, r1 := residuals(fitd)]
sw[, nb := .N, by=batter]; sw <- sw[nb >= 200]
sw[, tdev := r1 - mean(r1), by=batter]
sw[, adev := abs(tdev)]
sw[, bat_rv := delta_run_exp]              # offense positive
sw[, pit_rv := -delta_run_exp]
sw[, whiff := description %in% c("swinging_strike","swinging_strike_blocked","foul_tip")]
sw[, isbip := description == "hit_into_play" & bb_type != "" &
              is.finite(launch_speed) & is.finite(launch_speed_angle)]
sw[, brl := isbip & launch_speed_angle == 6]
sw[, hard := isbip & launch_speed >= 95]
sw[, cnt := paste0(balls,"-",strikes)]

# old magnitude-only disruption, for the head-to-head
sw[, b_adev := mean(adev), by=batter]
sw[, excess := adev - b_adev]

cat(sprintf("Swings on the axis: %s (FB %.1f%%, BR %.1f%%, OS %.1f%%)\n\n",
  format(nrow(sw), big.mark=","), 100*mean(sw$pgrp=="FB"),
  100*mean(sw$pgrp=="BR"), 100*mean(sw$pgrp=="OS")))

# =============================================================================
# 1. The run-value surface: depth x pitch group, within count
# =============================================================================
cat("############ 1. The weighting surface ############\n")
cat("bat_rv ~ ns(tdev,6) * pitch group + count fixed effects, fit on all swings.\n")
cat("This is what 'weighted' means: each inch of displacement is valued at the run value the\n")
cat("league actually earns there, for that pitch type, in that count.\n\n")
surf <- lm(bat_rv ~ ns(tdev, 6)*pgrp + factor(cnt), data=sw)
cat(sprintf("  Surface R2 = %.4f on %s swings\n", summary(surf)$r.squared,
            format(nrow(sw), big.mark=",")))
sw[, trv := predict(surf, newdata=sw)]

# Two things have to come out before this is a skill measure.
#   (a) the batter: some hitters are late against everybody
#   (b) the pitch group's own league mean: offspeed puts hitters in pitcher-friendly places
#       more often than fastballs do, so an uncentered score would just reward throwing
#       offspeed. Mix is a real lever but it is not "timing", so it is scored separately.
sw[, trv_b := trv - mean(trv), by=batter]        # batter removed, mix still in
sw[, trv_ex := trv_b - mean(trv_b), by=pg2]      # placement only
cat(sprintf("\n  League mean timing RV after batter adjustment: FB %+.5f, offspeed %+.5f\n",
  sw[pg2=="FB", mean(trv_b)], sw[pg2=="OFF", mean(trv_b)]))
cat("  That gap is pitch mix, not timing skill; it is removed from tscore and scored separately.\n")

# show the surface in count-neutral form
grid <- CJ(tdev = seq(-18, 18, 3), pgrp = c("FB","BR","OS"))
grid[, cnt := "0-0"]
grid[, rv := predict(surf, newdata=grid)]
base <- grid[, .(tdev, pgrp, rv)]
cat("\n  Fitted offensive run value per swing at 0-0, by depth and pitch group:\n")
cat("  (negative tdev = deep/late, positive = out front/early)\n\n")
print(dcast(base, tdev ~ pgrp, value.var="rv")[, lapply(.SD, function(x)
  if (is.numeric(x)) round(x,4) else x)], row.names=FALSE)
cat("\n  Best depth for the PITCHER, by pitch group:\n")
for (g in c("FB","BR","OS")) {
  b <- base[pgrp==g][which.min(rv)]
  w <- base[pgrp==g][which.max(rv)]
  cat(sprintf("     %-3s pitcher-best at %+3.0f in (rv %+.4f) | hitter-best at %+3.0f in (rv %+.4f)\n",
      g, b$tdev, b$rv, w$tdev, w$rv))
}

# =============================================================================
# 2. Pitcher metric and its two halves
# =============================================================================
agg <- function(d, minsw) {
  a <- d[, .(
    swings   = .N,
    fbsh     = 100*mean(pg2=="FB"),
    tscore   = -100*mean(trv_ex),
    se_ts    = 100*sd(trv_ex)/sqrt(.N),
    # the two halves, each expressed per 100 of the pitcher's TOTAL swings so they add up
    fb_part  = -100*sum(trv_ex[pg2=="FB"])/.N,
    off_part = -100*sum(trv_ex[pg2=="OFF"])/.N,
    # and each expressed per 100 of its own swings, for rate comparison
    fb_rate  = -100*mean(trv_ex[pg2=="FB"]),
    off_rate = -100*mean(trv_ex[pg2=="OFF"]),
    tdev_fb  = mean(tdev[pg2=="FB"]),
    tdev_off = mean(tdev[pg2=="OFF"]),
    mixscore = -100*mean(trv_b - trv_ex),   # runs from pitch mix alone, no placement
    total    = -100*mean(trv_b),            # placement + mix
    disrupt  = mean(excess),
    mtdev    = mean(tdev),
    whiffpct = 100*mean(whiff),
    barrel   = 100*sum(brl)/pmax(sum(isbip),1),
    hardhit  = 100*sum(hard)/pmax(sum(isbip),1),
    xw       = mean(estimated_woba_using_speedangle[isbip], na.rm=TRUE),
    wcon     = sum(woba_value[isbip], na.rm=TRUE)/pmax(sum(woba_denom[isbip], na.rm=TRUE),1),
    rv_swing = mean(pit_rv),
    se_rv    = sd(pit_rv)/sqrt(.N)
  ), by=pitcher]
  a[swings >= minsw]
}
shrink <- function(d, val, sev) {
  d <- copy(d); v <- d[[val]]; s <- d[[sev]]; ok <- is.finite(s)
  mu <- mean(v[ok]); tau2 <- pmax(var(v[ok]) - mean(s[ok]^2), 1e-10)
  d[, rel := fifelse(ok, tau2/(tau2 + s^2), 0)]
  d[[paste0(val,"_adj")]] <- mu + d$rel*(v - mu)
  cat(sprintf("  %-8s league mean %+.3f | true SD %.3f | median reliability %.2f\n",
    val, mu, sqrt(tau2), median(d$rel)))
  d
}
NAME <- unique(sw[, .(pitcher, player_name)], by="pitcher")

cat("\n############ 2. 2026 leaderboard (min 250 swings) ############\n")
P <- agg(sw[game_year==2026], 250)
P <- shrink(P, "tscore", "se_ts")
P <- shrink(P, "rv_swing", "se_rv")
P <- merge(P, NAME, by="pitcher")
cat("\ntscore = runs saved per 100 swings from WHERE he puts hitters on the depth axis,\n")
cat("net of those hitters' own norms. fb_part + off_part = tscore.\n\n")
fmt <- function(d, n=25, asc=FALSE) {
  d <- if (asc) d[order(tscore_adj)] else d[order(-tscore_adj)]
  d <- d[1:min(n,nrow(d))]
  d[, .(player_name, swings, fbsh=round(fbsh,0),
        tscore_adj=round(tscore_adj,2), fb_part=round(fb_part,2),
        off_part=round(off_part,2), mix=round(mixscore,2),
        tdev_fb=round(tdev_fb,1), tdev_off=round(tdev_off,1),
        disrupt=round(disrupt,2), whiff=round(whiffpct,1),
        barrel=round(barrel,1), wcon=round(wcon,3), rv=round(rv_swing,4))]
}
options(width = 200)
print(fmt(P, 25), row.names=FALSE)

cat("\n  -- worst 12 --\n")
print(fmt(P, 12, asc=TRUE), row.names=FALSE)

cat(sprintf("\n  Sanity: cor(tscore, fastball usage share) = %+.3f (was the mix confound; ~0 is the goal)\n",
            cor(P$tscore, P$fbsh)))
cat(sprintf("  Mix component: SD %.2f runs/100 sw, cor with placement tscore %+.3f\n",
            sd(P$mixscore), cor(P$mixscore, P$tscore)))
cat(sprintf("  A full-season starter sees ~2,200 swings, so 1 SD of tscore (%.2f) is ~%.1f implied runs.\n",
            sd(P$tscore), sd(P$tscore)*22))

cat("\n############ 2b. Starters only (min 900 swings) ############\n")
ST <- agg(sw[game_year==2026], 900)
ST <- shrink(ST, "tscore", "se_ts")
ST <- merge(ST, NAME, by="pitcher")
cat("\n")
print(fmt(ST, 20), row.names=FALSE)

# =============================================================================
# 3. Are the two halves the same skill?
# =============================================================================
cat("\n############ 3. Is fastball timing the same skill as offspeed timing? ############\n")
Q <- P[swings >= 400]
cat(sprintf("  n = %d pitchers with >=400 swings\n", nrow(Q)))
cat(sprintf("  cor(fastball rate, offspeed rate)        = %+.3f\n", cor(Q$fb_rate, Q$off_rate)))
cat(sprintf("  cor(mean tdev on FB, mean tdev on OFF)   = %+.3f\n", cor(Q$tdev_fb, Q$tdev_off)))
cat(sprintf("  cor(tscore, old magnitude disruption)    = %+.3f\n", cor(Q$tscore, Q$disrupt)))
cat(sprintf("  cor(tscore, mean tdev overall)           = %+.3f\n", cor(Q$tscore, Q$mtdev)))
cat(sprintf("\n  Contribution to tscore variance: FB half %.0f%%, offspeed half %.0f%%, 2cov %.0f%%\n",
  100*var(Q$fb_part)/var(Q$tscore), 100*var(Q$off_part)/var(Q$tscore),
  100*2*cov(Q$fb_part, Q$off_part)/var(Q$tscore)))
cat("\n  Pitchers grouped by which half they win on:\n")
Q[, style := fifelse(fb_rate > 0 & off_rate > 0, "both directions",
            fifelse(fb_rate > 0 & off_rate <= 0, "fastball only",
            fifelse(fb_rate <= 0 & off_rate > 0, "offspeed only", "neither")))]
print(Q[, .(pitchers=.N, tscore=round(mean(tscore),2), fb_rate=round(mean(fb_rate),2),
  off_rate=round(mean(off_rate),2), tdev_fb=round(mean(tdev_fb),1),
  tdev_off=round(mean(tdev_off),1), whiff=round(mean(whiffpct),1),
  barrel=round(mean(barrel),1), wcon=round(mean(wcon),3),
  rv=round(mean(rv_swing),4)), by=style][order(-tscore)], row.names=FALSE)

# =============================================================================
# 4. Repeatability
# =============================================================================
cat("\n############ 4. Repeatability ############\n")
sb <- function(r) 2*r/(1+r)
sw[, half := sample(rep_len(1:2, .N)), by=pitcher]
H1 <- agg(sw[game_year==2026 & half==1], 120)
H2 <- agg(sw[game_year==2026 & half==2], 120)
B <- merge(H1, H2, by="pitcher", suffixes=c("_a","_b"))
Y1 <- agg(sw[game_year==2025], 300); Y2 <- agg(sw[game_year==2026], 300)
C <- merge(Y1, Y2, by="pitcher", suffixes=c("_a","_b"))
MS <- c("tscore","fb_rate","off_rate","mixscore","tdev_fb","tdev_off","disrupt","mtdev",
        "whiffpct","hardhit","barrel","wcon","rv_swing")
cat(sprintf("\n  split-half within 2026 (n=%d) and 2025->2026 (n=%d)\n\n", nrow(B), nrow(C)))
R <- rbindlist(lapply(MS, function(v) {
  rb <- cor(B[[paste0(v,"_a")]], B[[paste0(v,"_b")]], use="complete.obs")
  rc <- cor(C[[paste0(v,"_a")]], C[[paste0(v,"_b")]], use="complete.obs")
  data.table(metric=v, split_half=round(rb,3), full_season=round(sb(rb),3),
             yoy=round(rc,3))
}))
print(R[order(-full_season)], row.names=FALSE)

# =============================================================================
# 5. Race against the old metric
# =============================================================================
cat("\n############ 5. Does the directional metric beat magnitude-only disruption? ############\n")
race <- function(M, lab) {
  cat(sprintf("\n  -- %s (n=%d) --\n", lab, nrow(M)))
  for (y in c("rv_swing","wcon","barrel","hardhit")) {
    cat(sprintf("     predicting later %-8s :", y))
    for (x in c("tscore","fb_rate","off_rate","disrupt")) {
      cat(sprintf("  %s %+.3f", x, cor(M[[paste0(x,"_a")]], M[[paste0(y,"_b")]],
                                       use="complete.obs")))
    }
    cat(sprintf("  | %s(own) %+.3f\n", y,
        cor(M[[paste0(y,"_a")]], M[[paste0(y,"_b")]], use="complete.obs")))
  }
  r2 <- function(f) summary(lm(as.formula(f), data=M))$r.squared
  for (y in c("barrel","hardhit","rv_swing")) {
    cat(sprintf("     R2 later %-8s: tscore %.3f | disrupt %.3f | ts+dis %.3f | halves %.3f | own %.3f | own+tscore %.3f | own+ts+dis %.3f\n",
      y, r2(sprintf("%s_b ~ tscore_a", y)), r2(sprintf("%s_b ~ disrupt_a", y)),
      r2(sprintf("%s_b ~ tscore_a + disrupt_a", y)),
      r2(sprintf("%s_b ~ fb_rate_a + off_rate_a", y)),
      r2(sprintf("%s_b ~ %s_a", y, y)),
      r2(sprintf("%s_b ~ %s_a + tscore_a", y, y)),
      r2(sprintf("%s_b ~ %s_a + tscore_a + disrupt_a", y, y))))
  }
}
race(B, "split-half within 2026")
race(C, "2025 -> 2026")

cat("\n  Same-season association (2026, min 250 swings):\n")
for (y in c("rv_swing","wcon","barrel","hardhit","whiffpct","xw")) {
  cat(sprintf("     %-9s vs tscore %+.3f | vs disrupt %+.3f\n", y,
    cor(P[[y]], P$tscore, use="complete.obs"),
    cor(P[[y]], P$disrupt, use="complete.obs")))
}

# =============================================================================
# 6. Individual pitches
# =============================================================================
cat("\n############ 6. Best individual pitches at hitting their group's target ############\n")
cat("Each pitch is scored against its OWN pitch-type league mean, so a slider is compared\n")
cat("to sliders. Target = deep for fastballs, out front for breaking/offspeed.\n")
PT <- sw[game_year==2026]
PT[, trv_pt := trv_b - mean(trv_b), by=pitch_type]
X <- PT[, .(n=.N, tscore=-100*mean(trv_pt), se=100*sd(trv_pt)/sqrt(.N),
            mtdev=mean(tdev), whiff=100*mean(whiff),
            barrel=100*sum(brl)/pmax(sum(isbip),1),
            wcon=sum(woba_value[isbip],na.rm=TRUE)/pmax(sum(woba_denom[isbip],na.rm=TRUE),1),
            rv=mean(pit_rv)), by=.(pitcher, pitch_type, pg2)][n >= 150]
X <- shrink(X, "tscore", "se")
X <- merge(X, NAME, by="pitcher")
show <- function(d, lab) {
  cat(sprintf("\n  -- %s --\n", lab))
  print(d[order(-tscore_adj)][1:12, .(player_name, pitch_type, n,
    tscore_adj=round(tscore_adj,2), mtdev=round(mtdev,1), whiff=round(whiff,1),
    barrel=round(barrel,1), wcon=round(wcon,3), rv=round(rv,4))], row.names=FALSE)
}
show(X[pg2=="FB"], "Fastballs: getting hitters DEEP")
show(X[pg2=="OFF"], "Breaking/offspeed: getting hitters OUT FRONT")
cat("\n  By pitch type, league-wide mean tdev and the direction that helps the pitcher:\n")
print(sw[game_year==2026, .(n=.N, mean_tdev=round(mean(tdev),2),
  pitcher_wants=fifelse(pg2[1]=="FB","deep (negative)","out front (positive)"),
  whiff=round(100*mean(whiff),1)), by=.(pitch_type,pg2)][n>=3000][order(pg2,-n)],
  row.names=FALSE)

fwrite(P[order(-tscore_adj)], file.path("data","statcast_2026","timing_directional_2026.csv"))
cat("\nWrote data/statcast_2026/timing_directional_2026.csv\n")
