#!/usr/bin/env Rscript

# Which pitchers are best at producing bad hitter timing?
#
# "Bad timing" is displacement away from the hitter's own contact-depth norm in EITHER
# direction, since both tails hurt the hitter (deep = weak contact, out front = whiffs).
# So the core metric is magnitude, not signed direction.
#
# Three complications handled explicitly:
#   (1) Some hitters are inherently easier to knock off balance. Raw mean |tdev| rewards a
#       pitcher for facing them, so disruption is measured net of each batter's own mean
#       |tdev| across all pitchers he faced. Units are inches of EXCESS displacement.
#   (2) The two tails are not worth the same amount. A second metric, timing-implied run
#       value, maps each swing's position on the axis to the league run value earned at that
#       position in that count, then averages. It answers "given only where this pitcher puts
#       hitters on the timing axis, what should that be worth?"
#   (3) Leaderboards are shrunk, and repeatability plus an out-of-sample run-value test are
#       reported so it is clear how much of the ranking is signal.

suppressPackageStartupMessages({ library(data.table); library(splines) })
set.seed(3)

cols <- c("game_year","game_type","player_name","pitcher","batter","stand","pitch_type",
          "description","bb_type","balls","strikes","plate_x","plate_z","sz_top","sz_bot",
          "release_speed","estimated_woba_using_speedangle","delta_run_exp","miss_distance",
          "intercept_ball_minus_batter_pos_y_inches")

dt <- rbindlist(lapply(c(2025, 2026), function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress = FALSE, select = cols)))
setnames(dt, "intercept_ball_minus_batter_pos_y_inches", "depth")
dt <- dt[game_type == "R" & pitch_type != "" & balls <= 3 & strikes <= 2 & is.finite(delta_run_exp)]
dt[, pit_rv := -delta_run_exp]

SWING <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip","hit_into_play")
sw <- dt[description %in% SWING & is.finite(depth) & is.finite(plate_x) &
         is.finite(plate_z) & is.finite(release_speed)]
FB<-c("FF","SI","FC"); BR<-c("SL","ST","CU","KC","SV","CS"); OS<-c("CH","FS","FO")
sw[, pgrp := fifelse(pitch_type %in% FB,"FB", fifelse(pitch_type %in% BR,"BR",
             fifelse(pitch_type %in% OS,"OS",NA_character_)))]
sw <- sw[!is.na(pgrp)]
sw[, px_bat := fifelse(stand=="R", -plate_x, plate_x)]
sw[, pz_rel := (plate_z - sz_bot)/pmax(sz_top - sz_bot, 0.1)]
fit <- lm(depth ~ ns(px_bat,5)*pgrp + ns(pz_rel,5) + ns(release_speed,4) + stand, data=sw)
sw[, r1 := residuals(fit)]
sw[, nb := .N, by=batter]; sw <- sw[nb >= 200]
sw[, tdev := r1 - mean(r1), by=batter]
sw[, whiff := description %in% c("swinging_strike","swinging_strike_blocked","foul_tip")]
sw[, bipf  := description == "hit_into_play" & bb_type != "" &
              is.finite(estimated_woba_using_speedangle)]
sw[, xw := estimated_woba_using_speedangle]

# ---- (1) disruption in inches, net of each batter's own disruptability ------
sw[, adev := abs(tdev)]
sw[, b_adev := mean(adev), by=batter]
sw[, excess := adev - b_adev]
cat(sprintf("Swings on axis: %s | league mean |tdev| = %.2f in | batter mean |tdev| ranges %.2f to %.2f\n",
  format(nrow(sw), big.mark=","), mean(sw$adev),
  min(sw[, mean(adev), by=batter]$V1), max(sw[, mean(adev), by=batter]$V1)))

# ---- (2) timing-implied run value ------------------------------------------
# League run value earned at each position on the axis, within count state, then applied back.
sw[, tbin := cut(tdev, c(-Inf, seq(-20, 20, 2), Inf))]
map <- sw[, .(bin_rv = mean(pit_rv), bin_n = .N), by=.(tbin, strikes)]
sw <- merge(sw, map, by=c("tbin","strikes"), all.x=TRUE)
cat(sprintf("Timing->run-value map: %d cells, min cell n = %d\n\n", nrow(map), min(map$bin_n)))

# ---- aggregate -------------------------------------------------------------
agg <- function(d, by) d[, .(
  swings   = .N,
  disrupt  = mean(excess),
  se_dis   = sd(excess)/sqrt(.N),
  trv      = mean(bin_rv),
  se_trv   = sd(bin_rv)/sqrt(.N),
  tdev     = mean(tdev),
  spread   = sd(tdev),
  late12   = 100*mean(tdev <= -12),
  early12  = 100*mean(tdev >=  12),
  ext      = 100*mean(abs(tdev) >= 12),
  whiffpct = 100*mean(whiff),
  miss     = mean(miss_distance[whiff], na.rm=TRUE),
  xwobacon = mean(xw[bipf], na.rm=TRUE),
  rv_swing = mean(pit_rv),
  se_rv    = sd(pit_rv)/sqrt(.N)
), by=by]

shrink <- function(d, val, sev) {
  d <- copy(d)
  v <- d[[val]]; s <- d[[sev]]
  ok <- is.finite(s)
  mu <- mean(v[ok]); tau2 <- pmax(var(v[ok]) - mean(s[ok]^2), 1e-10)
  d[, rel := fifelse(ok, tau2/(tau2 + s^2), 0)]
  d[[paste0(val,"_adj")]] <- mu + d$rel*(v - mu)
  cat(sprintf("  %-9s league mean %+.3f | true SD %.3f | median reliability %.2f\n",
    val, mu, sqrt(tau2), median(d$rel)))
  d
}

NAME <- unique(sw[, .(pitcher, player_name)], by="pitcher")
s26 <- sw[game_year == 2026]

cat("############ Shrinkage diagnostics, 2026 pitchers (min 250 swings) ############\n")
P26 <- agg(s26, "pitcher")[swings >= 250]
P26 <- shrink(P26, "disrupt", "se_dis")
P26 <- shrink(P26, "trv", "se_trv")
P26 <- shrink(P26, "rv_swing", "se_rv")
P26 <- merge(P26, NAME, by="pitcher")

fmt <- function(d, n=25, sortby="disrupt_adj") {
  d <- d[order(-get(sortby))][1:min(n, nrow(d))]
  d[, .(player_name,
        swings,
        disrupt = round(disrupt, 2),
        disrupt_adj = round(disrupt_adj, 2),
        trv = round(trv, 4),
        tdev = round(tdev, 1),
        spread = round(spread, 1),
        late12 = round(late12, 1),
        early12 = round(early12, 1),
        ext = round(ext, 1),
        whiff = round(whiffpct, 1),
        miss = round(miss, 2),
        xwobacon = round(xwobacon, 3),
        rv = round(rv_swing, 4))]
}

cat("\n############ 1. Best timing disruptors, 2026 (min 250 swings) ############\n")
cat("disrupt = inches of displacement above what the hitters he faced average against everyone.\n")
cat("ext = %% of swings pushed at least 12 inches off the hitter's own norm, either direction.\n\n")
print(fmt(P26, 25))

cat("\n  -- least disruptive 10 --\n")
print(P26[order(disrupt_adj)][1:10, .(player_name, swings, disrupt=round(disrupt,2),
  disrupt_adj=round(disrupt_adj,2), tdev=round(tdev,1), spread=round(spread,1),
  ext=round(ext,1), whiff=round(whiffpct,1), xwobacon=round(xwobacon,3), rv=round(rv_swing,4))])

cat("\n############ 1b. Starters only, 2026 (min 700 swings) ############\n")
cat("The list above is reliever-heavy because 250 swings is a full season for a bullpen arm.\n")
cat("Restricting to workhorse usage removes the platoon/short-burst advantage.\n\n")
ST <- agg(s26, "pitcher")[swings >= 700]
ST <- shrink(ST, "disrupt", "se_dis"); ST <- shrink(ST, "trv", "se_trv")
ST <- merge(ST, NAME, by="pitcher")
cat("\n")
print(fmt(ST, 20))

cat("\n############ 2. Same pitchers ranked on timing-implied run value ############\n")
cat("trv = league run value earned at each axis position in that count, averaged over his swings.\n")
cat("It credits the late tail and the far early tail, and penalises the +6 to +12 damage zone.\n\n")
print(fmt(P26, 25, "trv_adj"))

cat("\n############ 3. Two-season list, 2025+2026 pooled (min 500 swings) ############\n")
PB <- agg(sw, "pitcher")[swings >= 500]
PB <- shrink(PB, "disrupt", "se_dis")
PB <- shrink(PB, "trv", "se_trv")
PB <- merge(PB, NAME, by="pitcher")
cat("\n")
print(fmt(PB, 25))

cat("\n############ 4. Best individual pitches, 2026 (min 150 swings) ############\n")
PT <- agg(s26, c("pitcher","pitch_type"))[swings >= 150]
PT <- shrink(PT, "disrupt", "se_dis")
PT <- merge(PT, NAME, by="pitcher")
cat("\n")
print(PT[order(-disrupt_adj)][1:22, .(player_name, pitch_type, swings,
  disrupt=round(disrupt,2), disrupt_adj=round(disrupt_adj,2), tdev=round(tdev,1),
  spread=round(spread,1), late12=round(late12,1), early12=round(early12,1),
  whiff=round(whiffpct,1), miss=round(miss,2), xwobacon=round(xwobacon,3),
  rv=round(rv_swing,4))])

cat("\n############ 5. Is disruption repeatable, and does it pay? ############\n")
sw[, half := sample(rep_len(1:2, .N)), by=pitcher]
h1 <- agg(sw[game_year==2026 & half==1], "pitcher"); h2 <- agg(sw[game_year==2026 & half==2], "pitcher")
m <- merge(h1[swings>=120], h2[swings>=120], by="pitcher", suffixes=c("_a","_b"))
y1 <- agg(sw[game_year==2025], "pitcher"); y2 <- agg(sw[game_year==2026], "pitcher")
my <- merge(y1[swings>=300], y2[swings>=300], by="pitcher", suffixes=c("_a","_b"))
sb <- function(r) 2*r/(1+r)
cat(sprintf("  Repeatability of disruption:  split-half 2026 r = %+.3f (Spearman-Brown %+.3f, n=%d)\n",
  cor(m$disrupt_a, m$disrupt_b), sb(cor(m$disrupt_a, m$disrupt_b)), nrow(m)))
cat(sprintf("                                2025 -> 2026   r = %+.3f (n=%d)\n",
  cor(my$disrupt_a, my$disrupt_b), nrow(my)))
cat(sprintf("  Repeatability of timing-implied RV: split-half r = %+.3f | 2025->2026 r = %+.3f\n",
  cor(m$trv_a, m$trv_b), cor(my$trv_a, my$trv_b)))

cat("\n  Out-of-sample: does disruption in one half predict RUN VALUE in the other half?\n")
for (v in c("disrupt","trv","ext","whiffpct","spread")) {
  ct <- cor.test(m[[paste0(v,"_a")]], m$rv_swing_b)
  cat(sprintf("     %-9s -> other-half rv/swing  r = %+.3f (p = %.3g)\n", v, ct$estimate, ct$p.value))
}
cat("\n  And across seasons: 2025 profile -> 2026 run value?\n")
for (v in c("disrupt","trv","ext","whiffpct","spread")) {
  ct <- cor.test(my[[paste0(v,"_a")]], my$rv_swing_b)
  cat(sprintf("     %-9s -> 2026 rv/swing        r = %+.3f (p = %.3g)\n", v, ct$estimate, ct$p.value))
}

cat("\n############ 6. How disruptors do it: direction vs spread ############\n")
P26[, route := fifelse(late12 >= 10 & early12 >= 10, "both tails",
             fifelse(tdev <= -2, "pushes late",
             fifelse(tdev >=  2, "pushes early", "narrow / mixed")))]
print(P26[, .(pitchers=.N, disrupt=round(mean(disrupt),2), tdev=round(mean(tdev),1),
  spread=round(mean(spread),1), late12=round(mean(late12),1), early12=round(mean(early12),1),
  whiff=round(mean(whiffpct),1), miss=round(mean(miss),2),
  xwobacon=round(mean(xwobacon,na.rm=TRUE),3), trv=round(mean(trv),4),
  rv=round(mean(rv_swing),4)), by=route][order(-disrupt)])

fwrite(P26[order(-disrupt_adj)], file.path("data","statcast_2026","timing_disruptors_2026.csv"))
fwrite(PT[order(-disrupt_adj)], file.path("data","statcast_2026","timing_disruptors_2026_bypitch.csv"))
cat("\nWrote data/statcast_2026/timing_disruptors_2026.csv and _bypitch.csv\n")
