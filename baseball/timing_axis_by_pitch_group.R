#!/usr/bin/env Rscript

# The timing axis split by pitch group, scored in REALIZED wOBA and run value.
#
# Hypothesis under test: a fastball caught out front is worth a lot of runs, and a breaking
# ball caught deep is ALSO worth runs -- because staying back on a slider and driving it the
# other way beats the defensive alignment, which an EV/LA model cannot see.
#
# Everything is scored on realized outcomes:
#   wOBAcon  = sum(woba_value)/sum(woba_denom) over balls in play
#   bat_rv   = delta_run_exp, positive = good for the OFFENSE (the user's "worth runs")
# xwOBAcon is carried only as a reference column so the model's blind spot is visible.
#
# Spray angle is derived from hit coordinates and flipped so positive = pull side for both
# handednesses, which lets the "beat the shift" mechanism be checked rather than assumed.

suppressPackageStartupMessages({ library(data.table); library(splines) })

cols <- c("game_year","game_type","pitcher","batter","stand","p_throws","pitch_type",
          "description","bb_type","events","balls","strikes","plate_x","plate_z",
          "sz_top","sz_bot","release_speed","launch_speed","launch_angle",
          "launch_speed_angle","estimated_woba_using_speedangle","woba_value","woba_denom",
          "delta_run_exp","hc_x","hc_y","if_fielding_alignment",
          "intercept_ball_minus_batter_pos_y_inches")

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
sw[, pg2 := fifelse(pgrp=="FB", "Fastball", "Offspeed")]
sw[, px_bat := fifelse(stand=="R", -plate_x, plate_x)]
sw[, pz_rel := (plate_z - sz_bot)/pmax(sz_top - sz_bot, 0.1)]
fit <- lm(depth ~ ns(px_bat,5)*pgrp + ns(pz_rel,5) + ns(release_speed,4) + stand, data=sw)
sw[, r1 := residuals(fit)]
sw[, nb := .N, by=batter]; sw <- sw[nb >= 200]
sw[, tdev := r1 - mean(r1), by=batter]

bp <- sw[description == "hit_into_play" & bb_type != "" & events != "" &
         is.finite(launch_speed) & is.finite(launch_angle) & is.finite(launch_speed_angle)]
bp[, hr     := events == "home_run"]
bp[, hit1b3 := events %in% c("single","double","triple")]
bp[, bat_rv := delta_run_exp]          # positive = good for the offense
bp[, barrel := launch_speed_angle == 6]

# spray angle: positive = pull side for both hands
bp[, spray_raw := atan2(hc_x - 125.42, 198.27 - hc_y) * 180/pi]
bp[, pull_ang  := fifelse(stand == "R", -spray_raw, spray_raw)]
bp[, has_spray := is.finite(pull_ang) & abs(pull_ang) <= 90]

ZB <- c(-Inf,-12,-6,-3,3,6,12,Inf)
ZL <- c("Deep >12","Deep 6-12","Deep 3-6","On time +/-3","Front 3-6","Front 6-12","Front >12")
bp[, zone := cut(tdev, ZB, labels=ZL)]

cat(sprintf("Balls in play: %s | spray coords present on %.1f%%\n",
  format(nrow(bp), big.mark=","), 100*mean(bp$has_spray)))
cat(sprintf("League wOBAcon %.3f | xwOBAcon %.3f | bat_rv per BIP %+.4f\n\n",
  sum(bp$woba_value,na.rm=TRUE)/sum(bp$woba_denom,na.rm=TRUE),
  mean(bp$estimated_woba_using_speedangle,na.rm=TRUE), mean(bp$bat_rv)))
cat("Sign convention: bat_rv is from the OFFENSE's side. Positive = the hitter gained runs.\n\n")

W <- function(d) sum(d$woba_value,na.rm=TRUE)/sum(d$woba_denom,na.rm=TRUE)

# =============================================================================
cat("############ 1. Where hitters make contact against each pitch group ############\n")
cat("tdev is already residualized within pitch group, so these means show what is LEFT over.\n\n")
print(bp[, .(bip=.N, share=round(100*.N/nrow(bp),1), mean_tdev=round(mean(tdev),2),
             sd_tdev=round(sd(tdev),2), pct_deep12=round(100*mean(tdev <= -12),1),
             pct_front12=round(100*mean(tdev >= 12),1),
             raw_depth=round(mean(depth),1)), by=pgrp][order(-bip)], row.names=FALSE)

# =============================================================================
cat("\n############ 2. Realized wOBAcon and run value by pitch group x timing zone ############\n\n")
G <- bp[, .(bip=.N,
            wOBAcon = W(.SD),
            xwOBAcon = mean(estimated_woba_using_speedangle, na.rm=TRUE),
            BABIP = sum(hit1b3)/pmax(sum(!hr),1),
            HRpct = 100*mean(hr),
            GB = 100*mean(bb_type=="ground_ball"),
            barrel = 100*mean(barrel),
            bat_rv = mean(bat_rv)),
        by=.(pg2, zone), .SDcols=c("woba_value","woba_denom")][order(pg2, zone)]
for (g in c("Fastball","Offspeed")) {
  cat(sprintf("  -- %s (n = %s balls in play) --\n", g,
              format(sum(G[pg2==g]$bip), big.mark=",")))
  print(G[pg2==g, .(zone, bip, wOBAcon=round(wOBAcon,3), xwOBAcon=round(xwOBAcon,3),
                    gap=round(wOBAcon-xwOBAcon,3), BABIP=round(BABIP,3),
                    HRpct=round(HRpct,1), GB=round(GB,1), barrel=round(barrel,1),
                    bat_rv=round(bat_rv,4))], row.names=FALSE)
  cat("\n")
}

cat("  Same cells, three-way pitch group, wOBAcon only:\n\n")
G3 <- bp[, .(bip=.N, wOBAcon=W(.SD), bat_rv=mean(bat_rv)),
         by=.(pgrp, zone), .SDcols=c("woba_value","woba_denom")]
print(dcast(G3, zone ~ pgrp, value.var="wOBAcon")[
  , lapply(.SD, function(x) if (is.numeric(x)) round(x,3) else x)], row.names=FALSE)
cat("\n  and run value per ball in play (offense positive):\n\n")
print(dcast(G3, zone ~ pgrp, value.var="bat_rv")[
  , lapply(.SD, function(x) if (is.numeric(x)) round(x,4) else x)], row.names=FALSE)

# =============================================================================
cat("\n############ 3. The hypothesis, stated as two cells ############\n\n")
q <- function(g, z, lab) {
  d <- bp[pg2==g & zone %in% z]
  data.table(cell = lab, bip = nrow(d), wOBAcon = round(W(d),3),
    xwOBAcon = round(mean(d$estimated_woba_using_speedangle, na.rm=TRUE),3),
    gap = round(W(d) - mean(d$estimated_woba_using_speedangle, na.rm=TRUE),3),
    BABIP = round(sum(d$hit1b3)/pmax(sum(!d$hr),1),3),
    HRpct = round(100*mean(d$hr),1), GB = round(100*mean(d$bb_type=="ground_ball"),1),
    bat_rv = round(mean(d$bat_rv),4),
    pull = round(mean(d[has_spray==TRUE]$pull_ang),1))
}
print(rbindlist(list(
  q("Fastball", c("Front 6-12","Front >12"), "Fastball caught out front (>6 in)"),
  q("Offspeed", c("Front 6-12","Front >12"), "Offspeed caught out front (>6 in)"),
  q("Fastball", c("Deep >12","Deep 6-12"),   "Fastball caught deep (>6 in)"),
  q("Offspeed", c("Deep >12","Deep 6-12"),   "Offspeed caught deep (>6 in)"),
  q("Fastball", "On time +/-3",              "Fastball on time"),
  q("Offspeed", "On time +/-3",              "Offspeed on time"))), row.names=FALSE)

cat("\n  Deep contact, offspeed minus fastball:\n")
do <- bp[pg2=="Offspeed" & tdev <= -6]; df <- bp[pg2=="Fastball" & tdev <= -6]
cat(sprintf("     wOBAcon  %.3f vs %.3f  -> %+.3f\n", W(do), W(df), W(do)-W(df)))
cat(sprintf("     bat_rv   %+.4f vs %+.4f -> %+.4f\n",
            mean(do$bat_rv), mean(df$bat_rv), mean(do$bat_rv)-mean(df$bat_rv)))
tt <- t.test(do$bat_rv, df$bat_rv)
cat(sprintf("     run-value difference p = %.3g\n", tt$p.value))
cat("\n  Out-front contact, fastball minus offspeed:\n")
fo <- bp[pg2=="Offspeed" & tdev >= 6]; ff <- bp[pg2=="Fastball" & tdev >= 6]
cat(sprintf("     wOBAcon  %.3f vs %.3f  -> %+.3f\n", W(ff), W(fo), W(ff)-W(fo)))
cat(sprintf("     bat_rv   %+.4f vs %+.4f -> %+.4f\n",
            mean(ff$bat_rv), mean(fo$bat_rv), mean(ff$bat_rv)-mean(fo$bat_rv)))
cat(sprintf("     run-value difference p = %.3g\n", t.test(ff$bat_rv, fo$bat_rv)$p.value))

# =============================================================================
cat("\n############ 4. Is the interaction real? ############\n")
bp[, zc := relevel(factor(as.character(zone)), ref="On time +/-3")]
bp[, pgf := relevel(factor(pg2), ref="Fastball")]
m0 <- lm(bat_rv ~ zc + pgf, data=bp)
m1 <- lm(bat_rv ~ zc * pgf, data=bp)
an <- anova(m0, m1)
cat(sprintf("  Run value: adding pitch-group x zone interaction, F = %.2f, p = %.3g\n",
            an$F[2], an$`Pr(>F)`[2]))
cat("  Interaction terms (offspeed relative to fastball, within each zone):\n")
co <- summary(m1)$coefficients
ix <- grep(":pgfOffspeed", rownames(co))
print(round(co[ix, c(1,2,4)], 4))

# =============================================================================
cat("\n############ 5. Runs on the table ############\n")
cat("Run value above the league BIP average, per 100 balls in play, and the 2026 season total.\n\n")
lgrv <- mean(bp$bat_rv)
R <- bp[game_year==2026, .(bip=.N, rv=mean(bat_rv)), by=.(pg2, zone)]
R[, `:=`(per100 = 100*(rv - lgrv), season_runs = bip*(rv - lgrv))]
print(R[order(-abs(season_runs))][, .(pg2, zone, bip, rv=round(rv,4),
  per100=round(per100,2), season_runs=round(season_runs,0))], row.names=FALSE)
cat(sprintf("\n  Total 2026 offensive runs from out-front contact (>6 in) vs deep (>6 in):\n"))
for (g in c("Fastball","Offspeed")) {
  a <- bp[game_year==2026 & pg2==g & tdev >= 6]; b <- bp[game_year==2026 & pg2==g & tdev <= -6]
  cat(sprintf("     %-9s out front %+7.0f runs over %s BIP | deep %+7.0f runs over %s BIP\n",
    g, sum(a$bat_rv - lgrv), format(nrow(a), big.mark=","),
    sum(b$bat_rv - lgrv), format(nrow(b), big.mark=",")))
}

# =============================================================================
cat("\n############ 6. Does deep contact beat the alignment? ############\n")
cat("Spray angle is flipped so positive = pull side for both hands.\n\n")
SP <- bp[has_spray==TRUE, .(bip=.N, pull_ang=mean(pull_ang),
  oppo=100*mean(pull_ang < -15), mid=100*mean(abs(pull_ang) <= 15),
  pull=100*mean(pull_ang > 15), wOBAcon=W(.SD),
  xwOBAcon=mean(estimated_woba_using_speedangle,na.rm=TRUE)),
  by=.(pg2, zone), .SDcols=c("woba_value","woba_denom")][order(pg2, zone)]
for (g in c("Fastball","Offspeed")) {
  cat(sprintf("  -- %s --\n", g))
  print(SP[pg2==g, .(zone, bip, mean_pull_angle=round(pull_ang,1), oppo=round(oppo,1),
    mid=round(mid,1), pull=round(pull,1), gap=round(wOBAcon-xwOBAcon,3))], row.names=FALSE)
  cat("\n")
}
cat("  Correlation across the 14 cells between opposite-field rate and the wOBA-xwOBA gap:\n")
SP[, gap := wOBAcon - xwOBAcon]
cat(sprintf("     r = %+.3f\n", cor(SP$oppo, SP$gap)))

cat("\n  Ground balls only, by infield alignment -- the direct shift test:\n\n")
gbs <- bp[bb_type=="ground_ball" & if_fielding_alignment != ""]
gbs[, tzone := fifelse(tdev <= -6, "Deep >6", fifelse(tdev >= 6, "Front >6", "On time"))]
print(gbs[, .(bip=.N, wOBAcon=round(W(.SD),3), BABIP=round(sum(hit1b3)/pmax(sum(!hr),1),3),
  bat_rv=round(mean(bat_rv),4)), by=.(if_fielding_alignment, tzone),
  .SDcols=c("woba_value","woba_denom")][order(if_fielding_alignment, tzone)], row.names=FALSE)
cat("\n  Deep-minus-on-time BABIP on grounders, within each alignment:\n")
for (al in sort(unique(gbs$if_fielding_alignment))) {
  d <- gbs[if_fielding_alignment==al]
  a <- d[tzone=="Deep >6"]; b <- d[tzone=="On time"]
  if (nrow(a) > 200 && nrow(b) > 200)
    cat(sprintf("     %-16s %+.3f  (deep %.3f on %s BIP vs on-time %.3f)\n", al,
      sum(a$hit1b3)/sum(!a$hr) - sum(b$hit1b3)/sum(!b$hr),
      sum(a$hit1b3)/sum(!a$hr), format(nrow(a), big.mark=","),
      sum(b$hit1b3)/sum(!b$hr)))
}
