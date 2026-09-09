#!/usr/bin/env Rscript

# Final pass. Two things to settle:
#
# 1. Raw delta_run_exp per ball in play is not comparable across timing bins, because
#    the bins have different count mixes and the run-expectancy baseline moves with the
#    count (the same hard grounder is worth +0.034 to a pitcher at 0-1 strikes and
#    -0.011 at 2 strikes). So contact value is also reported count-neutral: run value
#    minus the mean run value of all balls in play from the same balls-strikes cell.
#
# 2. The per-swing objective: for each count state, does the timing band that maximises
#    total pitcher run value coincide with the band that maximises whiffs, or not?
#
# Also validates the sign of the timing axis against attack_direction and spray.

suppressPackageStartupMessages({ library(data.table); library(splines) })

SEASONS <- c(2025, 2026)
cols <- c("game_type","batter","stand","pitch_type","description","events","bb_type",
          "balls","strikes","outs_when_up","on_1b","plate_x","plate_z","sz_top","sz_bot",
          "release_speed","launch_speed","launch_angle","hc_x","hc_y","attack_direction",
          "estimated_woba_using_speedangle","woba_value","woba_denom","babip_value",
          "delta_run_exp","intercept_ball_minus_batter_pos_y_inches")

dt <- rbindlist(lapply(SEASONS, function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress = FALSE, select = cols)))
setnames(dt, "intercept_ball_minus_batter_pos_y_inches", "depth")

SWING <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip","hit_into_play")
sw <- dt[game_type=="R" & description %in% SWING & is.finite(depth) & is.finite(plate_x) &
         is.finite(plate_z) & is.finite(release_speed) & is.finite(delta_run_exp)]
FB<-c("FF","SI","FC"); BR<-c("SL","ST","CU","KC","SV","CS"); OS<-c("CH","FS","FO")
sw[, pgrp := fifelse(pitch_type %in% FB,"FB", fifelse(pitch_type %in% BR,"BR",
             fifelse(pitch_type %in% OS,"OS",NA_character_)))]
sw <- sw[!is.na(pgrp)]
sw[, px_bat := fifelse(stand=="R", -plate_x, plate_x)]
sw[, pz_rel := (plate_z - sz_bot)/pmax(sz_top - sz_bot, 0.1)]
sw[, pit_rv := -delta_run_exp]
sw[, whiff := description %in% c("swinging_strike","swinging_strike_blocked","foul_tip")]
sw[, bipf := description == "hit_into_play"]
one <- c("field_out","force_out","sac_fly","sac_bunt","fielders_choice_out","other_out","strikeout")
sw[, outs_made := fifelse(grepl("triple_play",events),3L, fifelse(grepl("double_play",events),2L,
                  fifelse(events %in% one,1L,0L)))]
sw[, cell := paste0(balls,"-",strikes)]
sw[, cs := fifelse(strikes==2, "2 strikes", "0-1 strikes")]

# one timing axis, fit on ALL swings so bins are comparable across outcome channels,
# then demeaned per hitter
fit <- lm(depth ~ ns(px_bat,5)*pgrp + ns(pz_rel,5) + ns(release_speed,4) + stand, data=sw)
sw[, r1 := residuals(fit)]; sw[, nb := .N, by=batter]; sw <- sw[nb>=200]
sw[, tdev := r1 - mean(r1), by=batter]
E <- c(-Inf,-9,-6,-3,0,3,6,9,Inf)
L <- c("<-9 very late","-9..-6 late","-6..-3 bit late","-3..0 on time","0..3 on time",
       "3..6 bit early","6..9 early","9+ very early")
sw[, bin := cut(tdev, E, labels=L)]

# ---- sign validation -------------------------------------------------------
cat("############ Sign check: does positive tdev really mean EARLY? ############\n")
bip <- sw[bipf==TRUE & bb_type!="" & is.finite(launch_speed)]
bip[, spray_pull := {
  ang <- atan2(hc_x - 125.42, 198.27 - hc_y) * 180/pi
  fifelse(stand=="R", -ang, ang)                       # positive = pulled
}]
print(bip[, .(n=.N, raw_depth_in=round(mean(depth),1),
              attack_dir=round(mean(attack_direction, na.rm=TRUE),1),
              spray_pull_deg=round(mean(spray_pull, na.rm=TRUE),1)), by=bin][order(bin)])
cat("  (attack_direction and spray both rise with tdev => positive tdev = out front / pulled = EARLY)\n")

# ---- count-neutral contact value ------------------------------------------
bip[, cellmean := mean(pit_rv), by=cell]
bip[, rv_adj := pit_rv - cellmean]
cat("\n############ Contact value per ball in play, count-neutral ############\n")
cat("(rv_adj > 0 = better for the pitcher than an average ball in play from the same count)\n")
print(bip[, .(n=.N,
   gb=round(100*mean(bb_type=="ground_ball"),1),
   soft_gb=round(100*mean(bb_type=="ground_ball" & launch_speed<85),1),
   ev=round(mean(launch_speed),1),
   xwobacon=round(mean(estimated_woba_using_speedangle,na.rm=TRUE),3),
   wobacon=round(sum(woba_value,na.rm=TRUE)/pmax(sum(woba_denom,na.rm=TRUE),1),3),
   babip=round(mean(babip_value,na.rm=TRUE),3),
   outs=round(mean(outs_made),3),
   rv_raw=round(mean(pit_rv),4),
   rv_adj=round(mean(rv_adj),4)), by=bin][order(bin)])

# ---- double plays, on the same timing axis ---------------------------------
bip[, dpopp := (on_1b > 0 & !is.na(on_1b)) & outs_when_up < 2]
cat("\n############ Double plays per ball in play (runner on 1B, <2 outs) ############\n")
print(bip[dpopp==TRUE, .(n=.N,
   gb=round(100*mean(bb_type=="ground_ball"),1),
   dp=round(100*mean(grepl("double_play", events)),2),
   dp_per_gb=round(100*sum(grepl("double_play", events))/pmax(sum(bb_type=="ground_ball"),1),1),
   outs=round(mean(outs_made),3),
   rv_adj=round(mean(rv_adj),4)), by=bin][order(bin)])

# ---- per-swing objective, by count ----------------------------------------
cat("\n############ Per-swing pitcher run value, by count state ############\n")
for (k in c("0-1 strikes","2 strikes")) {
  d <- sw[cs==k]
  tab <- d[, .(n=.N,
    whiff_pct = round(100*mean(whiff),1),
    foul_pct  = round(100*mean(description=="foul"),1),
    bip_pct   = round(100*mean(bipf),1),
    softgb_pct= round(100*mean(bipf & bb_type=="ground_ball" &
                               is.finite(launch_speed) & launch_speed<85),1),
    outs_sw   = round(mean(outs_made),3),
    rv_swing  = round(mean(pit_rv),4)), by=bin][order(bin)]
  cat(sprintf("\n-- %s (n=%s) --\n", k, format(nrow(d), big.mark=",")))
  print(tab)
  b <- tab[n>=5000]
  cat(sprintf("   best run value : %-16s (%.4f, whiff %.1f%%)\n",
              b$bin[which.max(b$rv_swing)], max(b$rv_swing),
              b$whiff_pct[which.max(b$rv_swing)]))
  cat(sprintf("   best whiff rate: %-16s (%.4f, whiff %.1f%%)\n",
              b$bin[which.max(b$whiff_pct)], b$rv_swing[which.max(b$whiff_pct)],
              max(b$whiff_pct)))
}

# ---- the headline trade ----------------------------------------------------
cat("\n############ The trade, stated directly ############\n")
for (k in c("0-1 strikes","2 strikes")) {
  d <- sw[cs==k]
  late  <- d[tdev < -6]; ontime <- d[abs(tdev) <= 3]; early <- d[tdev > 6]
  cat(sprintf("\n%s:\n", k))
  for (nm in c("late (tdev < -6in)","on time (|tdev| <= 3in)","early (tdev > +6in)")) {
    s <- switch(nm, "late (tdev < -6in)"=late, "on time (|tdev| <= 3in)"=ontime, early)
    cat(sprintf("   %-26s whiff %5.1f%%  soft-GB %4.1f%%  outs/swing %.3f  RV/swing %+.4f\n",
      nm, 100*mean(s$whiff),
      100*mean(s$bipf & s$bb_type=="ground_ball" & is.finite(s$launch_speed) & s$launch_speed<85),
      mean(s$outs_made), mean(s$pit_rv)))
  }
}

out <- file.path("data","statcast_model","swing_timing_verdict.csv")
fwrite(sw[, .(n=.N, whiff=mean(whiff), bip=mean(bipf),
   soft_gb=mean(bipf & bb_type=="ground_ball" & is.finite(launch_speed) & launch_speed<85),
   outs_sw=mean(outs_made), rv_swing=mean(pit_rv)), by=.(bin, cs)][order(cs,bin)], out)
cat(sprintf("\nWrote %s\n", out))
