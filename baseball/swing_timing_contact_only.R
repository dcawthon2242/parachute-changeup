#!/usr/bin/env Rscript

# Follow-up to swing_timing_soft_contact.R, fixing two problems:
#
# 1. On a whiff the "intercept point" is a closest-approach point, not a real contact
#    point, so extreme depth values are partly definitional. The soft-contact question
#    is therefore asked on BALLS IN PLAY ONLY, where depth is a true contact depth,
#    and the depth model is refit on contact alone.
# 2. Batted-ball rates are reported per BALL IN PLAY, not per swing, so they are not
#    contaminated by the contact rate varying across bins.
#
# Then the actual trade the question is about: what is a mistimed grounder worth to the
# pitcher relative to a whiff, count by count, once the extra outs are counted?

suppressPackageStartupMessages({ library(data.table); library(splines) })

SEASONS <- c(2025, 2026)
cols <- c("game_year","game_type","batter","stand","pitch_type","description","events",
          "bb_type","balls","strikes","outs_when_up","on_1b","plate_x","plate_z",
          "sz_top","sz_bot","release_speed","launch_speed","launch_angle",
          "estimated_woba_using_speedangle","delta_run_exp",
          "intercept_ball_minus_batter_pos_y_inches")

dt <- rbindlist(lapply(SEASONS, function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress = FALSE, select = cols)))
setnames(dt, "intercept_ball_minus_batter_pos_y_inches", "depth")

SWING <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip","hit_into_play")
sw <- dt[game_type == "R" & description %in% SWING & is.finite(depth) &
         is.finite(plate_x) & is.finite(plate_z) & is.finite(release_speed) &
         is.finite(delta_run_exp)]
FB <- c("FF","SI","FC"); BR <- c("SL","ST","CU","KC","SV","CS"); OS <- c("CH","FS","FO")
sw[, pgrp := fifelse(pitch_type %in% FB, "FB",
             fifelse(pitch_type %in% BR, "BR", fifelse(pitch_type %in% OS, "OS", NA_character_)))]
sw <- sw[!is.na(pgrp)]
sw[, px_bat := fifelse(stand == "R", -plate_x, plate_x)]
sw[, pz_rel := (plate_z - sz_bot) / pmax(sz_top - sz_bot, 0.1)]
sw[, pit_rv := -delta_run_exp]
one <- c("field_out","force_out","sac_fly","sac_bunt","fielders_choice_out","other_out","strikeout")
sw[, outs_made := fifelse(grepl("triple_play", events), 3L,
                  fifelse(grepl("double_play", events), 2L,
                  fifelse(events %in% one, 1L, 0L)))]
sw[, whiff := description %in% c("swinging_strike","swinging_strike_blocked","foul_tip")]

# =============================================================================
# PART 1 -- run value of each swing outcome, by count. Is a weak grounder really
# worth more to the pitcher than a swinging strike?
# =============================================================================
bip <- sw[description == "hit_into_play" & bb_type != "" & is.finite(launch_speed)]
sw[, outclass := fifelse(whiff, "whiff",
                fifelse(description == "foul", "foul",
                fifelse(bb_type == "ground_ball" & launch_speed <  85, "GB soft (<85)",
                fifelse(bb_type == "ground_ball" & launch_speed >= 85, "GB hard (85+)",
                fifelse(bb_type == "popup", "popup",
                fifelse(bb_type == "fly_ball", "fly ball",
                fifelse(bb_type == "line_drive", "line drive", NA_character_)))))))]

cat("############ PART 1: pitcher run value per swing outcome, by count ############\n")
cat("(positive = good for the pitcher; the units are runs)\n\n")
p1 <- sw[!is.na(outclass), .(n = .N, rv = round(mean(pit_rv), 4), outs = round(mean(outs_made), 3)),
         by = .(outclass, cnt = fifelse(strikes == 2, "2 strikes", "0-1 strikes"))]
print(dcast(p1, outclass ~ cnt, value.var = c("n","rv","outs")))

cat("\n  -> the comparison that matters:\n")
for (cs in c("0-1 strikes","2 strikes")) {
  d <- sw[!is.na(outclass) & fifelse(strikes==2,"2 strikes","0-1 strikes") == cs]
  w  <- d[outclass == "whiff", mean(pit_rv)]
  g  <- d[outclass == "GB soft (<85)", mean(pit_rv)]
  cat(sprintf("     %-12s  whiff = %+.4f runs | soft grounder = %+.4f runs | soft GB advantage = %+.4f\n",
              cs, w, g, g - w))
}

# =============================================================================
# PART 2 -- contact depth on BALLS IN PLAY ONLY
# =============================================================================
fitb <- lm(depth ~ ns(px_bat, 5) * pgrp + ns(pz_rel, 5) + ns(release_speed, 4) + stand,
           data = bip)
cat(sprintf("\nContact-only depth model: R2 = %.3f, residual SD = %.2f in (n=%s)\n",
            summary(fitb)$r.squared, sd(residuals(fitb)), format(nrow(bip), big.mark=",")))
bip[, r1 := residuals(fitb)]
bip[, nbat := .N, by = batter]
bip <- bip[nbat >= 100]
bip[, tdev := r1 - mean(r1), by = batter]

E <- c(-Inf,-9,-6,-3,0,3,6,9,Inf)
L <- c("<-9 (very late)","-9..-6 late","-6..-3","-3..0","0..3","3..6","6..9 early","9+ (very early)")
bip[, bin := cut(tdev, E, labels = L)]
bip[, dpopp := (on_1b > 0 & !is.na(on_1b)) & outs_when_up < 2]

cat("\n############ PART 2: batted-ball profile per BALL IN PLAY, by timing deviation ############\n")
cat("(negative tdev = LATE / deeper contact than this hitter's norm; positive = EARLY / out front)\n")
print(bip[, .(
  n        = .N,
  pct      = round(100*.N/nrow(bip),1),
  gb       = round(100*mean(bb_type=="ground_ball"),1),
  soft_gb  = round(100*mean(bb_type=="ground_ball" & launch_speed<85),1),
  ld       = round(100*mean(bb_type=="line_drive"),1),
  fb       = round(100*mean(bb_type=="fly_ball"),1),
  popup    = round(100*mean(bb_type=="popup"),1),
  barrel   = round(100*mean(launch_speed>=98 & launch_angle>=26 & launch_angle<=30),1),
  ev       = round(mean(launch_speed),1),
  hard95   = round(100*mean(launch_speed>=95),1),
  xwobacon = round(mean(estimated_woba_using_speedangle, na.rm=TRUE),3),
  outs     = round(mean(outs_made),3),
  rv_bip   = round(mean(pit_rv),4)
), by = bin][order(bin)])

cat("\n############ PART 3: double plays per ball in play (runner on 1B, <2 outs) ############\n")
print(bip[dpopp == TRUE, .(
  n       = .N,
  gb      = round(100*mean(bb_type=="ground_ball"),1),
  dp      = round(100*mean(grepl("double_play", events)),2),
  dp_pergb= round(100*sum(grepl("double_play", events)) /
                  pmax(sum(bb_type=="ground_ball"),1),2),
  outs    = round(mean(outs_made),3),
  rv_bip  = round(mean(pit_rv),4)
), by = bin][order(bin)])

# =============================================================================
# PART 4 -- put it together: expected pitcher run value per SWING as a function of
# timing, decomposed into the whiff channel and the contact channel.
# =============================================================================
cat("\n############ PART 4: per-swing decomposition, 0-1 strike counts ############\n")
sw2 <- sw[strikes < 2 & !is.na(outclass)]
fits <- lm(depth ~ ns(px_bat,5)*pgrp + ns(pz_rel,5) + ns(release_speed,4) + stand, data = sw2)
sw2[, r1 := residuals(fits)]; sw2[, nb := .N, by = batter]
sw2 <- sw2[nb >= 100]; sw2[, tdev := r1 - mean(r1), by = batter]
sw2[, bin := cut(tdev, E, labels = L)]
print(sw2[, .(
  n         = .N,
  whiff_pct = round(100*mean(whiff),1),
  rv_whiff  = round(mean(whiff) * mean(pit_rv[whiff]), 4),
  bip_pct   = round(100*mean(description=="hit_into_play"),1),
  rv_bip    = round(mean(description=="hit_into_play") *
                    mean(pit_rv[description=="hit_into_play"]), 4),
  rv_foul   = round(mean(description=="foul") * mean(pit_rv[description=="foul"]), 4),
  rv_total  = round(mean(pit_rv), 4)
), by = bin][order(bin)])

cat("\n############ PART 5: how much run value is on the table? ############\n")
b <- bip[, .(n=.N, rv=mean(pit_rv)), by=bin][order(bin)][n >= 3000]
best <- b[which.max(rv)]; worst <- b[which.min(rv)]
cat(sprintf("  Per ball in play, pitcher RV ranges from %+.4f (%s) to %+.4f (%s)\n",
            worst$rv, worst$bin, best$rv, best$bin))
cat(sprintf("  Spread = %.4f runs per ball in play. A starter allows ~18 BIP/start.\n",
            best$rv - worst$rv))
lateshift <- bip[tdev < -3, mean(pit_rv)] - bip[abs(tdev) <= 3, mean(pit_rv)]
cat(sprintf("  Moving contact from on-time (|tdev|<=3in) to late (tdev < -3in) is worth %+.4f runs/BIP.\n", lateshift))

out <- file.path("data","statcast_model","swing_timing_contact_only.csv")
fwrite(bip[, .(n=.N, gb=mean(bb_type=="ground_ball"),
               soft_gb=mean(bb_type=="ground_ball" & launch_speed<85),
               ev=mean(launch_speed), xwobacon=mean(estimated_woba_using_speedangle,na.rm=TRUE),
               outs=mean(outs_made), rv_bip=mean(pit_rv)), by=bin][order(bin)], out)
cat(sprintf("\nWrote %s\n", out))
