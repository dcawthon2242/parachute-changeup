#!/usr/bin/env Rscript

# Actual outcomes along the timing axis, not just the expected-value model.
#
# xwOBAcon is estimated from exit velocity and launch angle alone, so it treats every ball
# with the same EV/LA identically. BABIP is what actually happened, and it EXCLUDES home
# runs from both numerator and denominator. Those two facts make the deep and out-front
# tails look very different from each other depending which metric you use, which is worth
# seeing explicitly before drawing conclusions about "weak contact".

suppressPackageStartupMessages({ library(data.table); library(splines) })

cols <- c("game_year","game_type","pitcher","batter","stand","pitch_type","description",
          "bb_type","events","balls","strikes","plate_x","plate_z","sz_top","sz_bot",
          "release_speed","launch_speed","launch_angle","launch_speed_angle",
          "estimated_woba_using_speedangle","woba_value","woba_denom","babip_value",
          "delta_run_exp","intercept_ball_minus_batter_pos_y_inches")

dt <- rbindlist(lapply(c(2025, 2026), function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress = FALSE, select = cols)))
setnames(dt, "intercept_ball_minus_batter_pos_y_inches", "depth")
dt <- dt[game_type == "R" & pitch_type != "" & balls <= 3 & strikes <= 2 &
         is.finite(delta_run_exp)]
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

bp <- sw[description == "hit_into_play" & bb_type != "" & events != "" &
         is.finite(launch_speed) & is.finite(launch_angle) & is.finite(launch_speed_angle)]
bp[, hr     := events == "home_run"]
bp[, hit    := events %in% c("single","double","triple","home_run")]
bp[, hit1b3 := events %in% c("single","double","triple")]
bp[, hard   := launch_speed >= 95]
bp[, barrel := launch_speed_angle == 6]

cat(sprintf("Balls in play, 2025-2026: %s\n", format(nrow(bp), big.mark=",")))
cat(sprintf("League: xwOBAcon %.3f | wOBAcon %.3f | BABIP %.3f | HR/BIP %.1f%% | hit%% %.1f\n\n",
  mean(bp$estimated_woba_using_speedangle, na.rm=TRUE),
  sum(bp$woba_value, na.rm=TRUE)/sum(bp$woba_denom, na.rm=TRUE),
  sum(bp$hit1b3)/sum(!bp$hr), 100*mean(bp$hr), 100*mean(bp$hit)))
cat("BABIP = (singles + doubles + triples) / (balls in play - home runs). Home runs are\n")
cat("removed from BOTH sides, so BABIP is blind to exactly the outcome the out-front tail\n")
cat("of the timing axis produces most of.\n\n")

BINS <- c(-Inf,-15,-12,-9,-6,-3,0,3,6,9,12,15,Inf)
LB <- c("<= -15","-15 to -12","-12 to -9","-9 to -6","-6 to -3","-3 to 0",
        "0 to +3","+3 to +6","+6 to +9","+9 to +12","+12 to +15","> +15")
bp[, tb := cut(tdev, BINS, labels=LB)]

cat("############ Outcomes by timing displacement (inches from the hitter's own norm) ############\n\n")
S <- bp[, .(
  bip     = .N,
  ev      = mean(launch_speed),
  la      = mean(launch_angle),
  xwcon   = mean(estimated_woba_using_speedangle, na.rm=TRUE),
  wcon    = sum(woba_value, na.rm=TRUE)/sum(woba_denom, na.rm=TRUE),
  babip   = sum(hit1b3)/pmax(sum(!hr),1),
  hrpct   = 100*mean(hr),
  hitpct  = 100*mean(hit),
  gb      = 100*mean(bb_type=="ground_ball"),
  ld      = 100*mean(bb_type=="line_drive"),
  fb      = 100*mean(bb_type=="fly_ball"),
  pop     = 100*mean(bb_type=="popup"),
  barrel  = 100*mean(barrel),
  rv      = mean(pit_rv)
), by=tb][order(tb)]
print(S[, .(tb, bip, ev=round(ev,1), la=round(la,1), xwOBAcon=round(xwcon,3),
            wOBAcon=round(wcon,3), BABIP=round(babip,3), HRpct=round(hrpct,1),
            barrel=round(barrel,1))], row.names=FALSE)

cat("\n  Batted-ball mix by the same bins:\n\n")
print(S[, .(tb, bip, GB=round(gb,1), LD=round(ld,1), FB=round(fb,1), POP=round(pop,1),
            hitpct=round(hitpct,1), pitcher_rv=round(rv,4))], row.names=FALSE)

cat("\n############ The two tails side by side ############\n\n")
deep <- bp[tdev <= -12]; front <- bp[tdev >= 6 & tdev <= 12]; xfront <- bp[tdev > 15]
mid  <- bp[abs(tdev) <= 3]
cmp <- function(d, lab) data.table(
  zone = lab, bip = nrow(d),
  ev = round(mean(d$launch_speed),1), la = round(mean(d$launch_angle),1),
  xwOBAcon = round(mean(d$estimated_woba_using_speedangle, na.rm=TRUE),3),
  wOBAcon  = round(sum(d$woba_value,na.rm=TRUE)/sum(d$woba_denom,na.rm=TRUE),3),
  BABIP = round(sum(d$hit1b3)/pmax(sum(!d$hr),1),3),
  HRpct = round(100*mean(d$hr),1),
  GB = round(100*mean(d$bb_type=="ground_ball"),1),
  POP = round(100*mean(d$bb_type=="popup"),1),
  barrel = round(100*mean(d$barrel),1),
  pitcher_rv = round(mean(d$pit_rv),4))
print(rbindlist(list(
  cmp(deep,   "Deep / late (<= -12 in)"),
  cmp(mid,    "Well timed (|tdev| <= 3)"),
  cmp(front,  "Out front (+6 to +12 in)"),
  cmp(xfront, "Extreme out front (> +15)"))), row.names=FALSE)

cat("\n  Note how differently the two metrics rank the tails:\n")
cat(sprintf("    xwOBAcon: deep %.3f vs out front %.3f  -> out front is %.0f%% worse for the pitcher\n",
  mean(deep$estimated_woba_using_speedangle, na.rm=TRUE),
  mean(front$estimated_woba_using_speedangle, na.rm=TRUE),
  100*(mean(front$estimated_woba_using_speedangle,na.rm=TRUE)/
       mean(deep$estimated_woba_using_speedangle,na.rm=TRUE) - 1)))
cat(sprintf("    BABIP:    deep %.3f vs out front %.3f  -> difference of only %+.3f\n",
  sum(deep$hit1b3)/sum(!deep$hr), sum(front$hit1b3)/sum(!front$hr),
  sum(front$hit1b3)/sum(!front$hr) - sum(deep$hit1b3)/sum(!deep$hr)))
cat(sprintf("    HR rate:  deep %.1f%% vs out front %.1f%%  -> the gap BABIP throws away\n",
  100*mean(deep$hr), 100*mean(front$hr)))

cat("\n############ BABIP peaks where xwOBAcon does not ############\n")
cat("BABIP is maximized in the line-drive band, which sits closer to neutral timing than the\n")
cat("barrel band does, because barrels leave the park and stop counting.\n\n")
print(S[order(-babip)][1:5, .(tb, bip, BABIP=round(babip,3), xwOBAcon=round(xwcon,3),
        LD=round(ld,1), HRpct=round(hrpct,1))], row.names=FALSE)
cat("\n  And where xwOBAcon is highest:\n\n")
print(S[order(-xwcon)][1:5, .(tb, bip, xwOBAcon=round(xwcon,3), BABIP=round(babip,3),
        barrel=round(barrel,1), HRpct=round(hrpct,1))], row.names=FALSE)

cat("\n############ Pitcher-level: is BABIP allowed a repeatable skill at all? ############\n")
set.seed(31)
bp[, half := sample(rep_len(1:2, .N)), by=pitcher]
pag <- function(d) d[, .(bip=.N, babip=sum(hit1b3)/pmax(sum(!hr),1),
  xwcon=mean(estimated_woba_using_speedangle,na.rm=TRUE),
  wcon=sum(woba_value,na.rm=TRUE)/sum(woba_denom,na.rm=TRUE),
  barrel=100*mean(barrel), mtdev=mean(tdev)), by=pitcher]
h1 <- pag(bp[game_year==2026 & half==1])[bip>=50]
h2 <- pag(bp[game_year==2026 & half==2])[bip>=50]
m <- merge(h1, h2, by="pitcher", suffixes=c("_a","_b"))
sb <- function(r) 2*r/(1+r)
for (v in c("mtdev","barrel","xwcon","wcon","babip")) {
  r <- cor(m[[paste0(v,"_a")]], m[[paste0(v,"_b")]], use="complete.obs")
  cat(sprintf("  %-7s split-half r = %+.3f -> full season %+.3f\n", v, r, sb(r)))
}
y1 <- pag(bp[game_year==2025])[bip>=120]; y2 <- pag(bp[game_year==2026])[bip>=120]
my <- merge(y1, y2, by="pitcher", suffixes=c("_a","_b"))
cat(sprintf("\n  2025 -> 2026 (n = %d): BABIP r = %+.3f | xwOBAcon r = %+.3f | mean tdev r = %+.3f\n",
  nrow(my), cor(my$babip_a, my$babip_b), cor(my$xwcon_a, my$xwcon_b),
  cor(my$mtdev_a, my$mtdev_b)))
cat(sprintf("  cor(mean tdev, BABIP allowed) same season = %+.3f | with xwOBAcon = %+.3f\n",
  cor(y2$mtdev, y2$babip), cor(y2$mtdev, y2$xwcon)))
