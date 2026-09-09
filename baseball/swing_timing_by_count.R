#!/usr/bin/env Rscript

# Count-by-count: should a pitcher aim at SOFT CONTACT or at SWING-AND-MISS?
#
# The two objectives sit at opposite ends of the same timing axis, so the choice is
# directional: push the hitter LATE (deep contact) or push him EARLY (out front).
#
# Weak contact is now scored with xwOBAcon rather than soft-ground-ball rate, so popups
# and other weak air contact count as weak contact. Bins are 3 inches wide with extended
# tails, to locate the xwOBAcon minimum properly rather than lumping the whole tail.
#
# Timing axis: tdev = intercept depth - E[depth | location, velo, pitch group], demeaned
# per hitter. Positive = further out front than that hitter's norm = EARLY.

suppressPackageStartupMessages({ library(data.table); library(splines) })

SEASONS <- c(2025, 2026)
cols <- c("game_type","batter","stand","pitch_type","description","events","bb_type",
          "balls","strikes","plate_x","plate_z","sz_top","sz_bot","release_speed",
          "launch_speed","launch_angle","estimated_woba_using_speedangle","woba_value",
          "woba_denom","delta_run_exp","intercept_ball_minus_batter_pos_y_inches")

dt <- rbindlist(lapply(SEASONS, function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress = FALSE, select = cols)))
setnames(dt, "intercept_ball_minus_batter_pos_y_inches", "depth")

SWING <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip","hit_into_play")
sw <- dt[game_type=="R" & description %in% SWING & is.finite(depth) & is.finite(plate_x) &
         is.finite(plate_z) & is.finite(release_speed) & is.finite(delta_run_exp) &
         balls <= 3 & strikes <= 2]
FB<-c("FF","SI","FC"); BR<-c("SL","ST","CU","KC","SV","CS"); OS<-c("CH","FS","FO")
sw[, pgrp := fifelse(pitch_type %in% FB,"FB", fifelse(pitch_type %in% BR,"BR",
             fifelse(pitch_type %in% OS,"OS",NA_character_)))]
sw <- sw[!is.na(pgrp)]
sw[, px_bat := fifelse(stand=="R", -plate_x, plate_x)]
sw[, pz_rel := (plate_z - sz_bot)/pmax(sz_top - sz_bot, 0.1)]
sw[, pit_rv := -delta_run_exp]
sw[, whiff := description %in% c("swinging_strike","swinging_strike_blocked","foul_tip")]
sw[, bipf  := description == "hit_into_play" & bb_type != "" & is.finite(launch_speed)]
sw[, cnt   := paste0(balls,"-",strikes)]

fit <- lm(depth ~ ns(px_bat,5)*pgrp + ns(pz_rel,5) + ns(release_speed,4) + stand, data=sw)
sw[, r1 := residuals(fit)]; sw[, nb := .N, by=batter]; sw <- sw[nb>=200]
sw[, tdev := r1 - mean(r1), by=batter]

E <- c(-Inf,-15,-12,-9,-6,-3,0,3,6,9,12,15,Inf)
L <- c("<-15","-15..-12","-12..-9","-9..-6","-6..-3","-3..0",
       "0..3","3..6","6..9","9..12","12..15",">15")
sw[, bin := cut(tdev, E, labels=L)]
CNTS <- c("0-0","0-1","0-2","1-0","1-1","1-2","2-0","2-1","2-2","3-0","3-1","3-2")
sw[, cnt := factor(cnt, levels=CNTS)]

cat(sprintf("Competitive swings: %s | balls in play: %s | hitters: %d\n\n",
  format(nrow(sw), big.mark=","), format(sw[bipf==TRUE,.N], big.mark=","),
  uniqueN(sw$batter)))

# =============================================================================
# A. The weak-contact curve, now scored on xwOBAcon so popups count
# =============================================================================
bip <- sw[bipf==TRUE]
cat("############ A. Where is xwOBAcon minimised across the timing axis? ############\n")
cat("(all counts pooled; popup share shown to confirm weak air contact is being credited)\n")
A <- bip[, .(bip=.N,
  gb=round(100*mean(bb_type=="ground_ball"),1),
  popup=round(100*mean(bb_type=="popup"),1),
  fb=round(100*mean(bb_type=="fly_ball"),1),
  ld=round(100*mean(bb_type=="line_drive"),1),
  soft_gb=round(100*mean(bb_type=="ground_ball" & launch_speed<85),1),
  ev=round(mean(launch_speed),1),
  xwobacon=round(mean(estimated_woba_using_speedangle,na.rm=TRUE),3),
  wobacon=round(sum(woba_value,na.rm=TRUE)/pmax(sum(woba_denom,na.rm=TRUE),1),3)
), by=bin][order(bin)]
A <- merge(sw[, .(swings=.N, whiff_pct=round(100*mean(whiff),1),
                  rv_swing=round(mean(pit_rv),4)), by=bin], A, by="bin")[order(bin)]
print(A)
am <- A[bip>=1500]
cat(sprintf("\n  xwOBAcon minimum: %s (%.3f)   |   maximum: %s (%.3f)\n",
  am$bin[which.min(am$xwobacon)], min(am$xwobacon),
  am$bin[which.max(am$xwobacon)], max(am$xwobacon)))

W <- sw[, .(swings=.N, whiff=round(100*mean(whiff),1)), by=bin][order(bin)]
wm <- W[swings>=3000]
cat(sprintf("  whiff maximum  : %s (%.1f%%)\n",
  wm$bin[which.max(wm$whiff)], max(wm$whiff)))
cat("\n  => the two objectives point in OPPOSITE directions on this axis.\n")
cat("     SOFT CONTACT  = push the hitter LATE  (negative tdev, deep contact)\n")
cat("     SWING-AND-MISS = push the hitter EARLY (positive tdev, out front)\n")

# =============================================================================
# B. xwOBAcon by count x timing bin -- does the weak-contact target move?
# =============================================================================
cat("\n############ B. xwOBAcon by count x timing bin (cells with >=250 BIP) ############\n")
B <- bip[, .(n=.N, xw=mean(estimated_woba_using_speedangle,na.rm=TRUE)), by=.(cnt,bin)]
B[n<250, xw := NA_real_]
print(dcast(B, cnt ~ bin, value.var="xw")[, lapply(.SD, function(x)
  if (is.numeric(x)) round(x,3) else x)])

cat("\n  per-count xwOBAcon minimum (>=250 BIP in cell):\n")
Bm <- B[!is.na(xw)][order(cnt, xw)][, .SD[1], by=cnt]
print(Bm[, .(cnt, best_bin=bin, xwobacon=round(xw,3), bip=n)])

# =============================================================================
# C. whiff rate by count x timing bin
# =============================================================================
cat("\n############ C. Whiff % by count x timing bin (cells with >=400 swings) ############\n")
C <- sw[, .(n=.N, wh=100*mean(whiff)), by=.(cnt,bin)]
C[n<400, wh := NA_real_]
print(dcast(C, cnt ~ bin, value.var="wh")[, lapply(.SD, function(x)
  if (is.numeric(x)) round(x,1) else x)])

# =============================================================================
# D. pitcher run value per swing by count x timing bin
# =============================================================================
cat("\n############ D. Pitcher run value per swing, by count x timing bin ############\n")
D <- sw[, .(n=.N, rv=mean(pit_rv)), by=.(cnt,bin)]
D[n<400, rv := NA_real_]
print(dcast(D, cnt ~ bin, value.var="rv")[, lapply(.SD, function(x)
  if (is.numeric(x)) round(x,3) else x)])

# =============================================================================
# E. THE VERDICT. Five zones, because the two objectives do not sit symmetrically:
#    the whiff optimum is the EXTREME early tail, while the moderate early band is
#    the worst place on the axis for the pitcher. Lumping them together hides that.
# =============================================================================
sw[, zone := fifelse(tdev <= -12, "ext_late",
             fifelse(tdev <= -6,  "mod_late",
             fifelse(tdev <   6,  "middle",
             fifelse(tdev <  12,  "mod_early", "ext_early"))))]
ZL <- c("ext_late","mod_late","middle","mod_early","ext_early")
ZN <- c("late extreme (<=-12in)","late moderate (-12..-6)","middle (-6..+6)",
        "early moderate (+6..+12)","early extreme (>=+12in)")

cat("\n############ E. The five zones, pooled -- why the tails are not symmetric ############\n")
Z <- sw[, .(swings=.N, whiff=round(100*mean(whiff),1),
  bip=sum(bipf), xwobacon=round(mean(estimated_woba_using_speedangle[bipf],na.rm=TRUE),3),
  wobacon=round(sum(woba_value[bipf],na.rm=TRUE)/pmax(sum(woba_denom[bipf],na.rm=TRUE),1),3),
  rv=round(mean(pit_rv),4)), by=zone]
Z <- Z[match(ZL, zone)]; Z[, zone := ZN]
print(Z)

cat("\n############ E2. Per count: run value per swing in each zone ############\n")
Q <- sw[, .(n=.N, rv=mean(pit_rv), se=sd(pit_rv)/sqrt(.N),
            wh=100*mean(whiff),
            xw=mean(estimated_woba_using_speedangle[bipf],na.rm=TRUE)), by=.(cnt, zone)]
pz <- function(v, lab, dg=4) {
  m <- dcast(Q, cnt ~ zone, value.var=v)
  keep <- intersect(ZL, names(m))
  m <- m[, c("cnt", keep), with=FALSE]; setnames(m, c("cnt", keep))
  cat(sprintf("\n-- %s --\n", lab))
  print(m[, lapply(.SD, function(x) if (is.numeric(x)) round(x, dg) else x)])
}
pz("n",  "swings per cell", 0)
pz("rv", "pitcher run value per swing", 4)
pz("wh", "whiff % of swings", 1)
pz("xw", "xwOBAcon", 3)

cat("\n############ E3. Head to head: the two OPTIMA, count by count ############\n")
cat("SOFT CONTACT target  = late extreme  (tdev <= -12in), lowest xwOBAcon\n")
cat("SWING-AND-MISS target = early extreme (tdev >= +12in), highest whiff rate\n")
cat("Diff = soft-contact minus swing-and-miss run value per swing. Cells under 250 swings dropped.\n\n")
res <- rbindlist(lapply(CNTS, function(k){
  l <- Q[cnt==k & zone=="ext_late"]; e <- Q[cnt==k & zone=="ext_early"]
  if (!nrow(l) || !nrow(e) || l$n < 250 || e$n < 250) return(NULL)
  d <- l$rv - e$rv; sed <- sqrt(l$se^2 + e$se^2)
  data.table(cnt=k, n_soft=l$n, n_miss=e$n,
    rv_soft=round(l$rv,4), rv_miss=round(e$rv,4),
    wh_soft=round(l$wh,1), wh_miss=round(e$wh,1),
    xw_soft=round(l$xw,3), xw_miss=round(e$xw,3),
    diff=round(d,4), se=round(sed,4), z=round(d/sed,2),
    aim=fifelse(abs(d/sed) < 1.64, "either (tie)",
         fifelse(d > 0, "SOFT CONTACT", "SWING-AND-MISS")))
}))
print(res)

cat("\n############ F. Collapsed by strike count ############\n")
sw[, scls := paste0(strikes, " strike", fifelse(strikes==1,"","s"))]
V2 <- sw[, .(n=.N, rv=mean(pit_rv), se=sd(pit_rv)/sqrt(.N), wh=100*mean(whiff),
    xw=mean(estimated_woba_using_speedangle[bipf], na.rm=TRUE)), by=.(scls, zone)]
for (s in sort(unique(sw$scls))) {
  l <- V2[scls==s & zone=="ext_late"]; e <- V2[scls==s & zone=="ext_early"]
  m <- V2[scls==s & zone=="middle"]
  ml <- V2[scls==s & zone=="mod_late"]; me <- V2[scls==s & zone=="mod_early"]
  d <- l$rv - e$rv; sed <- sqrt(l$se^2 + e$se^2)
  cat(sprintf("\n  %s\n", s))
  cat(sprintf("    late extreme  %+.4f (whiff %4.1f%%, xwOBAcon %.3f, n=%s)\n", l$rv, l$wh, l$xw, format(l$n,big.mark=",")))
  cat(sprintf("    late moderate %+.4f (whiff %4.1f%%, xwOBAcon %.3f, n=%s)\n", ml$rv, ml$wh, ml$xw, format(ml$n,big.mark=",")))
  cat(sprintf("    middle        %+.4f (whiff %4.1f%%, xwOBAcon %.3f, n=%s)\n", m$rv, m$wh, m$xw, format(m$n,big.mark=",")))
  cat(sprintf("    early moderate%+.4f (whiff %4.1f%%, xwOBAcon %.3f, n=%s)\n", me$rv, me$wh, me$xw, format(me$n,big.mark=",")))
  cat(sprintf("    early extreme %+.4f (whiff %4.1f%%, xwOBAcon %.3f, n=%s)\n", e$rv, e$wh, e$xw, format(e$n,big.mark=",")))
  cat(sprintf("    soft-contact minus swing-and-miss = %+.4f (z=%+.1f) -> %s\n", d, d/sed,
    if (abs(d/sed) < 1.64) "either" else if (d>0) "SOFT CONTACT" else "SWING-AND-MISS"))
}

cat("\n############ G. Why the early route decays: whiff rate in each extreme zone by balls ############\n")
G <- sw[zone %in% c("ext_late","ext_early"),
        .(swings=.N, whiff=round(100*mean(whiff),1),
          xwobacon=round(mean(estimated_woba_using_speedangle[bipf],na.rm=TRUE),3),
          rv=round(mean(pit_rv),4)), by=.(balls, zone)]
print(dcast(G, balls ~ zone, value.var=c("swings","whiff","xwobacon","rv")))

fwrite(res, file.path("data","statcast_model","swing_timing_by_count.csv"))
fwrite(dcast(B, cnt ~ bin, value.var="xw"),
       file.path("data","statcast_model","swing_timing_xwobacon_matrix.csv"))
cat("\nWrote data/statcast_model/swing_timing_by_count.csv and _xwobacon_matrix.csv\n")
