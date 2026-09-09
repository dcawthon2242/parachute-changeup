#!/usr/bin/env Rscript

# Does a pitcher's timing-disruption skill do anything for barrel suppression?
#
# The two prior results this joins:
#   (a) barrel% = hardhit% x conversion. 71% of between-pitcher barrel variance lives in
#       conversion (the launch-angle factor), and launch angle is the only batted-ball
#       quantity where the pitcher explains as much variance as the hitter.
#   (b) induced contact depth is the single most repeatable pitcher trait measured
#       (r = 0.85 year over year), while conversion is the least (r = 0.33).
#
# If disruption acts on the ANGLE channel, it is a repeatable leading indicator for the
# half of barrel rate that does not predict itself. That is the only way this is "usable".
#
# The forecast tests are built to have no mechanical overlap:
#   Test A  disruption measured ONLY on swings that did not become balls in play,
#           predicting contact quality on the balls in play. Same season, disjoint pitches.
#   Test B  split-half: first half disruption -> second half barrel rate.
#   Test C  2025 disruption -> 2026 barrel rate.

suppressPackageStartupMessages({ library(data.table); library(splines) })
set.seed(11)

cols <- c("game_year","game_type","player_name","pitcher","batter","stand","pitch_type",
          "description","bb_type","balls","strikes","plate_x","plate_z","sz_top","sz_bot",
          "release_speed","launch_speed","launch_angle","launch_speed_angle","delta_run_exp",
          "estimated_woba_using_speedangle","intercept_ball_minus_batter_pos_y_inches")

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
sw[, adev := abs(tdev)]
sw[, b_adev := mean(adev), by=batter]
sw[, excess := adev - b_adev]

sw[, isbip := description == "hit_into_play" & bb_type != "" &
              is.finite(launch_speed) & is.finite(launch_angle) &
              is.finite(launch_speed_angle)]
sw[, hard   := isbip & launch_speed >= 95]
sw[, barrel := isbip & launch_speed_angle == 6]

cat(sprintf("Swings on the timing axis: %s | of which balls in play: %s\n",
  format(nrow(sw), big.mark=","), format(sum(sw$isbip), big.mark=",")))
cat(sprintf("League mean |tdev| = %.2f in | BIP hard-hit %.1f%% | barrel %.1f%% | conversion %.1f%%\n\n",
  mean(sw$adev), 100*sum(sw$hard)/sum(sw$isbip), 100*sum(sw$barrel)/sum(sw$isbip),
  100*sum(sw$barrel)/sum(sw$hard)))

# =============================================================================
# 1. Which channel does mistiming actually attack?
# =============================================================================
cat("############ 1. The two channels along the timing axis ############\n")
cat("Balls in play only, binned by displacement from the hitter's own contact-depth norm.\n")
cat("Negative = deep/late contact, positive = out front/early.\n\n")
bp <- sw[isbip == TRUE]
bp[, sbin := cut(tdev, c(-Inf,-15,-12,-9,-6,-3,0,3,6,9,12,15,Inf))]
S <- bp[, .(bip=.N, ev=mean(launch_speed), la=mean(launch_angle),
            hard=100*mean(hard), barrel=100*mean(barrel),
            conv=100*sum(barrel)/pmax(sum(hard),1),
            xw=mean(estimated_woba_using_speedangle, na.rm=TRUE)), by=sbin][order(sbin)]
print(S[, .(sbin, bip, ev=round(ev,1), la=round(la,1), hard=round(hard,1),
            conv=round(conv,1), barrel=round(barrel,1), xwobacon=round(xw,3))],
      row.names=FALSE)

cat("\n  Magnitude only -- what does being knocked off the norm cost each channel?\n")
bp[, mbin := cut(adev, c(-Inf,3,6,9,12,15,Inf),
                 labels=c("0-3","3-6","6-9","9-12","12-15","15+"))]
M <- bp[, .(bip=.N, hard=100*mean(hard), conv=100*sum(barrel)/pmax(sum(hard),1),
            barrel=100*mean(barrel)), by=mbin][order(mbin)]
base <- M[mbin=="0-3"]
M[, `:=`(hard_rel = round(100*hard/base$hard - 100, 1),
         conv_rel = round(100*conv/base$conv - 100, 1),
         brl_rel  = round(100*barrel/base$barrel - 100, 1))]
print(M[, .(displacement_in = mbin, bip, hard=round(hard,1), conv=round(conv,1),
            barrel=round(barrel,1), hard_pct_change=hard_rel, conv_pct_change=conv_rel,
            barrel_pct_change=brl_rel)], row.names=FALSE)
cat("\n  (pct_change is relative to the 0-3 inch reference row, i.e. well-timed swings)\n")

# =============================================================================
# 2. Pitcher level: does the skill map onto the two factors?
# =============================================================================
cat("\n############ 2. Pitcher-level association ############\n")
agg <- function(d, minsw, minbip) {
  a <- d[, .(
    swings   = .N,
    bip      = sum(isbip),
    disrupt  = mean(excess),
    tdev     = mean(tdev),
    spread   = sd(tdev),
    late12   = 100*mean(tdev <= -12),
    early12  = 100*mean(tdev >=  12),
    ext      = 100*mean(adev >= 12),
    # contact outcomes, computed on the balls in play only
    hardhit  = 100*sum(hard)/pmax(sum(isbip),1),
    conv     = 100*sum(barrel)/pmax(sum(hard),1),
    barrelpc = 100*sum(barrel)/pmax(sum(isbip),1),
    ev       = mean(launch_speed[isbip], na.rm=TRUE),
    la       = mean(launch_angle[isbip], na.rm=TRUE),
    gb       = 100*mean(bb_type[isbip]=="ground_ball"),
    xw       = mean(estimated_woba_using_speedangle[isbip], na.rm=TRUE),
    rv       = mean(pit_rv)
  ), by=pitcher]
  a[swings >= minsw & bip >= minbip]
}
P <- agg(sw[game_year==2026], 250, 100)
cat(sprintf("2026 pitchers with >=250 swings and >=100 balls in play: %d\n\n", nrow(P)))

TARG <- c("hardhit","conv","barrelpc","ev","la","gb","xw")
TL <- c(hardhit="Hard-hit % (EV channel)", conv="Conversion % (angle channel)",
        barrelpc="Barrel %", ev="Mean exit velocity", la="Mean launch angle",
        gb="Ground-ball %", xw="xwOBAcon")
PRED <- c("disrupt","ext","spread","tdev","late12","early12")
PLB <- c(disrupt="Disruption (excess in)", ext="% pushed >=12 in", spread="SD of tdev",
         tdev="Mean tdev (signed)", late12="% late >=12 in", early12="% early >=12 in")
tab <- rbindlist(lapply(TARG, function(y) {
  r <- sapply(PRED, function(x) cor(P[[x]], P[[y]], use="complete.obs"))
  as.data.table(c(list(outcome = TL[[y]]), setNames(as.list(round(r,3)), PLB[PRED])))
}))
print(tab, row.names=FALSE)

cat("\n  Same table, but with hard-hit%% partialled out of the angle channel:\n")
P[, conv_res := residuals(lm(conv ~ hardhit, data=P))]
for (x in PRED) cat(sprintf("     %-24s vs conversion net of hard-hit: r = %+.3f\n",
  PLB[[x]], cor(P[[x]], P$conv_res)))

# =============================================================================
# 3. Test A -- orthogonal within season, zero shared pitches
# =============================================================================
cat("\n############ 3. Test A: disruption on non-contact swings -> contact quality ############\n")
cat("Disruption is measured only on whiffs and fouls; the outcomes are measured only on the\n")
cat("balls in play. The two sides share no pitch, so any link is not mechanical.\n\n")
D <- sw[game_year==2026 & isbip==FALSE, .(nsw_off = .N, disrupt_off = mean(excess),
        tdev_off = mean(tdev), ext_off = 100*mean(adev >= 12)), by=pitcher][nsw_off >= 150]
A <- merge(P, D, by="pitcher")
cat(sprintf("  n = %d pitchers\n", nrow(A)))
for (y in TARG) {
  ct <- cor.test(A$disrupt_off, A[[y]])
  cat(sprintf("     disruption (non-contact swings) -> %-28s r = %+.3f (p = %.3g)\n",
              TL[[y]], ct$estimate, ct$p.value))
}
cat("\n")
for (y in c("conv","barrelpc")) {
  ct <- cor.test(A$tdev_off, A[[y]])
  cat(sprintf("     mean tdev  (non-contact swings) -> %-28s r = %+.3f (p = %.3g)\n",
              TL[[y]], ct$estimate, ct$p.value))
}

# =============================================================================
# 4. Tests B and C -- does it forecast better than the rate itself?
# =============================================================================
cat("\n############ 4. Tests B and C: forecasting ############\n")
sb <- function(r) 2*r/(1+r)
sw[, half := sample(rep_len(1:2, .N)), by=pitcher]
H1 <- agg(sw[game_year==2026 & half==1], 120, 50)
H2 <- agg(sw[game_year==2026 & half==2], 120, 50)
B  <- merge(H1, H2, by="pitcher", suffixes=c("_a","_b"))
Y1 <- agg(sw[game_year==2025], 250, 100)
Y2 <- agg(sw[game_year==2026], 250, 100)
C  <- merge(Y1, Y2, by="pitcher", suffixes=c("_a","_b"))

cat("\n  Reliability of each side (split-half within 2026, Spearman-Brown corrected):\n")
for (v in c("disrupt","ext","tdev","la","gb","hardhit","conv","barrelpc")) {
  r <- cor(B[[paste0(v,"_a")]], B[[paste0(v,"_b")]], use="complete.obs")
  cat(sprintf("     %-9s r = %+.3f -> full season %+.3f\n", v, r, sb(r)))
}

fore <- function(M, lab) {
  cat(sprintf("\n  -- %s (n = %d) --\n", lab, nrow(M)))
  cat("     single predictors, correlation with the later barrel %:\n")
  for (v in c("barrelpc","conv","hardhit","disrupt","ext","la","gb","tdev","spread")) {
    cat(sprintf("        %-9s r = %+.3f\n", v, cor(M[[paste0(v,"_a")]], M$barrelpc_b,
                                                  use="complete.obs")))
  }
  r2 <- function(f) summary(lm(as.formula(f), data=M))$r.squared
  cat("     R2 predicting the later barrel %:\n")
  cat(sprintf("        barrel%% alone                        %.3f\n", r2("barrelpc_b ~ barrelpc_a")))
  cat(sprintf("        hard-hit%% alone                      %.3f\n", r2("barrelpc_b ~ hardhit_a")))
  cat(sprintf("        disruption alone                     %.3f\n", r2("barrelpc_b ~ disrupt_a")))
  cat(sprintf("        hard-hit%% + disruption               %.3f\n", r2("barrelpc_b ~ hardhit_a + disrupt_a")))
  cat(sprintf("        hard-hit%% + launch angle + GB%%       %.3f\n", r2("barrelpc_b ~ hardhit_a + la_a + gb_a")))
  cat(sprintf("        the above + disruption               %.3f\n", r2("barrelpc_b ~ hardhit_a + la_a + gb_a + disrupt_a")))
  cat("     and predicting the later CONVERSION (the angle factor):\n")
  cat(sprintf("        conversion alone                     %.3f\n", r2("conv_b ~ conv_a")))
  cat(sprintf("        disruption alone                     %.3f\n", r2("conv_b ~ disrupt_a")))
  cat(sprintf("        launch angle + GB%%                   %.3f\n", r2("conv_b ~ la_a + gb_a")))
  cat(sprintf("        launch angle + GB%% + disruption      %.3f\n", r2("conv_b ~ la_a + gb_a + disrupt_a")))
}
fore(B, "Test B: split-half within 2026")
fore(C, "Test C: 2025 -> 2026")

# =============================================================================
# 5. Where the timing profile disagrees with the barrel rate
# =============================================================================
cat("\n############ 5. Timing-implied barrel rate vs actual, 2026 ############\n")
cat("Expected barrel%% is fitted from the timing profile plus hard-hit%% only -- it never sees\n")
cat("the pitcher's own barrel or conversion rate. Residual = actual minus expected.\n\n")
mod <- lm(barrelpc ~ hardhit + disrupt + tdev + ext, data=P)
cat(sprintf("  Model R2 (in sample, 2026) = %.3f\n", summary(mod)$r.squared))
print(round(summary(mod)$coefficients, 4))
P[, exp_barrel := fitted(mod)]
P[, resid_barrel := barrelpc - exp_barrel]
NAME <- unique(sw[, .(pitcher, player_name)], by="pitcher")
P <- merge(P, NAME, by="pitcher")
shw <- function(d, n=12) d[1:min(n,nrow(d)), .(player_name, swings, bip,
  disrupt=round(disrupt,2), tdev=round(tdev,1), ext=round(ext,1),
  hardhit=round(hardhit,1), conv=round(conv,1), barrel=round(barrelpc,1),
  expected=round(exp_barrel,1), resid=round(resid_barrel,1),
  la=round(la,1), gb=round(gb,1), xwobacon=round(xw,3))]
cat("\n  -- barrel rate WORSE than the timing profile implies (regression candidates) --\n")
print(shw(P[order(-resid_barrel)]), row.names=FALSE)
cat("\n  -- barrel rate BETTER than the timing profile implies (may not hold) --\n")
print(shw(P[order(resid_barrel)]), row.names=FALSE)

cat("\n  Does the residual carry to the next season, or does it wash out?\n")
CY <- merge(P[, .(pitcher, resid_barrel, barrelpc, disrupt)],
            Y2[, .(pitcher)], by="pitcher")
R25 <- agg(sw[game_year==2025], 250, 100)
R25[, exp25 := predict(mod, newdata=R25)]
R25[, res25 := barrelpc - exp25]
CC <- merge(R25[, .(pitcher, res25, barrel25=barrelpc, disrupt25=disrupt)],
            Y2[, .(pitcher, barrel26=barrelpc, conv26=conv)], by="pitcher")
cat(sprintf("     n = %d pitchers in both seasons\n", nrow(CC)))
cat(sprintf("     2025 residual        -> 2026 barrel%%  r = %+.3f\n", cor(CC$res25, CC$barrel26)))
cat(sprintf("     2025 barrel%%         -> 2026 barrel%%  r = %+.3f\n", cor(CC$barrel25, CC$barrel26)))
cat(sprintf("     2025 disruption      -> 2026 barrel%%  r = %+.3f\n", cor(CC$disrupt25, CC$barrel26)))

# =============================================================================
# 6. Bridge: do the timing levers line up with the four contact-manager types?
# =============================================================================
cat("\n############ 6. Timing levers vs the four contact-manager types ############\n")
cat(sprintf("  The two levers are close to independent: cor(disruption, mean tdev) = %+.3f\n",
            cor(P$disrupt, P$tdev)))
cat("  So a pitcher can push hitters off their norm, aim them late, both, or neither.\n\n")

hmed <- median(P$hardhit); cmed <- median(P$conv)
P[, ctype := fifelse(hardhit <= hmed & conv <= cmed, "Weak contact",
             fifelse(hardhit >  hmed & conv <= cmed, "Angle manager",
             fifelse(hardhit <= hmed & conv >  cmed, "Velocity suppressor",
                     "Vulnerable")))]
print(P[, .(pitchers=.N, disrupt=round(mean(disrupt),2), tdev=round(mean(tdev),2),
            late12=round(mean(late12),1), early12=round(mean(early12),1),
            ext=round(mean(ext),1), hardhit=round(mean(hardhit),1),
            conv=round(mean(conv),1), barrel=round(mean(barrelpc),1),
            la=round(mean(la),1), gb=round(mean(gb),1)),
        by=ctype][order(barrel)], row.names=FALSE)

cat("\n  Route taken, cross-tabbed against contact-manager type:\n")
P[, route := fifelse(late12 >= 10 & early12 >= 10, "both tails",
            fifelse(tdev <= -1, "pushes late",
            fifelse(tdev >=  1, "pushes early", "neutral")))]
print(table(P$ctype, P$route))
cat("\n  Mean tdev by type is the discriminator that matters:\n")
for (ty in c("Weak contact","Angle manager","Velocity suppressor","Vulnerable")) {
  d <- P[ctype==ty]
  cat(sprintf("     %-20s tdev %+.2f in | disruption %+.2f in | LA %.1f deg | GB %.1f%%\n",
              ty, mean(d$tdev), mean(d$disrupt), mean(d$la), mean(d$gb)))
}

fwrite(P[order(resid_barrel)], file.path("data","statcast_2026","timing_vs_barrel_2026.csv"))
cat("\nWrote data/statcast_2026/timing_vs_barrel_2026.csv\n")
