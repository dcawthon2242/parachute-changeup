#!/usr/bin/env Rscript

# How much of a pitcher's TRUE contact-quality skill does his repeatable timing ability
# explain?
#
# Every correlation reported so far is between two noisily measured quantities, and barrel
# rate and conversion are measured very badly (full-season reliability 0.27 and 0.24). That
# drags observed correlations toward zero and makes timing look weaker than it is. The fix
# is to estimate each metric's reliability by split-half, then disattenuate:
#
#     r_true = r_obs / sqrt(rel_x * rel_y)
#
# The question this answers: of the contact-quality skill that actually persists, what share
# is traceable to timing? That is different from "what share of observed variance", which is
# mostly noise.

suppressPackageStartupMessages({ library(data.table); library(splines) })
set.seed(23)

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
sw[, isbip := description == "hit_into_play" & bb_type != "" & is.finite(launch_speed) &
              is.finite(launch_angle) & is.finite(launch_speed_angle)]
sw[, hard := isbip & launch_speed >= 95]
sw[, brl  := isbip & launch_speed_angle == 6]

# league rates within timing bin, for the timing-mix component
sw[, tb := cut(tdev, c(-Inf,-15,-12,-9,-6,-3,0,3,6,9,12,15,Inf))]
LG <- sw[isbip == TRUE, .(lg_hard = mean(hard), lg_brl = mean(brl)), by=tb]
sw <- merge(sw, LG, by="tb", all.x=TRUE)

agg <- function(d) d[, .(
  swings   = .N,
  bip      = sum(isbip),
  # repeatable timing abilities, measured on ALL swings
  disrupt  = mean(excess),
  mtdev    = mean(tdev),
  ext      = 100*mean(adev >= 12),
  # timing-mix implied contact quality
  t_hard   = 100*mean(lg_hard[isbip], na.rm=TRUE),
  t_brl    = 100*mean(lg_brl[isbip],  na.rm=TRUE),
  # actual contact quality
  hardhit  = 100*sum(hard)/pmax(sum(isbip),1),
  brlpc    = 100*sum(brl)/pmax(sum(isbip),1),
  conv     = 100*sum(brl)/pmax(sum(hard),1),
  la       = mean(launch_angle[isbip], na.rm=TRUE),
  gb       = 100*mean(bb_type[isbip]=="ground_ball"),
  xw       = mean(estimated_woba_using_speedangle[isbip], na.rm=TRUE)
), by=pitcher]

M <- c("disrupt","mtdev","ext","t_hard","t_brl","hardhit","brlpc","conv","la","gb","xw")
LAB <- c(disrupt="Disruption (magnitude)", mtdev="Mean tdev (direction)",
         ext="% pushed >=12 in", t_hard="Timing-implied hard-hit%",
         t_brl="Timing-implied barrel%", hardhit="Hard-hit% allowed",
         brlpc="Barrel% allowed", conv="Conversion%", la="Mean launch angle",
         gb="Ground-ball%", xw="xwOBAcon allowed")

# ---- reliability by split-half, projected to full season -------------------
sb <- function(r) 2*r/(1+r)
sw[, half := sample(rep_len(1:2, .N)), by=pitcher]
H1 <- agg(sw[game_year==2026 & half==1])[swings>=150 & bip>=50]
H2 <- agg(sw[game_year==2026 & half==2])[swings>=150 & bip>=50]
HH <- merge(H1, H2, by="pitcher", suffixes=c("_a","_b"))
REL <- sapply(M, function(v) {
  r <- cor(HH[[paste0(v,"_a")]], HH[[paste0(v,"_b")]], use="complete.obs"); sb(r)
})
cat("############ Reliability, full-season equivalent (split-half, Spearman-Brown) ############\n")
cat(sprintf("  n = %d pitchers\n\n", nrow(HH)))
for (v in M) cat(sprintf("  %-26s %.3f\n", LAB[[v]], REL[[v]]))

# ---- observed vs disattenuated --------------------------------------------
P <- agg(sw[game_year==2026])[swings>=250 & bip>=120]
cat(sprintf("\n############ Observed vs disattenuated correlations, 2026 (n = %d) ############\n",
            nrow(P)))
cat("Disattenuated = observed / sqrt(rel_x * rel_y). This is the correlation between the two\n")
cat("underlying skills, with sampling noise in both removed.\n\n")
pairs <- list(
  c("disrupt","hardhit"), c("disrupt","brlpc"), c("disrupt","conv"), c("disrupt","xw"),
  c("mtdev","la"), c("mtdev","conv"), c("mtdev","brlpc"), c("mtdev","gb"),
  c("t_hard","hardhit"), c("t_brl","brlpc")
)
out <- rbindlist(lapply(pairs, function(p) {
  ro <- cor(P[[p[1]]], P[[p[2]]], use="complete.obs")
  rd <- ro/sqrt(REL[[p[1]]]*REL[[p[2]]])
  data.table(predictor = LAB[[p[1]]], outcome = LAB[[p[2]]],
             observed_r = round(ro,3), disattenuated_r = round(max(-1,min(1,rd)),3),
             true_R2 = round(min(1, rd^2), 3))
}))
print(out, row.names=FALSE)

# ---- share of TRUE contact-quality variance explained by the timing mix ----
cat("\n############ Share of TRUE skill variance traceable to the timing mix ############\n")
cat("Observed-variance shares understate this badly, because most observed variance in\n")
cat("barrel% and hard-hit% is sampling noise rather than skill.\n\n")
for (p in list(c("t_hard","hardhit"), c("t_brl","brlpc"))) {
  ro <- cor(P[[p[1]]], P[[p[2]]], use="complete.obs")
  obs_share <- var(P[[p[1]]])/var(P[[p[2]]])
  rd <- ro/sqrt(REL[[p[1]]]*REL[[p[2]]])
  cat(sprintf("  %s -> %s\n", LAB[[p[1]]], LAB[[p[2]]]))
  cat(sprintf("     share of OBSERVED variance         %.0f%%\n", 100*obs_share))
  cat(sprintf("     reliability of the timing mix      %.3f\n", REL[[p[1]]]))
  cat(sprintf("     reliability of the outcome         %.3f\n", REL[[p[2]]]))
  cat(sprintf("     disattenuated r                    %+.3f\n", rd))
  cat(sprintf("     share of TRUE skill variance       %.0f%%\n\n", 100*min(1, rd^2)))
}

# ---- the practical question: forecasting value at equal sample -------------
cat("############ Forecast value: 2025 -> 2026 barrel%, and why timing competes ############\n")
Y1 <- agg(sw[game_year==2025])[swings>=250 & bip>=120]
Y2 <- agg(sw[game_year==2026])[swings>=250 & bip>=120]
C  <- merge(Y1, Y2, by="pitcher", suffixes=c("_a","_b"))
cat(sprintf("  n = %d pitchers in both seasons\n\n", nrow(C)))
cat("  Each 2025 metric's own reliability, its raw correlation with 2026 barrel%,\n")
cat("  and the ratio -- how much predictive punch it delivers per unit of measurement quality.\n\n")
fr <- rbindlist(lapply(c("brlpc","conv","hardhit","la","gb","disrupt","mtdev","ext","t_brl"),
  function(v) {
    r <- cor(C[[paste0(v,"_a")]], C$brlpc_b, use="complete.obs")
    data.table(metric = LAB[[v]], reliability = round(REL[[v]],3),
               r_with_2026_barrel = round(r,3),
               per_reliability = round(abs(r)/sqrt(REL[[v]]),3))
  }))
print(fr[order(-abs(r_with_2026_barrel))], row.names=FALSE)

r2 <- function(f) summary(lm(as.formula(f), data=C))$r.squared
cat("\n  R2 predicting 2026 barrel%:\n")
cat(sprintf("     2025 barrel%% alone                              %.3f\n", r2("brlpc_b ~ brlpc_a")))
cat(sprintf("     2025 timing only (disrupt + mtdev + ext)        %.3f\n",
            r2("brlpc_b ~ disrupt_a + mtdev_a + ext_a")))
cat(sprintf("     2025 timing + hard-hit%%                         %.3f\n",
            r2("brlpc_b ~ disrupt_a + mtdev_a + ext_a + hardhit_a")))
cat(sprintf("     2025 launch angle + GB%% + hard-hit%%             %.3f\n",
            r2("brlpc_b ~ la_a + gb_a + hardhit_a")))
cat(sprintf("     everything                                     %.3f\n",
            r2("brlpc_b ~ la_a + gb_a + hardhit_a + disrupt_a + mtdev_a + ext_a")))
cat(sprintf("\n  Ceiling: with 2026 barrel%% reliability at %.3f, the maximum achievable R2\n",
            REL[["brlpc"]]))
cat(sprintf("  for ANY predictor of the observed 2026 rate is about %.3f.\n", REL[["brlpc"]]))
cat(sprintf("  Timing alone reaches %.0f%% of that ceiling; the full model reaches %.0f%%.\n",
  100*r2("brlpc_b ~ disrupt_a + mtdev_a + ext_a")/REL[["brlpc"]],
  100*r2("brlpc_b ~ la_a + gb_a + hardhit_a + disrupt_a + mtdev_a + ext_a")/REL[["brlpc"]]))
