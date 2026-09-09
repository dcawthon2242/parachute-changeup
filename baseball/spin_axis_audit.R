#!/usr/bin/env Rscript

# IS "SPIN SIMILARITY" MEASURING SPIN AT ALL?
#
# Three challenges, in order of how badly each would invalidate the cue:
#
#   1. Is Statcast's spin_axis an independent measurement, or is it just the movement
#      vector rewritten as an angle? If the latter, "spin similarity" is not a spin cue -
#      it is a movement-direction cue wearing a spin label, and nothing a hitter does with
#      seeing the ball rotate is involved.
#   2. If it is movement direction, does axis similarity simply track how little the two
#      pitches separate - i.e. is it a tunneling proxy?
#   3. Does the +0.42 on sliders survive once movement separation and tunnel tightness are
#      controlled for, or does it vanish?

suppressPackageStartupMessages({ library(data.table) })
options(width = 200)
MDIR <- "data/statcast_model"

COLS <- c("pitch_type","player_name","pitcher","p_throws","release_speed","pfx_x","pfx_z",
          "ax","az","spin_axis","release_spin_rate","game_year")
d <- fread("data/statcast_2025/statcast_2025_all.csv", select = COLS, showProgress = FALSE)
d <- d[is.finite(spin_axis) & is.finite(pfx_x) & is.finite(ax)]

## ---- 1. is spin_axis just the movement angle? ---------------------------------
# Statcast reports the axis of the deflection. Pure backspin = movement straight up,
# which is reported as 180 deg, so the mapping to test is atan2(up, side) + 90.
ang <- function(x, z) (atan2(z, x)*180/pi + 90) %% 360
d[, `:=`(pred_pfx = ang(pfx_x, pfx_z),
         pred_acc = ang(ax, az + 32.174))]          # az carries gravity; Magnus is what's left
dev <- function(a, b) { e <- abs(a - b) %% 360; pmin(e, 360 - e) }
d[, `:=`(err_pfx = dev(spin_axis, pred_pfx), err_acc = dev(spin_axis, pred_acc))]

cat("=== 1. CAN spin_axis BE RECONSTRUCTED FROM MOVEMENT ALONE? ===\n")
cat("   deviation of the reported axis from the angle implied by the movement vector\n\n")
print(d[, .(n = .N,
            median_err_from_pfx = round(median(err_pfx),2),
            p95_err_from_pfx    = round(quantile(err_pfx,.95),2),
            median_err_from_acc = round(median(err_acc),2),
            p95_err_from_acc    = round(quantile(err_acc,.95),2)),
        by = pitch_type][order(median_err_from_acc)][pitch_type %in%
          c("FF","SI","FC","SL","ST","CU","KC","CH","FS")], row.names = FALSE)
cat("\n   ALL PITCHES: median", round(median(d$err_acc),2), "deg,",
    round(100*mean(d$err_acc < 2),1), "% within 2 deg of the movement-implied angle\n")

## ---- 2. does axis similarity just mean "these two barely separate"? -----------
d[, `:=`(hb = fifelse(p_throws=="R",-1,1)*pfx_x*12, ivb = pfx_z*12,
         ax_m = fifelse(p_throws=="R", spin_axis, (360-spin_axis)%%360))]
agg <- d[, .(n=.N, velo=mean(release_speed), hb=mean(hb), ivb=mean(ivb),
             sx=mean(sin(ax_m*pi/180)), cx=mean(cos(ax_m*pi/180))),
         by=.(pitcher, pitch_type)][n >= 200]
agg[, axis := (atan2(sx,cx)*180/pi) %% 360]
fb <- agg[pitch_type %in% c("FF","SI")][order(pitcher,-n)][, .SD[1], by=pitcher]
sl <- agg[pitch_type == "SL"]
p <- merge(sl, fb[, .(pitcher, fb_axis=axis, fb_hb=hb, fb_ivb=ivb, fb_velo=velo)], by="pitcher")
p[, gap := abs(axis-fb_axis)][gap > 180, gap := 360-gap]
p[, `:=`(mv_sep = sqrt((hb-fb_hb)^2 + (ivb-fb_ivb)^2),   # inches of total separation
         sl_mv  = sqrt(hb^2 + ivb^2),
         velo_gap = fb_velo - velo)]

sp <- function(x,y){ ct <- suppressWarnings(cor.test(x,y,method="spearman",exact=FALSE))
                     sprintf("%+.3f (p=%.2g)", ct$estimate, ct$p.value) }
cat("\n=== 2. WHAT IS AXIS SIMILARITY ACTUALLY TRACKING? (sliders, n =", nrow(p), ") ===\n")
cat("   axis similarity vs movement separation from the FB :", sp(-p$gap, p$mv_sep), "\n")
cat("   axis similarity vs the slider's own total movement :", sp(-p$gap, p$sl_mv), "\n")
cat("   axis similarity vs velocity gap from the FB        :", sp(-p$gap, p$velo_gap), "\n")

## ---- 3. does the effect survive controlling for separation and tunnel? --------
oof <- readRDS(file.path(MDIR,"oof_whiff_resid.rds")); setDT(oof)
oof[, res_loc := NA_real_]
for (g in c("breaking","offspeed")) { i <- which(oof$grp==g)
  oof$res_loc[i] <- residuals(lm(wres ~ poly(plate_x,3)*poly(plate_z,3)+below_zone+VAA+HAA,
                                 data=oof[i])) }
md <- readRDS(file.path(MDIR,"miss_grade_data_activespin.rds"))
k <- c("game_pk","at_bat_number","pitch_number","season"); setkeyv(md,k); setkeyv(oof,k)
b <- md[, .(game_pk,at_bat_number,pitch_number,season,axis_diff,path_ratio,
            ax_diff,az_diff,speed_diff)][oof]
b <- b[pitch_type=="SL" & is.finite(res_loc) & is.finite(axis_diff)]
A <- b[, .(n=.N, axis_sim=-mean(axis_diff), path_ratio=mean(path_ratio,na.rm=TRUE),
           mv_sep=mean(sqrt(ax_diff^2+az_diff^2),na.rm=TRUE),
           velo_gap=-mean(speed_diff,na.rm=TRUE), resid=100*mean(res_loc)),
       by=pitcher][n>=200]
A <- A[is.finite(path_ratio) & is.finite(mv_sep)]

pcor <- function(y, x, z) {   # spearman partial: rank everything, regress out z
  R <- as.data.table(lapply(c(list(y=y,x=x), z), rank))
  ry <- residuals(lm(y ~ ., data=R[, !"x"])); rx <- residuals(lm(x ~ ., data=R[, !"y"]))
  ct <- cor.test(rx, ry); sprintf("%+.3f (p=%.2g)", ct$estimate, ct$p.value)
}
cat("\n=== 3. DOES IT SURVIVE CONTROLS? (sliders, n =", nrow(A), " pitchers) ===\n")
cat("   axis similarity vs location-adjusted overperformance\n")
cat("     raw                                    :", sp(A$axis_sim, A$resid), "\n")
cat("     partialling out movement separation    :", pcor(A$resid, A$axis_sim, list(m=A$mv_sep)), "\n")
cat("     partialling out path_ratio (tunnel)    :", pcor(A$resid, A$axis_sim, list(t=A$path_ratio)), "\n")
cat("     partialling out velo gap               :", pcor(A$resid, A$axis_sim, list(v=A$velo_gap)), "\n")
cat("     partialling out all three              :",
    pcor(A$resid, A$axis_sim, list(m=A$mv_sep, t=A$path_ratio, v=A$velo_gap)), "\n")
cat("\n   for reference, the competing cues on their own:\n")
cat("     movement separation vs overperformance :", sp(-A$mv_sep, A$resid), "\n")
cat("     tunnel path_ratio vs overperformance   :", sp(-A$path_ratio, A$resid), "\n")
cat("     axis similarity vs movement separation :", sp(A$axis_sim, -A$mv_sep), "\n")
cat("     axis similarity vs path_ratio          :", sp(A$axis_sim, -A$path_ratio), "\n")
