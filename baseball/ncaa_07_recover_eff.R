#!/usr/bin/env Rscript

# CAN THE MEASURED-VERSUS-INFERRED AXIS DISCREPANCY RECOVER SPIN EFFICIENCY?
# AND DOES RESTRICTING TO HIGH-EFFICIENCY CHANGEUPS MAKE THE TRACKMAN AXIS USABLE?
#
# The premise is physically sound and worth testing properly: as active spin approaches 100
# percent, the ball's rotation and its break must line up, so measured and movement-derived axes
# should converge. If they do, then screening on high inferred efficiency also screens for the
# regime where TrackMan's axis is trustworthy, and the college analysis is back on.
#
# The earlier crude test used atan2(ax, az + g), which is wrong in a way that matters here: the
# total acceleration contains drag as well as gravity and Magnus, and a changeup and a
# four-seamer shed speed differently, so drag contaminates any comparison between them. This
# script does the decomposition properly. Magnus force is by definition perpendicular to
# velocity and drag is antiparallel to it, so subtracting gravity and then splitting the
# remainder into components parallel and perpendicular to the velocity vector isolates Magnus
# exactly. Statcast gives vx0/vy0/vz0 alongside ax/ay/az, so this is available per pitch.
#
# Four questions, in order:
#   1. does the axis discrepancy track measured active spin at all?
#   2. does the discrepancy shrink at high efficiency, as the premise predicts?
#   3. at high efficiency, does the movement-based CH-to-FF gap agree with the measured one?
#   4. does the MLB parachute result survive being rebuilt with only NCAA-available inputs?
#
# Question 4 is the decisive one. It runs the exact recipe college data would force - inferred
# efficiency, movement-based axis - on MLB, where the answer is known.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
set.seed(9); options(width = 205)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")

M <- readRDS(file.path(MDIR, "parachute_rv.rds"))[pitch_type %in% c("CH","FF")]
M <- M[is.finite(vx0) & is.finite(vy0) & is.finite(vz0) & is.finite(ax) & is.finite(az) &
       is.finite(sax) & is.finite(release_spin_rate) & release_spin_rate > 500]
M[, spin_axis := (atan2(sax, cax)*180/pi) %% 360]

# Gravity out, then split what remains into the drag component (along velocity) and the Magnus
# component (across it).
M[, `:=`(wx = ax, wy = ay, wz = az + 32.174)]
M[, vmag := sqrt(vx0^2 + vy0^2 + vz0^2)]
M[, dotp := (wx*vx0 + wy*vy0 + wz*vz0) / vmag^2]
M[, `:=`(mx = wx - dotp*vx0, my = wy - dotp*vy0, mz = wz - dotp*vz0)]
M[, magnus := sqrt(mx^2 + my^2 + mz^2)]
M[, mv_axis := (atan2(mx, mz)*180/pi) %% 360]

# Ratio of realised Magnus to what the measured total spin could produce. This is the magnitude
# route to efficiency; the axis route is what is being tested against it.
M[, ratio := magnus / (release_spin_rate * release_speed / 1000)]

P <- M[, .(n = .N, sa_s = mean(sin(spin_axis*pi/180)), sa_c = mean(cos(spin_axis*pi/180)),
           mv_s = mean(sin(mv_axis*pi/180)), mv_c = mean(cos(mv_axis*pi/180)),
           ratio = mean(ratio), magnus = mean(magnus), spin = mean(release_spin_rate),
           velo = mean(release_speed)), by = .(pitcher, season, pitch_type)][n >= 40]
P[, `:=`(sa = (atan2(sa_s, sa_c)*180/pi) %% 360, mv = (atan2(mv_s, mv_c)*180/pi) %% 360)]
AS <- readRDS(file.path(MDIR, "active_spin_long.rds"))
P <- merge(P, AS[, .(pitcher, season, pitch_type, eff = active_spin)],
           by = c("pitcher","season","pitch_type"))[is.finite(eff)]

# The two axes use different conventions, so the best rigid alignment is removed before any
# discrepancy is measured. Both a rotation and a possible reflection have to be tried: script 01
# found the reflected orientation fits, and removing only a rotation would leave an enormous
# residual that says nothing about the physics.
P[, `:=`(off = 0, flip = 1)]
for (pt in c("CH","FF")) {
  best <- NULL
  for (s in c(1, -1)) {
    d <- P[pitch_type == pt, ((s*sa - mv) %% 360) * pi/180]
    Rb <- Mod(mean(exp(1i*d)))
    if (is.null(best) || Rb > best$R) best <- list(R = Rb, s = s, o = Arg(mean(exp(1i*d)))*180/pi)
  }
  cat(sprintf("  alignment for %s: %s, offset %.1f deg, resultant %.3f\n", pt,
              if (best$s < 0) "reflected" else "direct", best$o, best$R))
  P[pitch_type == pt, `:=`(off = best$o, flip = best$s)]
}
P[, disc := abs(((flip*sa - mv - off + 180) %% 360) - 180)]
cat(sprintf("pitcher-season-pitch-types: %d (CH %d, FF %d)\n\n", nrow(P),
            sum(P$pitch_type=="CH"), sum(P$pitch_type=="FF")))

cat("=== Q1. does the axis discrepancy track measured active spin? ===\n")
for (pt in c("CH","FF")) { D <- P[pitch_type == pt]
  r1 <- cor.test(D$disc, D$eff); r2 <- cor.test(D$ratio, D$eff)
  cat(sprintf("  %s  axis discrepancy vs efficiency: r = %+.3f (p = %.2g)  |  magnitude ratio: r = %+.3f\n",
              pt, r1$estimate, r1$p.value, r2$estimate)) }
cat("  It does, and strongly. The first run of this script said otherwise because it removed only\n")
cat("  a rotation when aligning the two conventions; the fit is reflected, and forcing the direct\n")
cat("  orientation buried the signal under a systematic error larger than the effect.\n")

cat("\n=== Q2. does the discrepancy shrink at high efficiency, as the premise predicts? ===\n")
P[, ebin := cut(eff, c(0,.70,.80,.85,.90,.95,1.01), right = FALSE,
                labels = c("<.70",".70-.80",".80-.85",".85-.90",".90-.95",".95+"))]
print(P[pitch_type == "CH", .(seasons = .N, median_discrepancy = round(median(disc),1),
        p75 = round(quantile(disc,.75),1)), by = ebin][order(ebin)], row.names = FALSE)

cat("\n  where does the relationship hold? (the bands above are not monotone at the bottom)\n")
for (pt in c("CH","FF")) for (rg in list(c(.85,2), c(.75,2), c(0,.75))) {
  D <- P[pitch_type == pt & eff >= rg[1] & eff < rg[2]]
  ct <- cor.test(D$disc, D$eff)
  cat(sprintf("    %s  efficiency in [%.2f,%.2f): n = %4d  axis r = %+.3f (p = %.2g)  break-magnitude r = %+.3f\n",
              pt, rg[1], min(rg[2],1), nrow(D), ct$estimate, ct$p.value, cor(D$ratio, D$eff))) }
cat("  The axis route is a high-efficiency instrument. Below about .75 it carries nothing on\n")
cat("  changeups, which is why Avila's true .24 comes back as .56 - his discrepancy is only 9\n")
cat("  degrees. If the goal is specifically to isolate the extreme low tail, use break magnitude.\n")

cat("\n=== Q3. at high efficiency, do the two gap definitions agree? ===\n")
W <- dcast(P, pitcher + season ~ pitch_type, value.var = c("sa","mv","eff","ratio","n"))
W <- W[is.finite(sa_CH) & is.finite(sa_FF) & n_CH >= 50]
gp <- function(a,b) pmin(abs(a-b), 360-abs(a-b))
W[, `:=`(g_meas = gp(sa_CH, sa_FF), g_move = gp(mv_CH, mv_FF), emin = pmin(eff_CH, eff_FF))]
print(rbindlist(lapply(list(c(0,2), c(.85,2), c(.90,2), c(.95,2)), function(k) {
  D <- W[emin >= k[1]]
  data.table(efficiency_floor = k[1], seasons = nrow(D),
             mean_measured = round(mean(D$g_meas),1), mean_movement = round(mean(D$g_move),1),
             r = round(cor(D$g_meas, D$g_move),3),
             median_abs_diff = round(median(abs(D$g_meas - D$g_move)),1),
             pct_within10_meas = round(100*mean(D$g_meas <= 10),1),
             pct_within10_move = round(100*mean(D$g_move <= 10),1)) })), row.names = FALSE)

## ---- Q4: turn the discrepancy into a usable efficiency estimate ---------------------------------
#
# Two independent signals are now on the table: the angle between the measured axis and the
# Magnus direction, and the size of the break relative to what the measured spin rate could
# produce. Fitted together, out of fold, they say how much of measured active spin is
# recoverable from the pitch record alone. Five-fold cross-validation matters here because the
# alternative - reporting the in-sample fit of a four-term regression - would overstate it.
cat("\n=== Q4. combining the two signals into an efficiency estimate (5-fold out of fold) ===\n")
est <- rbindlist(lapply(c("CH","FF"), function(pt) {
  D <- copy(P[pitch_type == pt]); D[, fold := sample(rep_len(1:5, .N))]
  for (k in 1:5) {
    fit <- lm(eff ~ disc + ratio + spin + velo, D[fold != k])
    D[fold == k, pred := predict(fit, .SD)]
  }
  D[, `:=`(pred_disc = NA_real_, pred_ratio = NA_real_)]
  for (k in 1:5) {
    D[fold == k, pred_disc  := predict(lm(eff ~ disc,  D[fold != k]), .SD)]
    D[fold == k, pred_ratio := predict(lm(eff ~ ratio, D[fold != k]), .SD)]
  }
  D[] }))
print(est[, .(seasons = .N,
              r_axis_only = round(cor(pred_disc, eff),3), r_break_only = round(cor(pred_ratio, eff),3),
              r_combined  = round(cor(pred, eff),3),      R2 = round(cor(pred, eff)^2,3),
              rmse = round(sqrt(mean((pred-eff)^2)),3),
              sd_of_truth = round(sd(eff),3)), by = pitch_type], row.names = FALSE)

# A correlation is not the same as a usable gate. The question that matters for bin-building is
# whether the estimate can pick out the pitches below a .85 efficiency floor.
cat("\n  recovering the 'below .85 active spin' flag, threshold set to match the true base rate:\n")
for (pt in c("CH","FF")) { D <- est[pitch_type == pt]
  thr <- quantile(D$pred, mean(D$eff < .85)); flag <- D$pred < thr; truth <- D$eff < .85
  cat(sprintf("    %s  base rate %.1f%%  |  precision %.1f%%  recall %.1f%%  |  median true efficiency: flagged %.2f vs rest %.2f\n",
              pt, 100*mean(truth), 100*mean(truth[flag]), 100*mean(flag[truth]),
              median(D$eff[flag]), median(D$eff[!flag]))) }

cat("\n  lowest-efficiency changeups the estimate identifies, against the measured truth.\n")
cat("  The ranking is sound but a linear fit compresses the tail: Avila's true .24 comes back as .56.\n")
NM <- unique(readRDS(file.path(MDIR, "parachute_rv.rds"))[, .(pitcher, player_name)])
L <- merge(est[pitch_type == "CH"], NM, by = "pitcher", all.x = TRUE)[order(pred)][1:15]
print(L[, .(Pitcher = trimws(paste(sub(".*,\\s*","",player_name), sub(",.*","",player_name))),
            Season = season, CH = n, Discrepancy = round(disc,1),
            Estimated = round(pred,2), Measured = round(eff,2))], row.names = FALSE)

## ---- Q5: rebuild the MLB result using only NCAA-available inputs -------------------------------
cat("\n=== Q5. does the MLB parachute result survive an NCAA-only recipe? ===\n")
# The sample and the measured columns are rebuilt exactly as the Core bin defines them, so the
# first row below has to reproduce the known result before any substitution can be read.
R  <- readRDS(file.path(MDIR, "whiff_tjstuff.rds"))
Q  <- R[, .(nsw = .N, axis = mean(axis_diff), arm = mean(arm_angle),
            velo_sep = -mean(speed_diff), w = 100*mean(r4)), by = .(pitcher, season)]
Q <- merge(Q, AS[pitch_type == "CH", .(pitcher, season, as_ch = active_spin)], by = c("pitcher","season"))
Q <- merge(Q, AS[pitch_type == "FF", .(pitcher, season, as_fb = active_spin)], by = c("pitcher","season"))
Q <- Q[nsw >= 60]

# Movement-axis gap built the same way the measured one is: each changeup against the pitcher's
# mean four-seam axis, then averaged - not the gap between two means, which is a smaller number.
FB <- M[pitch_type == "FF", .(n = .N, s = mean(sin(mv_axis*pi/180)), c = mean(cos(mv_axis*pi/180))),
        by = .(pitcher, season)][n >= 40]
FB[, fbmv := (atan2(s, c)*180/pi) %% 360]
CH <- merge(M[pitch_type == "CH"], FB[, .(pitcher, season, fbmv)], by = c("pitcher","season"))
CH[, gi := abs(((mv_axis - fbmv + 180) %% 360) - 180)]
Q <- merge(Q, CH[, .(g_move = mean(gi), ratio_CH = mean(ratio)), by = .(pitcher, season)],
           by = c("pitcher","season"))
Q <- merge(Q, M[pitch_type == "FF", .(n = .N, ratio_FF = mean(ratio)),
                by = .(pitcher, season)][n >= 40, -"n"], by = c("pitcher","season"))
# Inferred efficiency is a rank within pitch type, since the magnitude ratio has no natural
# scale; the measured-efficiency gap is converted to the same rank space for a like-for-like cut.
Q[, `:=`(inf_CH = frank(ratio_CH)/.N, inf_FF = frank(ratio_FF)/.N,
         r_ch = frank(as_ch)/.N, r_fb = frank(as_fb)/.N)]
Q[, `:=`(meas_gap = abs(as_ch - as_fb), inf_gap = abs(inf_CH - inf_FF))]
# Thresholds for the substituted variables are set to admit the same share of seasons as the
# measured cuts they replace, so bin size is held fixed and only the input quality varies.
tm <- quantile(Q$g_move, mean(Q$axis <= 10)); ti <- quantile(Q$inf_gap, mean(Q$meas_gap <= .10))
cat(sprintf("  %d seasons; movement-gap cut %.1f deg matches axis<=10 (%.1f%%); inferred-gap cut %.2f matches |as_gap|<=.10 (%.1f%%)\n",
            nrow(Q), tm, 100*mean(Q$axis <= 10), ti, 100*mean(Q$meas_gap <= .10)))
cat("  Caveat: requiring 40+ four-seamers in the trajectory file cuts the sample from 1,086 seasons\n")
cat("  to 653, so the first row is not the published Core bin and none of this is powered. Read the\n")
cat("  overlap column, not the p-values: substituting the axis keeps only 3 of the 12 true members.\n")

recipes <- list(
  list(l = "MLB truth: measured axis + measured efficiency", i = quote(axis <= 10 & meas_gap <= .10 & arm >= 44)),
  list(l = "measured axis, INFERRED efficiency",             i = quote(axis <= 10 & inf_gap <= ti   & arm >= 44)),
  list(l = "MOVEMENT axis, measured efficiency",             i = quote(g_move <= tm & meas_gap <= .10 & arm >= 44)),
  list(l = "NCAA recipe: both substituted",                  i = quote(g_move <= tm & inf_gap <= ti & arm >= 44)))
print(rbindlist(lapply(recipes, function(rc) { i <- Q[, eval(rc$i)]
  if (sum(i) < 5) return(data.table(recipe = rc$l, n = sum(i), whiff = NA_real_, p = NA_real_,
                                    velo_r = NA_real_, velo_p = NA_real_, overlap = NA_integer_))
  t <- t.test(Q$w[i], Q$w[!i]); ct <- cor.test(Q$velo_sep[i], Q$w[i])
  data.table(recipe = rc$l, n = sum(i), whiff = round(diff(rev(t$estimate)),2),
             p = round(t$p.value,4), velo_r = round(ct$estimate,3), velo_p = round(ct$p.value,4),
             overlap = sum(i & Q[, eval(recipes[[1]]$i)])) })), row.names = FALSE)

## ---- figure ------------------------------------------------------------------------------------
PL <- copy(P); PL[, pt := factor(pitch_type, c("CH","FF"), c("Changeups","Four-seamers"))]
# Binned medians go on top of the smoother because the smoother alone misreads the sparse left
# tail on changeups, where the relationship genuinely reverses rather than merely flattening.
MB <- PL[, .(m = median(disc), n = .N), by = .(pt, b = round(pmin(pmax(eff,.4),.99)/.025)*.025)][n >= 12]
rl <- PL[, .(lab = sprintf("above .75:  r = %+.3f\nbelow .75:  r = %+.3f",
                           cor(disc[eff>=.75], eff[eff>=.75]), cor(disc[eff<.75], eff[eff<.75]))), by = pt]
gA <- ggplot(PL, aes(eff, disc)) +
  geom_point(alpha = .13, size = .8, colour = "grey35") +
  geom_vline(xintercept = .75, linetype = 2, colour = "#b4632a", linewidth = .6) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), colour = "#1d7870",
              fill = "#1d787033", linewidth = .9) +
  geom_point(data = MB, aes(b, m), colour = "#0d3f3a", size = 1.9, shape = 18) +
  geom_text(data = rl, aes(x = .40, y = 62, label = lab), hjust = 0, vjust = 1, size = 3.2,
            lineheight = 1.15, fontface = "bold", colour = "#1d7870") +
  facet_wrap(~pt) + coord_cartesian(ylim = c(0, 68)) +
  labs(subtitle = "Angle between the measured spin axis and the Magnus direction, against measured active spin. Diamonds are binned medians; dashed line marks .75",
       x = "Measured active spin (Savant leaderboards)", y = "Axis discrepancy (degrees)") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 9.5), panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8.3))

B <- P[pitch_type == "CH", .(med = median(disc), lo = quantile(disc,.25), hi = quantile(disc,.75),
                             n = .N), by = ebin][order(ebin)]
gB <- ggplot(B, aes(ebin, med)) +
  geom_linerange(aes(ymin = lo, ymax = hi), colour = "#1d787066", linewidth = 3.5) +
  geom_line(aes(group = 1), colour = "#1d7870", linewidth = .8) +
  geom_point(size = 3, colour = "#1d7870") +
  geom_text(aes(label = sprintf("%.1f", med)), vjust = -1.4, size = 3.1, fontface = "bold") +
  geom_text(aes(y = 0, label = paste0("n=", n)), vjust = 1.4, size = 2.7, colour = "grey45") +
  coord_cartesian(ylim = c(-2, 34), clip = "off") +
  labs(subtitle = "Changeups only: median discrepancy by efficiency band, bars are the interquartile range",
       x = "Measured active spin", y = "Axis discrepancy (degrees)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 8.6),
        plot.margin = margin(5, 5, 14, 5))

suppressPackageStartupMessages(library(patchwork))
gg <- (gA / gB) + plot_layout(heights = c(1, .85)) +
  plot_annotation(
    title = "Yes, and it is the better of the two signals - but on changeups it works above .75 active spin and reverses below",
    subtitle = paste0("Statcast measures a spin axis from the ball's rotation and separately gives enough of the trajectory to compute where the ball actually broke. Removing\n",
                      "gravity and splitting what remains into components along and across the velocity vector isolates the Magnus force exactly, which is cleaner than using pfx\n",
                      "because it also removes drag - and drag differs between a changeup and a four-seamer, so it would otherwise contaminate a comparison between them. The\n",
                      "angle between the measured axis and that Magnus direction is near zero when a pitch is almost all transverse spin and opens up as gyro spin grows: on\n",
                      "changeups the median runs 5.7 degrees above .95 active spin and 16.4 in the .85-to-.90 band. Across that upper range it beats the break-magnitude route\n",
                      "outright, at r = -.63 against +.43. Below .75 it stops working and the sign flips to +.29. That is the region where seams, not spin, generate much of the\n",
                      "movement, so the axis stops predicting the break - and it is exactly where the kick-change population lives, which is why Pedro Avila's true .24 comes back\n",
                      "as .56 off a discrepancy of only 9 degrees. Use break magnitude for that tail. On four-seamers the sign holds throughout. Two cautions. The alignment\n",
                      "between the two conventions is REFLECTED, not a simple rotation; an earlier version of this test removed only a rotation and the resulting systematic error\n",
                      "was larger than the effect, which made the correlation look like zero. And this route is specific to Hawk-Eye - TrackMan derives its axis from movement, so\n",
                      "the discrepancy there is near zero by construction and carries no signal at all."),
    caption = "Source: MLB Statcast 2020-2026 - pitcher-seasons with 40+ pitches of the type - Magnus isolated from vx0/vy0/vz0 and ax/ay/az - truth is Savant measured active spin",
    theme = theme(plot.title = element_text(face = "bold", size = 11.5),
                  plot.subtitle = element_text(size = 8.2),
                  plot.caption = element_text(size = 7.5, colour = "grey40")))
ggsave(file.path(AST, "fig41_axis_vs_efficiency.png"), gg, width = 11.5, height = 9.4, dpi = 150)
cat("\nwrote fig41_axis_vs_efficiency.png\n")
