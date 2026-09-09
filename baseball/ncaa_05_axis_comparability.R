#!/usr/bin/env Rscript

# IS THE D1 AXIS GAP WIDER THAN MLB'S BECAUSE D1 CHANGEUPS ARE DIFFERENT, OR BECAUSE THE TWO
# SYSTEMS MEASURE DIFFERENT THINGS?
#
# D1 changeups average a 31.5-degree gap to the four-seamer against 22.4 in MLB, which shrank
# the matched bin to 18 seasons and left the replication underpowered. Two explanations, with
# opposite consequences:
#
#   real       college changeups genuinely sit further off the fastball, and the archetype is
#              simply rarer at that level. The replication is then honest but underpowered.
#   artefact   TrackMan's axis is essentially the break direction (script 01), while Statcast's
#              is imaged rotation. If those two quantities have different dispersions, the same
#              10-degree threshold selects different things and the comparison is invalid.
#
# The test recomputes the MLB gap from MOVEMENT DIRECTION - the quantity TrackMan is really
# reporting - and asks whether MLB then looks like D1. If it does, the fix is to use the
# movement-based definition on both sides.

suppressPackageStartupMessages({ library(data.table) })
options(width = 200)
MDIR <- "data/statcast_model"

M <- readRDS(file.path(MDIR, "parachute_rv.rds"))[pitch_type %in% c("CH","FF")]
M[, spin_axis := (atan2(sax, cax)*180/pi) %% 360]
M[, mv_axis := (atan2(ax, az + 32.174)*180/pi) %% 360]
M <- M[is.finite(spin_axis) & is.finite(mv_axis)]

P <- M[, .(n = .N, sa = mean(spin_axis), mv_s = mean(sin(mv_axis*pi/180)),
           mv_c = mean(cos(mv_axis*pi/180))), by = .(pitcher, season, pitch_type)][n >= 30]
P[, mv := (atan2(mv_s, mv_c)*180/pi) %% 360]
W <- dcast(P, pitcher + season ~ pitch_type, value.var = c("sa","mv","n"))
W <- W[is.finite(sa_CH) & is.finite(sa_FF) & n_CH >= 50]
gap <- function(a, b) pmin(abs(a-b), 360-abs(a-b))
W[, `:=`(measured_gap = gap(sa_CH, sa_FF), movement_gap = gap(mv_CH, mv_FF))]

cat("=== MLB changeup-to-four-seam axis gap, two definitions ===\n")
cat(sprintf("  measured spin axis (Hawk-Eye):  mean %.1f  median %.1f  sd %.1f  share <= 10 deg %.1f%%\n",
    mean(W$measured_gap), median(W$measured_gap), sd(W$measured_gap), 100*mean(W$measured_gap <= 10)))
cat(sprintf("  movement direction:             mean %.1f  median %.1f  sd %.1f  share <= 10 deg %.1f%%\n",
    mean(W$movement_gap), median(W$movement_gap), sd(W$movement_gap), 100*mean(W$movement_gap <= 10)))
cat(sprintf("  correlation between the two definitions: r = %+.3f over %d pitcher-seasons\n",
            cor(W$measured_gap, W$movement_gap), nrow(W)))

cat("\n=== and D1, whose axis is effectively the movement direction ===\n")
S <- fread(file.path(MDIR, "article_assets", "ncaa_parachute_seasons.csv"))
cat(sprintf("  TrackMan SpinAxis gap:          mean %.1f  median %.1f  sd %.1f  share <= 10 deg %.1f%%\n",
    mean(S$axis), median(S$axis), sd(S$axis), 100*mean(S$axis <= 10)))

cat("\n=== reading ===\n")
d1 <- 100*mean(S$axis <= 10); mlb_m <- 100*mean(W$measured_gap <= 10); mlb_v <- 100*mean(W$movement_gap <= 10)
cat(sprintf("  share of seasons inside 10 degrees:  MLB measured %.1f%%  |  MLB movement %.1f%%  |  D1 %.1f%%\n",
            mlb_m, mlb_v, d1))
if (abs(mlb_v - d1) < abs(mlb_m - d1)) {
  cat("  D1 resembles the MLB MOVEMENT definition more than the measured one, which is what the\n")
  cat("  script-01 finding predicted. The two populations are closer than the raw numbers suggest,\n")
  cat("  and the honest comparison is movement-gap against movement-gap on both sides.\n")
} else {
  cat("  D1 does not line up with either MLB definition, so the difference is not purely one of\n")
  cat("  measurement and college changeups do appear to sit further off the fastball.\n")
}

# The consequence that matters: if the MLB parachute result was built on the measured axis and
# the two definitions disagree, then the movement definition would not have found it either -
# in which case D1 was never able to test the same hypothesis.
cat("\n=== would the MLB result have survived a movement-based axis gap? ===\n")
AS <- readRDS(file.path(MDIR, "active_spin_long.rds"))
R <- readRDS(file.path(MDIR, "whiff_tjstuff.rds"))
Q <- R[, .(nsw = .N, arm = mean(arm_angle), velo_sep = -mean(speed_diff), w = 100*mean(r4)),
       by = .(pitcher, season)]
Q <- merge(Q, W[, .(pitcher, season, measured_gap, movement_gap)], by = c("pitcher","season"))
Q <- merge(Q, AS[pitch_type=="CH", .(pitcher, season, as_ch = active_spin)], by = c("pitcher","season"))
Q <- merge(Q, AS[pitch_type=="FF", .(pitcher, season, as_fb = active_spin)], by = c("pitcher","season"))
Q <- Q[nsw >= 60][, as_gap := round(as_ch - as_fb, 4)]
for (g in c("measured_gap","movement_gap")) {
  i <- Q[[g]] <= 10 & abs(Q$as_gap) <= .10 & Q$arm >= 44
  if (sum(i) < 4) { cat(sprintf("  %-14s only %d seasons qualify\n", g, sum(i))); next }
  t <- t.test(Q$w[i], Q$w[!i]); ct <- cor.test(Q$velo_sep[i], Q$w[i])
  cat(sprintf("  %-14s n=%2d  whiff above model %+.2f pp (p=%.3f)  velo r=%+.3f (p=%.3f)\n",
      g, sum(i), diff(rev(t$estimate)), t$p.value, ct$estimate, ct$p.value)) }
