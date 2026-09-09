#!/usr/bin/env Rscript

# TRACKMAN vs HAWK-EYE SPIN DATA, SIDE BY SIDE, ON ONE PITCHER WHO SPANS BOTH.
#
# The question is whether 2017-2019 spin data can be converted onto the same footing as
# 2020+ so the parachute bin can use six seasons instead of three. Kershaw is the test case
# because he threw the same four pitches through both regimes with a stable delivery.
#
# The thing to check is not whether the NUMBERS are comparable - a scale or offset would be
# easy to correct. It is whether they are the same KIND of measurement. Hawk-Eye observes the
# ball rotating and reports the axis it sees. Trackman did not observe rotation; Savant
# inferred an axis from the ball's movement. If that is right, then in the Trackman era the
# spin axis is a deterministic function of the movement, and the "axis gap" this whole
# analysis rests on would be a restatement of the movement gap rather than an independent cue.
#
# That is testable directly: reconstruct the axis implied by the observed break and see how
# tightly it reproduces the published spin_axis in each era. A near-perfect fit means the
# column carries no information beyond movement.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
options(width = 205)
AST <- "data/statcast_model/article_assets"
d <- readRDS("data/kershaw_era_compare.rds")
d[, era := fifelse(season <= 2019L, "Trackman (2018-19)", "Hawk-Eye (2021-22)")]
d <- d[!is.na(spin_axis) & !is.na(pfx_x) & !is.na(pfx_z) & pitch_type %in% c("FF","SL","CU","CH")]

## ---- 1. raw formatting, one pitch from each era ----------------------------------
cat("=== how a single pitch is reported in each era (Kershaw four-seam) ===\n")
SHOW <- c("season","pitch_type","release_speed","release_spin_rate","spin_axis","spin_dir",
          "spin_rate_deprecated","pfx_x","pfx_z","ax","az","arm_angle")
print(d[pitch_type == "FF"][order(season)][, .SD[1], by = season][, ..SHOW], row.names = FALSE)

cat("\n=== column population by era ===\n")
pop <- d[, lapply(.SD, function(x) sprintf("%.0f%%", 100*mean(!is.na(x)))),
         by = era, .SDcols = c("spin_axis","release_spin_rate","spin_dir",
                               "spin_rate_deprecated","arm_angle","pfx_x")]
print(pop, row.names = FALSE)

cat("\n=== season averages by pitch type ===\n")
print(d[, .(n = .N, velo = round(mean(release_speed, na.rm=TRUE),1),
            spin = round(mean(release_spin_rate, na.rm=TRUE)),
            spin_axis = round(mean(spin_axis, na.rm=TRUE)),
            pfx_x = round(mean(pfx_x),2), pfx_z = round(mean(pfx_z),2)),
        by = .(season, pitch_type)][order(pitch_type, season)], row.names = FALSE)

## ---- 2. is spin_axis just the movement direction? --------------------------------
# Savant's spin_axis is a clock face where 180 = pure backspin. The Magnus deflection points
# 90 degrees off the spin axis, so the axis implied purely by the observed break is:
d[, mv_axis := (atan2(-pfx_x, -pfx_z)*180/pi) %% 360]
d[, gap := { z <- abs(spin_axis - mv_axis) %% 360; pmin(z, 360 - z) }]
cat("\n=== published spin_axis vs the axis implied by the observed break ===\n")
print(d[, .(n = .N, mean_abs_gap = round(mean(gap),2), median_gap = round(median(gap),2),
            sd_gap = round(sd(gap),2), pct_within_1deg = sprintf("%.1f%%", 100*mean(gap < 1)),
            pct_within_5deg = sprintf("%.1f%%", 100*mean(gap < 5))),
        by = .(era, pitch_type)][order(pitch_type, era)], row.names = FALSE)
cat("\n  overall:\n")
print(d[, .(n = .N, mean_abs_gap = round(mean(gap),2),
            pct_within_1deg = sprintf("%.1f%%", 100*mean(gap < 1)),
            R2_of_movement_axis = round(1 - sum(gap^2)/sum((d$spin_axis-mean(d$spin_axis))^2), 4)),
        by = era], row.names = FALSE)

## ---- 3. what that means for spin efficiency --------------------------------------
cat("\n=== the consequence for spin efficiency ===\n")
cat("  Active-spin leaderboards begin in 2020; Savant publishes nothing for 2017-2019.\n")
cat("  Kershaw measured active spin: 2020 FF 89.2 CU 79.9 SL 48.7 | 2021 FF 83.5 CU 78.5 SL 46.0 CH 99.2\n")
cat("  2022 FF 86.6 CU 83.7 SL 43.2 CH 99.1. No equivalent exists pre-2020 to convert FROM.\n")
cat("  Efficiency can only be inferred pre-2020 by comparing observed break to the break a\n")
cat("  100%-efficient ball would produce - but that inference uses movement, and an earlier\n")
cat("  pass in this project found inferred spin_eff unreliable enough to discard.\n")

## ---- 4. the question that actually decides it: is the FB-relative gap stable? ----
# The bin never uses the absolute axis, only the gap to the primary fastball, and a
# difference cancels any offset common to both pitches. So a drift in the absolute number is
# survivable; a drift in the GAP is not.
FBx <- d[pitch_type == "FF", .(fb_axis = mean(spin_axis)), by = season]
S <- merge(d[, .(n = .N, axis = mean(spin_axis), velo = mean(release_speed, na.rm=TRUE),
                 pfx_x = mean(pfx_x), pfx_z = mean(pfx_z)), by = .(season, pitch_type)],
           FBx, by = "season")
S[, gap := { z <- abs(axis - fb_axis) %% 360; pmin(z, 360-z) }]
cat("\n=== the quantity the bin actually uses: gap to the primary fastball ===\n")
print(dcast(S[pitch_type != "FF"], pitch_type ~ season, value.var = "gap")[
        , lapply(.SD, function(x) if (is.numeric(x)) round(x,1) else x)], row.names = FALSE)
cat("  sample sizes:\n")
print(dcast(S[pitch_type != "FF"], pitch_type ~ season, value.var = "n"), row.names = FALSE)
cat("\n  velocity and break over the same span, to show the pitches themselves barely moved:\n")
print(dcast(S, pitch_type ~ season, value.var = "velo")[
        , lapply(.SD, function(x) if (is.numeric(x)) round(x,1) else x)], row.names = FALSE)

## ---- figure ------------------------------------------------------------------------
P <- S[pitch_type != "FF" & n >= 100]
# Seasons are mapped to an index so the era band can be drawn between 2019 and 2021; a
# discrete x scale will not accept a rectangle at a fractional position.
YR <- sort(unique(P$season)); P[, xi := match(season, YR)]
gg <- ggplot(P, aes(xi, gap, colour = pitch_type, group = pitch_type)) +
  annotate("rect", xmin = 0.4, xmax = 2.5, ymin = -Inf, ymax = Inf, alpha = .09, fill = "#c0392b") +
  annotate("text", x = 1.45, y = 176, label = "Trackman", size = 3.6, fontface = "bold",
           colour = "#c0392b") +
  annotate("text", x = 3.5, y = 176, label = "Hawk-Eye", size = 3.6, fontface = "bold",
           colour = "#1d7870") +
  scale_x_continuous(breaks = seq_along(YR), labels = YR, expand = expansion(add = .45)) +
  geom_line(linewidth = 1.05) +
  geom_point(size = 3.2) +
  geom_text(aes(label = sprintf("%.0f", gap)), vjust = -1.25, size = 3.1, fontface = "bold",
            show.legend = FALSE) +
  scale_colour_manual(values = c(SL = "#1d3557", CU = "#2a9d8f"), name = NULL,
                      labels = c(CU = "Curveball", SL = "Slider")) +
  coord_cartesian(ylim = c(0, 185)) +
  labs(title = "Kershaw's spin axis relative to his own fastball, across the tracking change",
       subtitle = paste0("Mean spin-axis gap to his four-seam, the exact quantity the parachute bin is built on. His curveball reads 120 degrees off the fastball in 2018 and 166 in 2022,\n",
                         "a 46-degree move, while the pitch itself went from 72.9 to 73.2 mph with essentially unchanged break. Twenty-five of those degrees arrive between 2018 and\n",
                         "2019, inside the Trackman era, so this is not a clean one-time offset that could be subtracted out. The changeup is excluded here: Kershaw threw 11 to 16 of\n",
                         "them per season, far too few to mean anything, which is its own reason he cannot settle this question alone."),
       x = NULL, y = "Spin-axis gap vs his own four-seam (degrees)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 8.2),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "top")
ggsave(file.path(AST, "fig27_trackman_vs_hawkeye.png"), gg, width = 11, height = 6.4, dpi = 150)
cat("\nwrote fig27_trackman_vs_hawkeye.png\n")
