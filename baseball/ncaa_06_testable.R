#!/usr/bin/env Rscript

# WHAT THE D1 DATA CAN ACTUALLY TEST.
#
# The spin-axis half of the parachute hypothesis is not testable here: TrackMan's axis is
# effectively the break direction, and under a movement-based definition only one MLB season
# would ever have qualified for the Core bin, so the D1 null on that cue is uninformative.
#
# Two claims survive that, because neither needs an imaged spin axis:
#
#   VELOCITY   does separation from the fastball predict whiff the model did not expect? In MLB
#              this was flat league-wide and positive only inside the bin. D1 has 1,556 seasons,
#              which is real power on the league-wide half of that statement.
#   EFFICIENCY inside the MLB Core bin, active-spin LEVEL predicted overperformance at r = +.50,
#              better than velocity did. Inferred efficiency validates well enough to carry this
#              (r = +.92 on four-seamers, +.75 on changeups), and 1,556 seasons is far more power
#              than the 21 the MLB finding rested on.
#
# Both are pre-specified from the MLB work rather than searched for here.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
set.seed(4); options(width = 205)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")
S <- fread(file.path(AST, "ncaa_parachute_seasons.csv"))
cat(sprintf("D1 pitcher-seasons: %d\n\n", nrow(S)))

ct <- function(x, y, lab, D = S) { r <- cor.test(D[[x]], D[[y]])
  data.table(test = lab, n = nrow(D), r = round(r$estimate,3),
             ci = sprintf("[%+.3f, %+.3f]", r$conf.int[1], r$conf.int[2]),
             slope = round(coef(lm(D[[y]] ~ D[[x]]))[2],3), p = signif(r$p.value,3)) }

cat("=== the two testable claims, at full D1 sample ===\n")
S[, eff_min := pmin(eff_ch, eff_ff)]
print(rbindlist(list(
  ct("velo_sep","w",  "velocity separation vs whiff residual"),
  ct("eff_ch","w",    "changeup inferred efficiency vs residual"),
  ct("eff_min","w",   "lower of the two efficiencies vs residual"),
  ct("eff_ff","w",    "four-seam inferred efficiency vs residual"),
  ct("axis","w",      "axis gap vs residual (movement-based, weak)"),
  ct("slot","w",      "arm-slot proxy vs residual"))), row.names = FALSE)

cat("\n=== the MLB comparison, same quantities ===\n")
R <- readRDS(file.path(MDIR, "whiff_tjstuff.rds")); AS <- readRDS(file.path(MDIR, "active_spin_long.rds"))
M <- R[, .(nsw = .N, axis = mean(axis_diff), arm = mean(arm_angle), velo_sep = -mean(speed_diff),
           w = 100*mean(r4)), by = .(pitcher, season)]
M <- merge(M, AS[pitch_type=="CH", .(pitcher, season, eff_ch = active_spin)], by = c("pitcher","season"))
M <- merge(M, AS[pitch_type=="FF", .(pitcher, season, eff_ff = active_spin)], by = c("pitcher","season"))
M <- M[nsw >= 60][, `:=`(eff_min = pmin(eff_ch, eff_ff), as_gap = round(eff_ch - eff_ff, 4))]
M[, core := axis <= 10 & abs(as_gap) <= .10 & arm >= 44]
print(rbindlist(list(
  ct("velo_sep","w","velocity separation, all MLB", M),
  ct("eff_ch","w",  "changeup active spin, all MLB", M),
  ct("eff_min","w", "lower active spin, all MLB", M),
  ct("velo_sep","w","velocity separation, MLB Core bin", M[core == TRUE]),
  ct("eff_min","w", "lower active spin, MLB Core bin", M[core == TRUE]))), row.names = FALSE)

# Efficiency is not equally informative everywhere. In MLB it mattered inside a matched-axis
# bin; D1 has the sample to ask whether it matters as a main effect and whether it strengthens
# where the changeup is closest to the fastball.
cat("\n=== does efficiency matter more where the changeup sits closest to the fastball? ===\n")
S[, ax_q := cut(axis, quantile(axis, 0:4/4), include.lowest = TRUE, labels = c("Q1 tightest","Q2","Q3","Q4 widest"))]
print(S[, { r <- cor.test(eff_min, w)
  .(seasons = .N, axis_range = sprintf("%.0f-%.0f deg", min(axis), max(axis)),
    r_eff = round(r$estimate,3), p_eff = signif(r$p.value,3),
    r_velo = round(cor(velo_sep, w),3)) }, by = ax_q][order(ax_q)], row.names = FALSE)

cat("\n=== and by season, to check nothing is being carried by one year ===\n")
print(S[, { r <- cor.test(eff_min, w)
  .(seasons = .N, r_eff = round(r$estimate,3), p = signif(r$p.value,3),
    r_velo = round(cor(velo_sep, w),3)) }, by = season][order(season)], row.names = FALSE)

## ---- figure ------------------------------------------------------------------------------------
mk <- function(D, xv, lab, src) { r <- cor.test(D[[xv]], D$w)
  data.table(source = src, cue = lab, r = r$estimate, lo = r$conf.int[1], hi = r$conf.int[2],
             n = nrow(D), p = r$p.value) }
P <- rbindlist(list(
  mk(S, "velo_sep","Velocity separation","NCAA D1 2023-25 (1,556 seasons)"),
  mk(S, "eff_min", "Spin efficiency (lower of the two)","NCAA D1 2023-25 (1,556 seasons)"),
  mk(S, "axis",    "Axis gap to the fastball","NCAA D1 2023-25 (1,556 seasons)"),
  mk(S, "slot",    "Arm slot","NCAA D1 2023-25 (1,556 seasons)"),
  mk(M, "velo_sep","Velocity separation","MLB 2020-26, all (1,086 seasons)"),
  mk(M, "eff_min", "Spin efficiency (lower of the two)","MLB 2020-26, all (1,086 seasons)"),
  mk(M, "axis",    "Axis gap to the fastball","MLB 2020-26, all (1,086 seasons)"),
  mk(M, "arm",     "Arm slot","MLB 2020-26, all (1,086 seasons)"),
  mk(M[core == TRUE], "velo_sep","Velocity separation","MLB Core bin (21 seasons)"),
  mk(M[core == TRUE], "eff_min", "Spin efficiency (lower of the two)","MLB Core bin (21 seasons)")))
P[, cue := factor(cue, levels = rev(c("Velocity separation","Spin efficiency (lower of the two)",
                                      "Axis gap to the fastball","Arm slot")))]
P[, source := factor(source, levels = c("MLB Core bin (21 seasons)","MLB 2020-26, all (1,086 seasons)",
                                        "NCAA D1 2023-25 (1,556 seasons)"))]
P[, sig := fifelse(p < .05, "p < .05", "not significant")]
gg <- ggplot(P, aes(r, cue, colour = sig)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = .2, linewidth = .5) +
  geom_point(size = 2.8) +
  geom_text(aes(label = sprintf("%+.2f", r)), vjust = -1.1, size = 2.9, show.legend = FALSE) +
  facet_wrap(~source, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("p < .05" = "#1d7870", "not significant" = "grey55"), name = NULL) +
  labs(title = "What ports from MLB to college, and what cannot be tested there at all",
       subtitle = paste0("Correlation between each cue and whiff above an out-of-fold pitch-grade model, with 95 percent intervals. The D1 model is built on the same tjStuff+ v3.0\n",
                         "features plus location and reaches R-squared .198. Two cautions govern the bottom panel. The axis gap is not comparable: TrackMan reports the break\n",
                         "direction rather than imaged rotation, and recomputing the MLB gap that way leaves only 0.9 percent of seasons inside ten degrees against 10.2 percent\n",
                         "measured - so the matched-axis bin that the MLB result was built on cannot be constructed from college data, and its absence here is not a failed\n",
                         "replication. Spin efficiency is inferred from spin rate and break magnitude rather than measured, validated at r = +.92 on four-seamers and +.75 on\n",
                         "changeups. What does port is the league-wide null on velocity separation, which is near zero in both populations - the MLB claim was never that\n",
                         "separation helps everyone, only that it helps inside a bin that college data cannot build."),
       x = "Correlation with whiff above model", y = NULL,
       caption = "Sources: NCAA D1 TrackMan 2023-2025 (4.2M pitches) and MLB Statcast 2020-2026 - four-seam anchor - 40+ changeup swings (D1) / 60+ (MLB)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 8.2),
        strip.text = element_text(face = "bold", size = 9.5), panel.grid.minor = element_blank(),
        legend.position = "top", plot.caption = element_text(size = 7.5, colour = "grey40"))
ggsave(file.path(AST, "fig40_ncaa_vs_mlb.png"), gg, width = 11, height = 8, dpi = 150)
cat("\nwrote fig40_ncaa_vs_mlb.png\n")
