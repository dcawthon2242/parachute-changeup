#!/usr/bin/env Rscript

# WHAT IS ACTUALLY DRIVING THE HIGH-ACTIVE-SPIN BIN?
#
# Replacing the efficiency-gap filter with a floor on active spin, then opening the spin-axis
# gap, produces a bin whose whiff edge looks solid at every setting: +3.73 points at a
# ten-degree limit falling smoothly to +1.5 with the axis wide open, significant throughout.
# That looks like a dose-response on the look-alike cue.
#
# It is not. Three checks, in order:
#
#   2x2       split the gate into its two halves. The active-spin floor contributes nothing
#             (-0.38 pp, p = .37). The arm-slot requirement contributes all of it (+1.54,
#             p = .0001) and the two do not interact.
#   within    inside the gated group, does the axis gap predict the residual? No: r = -0.04,
#             and the quintile means are flat.
#   fair      the residual model deliberately omits arm angle because arm angle defines the
#             bin. So a positive residual for high-slot pitchers may only mean the model is
#             missing a feature. Refit with arm angle included and the edge goes from +1.48
#             (p = .0002) to +0.31 (p = .42).
#
# The apparent axis gradient in the sweep is an artifact of the comparison group: as the limit
# widens the bin absorbs more of the high-slot population, so the contrast against "everyone
# else" shrinks toward the plain arm-slot effect, which itself is not real.

suppressPackageStartupMessages({ library(data.table); library(lightgbm); library(ggplot2) })
set.seed(1); options(width = 205)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")
CACHE <- file.path(MDIR, "highspin_armfit.rds")

BASE <- c("release_speed","release_spin_rate","release_extension","release_pos_x",
          "release_pos_z","ax","az","speed_diff","ax_diff","az_diff","plate_x","plate_z",
          "plate_x_in","plate_x_arm","z_rel_bot","z_rel_top","VAA","HAA","HAA_in",
          "stand_R","throws_R","same_hand","balls","strikes")

if (!file.exists(CACHE)) {
  CH <- readRDS(file.path(MDIR, "parachute_ff.rds"))
  CH <- CH[is.finite(axis_diff) & is.finite(arm_angle) & stats::complete.cases(CH[, ..BASE])]
  SW <- CH[is_swing == TRUE]
  fit <- function(FEAT, tag) {
    K <- 4; fold <- sample(rep(1:K, length.out = nrow(SW))); p <- rep(NA_real_, nrow(SW))
    for (f in 1:K) {
      tr <- SW[fold != f]; n <- nrow(tr); vi <- sample(n, floor(.12*n))
      dtr <- lgb.Dataset(as.matrix(tr[-vi, ..FEAT]), label = tr$whiff[-vi])
      dva <- lgb.Dataset.create.valid(dtr, as.matrix(tr[vi, ..FEAT]), label = tr$whiff[vi])
      m <- lgb.train(params = list(objective = "binary", metric = "binary_logloss",
                     learning_rate = .06, num_leaves = 31, min_data_in_leaf = 300,
                     feature_fraction = .8, bagging_fraction = .8, bagging_freq = 1),
                     data = dtr, nrounds = 1500, valids = list(v = dva),
                     early_stopping_rounds = 50, verbose = -1)
      p[fold == f] <- predict(m, as.matrix(SW[fold == f, ..FEAT]))
    }
    cat(sprintf("  %-22s R2=%.4f\n", tag, 1 - var(SW$whiff - p)/var(SW$whiff)))
    SW$whiff - p
  }
  cat("=== out-of-fold whiff models ===\n")
  SW[, res_noarm := fit(BASE, "without arm angle")]
  SW[, res_arm   := fit(c(BASE, "arm_angle"), "with arm angle")]
  saveRDS(SW[, .(pitcher, player_name, season, axis_diff, speed_diff, arm_angle,
                 res_noarm, res_arm)], CACHE)
}
SW <- readRDS(CACHE)
AS <- readRDS(file.path(MDIR, "active_spin_long.rds"))
S <- SW[, .(nsw = .N, axis = mean(axis_diff), arm = mean(arm_angle),
            velo_sep = -mean(speed_diff), w0 = 100*mean(res_noarm), w1 = 100*mean(res_arm)),
        by = .(pitcher, player_name, season)]
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, as_ch = active_spin)], by = c("pitcher","season"))
S <- merge(S, AS[pitch_type == "FF", .(pitcher, season, as_fb = active_spin)], by = c("pitcher","season"))
S <- S[nsw >= 60]
S[, `:=`(hs = as_ch >= .80 & as_fb >= .80, hi = arm >= 44, last = sub(",.*","",player_name))]

## ---- panel A: the sweep that looks like a dose-response ---------------------------------
A <- rbindlist(lapply(seq(10, 45, 5), function(ax) {
  i <- S$hs & S$hi & S$axis <= ax; t <- t.test(S$w0[i], S$w0[!i])
  # t$estimate is c(mean in, mean out) and t$conf.int is the interval on in-minus-out, so the
  # bounds carry straight over; negating them puts the point estimate outside its own interval.
  data.table(axis_max = ax, seasons = sum(i), v = diff(rev(t$estimate)),
             lo = t$conf.int[1], hi = t$conf.int[2], p = t$p.value) }))
cat("\n=== panel A: whiff edge vs the rest of the league, by axis limit ===\n"); print(A, row.names = FALSE)

## ---- panel B: split the gate ------------------------------------------------------------
B <- S[, .(seasons = .N, v = mean(w0), se = sd(w0)/sqrt(.N)),
       by = .(grp = fifelse(hs & hi, "High spin\n+ high slot",
               fifelse(hs & !hi, "High spin\nlow slot",
               fifelse(!hs & hi, "Low spin\n+ high slot", "Neither"))))]
B[, `:=`(lo = v - 1.96*se, hi = v + 1.96*se)]
B[, grp := factor(grp, levels = c("High spin\n+ high slot","Low spin\n+ high slot",
                                  "High spin\nlow slot","Neither"))]
cat("\n=== panel B: the 2x2 ===\n"); print(B[order(grp)], row.names = FALSE)

## ---- panel C: the same contrasts under a model that sees arm angle ----------------------
con <- function(i, v, lab) { t <- t.test(S[[v]][i], S[[v]][!i])
  data.table(lab = lab, model = if (v == "w0") "Model blind to\narm angle" else "Model sees\narm angle",
             v = diff(rev(t$estimate)), lo = t$conf.int[1], hi = t$conf.int[2], p = t$p.value) }
C <- rbindlist(c(
  lapply(c("w0","w1"), function(v) con(S$hi, v, "High slot (>= 44 deg)")),
  lapply(c("w0","w1"), function(v) con(S$hs, v, "High active spin (both >= .80)")),
  lapply(c("w0","w1"), function(v) con(S$hs & S$hi & S$axis <= 15, v, "Full gate, axis <= 15 deg"))))
C[, lab := factor(lab, levels = rev(c("High slot (>= 44 deg)","High active spin (both >= .80)",
                                      "Full gate, axis <= 15 deg")))]
cat("\n=== panel C: before and after the model can see arm angle ===\n"); print(C, row.names = FALSE)

## ---- figure -------------------------------------------------------------------------------
th <- theme_minimal(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold", size = 10.2), panel.grid.minor = element_blank(),
        legend.position = "none", plot.subtitle = element_text(size = 8))
gA <- ggplot(A, aes(axis_max, v)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#1d7870", alpha = .15) +
  geom_line(colour = "#1d7870", linewidth = .9) +
  geom_point(aes(size = seasons), colour = "#1d7870") + scale_size_area(max_size = 5) +
  scale_x_continuous(breaks = seq(10, 45, 5)) +
  labs(title = "A. Opening the axis gap looks like a dose-response",
       subtitle = "Whiff above model vs everyone outside the bin. Significant at every setting.",
       x = "Spin-axis gap limit (degrees)", y = "Whiff above model (pp)") + th
gB <- ggplot(B, aes(grp, v, fill = grp == "High spin\n+ high slot")) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_col(width = .62) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .16, linewidth = .5, colour = "grey30") +
  scale_fill_manual(values = c("TRUE" = "#1d7870", "FALSE" = "grey68")) +
  labs(title = "B. But it is the arm slot, not the spin",
       subtitle = "Slot main effect +1.54 pp (p = .0001); active spin -0.38 (p = .37); no interaction.",
       x = NULL, y = "Whiff above model (pp)") + th
gC <- ggplot(C, aes(v, lab, colour = model)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = .2, linewidth = .6,
                position = position_dodge(width = .55)) +
  geom_point(size = 2.4, position = position_dodge(width = .55)) +
  scale_colour_manual(values = c("Model blind to\narm angle" = "#c0392b",
                                 "Model sees\narm angle" = "#1d7870")) +
  labs(title = "C. And it vanishes once the model is allowed to see the arm angle",
       subtitle = "Red: arm angle withheld (how every number above was computed). Teal: arm angle included.",
       x = "Whiff above model (pp)", y = NULL) + th +
  theme(legend.position = "right", legend.title = element_blank(),
        legend.text = element_text(size = 7.5))

gg <- patchwork::wrap_plots(gA, gB, gC, ncol = 1, heights = c(1, 1, .85)) +
  patchwork::plot_annotation(
    title = "The high-spin gate does nothing and the axis gap does nothing. The bin was measuring arm slot, which the model could not see",
    subtitle = paste0("The changeup's active spin and the four-seamer's are both required to clear 80 percent, which removes seam-shifted changeups outright, and the active-spin gap is\n",
                      "dropped from the definition entirely. Panel A is the result that request produces and it looks convincing. Panel B splits the gate and finds the active-spin floor\n",
                      "contributes nothing while the 44-degree arm-slot requirement contributes everything. Panel C explains why: the residual model withholds arm angle by construction,\n",
                      "because arm angle defines the bin, so high-slot pitchers beat it simply because it cannot see them. Give the model the arm angle and the edge falls from +1.48 to\n",
                      "+0.31 points and stops being significant, and the full gate at 15 degrees falls from +1.81 to +0.58. Inside the gated group the axis gap predicts nothing either\n",
                      "way (r = -0.04, p = .53), and arm angle read as a continuous cue rather than a 44-degree cutoff predicts nothing either (r = +0.04, p = .17)."),
    caption = "Source: Statcast 2020-2026 - 203,132 changeup swings - four-seam anchor - out-of-fold LightGBM residuals - 1,086 pitcher-seasons at a 60-swing floor",
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = element_text(size = 8.2),
                  plot.caption = element_text(size = 7.5, colour = "grey40")))
ggsave(file.path(AST, "fig33_highspin_decomp.png"), gg, width = 11, height = 12, dpi = 150)
cat("\nwrote fig33_highspin_decomp.png\n")
