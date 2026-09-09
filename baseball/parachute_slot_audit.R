#!/usr/bin/env Rscript

# AUDITING THE SLOT INTERACTION BEFORE BELIEVING IT.
#
# The high-slot split produced the first surviving result in this line of work: among
# top-quartile arm angles, a smaller changeup spin-axis gap predicts whiffs above a
# shape-and-location model at r = -0.191 (p = .01), while low-slot pitchers show nothing, and
# the axis x arm interaction lands at p = .002. That is exactly the shape of the result that
# was retracted earlier in this project, so it gets the same four tests that killed the last one.
#
#   A  MULTIPLICITY. Nine stratum tests were run. Report which survive Bonferroni.
#   B  CONFOUNDS. High-slot pitchers differ from low-slot pitchers in ways the outcome model
#      does not absorb. Re-run the interaction controlling for spin rate, drop, velocity
#      separation, release height and handedness, and check it is not a release-height effect
#      wearing an arm-angle costume.
#   C  WITHIN PITCHER. The decisive one. Between-pitcher correlations are what failed before.
#      If the mechanism is real, a high-slot pitcher who closes his axis gap should gain more
#      than a low-slot pitcher who closes his by the same amount.
#   D  HOLD-OUT. Fit the interaction on 2023H2-2025 and test it on 2026.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
set.seed(1); options(width = 210)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")

L <- readRDS(file.path(MDIR, "parachute_within.rds"))
S <- L$CH[, .(np = .N, nsw = sum(is_swing), axis = mean(axis_diff), kill = -mean(az_diff),
              hmove = mean(ax_diff), velo_sep = -mean(speed_diff),
              spin = mean(release_spin_rate, na.rm = TRUE), rv100 = 100*mean(rv_res)),
          by = .(pitcher, player_name, season)]
S <- merge(S, L$SWD[, .(wh = 100*mean(wh_res)), by = .(pitcher, season)], by = c("pitcher","season"))
S <- merge(S, L$BIP[, .(gb = 100*mean(gb_res)), by = .(pitcher, season)], by = c("pitcher","season"))
aa <- fread(file.path(MDIR, "arm_angle_tunnel.csv"))[pitch_type == "CH",
        .(pitcher, season, arm = mean_arm, relz = mean_relz, relx = mean_relx)]
S <- merge(S, unique(aa, by = c("pitcher","season")), by = c("pitcher","season"))
S <- S[nsw >= 60 & is.finite(axis) & is.finite(arm)]
qs <- quantile(S$arm, c(.25,.75)); S[, hi_slot := arm >= qs[2]]

cat(sprintf("n = %d pitcher-seasons with arm angle and >= 60 changeup swings\n\n", nrow(S)))

## ---- A. multiplicity --------------------------------------------------------------
cat("=== A. multiplicity: 9 stratum tests were run, 3 outcomes x 3 slot strata ===\n")
S[, slot3 := cut(arm, breaks = c(-Inf, qs[1], qs[2], Inf),
                 labels = c("low","mid","HIGH"))]
gr <- CJ(v = c("rv100","wh","gb"), s = c("low","mid","HIGH"))
gr[, c("r","p") := { z <- mapply(function(vv, ss) {
    k <- suppressWarnings(cor.test(S[slot3 == ss][[vv]], S[slot3 == ss]$axis, method = "spearman"))
    c(unname(k$estimate), k$p.value) }, v, s); list(z[1,], z[2,]) }]
gr[, `:=`(bonf = p < .05/9, raw = p < .05)]
print(gr[order(p)][, .(outcome = v, slot = s, r = round(r,3), p = round(p,4), raw, bonf)],
      row.names = FALSE)
cat(sprintf("  -> %d of 9 clear p<.05; %d clear Bonferroni (p < %.4f)\n", sum(gr$raw), sum(gr$bonf), .05/9))

## ---- B. is it arm angle, or release height / spin / shape? ------------------------
cat("\n=== B. the interaction with controls stacked on ===\n")
mods <- list(
  "axis x arm, bare"                 = "wh ~ axis * arm",
  "+ drop, velo sep, spin"           = "wh ~ axis * arm + kill + velo_sep + spin",
  "+ release height and side"        = "wh ~ axis * arm + kill + velo_sep + spin + relz + relx",
  "release HEIGHT swapped for arm"   = "wh ~ axis * relz + kill + velo_sep + spin",
  "+ n as weight (see below)"        = "wh ~ axis * arm + kill + velo_sep + spin + relz + relx")
for (i in seq_along(mods)) {
  m <- if (i == 5) lm(as.formula(mods[[i]]), data = S, weights = np) else
                   lm(as.formula(mods[[i]]), data = S)
  ky <- grep(":", rownames(summary(m)$coefficients), value = TRUE)[1]
  co <- summary(m)$coefficients[ky, ]
  cat(sprintf("  %-32s %-14s beta = %+.5f   p = %.4f\n", names(mods)[i], ky, co[1], co[4]))
}

## ---- C. within pitcher: the test that killed the last one -------------------------
PREV <- S[, .(pitcher, season = season + 1L, axis_p = axis, kill_p = kill, velo_p = velo_sep,
              spin_p = spin, wh_p = wh, rv_p = rv100, gb_p = gb, arm_p = arm)]
P <- merge(S, PREV, by = c("pitcher","season"))
P[, `:=`(d_axis = axis_p - axis, d_wh = wh - wh_p, d_rv = rv100 - rv_p, d_gb = gb - gb_p,
         d_kill = kill - kill_p, d_velo = velo_sep - velo_p, d_spin = spin - spin_p,
         arm_avg = (arm + arm_p)/2)]
P[, hi := arm_avg >= qs[2]]
cat(sprintf("\n=== C. within pitcher: %d season pairs, %d at high slot ===\n", nrow(P), sum(P$hi)))
cat("    d_axis > 0 means the gap CLOSED toward the fastball; positive d_wh means improvement.\n")
for (g in c(TRUE, FALSE)) {
  z <- suppressWarnings(cor.test(P[hi == g]$d_axis, P[hi == g]$d_wh, method = "spearman"))
  cat(sprintf("    %-10s n=%3d   r = %+.3f   p = %.3f\n",
              if (g) "HIGH slot" else "rest", sum(P$hi == g), z$estimate, z$p.value))
}
mi <- summary(lm(d_wh ~ d_axis * arm_avg + d_kill + d_velo + d_spin, data = P))$coefficients
cat(sprintf("    within-pitcher d_axis x arm interaction: beta = %+.5f   p = %.3f\n",
            mi["d_axis:arm_avg","Estimate"], mi["d_axis:arm_avg","Pr(>|t|)"]))

## ---- D. hold out 2026 --------------------------------------------------------------
cat("\n=== D. fit on 2023H2-2025, test on 2026 ===\n")
TR <- S[season <= 2025]; TE <- S[season == 2026]
for (g in list(c(TRUE,"HIGH slot"), c(FALSE,"rest"))) {
  a <- suppressWarnings(cor.test(TR[hi_slot == as.logical(g[1])]$axis,
                                 TR[hi_slot == as.logical(g[1])]$wh, method = "spearman"))
  b <- suppressWarnings(cor.test(TE[hi_slot == as.logical(g[1])]$axis,
                                 TE[hi_slot == as.logical(g[1])]$wh, method = "spearman"))
  cat(sprintf("  %-10s train r = %+.3f (p=%.3f, n=%d)    2026 held out r = %+.3f (p=%.3f, n=%d)\n",
      g[2], a$estimate, a$p.value, sum(TR$hi_slot == as.logical(g[1])),
      b$estimate, b$p.value, sum(TE$hi_slot == as.logical(g[1]))))
}

## ---- the redefined bin: high slot + matched spin, no drop requirement --------------
# Cease sits at 59.6 deg with a 5.2 deg gap and +18.2 whiff points over model but fails the
# frozen bin on IVB kill. If the mechanism is slot plus spin match, the kill filter is the
# part that was wrong. This is a redefinition after seeing data and is flagged as exploratory.
S[, para2 := arm >= qs[2] & axis <= 15]
cat(sprintf("\n=== EXPLORATORY: high slot (>= %.0f deg) + axis gap <= 15, no drop filter: n = %d ===\n",
            qs[2], sum(S$para2)))
for (v in c("wh","rv100","gb")) { tt <- t.test(S[para2 == TRUE][[v]], S[para2 == FALSE][[v]])
  cat(sprintf("  %-6s  %+.3f vs %+.3f   diff %+.3f   p = %.3f\n", v, mean(S[para2==TRUE][[v]]),
      mean(S[para2==FALSE][[v]]), diff(rev(tt$estimate)), tt$p.value)) }
print(S[para2 == TRUE][order(-wh), .(player_name, season, np, arm = round(arm,1),
      axis = round(axis,1), kill = round(kill,1), velo_sep = round(velo_sep,1),
      whiff_over = round(wh,1), rv100 = round(rv100,2))], row.names = FALSE)
fwrite(S[order(-para2, -wh)], file.path(AST, "ext_parachute_slot_audit.csv"))

## ---- figures -------------------------------------------------------------------------
# Labelling all 36 members of the redefined bin buries the plot, so this names the two
# canonical high-slot cases, the best and worst performers inside the corner, and the
# low-slot pitchers the old drop-based bin wrongly admitted.
KEEP <- c("Cease, Dylan","Feltner, Ryan","Blackburn, Paul","Bibee, Tanner","Kikuchi, Yusei",
          "Banks, Tanner","Jones, Jared","Rogers, Trevor","Uceta, Edwin","Kelly, Merrill",
          "Rodón, Carlos","Hendricks, Kyle")
S[, lab := fifelse(player_name %chin% KEEP, sub(",.*","", player_name), NA_character_)]
S[!is.na(lab), lab := fifelse(seq_len(.N) == which.max(np), lab, NA_character_), by = lab]
g1 <- ggplot(S, aes(arm, axis)) +
  annotate("rect", xmin = qs[2], xmax = Inf, ymin = -Inf, ymax = 15, alpha = .11, fill = "#2a9d8f") +
  geom_hline(yintercept = 15, linetype = "dashed", colour = "#c0392b", linewidth = .45) +
  geom_vline(xintercept = qs[2], linetype = "dashed", colour = "#2a9d8f", linewidth = .45) +
  geom_point(aes(size = np, colour = wh), alpha = .8) +
  geom_text_repel(aes(label = lab), size = 3, max.overlaps = Inf, seed = 7, box.padding = .55,
                  min.segment.length = 0, segment.size = .3, segment.colour = "grey55",
                  colour = "grey10", fontface = "bold") +
  scale_colour_gradient2(low = "#c0392b", mid = "grey85", high = "#1d7870", midpoint = 0,
                         name = "Whiff% over model") +
  scale_size_continuous(range = c(.9, 4.5), guide = "none") +
  annotate("text", x = 74, y = 17.2, hjust = 1, size = 3.3, colour = "#1d7870",
           fontface = "bold", label = "high slot AND fastball-matched spin") +
  labs(title = "The parachute bin was catching the wrong pitchers",
       subtitle = paste0("One point per pitcher-season, min 60 changeup swings, 2023H2-2026, coloured by whiff rate above a shape-and-location model that never sees the axis gap.\n",
                         "The old bin also required a large drop, which admitted low-slot pitchers like Rogers at 22 degrees and Uceta at 16, and excluded Cease, whose changeup matches\n",
                         "his fastball's axis to 5 degrees off a 60-degree slot but does not drop enough to qualify. Arm angle is measured from 2023, so Hellickson cannot appear here."),
       x = "Arm angle (degrees; 90 = overhead, 45 = three-quarters, 0 = sidearm)",
       y = "Mean spin-axis gap vs primary fastball (degrees)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 8.2), panel.grid.minor = element_blank(),
        legend.position = "right")
ggsave(file.path(AST, "fig18_parachute_slot.png"), g1, width = 12.5, height = 7.2, dpi = 150)

ann <- S[, { z <- suppressWarnings(cor.test(axis, wh, method = "spearman"))
             .(lab = sprintf("r = %+.3f   p = %.3f   n = %d", z$estimate, z$p.value, .N)) }, by = slot3]
LV <- c(low = "Low slot (bottom 25%)", mid = "Middle 50%", HIGH = "High slot (top 25%)")
S[, sl := factor(LV[as.character(slot3)], levels = LV)]; ann[, sl := factor(LV[as.character(slot3)], levels = LV)]
g2 <- ggplot(S, aes(axis, wh)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = .35) +
  geom_point(alpha = .4, size = 1.4, colour = "#1d3557") +
  geom_smooth(method = "lm", formula = y ~ x, colour = "#c0392b", fill = "#c0392b",
              alpha = .13, linewidth = .85) +
  geom_text(data = ann, aes(x = Inf, y = Inf, label = lab), hjust = 1.05, vjust = 1.6,
            size = 3.1, fontface = "bold", inherit.aes = FALSE) +
  facet_wrap(~ sl, nrow = 1) +
  labs(title = "Spin match only tracks whiffs when the pitcher throws from over the top",
       subtitle = paste0("Between-pitcher, one point per pitcher-season. A downward slope means a SMALLER axis gap goes with more whiffs than shape and location predict. The slope is\n",
                         "flat for the bottom 75 percent of arm angles and clearly negative for the top quartile; the axis-by-arm-angle interaction is p = .002 across the full sample.\n",
                         "This is a between-pitcher result of exactly the kind that failed audit earlier in this project, so see the within-pitcher and hold-out panels before relying on it."),
       x = "Mean spin-axis gap vs primary fastball (degrees)", y = "Whiff% over model") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 8.2), panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 10), panel.spacing.x = unit(15,"pt"))
ggsave(file.path(AST, "fig19_slot_interaction.png"), g2, width = 13, height = 5.6, dpi = 150)
cat("\nwrote fig18_parachute_slot.png, fig19_slot_interaction.png, ext_parachute_slot_audit.csv\n")
