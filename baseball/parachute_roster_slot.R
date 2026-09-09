#!/usr/bin/env Rscript

# WHO HAS ACTUALLY BEEN IN THE PARACHUTE BIN, AND DOES THE HIGH-SLOT CRITERION DESCRIBE THEM?
#
# The bin so far is purely kinematic - mean spin-axis gap <= 15 deg from the primary fastball
# and IVB kill >= 13.9 - and says nothing about delivery. The claim on the table is that a
# real parachute changeup requires an over-the-top slot. That is a mechanism claim, and it is
# checkable two ways:
#
#   1. Descriptively. Print the roster with each pitcher's Hawk-Eye arm angle next to the
#      league distribution, and see whether the bin is already a high-slot group or whether
#      it is full of low-slot pitchers who happen to satisfy the kinematics.
#
#   2. As an interaction, which is the only properly powered version. Adding a third filter to
#      an 18-season bin leaves nothing to test. Instead ask, across ALL changeup throwers,
#      whether the axis-gap effect strengthens as the slot rises. If the parachute needs a high
#      slot, axis gap should predict performance among over-the-top pitchers and not among
#      low-slot ones, and the interaction term carries that.
#
# Statcast arm angle: 90 = directly overhead, ~45 = three-quarters, 0 = sidearm, negative =
# submarine. Data begins in 2023, so Hellickson (last MLB pitch 2020) cannot be measured here.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
options(width = 210)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")

L <- readRDS(file.path(MDIR, "parachute_within.rds"))
CH <- L$CH; SWD <- L$SWD; BIP <- L$BIP
S <- CH[, .(np = .N, nsw = sum(is_swing), axis = mean(axis_diff), kill = -mean(az_diff),
            velo_sep = -mean(speed_diff), spin = mean(release_spin_rate, na.rm = TRUE),
            rv100 = 100*mean(rv_res)), by = .(pitcher, player_name, season)]
S <- merge(S, SWD[, .(wh = 100*mean(wh_res)), by = .(pitcher, season)], by = c("pitcher","season"))
S <- merge(S, BIP[, .(gb = 100*mean(gb_res)), by = .(pitcher, season)], by = c("pitcher","season"))

aa <- fread(file.path(MDIR, "arm_angle_tunnel.csv"))[pitch_type == "CH",
        .(pitcher, season, arm = mean_arm)]
S <- merge(S, unique(aa, by = c("pitcher","season")), by = c("pitcher","season"), all.x = TRUE)
S <- S[nsw >= 60 & is.finite(axis) & is.finite(kill)]
S[, para := axis <= 15 & kill >= 13.9]

## ---- 1. the roster ---------------------------------------------------------------
qs <- quantile(S$arm, c(.1,.25,.5,.75,.9), na.rm = TRUE)
cat(sprintf("Arm angle among all %d qualifying changeup pitcher-seasons (deg, 90 = overhead):\n", nrow(S)))
cat(sprintf("  10th %.0f   25th %.0f   median %.0f   75th %.0f   90th %.0f\n\n",
            qs[1], qs[2], qs[3], qs[4], qs[5]))
HI <- unname(qs[4])   # "high slot" = top quartile of actual MLB changeup throwers
cat(sprintf("=== the frozen parachute bin (axis <= 15 deg, IVB kill >= 13.9), n = %d ===\n", sum(S$para)))
R <- S[para == TRUE][order(-arm)][, .(player_name, season, np, arm = round(arm,1),
        axis = round(axis,1), kill = round(kill,1), velo_sep = round(velo_sep,1),
        rv100 = round(rv100,2), whiff_over = round(wh,1), gb_over = round(gb,1),
        slot = fifelse(arm >= HI, "HIGH", fifelse(arm >= qs[3], "above med", "low/avg")))]
print(R, row.names = FALSE)
cat(sprintf("\n  %d of %d bin members are top-quartile slot (>= %.0f deg). Bin mean arm %.1f vs league %.1f.\n",
            sum(R$slot == "HIGH"), nrow(R), HI, mean(S[para == TRUE]$arm, na.rm=TRUE),
            mean(S$arm, na.rm = TRUE)))
t <- t.test(S[para == TRUE]$arm, S[para == FALSE]$arm)
cat(sprintf("  parachute vs rest on arm angle: %+.2f deg, p = %.3f\n", diff(rev(t$estimate)), t$p.value))

## ---- the two pitchers named ------------------------------------------------------
cat("\n=== the two names raised ===\n")
for (nm in c("Cease","Vesia","Hellickson")) {
  z <- S[grepl(nm, player_name)]
  if (!nrow(z)) { cat(sprintf("  %-11s no qualifying changeup seasons in the 2023-2026 window\n", nm)); next }
  print(z[order(season), .(player_name, season, np, arm = round(arm,1), axis = round(axis,1),
          kill = round(kill,1), velo_sep = round(velo_sep,1), rv100 = round(rv100,2),
          whiff_over = round(wh,1), in_bin = para)], row.names = FALSE)
}

## ---- 2. the interaction, on the full population ----------------------------------
A <- S[is.finite(arm)]
cat(sprintf("\n=== does the axis-gap effect depend on slot? (n = %d pitcher-seasons with arm angle) ===\n", nrow(A)))
for (v in c("rv100","wh","gb")) {
  f <- summary(lm(as.formula(sprintf("%s ~ axis * arm + kill + velo_sep + spin", v)), data = A))
  co <- f$coefficients
  cat(sprintf("  %-6s  axis:arm interaction  beta = %+.5f   p = %.3f      (axis main %+.4f, p=%.3f)\n",
              v, co["axis:arm","Estimate"], co["axis:arm","Pr(>|t|)"],
              co["axis","Estimate"], co["axis","Pr(>|t|)"]))
}
# Same question without assuming linearity: split the population at the slot quartiles and
# correlate axis gap with performance inside each stratum.
A[, slot3 := cut(arm, breaks = c(-Inf, qs[2], qs[4], Inf),
                 labels = c("low slot (bottom 25%)","middle 50%","HIGH slot (top 25%)"))]
sp <- function(x,y){ z <- suppressWarnings(cor.test(x,y,method="spearman"))
                     sprintf("%+.3f (p=%.2f)", z$estimate, z$p.value) }
cat("\n  Spearman of axis gap with each outcome, by slot stratum (negative = smaller gap is better):\n")
print(A[, .(n = .N, rv = sp(axis, rv100), whiff = sp(axis, wh), gb = sp(axis, gb)),
        by = slot3][order(slot3)], row.names = FALSE)
cat("\n  Within the HIGH-slot stratum only, parachute vs rest:\n")
H <- A[slot3 == "HIGH slot (top 25%)"]
for (v in c("rv100","wh","gb")) { tt <- t.test(H[para == TRUE][[v]], H[para == FALSE][[v]])
  cat(sprintf("    %-6s  parachute %+.3f (n=%d)  rest %+.3f (n=%d)   p = %.3f\n", v,
      mean(H[para==TRUE][[v]]), sum(H$para), mean(H[para==FALSE][[v]]), sum(!H$para), tt$p.value)) }

fwrite(S[order(-para, -arm)], file.path(AST, "ext_parachute_roster_slot.csv"))
cat("\nwrote ext_parachute_roster_slot.csv\n")

## ---- figure ----------------------------------------------------------------------
S[, lab := fifelse(para, sub(",.*","", player_name), NA_character_)]
gg <- ggplot(S[is.finite(arm)], aes(arm, axis)) +
  annotate("rect", xmin = HI, xmax = Inf, ymin = -Inf, ymax = 15, alpha = .1, fill = "#2a9d8f") +
  geom_hline(yintercept = 15, linetype = "dashed", colour = "#c0392b", linewidth = .45) +
  geom_vline(xintercept = HI, linetype = "dashed", colour = "#2a9d8f", linewidth = .45) +
  geom_point(aes(size = np, colour = para), alpha = .55) +
  ggrepel::geom_text_repel(aes(label = lab), size = 2.7, max.overlaps = 30, seed = 1,
                           segment.size = .25, colour = "grey20") +
  scale_colour_manual(values = c(`FALSE` = "#9aa5b1", `TRUE` = "#1d3557"), guide = "none") +
  scale_size_continuous(range = c(.8, 4), guide = "none") +
  annotate("text", x = HI + 1, y = 3, hjust = 0, size = 3.2, colour = "#1d7870",
           fontface = "bold", label = "high slot AND fastball-matched spin") +
  labs(title = "Arm slot against spin-axis match: the parachute bin is not a high-slot group",
       subtitle = paste0("One point per pitcher-season, min 60 changeup swings, 2023H2-2026. Dashed lines are the 15-degree axis-gap cutoff that defines the bin and the top-quartile\n",
                         "arm angle among actual changeup throwers. Statcast measures arm angle from 2023, so Hellickson cannot appear. If the parachute changeup required an\n",
                         "over-the-top delivery, the labelled points would sit in the shaded corner."),
       x = "Arm angle (degrees; 90 = directly overhead, 45 = three-quarters, 0 = sidearm)",
       y = "Mean spin-axis gap vs primary fastball (degrees)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12.5),
        plot.subtitle = element_text(size = 8.2), panel.grid.minor = element_blank())
ggsave(file.path(AST, "fig18_parachute_slot.png"), gg, width = 12, height = 7, dpi = 150)
cat("wrote fig18_parachute_slot.png\n")
