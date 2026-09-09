#!/usr/bin/env Rscript

# THE PARACHUTE BIN WITH THE EFFICIENCY FILTER APPLIED.
#
# Definition is now all three conditions: arm angle >= 44 deg, spin-axis gap <= 10 deg from
# the primary fastball, and measured active spin within 10 points of the fastball. The third
# condition is what makes "same spin" mean the whole spin vector rather than its shadow on the
# clock face, and it removes Rodon (-21 points of active spin vs his fastball), Blackburn's
# 2023 and 2024 seasons (-12 and -14) and Feltner (-12).
#
# This filter removes two of the four best whiff performers in the old bin, so it is a hostile
# test of the whiff result rather than a favourable one. The efficiency threshold is swept
# rather than fixed, because Feltner sits at -0.12 and a 10-point cutoff is the only thing
# keeping him out.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
set.seed(1); options(width = 205)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")
A <- fread(file.path(AST, "ext_parachute_spin_efficiency.csv"))
A <- A[is.finite(as_gap)]
AXC <- 10; ARC <- 44; EFC <- 0.10
A[, para := axis <= AXC & arm >= ARC & abs(as_gap) <= EFC]

cat(sprintf("population %d pitcher-seasons with measured active spin on both pitches\n", nrow(A)))
cat(sprintf("=== ROSTER: arm >= %d, axis gap <= %d, |active-spin gap| <= %.2f  ->  n = %d ===\n",
            ARC, AXC, EFC, sum(A$para)))
print(A[para == TRUE][order(-wh), .(player_name, season, pitches = np, arm = round(arm,1),
        axis = round(axis,1), as_ch = round(as_ch,2), as_fb = round(as_fb,2),
        as_gap = round(as_gap,2), velo_sep = round(velo_sep,1), ivb_kill = round(kill,1),
        whiff_over = round(wh,1), gb_over = round(gb,1), rv100 = round(rv100,2))], row.names = FALSE)
cat(sprintf("\n  NOTE: Blackburn 2026 passes the filter at %+.2f. Only his 2023 and 2024 seasons fail.\n",
            A[para == TRUE & grepl("Blackburn", player_name)]$as_gap[1]))

cat("\n=== bin vs everyone else ===\n")
for (v in c("wh","gb","rv100")) {
  t <- t.test(A[para == TRUE][[v]], A[para == FALSE][[v]]); ci <- t$conf.int
  cat(sprintf("  %-6s  bin %+.3f   rest %+.3f   diff %+.3f   95%% CI [%+.3f, %+.3f]   p = %.4f\n",
      v, mean(A[para==TRUE][[v]]), mean(A[para==FALSE][[v]]), diff(rev(t$estimate)),
      ci[1], ci[2], t$p.value)) }

## ---- where does Feltner enter, and does it matter? --------------------------------
cat("\n=== sweeping the efficiency threshold (Feltner enters at 0.125) ===\n")
sw <- rbindlist(lapply(c(0.04, 0.06, 0.08, 0.10, 0.125, 0.15, 0.20, 1.00), function(e) {
  i <- A$axis <= AXC & A$arm >= ARC & abs(A$as_gap) <= e
  if (sum(i) < 3) return(NULL)
  tw <- t.test(A$wh[i], A$wh[!i]); tg <- t.test(A$gb[i], A$gb[!i])
  data.table(eff = e, n = sum(i), whiff = diff(rev(tw$estimate)), p_wh = tw$p.value,
             gb = diff(rev(tg$estimate)), p_gb = tg$p.value,
             who = paste(sort(unique(sub(",.*","", A$player_name[i]))), collapse = " ")) }))
print(sw[, .(eff, n, whiff = round(whiff,2), p_wh = round(p_wh,4),
             gb = round(gb,2), p_gb = round(p_gb,4))], row.names = FALSE)
cat("\n  membership at each threshold:\n")
for (i in seq_len(nrow(sw))) cat(sprintf("   %.3f (n=%2d): %s\n", sw$eff[i], sw$n[i], sw$who[i]))

## ---- is the surviving effect a gradient or just a bin contrast? -------------------
cat("\n=== continuous test inside the efficiency-matched population ===\n")
E <- A[abs(as_gap) <= EFC]
for (v in c("wh","gb")) {
  co <- summary(lm(as.formula(sprintf("%s ~ axis * arm + kill + velo_sep + spin", v)), data = E))$coefficients
  cat(sprintf("  %-3s  axis:arm  beta = %+.5f   p = %.4f   (n = %d efficiency-matched seasons)\n",
      v, co["axis:arm","Estimate"], co["axis:arm","Pr(>|t|)"], nrow(E))) }
sp <- function(x,y){ z <- suppressWarnings(cor.test(x,y,method="spearman"))
                     sprintf("%+.3f (p=%.3f)", z$estimate, z$p.value) }
H <- E[arm >= ARC]
cat(sprintf("  within efficiency-matched HIGH-slot seasons (n=%d): axis vs whiff %s, vs grounders %s\n",
            nrow(H), sp(H$axis, H$wh), sp(H$axis, H$gb)))

## ---- robustness ---------------------------------------------------------------------
cat("\n=== robustness of the filtered bin ===\n")
B <- A[para == TRUE]
for (v in c("wh","gb")) {
  lo <- sapply(unique(B$player_name), function(q) {
    i <- A$para & A$player_name != q; mean(A[[v]][i]) - mean(A[[v]][!A$para]) })
  cat(sprintf("  %-3s  edge %+.2f pp   leave-one-pitcher-out %+.2f to %+.2f   (most influential %s)\n",
      v, mean(B[[v]]) - mean(A[para==FALSE][[v]]), min(lo), max(lo), names(which.min(lo)))) }
mem <- A[season <= 2025 & para == TRUE, unique(pitcher)]
T <- A[season == 2026]; T[, pr := pitcher %in% mem]
if (sum(T$pr) >= 3) for (v in c("wh","gb","rv100")) {
  t <- t.test(T[pr==TRUE][[v]], T[pr==FALSE][[v]])
  cat(sprintf("  held out 2026 %-6s prior-bin %+.2f (n=%d)  rest %+.2f   p = %.3f\n",
      v, mean(T[pr==TRUE][[v]]), sum(T$pr), mean(T[pr==FALSE][[v]]), t$p.value)) }

fwrite(A[order(-para, -wh)], file.path(AST, "ext_parachute_filtered.csv"))

## ---- figure ---------------------------------------------------------------------------
A[, status := fifelse(para, "In the bin",
              fifelse(axis <= AXC & arm >= ARC, "Removed by the efficiency filter", "Everyone else"))]
A[, lab := fifelse(status != "Everyone else", sub(",.*","", player_name), NA_character_)]
A[!is.na(lab), lab := fifelse(seq_len(.N) == which.max(np), lab, NA_character_), by = .(lab, status)]
PAL <- c(`In the bin` = "#1d7870", `Removed by the efficiency filter` = "#c0392b",
         `Everyone else` = "grey78")
gg <- ggplot(A, aes(axis, as_gap)) +
  annotate("rect", xmin = -Inf, xmax = AXC, ymin = -EFC, ymax = EFC, alpha = .13, fill = "#2a9d8f") +
  geom_hline(yintercept = 0, colour = "grey45", linewidth = .4) +
  geom_vline(xintercept = AXC, linetype = "dashed", colour = "#c0392b", linewidth = .45) +
  geom_point(aes(size = np, colour = status), alpha = .8) +
  geom_text_repel(aes(label = lab, colour = status), size = 3, max.overlaps = Inf, seed = 3,
                  box.padding = .5, min.segment.length = 0, segment.size = .3,
                  segment.colour = "grey55", fontface = "bold", show.legend = FALSE) +
  scale_colour_manual(values = PAL, name = NULL,
                      breaks = c("In the bin","Removed by the efficiency filter")) +
  scale_size_continuous(range = c(.9, 4.3), guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  coord_cartesian(xlim = c(0, 40), ylim = c(-.30, .22)) +
  labs(title = "The parachute bin after requiring the spin to match in three dimensions, not just on the clock face",
       subtitle = paste0("One point per pitcher-season with measured active spin on both the changeup and the primary fastball, min 60 changeup swings. The shaded box is the bin: axis\n",
                         "within 10 degrees, active spin within 10 points, arm angle at or above 44 degrees (points outside the box on arm angle are not shown as bin members).\n",
                         "The filter removes Rodon, Feltner and Blackburn's 2023-24 seasons, two of which were the strongest whiff performers in the unfiltered version - so this is a\n",
                         "hostile test. Blackburn's 2026 changeup passes at -0.03 and stays."),
       x = "Spin-axis gap vs primary fastball (degrees)",
       y = "Active-spin gap vs primary fastball\n(negative = changeup is more gyro)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12.4), plot.subtitle = element_text(size = 8.1),
        panel.grid.minor = element_blank(), legend.position = "top")
ggsave(file.path(AST, "fig23_parachute_filtered.png"), gg, width = 12.5, height = 7.2, dpi = 150)

sw2 <- melt(sw[eff < 1], id.vars = c("eff","n"), measure.vars = list(c("whiff","gb"), c("p_wh","p_gb")),
            value.name = c("d","p"), variable.name = "out")
sw2[, out := factor(out, labels = c("Whiff% over model","Ground-ball% over model"))]
g2 <- ggplot(sw2, aes(factor(eff), d, fill = p < .05)) +
  geom_hline(yintercept = 0, linewidth = .4) +
  geom_col(width = .66, colour = "black", linewidth = .25) +
  geom_text(aes(label = sprintf("%+.1f\np=%.3f\nn=%d", d, p, n),
                vjust = fifelse(d >= 0, -0.18, 1.1)), size = 2.7, lineheight = .95) +
  geom_vline(xintercept = 4.5, linetype = "dashed", colour = "#c0392b", linewidth = .5) +
  annotate("text", x = 4.6, y = 11.5, hjust = 0, size = 3, colour = "#c0392b", fontface = "bold",
           label = "Feltner enters here") +
  scale_fill_manual(values = c(`FALSE` = "#c9ced6", `TRUE` = "#1d7870"), guide = "none") +
  facet_wrap(~ out, nrow = 1) + coord_cartesian(ylim = c(-1, 13)) +
  labs(title = "The grounder effect does not care where the efficiency line is drawn. The whiff effect never clears significance anywhere",
       subtitle = paste0("Each bar loosens the active-spin requirement by one step, holding arm angle at 44 degrees and the axis gap at 10. Teal bars are p < .05. The whiff edge wanders\n",
                         "between +2.3 and +4.5 with p ranging from .07 to .41 and is never significant at any cutoff; it moves non-monotonically as single pitchers enter, which is what a\n",
                         "nine-to-fifteen-season sample does. The ground-ball edge is p < .005 at all eight thresholds and is largest when the efficiency requirement is strictest."),
       x = "Maximum allowed active-spin gap vs the fastball", y = "Bin minus rest (percentage points)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12.2), plot.subtitle = element_text(size = 8.1),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        strip.text = element_text(face = "bold", size = 10), panel.spacing.x = unit(14, "pt"))
ggsave(file.path(AST, "fig24_efficiency_sweep.png"), g2, width = 13.2, height = 6, dpi = 150)
cat("\nwrote fig23_parachute_filtered.png, fig24_efficiency_sweep.png, ext_parachute_filtered.csv\n")
