#!/usr/bin/env Rscript

# THE FIGURE, AND WHERE MARTINEZ ACTUALLY SITS.
#
# Panel A is the claim: bin pitcher-seasons by velocity separation and plot how far actual whiff
# rate lands from each model's grade. A fastball-blind model should sit on zero everywhere if it is
# well calibrated. It does not - it runs three points low on the big separators and three points
# high on the pitchers who barely change speeds.
#
# Panel B is the caveat. Martinez is the case that prompted this, and separation turns out to
# explain only a fifth of his edge. He beats the separation-aware model nearly as badly as he beats
# the blind one, so "the model does not know about his changeup gap" is not the whole story for him
# even though it is a real story leaguewide.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
options(width = 205)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")
S <- readRDS(file.path(MDIR, "velo_sep_seasons.rds")); setDT(S)

## ---- panel A -------------------------------------------------------------------------------------
S[, dec := cut(sep, quantile(sep, 0:10/10), labels = FALSE, include.lowest = TRUE)]
D <- S[, .(sep = mean(sep), n = .N,
           blind = weighted.mean(wh_blind, nsw), aware = weighted.mean(wh_aware, nsw),
           se_b = sd(wh_blind)/sqrt(.N), se_a = sd(wh_aware)/sqrt(.N)), by = dec]
M <- rbind(D[, .(sep, dec, est = blind, se = se_b, model = "Fastball-blind (a stuff model)")],
           D[, .(sep, dec, est = aware, se = se_a, model = "Separation priced in")])

pA <- ggplot(M, aes(sep, est, colour = model, fill = model)) +
  geom_hline(yintercept = 0, colour = "grey35", linewidth = .4) +
  geom_ribbon(aes(ymin = est - se, ymax = est + se), alpha = .16, colour = NA) +
  geom_line(linewidth = .95) + geom_point(size = 2.1) +
  scale_colour_manual(values = c("Fastball-blind (a stuff model)" = "#c0392b",
                                 "Separation priced in" = "#1d7870"), name = NULL) +
  scale_fill_manual(values = c("Fastball-blind (a stuff model)" = "#c0392b",
                               "Separation priced in" = "#1d7870"), name = NULL) +
  annotate("text", x = 4.6, y = -2.5, hjust = 0, size = 3.1, colour = "#c0392b", fontface = "bold",
           label = "graded too high") +
  annotate("text", x = 11.4, y = 2.9, hjust = 1, size = 3.1, colour = "#c0392b", fontface = "bold",
           label = "graded too low") +
  labs(title = "A stuff model that never sees the fastball misgrades changeups by six points of whiff, end to end",
       subtitle = paste0("Pitcher-seasons with 40+ changeup swings, 2020-2026, in velocity-separation deciles. The y axis is actual whiff rate minus the model's out-of-fold grade,\n",
                         "so zero is a calibrated model. Both models see the changeup's velocity, spin, movement, extension and release point; only the green one also sees how\n",
                         "far off the fastball it is. Pricing in separation absorbs 77 percent of the slope. Shaded band is one standard error."),
       x = "Velocity separation from the four-seamer (mph)", y = "Whiff% above the model's grade") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12.5), plot.subtitle = element_text(size = 8.2),
        panel.grid.minor = element_blank(), legend.position = "top")
ggsave(file.path(AST, "fig20_velo_sep_calibration.png"), pA, width = 12, height = 6.4, dpi = 150)

## ---- panel B -------------------------------------------------------------------------------------
B <- S[nsw >= 150]
# Labelling every extreme buries the point. This names Martinez, the handful of seasons at each
# vertical extreme, and the widest separators, and leaves the rest of the cloud unlabelled.
# Keyed on the row, not the printed name: three different Martinezes throw changeups here and
# matching on "Martinez '22" labels the wrong dot.
B[, `:=`(row = .I, nm = paste0(sub(",.*","",player_name), " '", substr(season,3,4)))]
B[grepl("Martinez, Nick", player_name), nm := paste0("N. ", nm)]
key <- unique(c(B[grepl("Martinez, Nick", player_name)]$row,
                B[order(-wh_blind)][1:6]$row, B[order(wh_blind)][1:4]$row,
                B[order(-sep)][1:5]$row, B[order(sep)][1:2]$row))
B[, lab := fifelse(row %in% key, nm, NA_character_)]
pB <- ggplot(B, aes(sep, wh_blind)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = .35) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "#c0392b", fill = "#c0392b",
              alpha = .12, linewidth = .8) +
  geom_point(aes(size = nsw, colour = grepl("Martinez, Nick", player_name)), alpha = .75) +
  geom_text_repel(aes(label = lab), size = 2.7, max.overlaps = Inf, seed = 4, box.padding = .4,
                  min.segment.length = 0, segment.size = .25, segment.colour = "grey60",
                  colour = "grey15") +
  scale_colour_manual(values = c(`TRUE` = "#c0392b", `FALSE` = "#33608f"), guide = "none") +
  scale_size_continuous(range = c(1, 4.2), guide = "none") +
  labs(title = "Martinez is the archetype, but separation explains only a fifth of his edge",
       subtitle = paste0("Pitcher-seasons with 150+ changeup swings. Martinez in red. He sits at the far right of the separation distribution and well above the fit line, and he stays\n",
                         "roughly as far above a model that DOES price in separation - so the missing fastball context is a real leaguewide bias but not a complete account of him."),
       x = "Velocity separation from the four-seamer (mph)",
       y = "Whiff% above the fastball-blind grade") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12.5), plot.subtitle = element_text(size = 8.2),
        panel.grid.minor = element_blank())
ggsave(file.path(AST, "fig21_martinez.png"), pB, width = 12, height = 6.6, dpi = 150)

## ---- Martinez arithmetic ----------------------------------------------------------------------------
NM <- S[grepl("Martinez, Nick", player_name)]
lg <- mean(S$sep); sl <- coef(lm(wh_blind ~ sep, S))["sep"]
cat("=== Nick Martinez, five qualifying seasons ===\n")
cat(sprintf("  mean separation %.1f mph vs a league mean of %.1f, so %+.1f mph above average\n",
            mean(NM$sep), lg, mean(NM$sep) - lg))
cat(sprintf("  the leaguewide slope predicts %+.2f whiff points of underrating from separation alone\n",
            (mean(NM$sep) - lg) * sl))
cat(sprintf("  he actually beats the blind grade by %+.2f, and the separation-aware grade by %+.2f\n",
            mean(NM$wh_blind), mean(NM$wh_aware)))
cat(sprintf("  so pricing in separation recovers %.2f of his %.2f points, or %.0f%%\n",
            mean(NM$wh_blind) - mean(NM$wh_aware), mean(NM$wh_blind),
            100*(mean(NM$wh_blind) - mean(NM$wh_aware))/mean(NM$wh_blind)))
cat(sprintf("  run value: %+.2f above the blind grade per 100, %+.2f above the aware grade\n",
            mean(NM$rv_blind), mean(NM$rv_aware)))

cat("\n=== who else is systematically underrated by a fastball-blind model? ===\n")
cat("   arms with 3+ qualifying seasons, ranked by mean whiff above the blind grade\n\n")
A <- S[, .(seasons = .N, swings = sum(nsw), sep = round(mean(sep),1),
           blind = round(mean(wh_blind),2), aware = round(mean(wh_aware),2),
           from_sep = round(mean(wh_blind) - mean(wh_aware),2)),
       by = .(player_name = sub(",.*","",player_name), id)][seasons >= 3]
print(head(A[order(-blind)][, !"id"], 12), row.names = FALSE)
cat("\n   and the arms a blind model most OVERRATES:\n\n")
print(head(A[order(blind)][, !"id"], 8), row.names = FALSE)
fwrite(A[order(-blind)], file.path(AST, "ext_velo_sep_arms.csv"))
cat(sprintf("\nwrote fig20_velo_sep_calibration.png, fig21_martinez.png, ext_velo_sep_arms.csv\n"))
