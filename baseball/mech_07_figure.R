#!/usr/bin/env Rscript

# THE DECISION FIGURE.
#
# One panel for the estimate that survived and one for the tests that were meant to explain it. The
# numbers are transcribed from the scripts that produced them rather than recomputed, so that this
# file cannot quietly disagree with the analysis it is drawing.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
MDIR <- file.path("data","statcast_model"); AST <- file.path(MDIR,"article_assets")
dir.create(AST, showWarnings = FALSE, recursive = TRUE)
theme_set(theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(),
        plot.title.position = "plot"))
ACC <- "#1b6ca8"; POS <- "#2a9d8f"; NEG <- "#e76f51"; GREY <- "#9aa0a6"

E <- data.table(
  panel = c(rep("The estimate", 3), rep("Tests that would explain it", 6)),
  label = c("D1 (14 seasons, 13 arms)", "MLB (20 seasons, 14 arms)", "Pooled, one row per arm",
            "M1  four-seam usage dose", "M2  after a four-seamer", "M3  well-tunnelled pitches",
            "F1  spin-rate gap control", "F2  held-out 2024-2026", "F3  within pitcher"),
  est = c(1.15, 3.07, 1.78,  -0.75, 2.34, 0.01, -3.85, -0.31, 2.50),
  se  = c(1.95, 1.38, 0.97,   1.04, 2.62, 1.14,  0.73,  2.94, 2.26),
  kind = c("league","league","pooled", rep("mech",3), rep("falsify",3)))
E[, label := factor(label, levels = rev(label))]
E[, panel := factor(panel, levels = c("The estimate","Tests that would explain it"))]

p <- ggplot(E, aes(est, label, colour = kind)) +
  geom_vline(xintercept = 0, colour = GREY, linewidth = .4) +
  geom_errorbarh(aes(xmin = est - 1.96*se, xmax = est + 1.96*se), height = .22, linewidth = .7) +
  geom_point(aes(size = kind == "pooled")) +
  scale_size_manual(values = c("FALSE" = 2.6, "TRUE" = 4.2), guide = "none") +
  scale_colour_manual(values = c(league = GREY, pooled = ACC, mech = POS, falsify = NEG),
                      guide = "none") +
  facet_grid(panel ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(title = "Four-seam-primary rerun: the pooled estimate holds, the out-of-sample test does not",
       subtitle = paste("Population is pitcher-seasons with four-seam usage at least sinker usage. Bin gates unchanged:",
                        "active spin >= .85 on both pitches, axis gap <= 10 degrees, arm slot in the league top third, 40+ swings.",
                        "The blue pooled point is one row per arm, the number the pre-commit says decides. Bars are 95% intervals.", sep = "\n"),
       x = "Whiff points above model", y = NULL,
       caption = paste("M1 is the interaction per +10 points of four-seam usage. M2 and M3 are interactions, not levels.",
                       "F1 is a negative control and should sit at zero. F2 is the held-out MLB era: the same spec reads +4.83 in 2020-2023.",
                       sep = "\n")) +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 90, face = "bold"),
        plot.caption = element_text(colour = GREY, hjust = 0),
        plot.subtitle = element_text(colour = "grey30", size = 10.5))
ggsave(file.path(AST, "fig_parachute_verdict.png"), p, width = 10.5, height = 7, dpi = 150,
       bg = "white")
fwrite(E, file.path(AST, "ext_parachute_verdict.csv"))
cat(sprintf("wrote %s\n", file.path(AST, "fig_parachute_verdict.png")))
