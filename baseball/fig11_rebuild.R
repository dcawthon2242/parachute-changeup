#!/usr/bin/env Rscript

# FIG 11, REBUILT AT THE PITCHER LEVEL - the only unit these cues support.
#
# Two corrections to the original.
#
# 1. THE UNIT. The original correlated each cue against per-pitch residuals over ~300k
#    pitches. But path_ratio and arm_diff never vary within a pitcher x pitch-type x
#    season - path_ratio varies in 0 of 5,301 groups - because both are stored as
#    pitcher-season aggregates and broadcast onto every pitch. Correlating 296,774 copies
#    of ~5,300 distinct numbers against per-pitch outcomes inflates significance by an
#    enormous factor while estimating nothing but a badly weighted between-pitcher
#    relationship. Every claim in the article is a between-pitcher claim anyway.
#
# 2. THE RESIDUAL. Three stages shown, because which one you pick decides the answer:
#    shape only; shape plus the published post-hoc cubic location fit; and location with
#    handedness carried inside the model.
#
# All cues oriented so higher = more like the primary fastball.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")
MINP <- 200

d <- rbindlist(lapply(c("breaking","offspeed"),
  function(g) readRDS(file.path(MDIR, sprintf("resid_old_vs_new_%s.rds", g)))))
aa <- fread(file.path(MDIR, "arm_angle_tunnel.csv"))[, .(pitcher, season, pitch_type, arm_diff)]
d <- merge(d, aa, by = c("pitcher","season","pitch_type"), all.x = TRUE)
d[, group := fifelse(grp == "breaking", "Breaking", "Offspeed")]
d[, `:=`(`Arm Angle` = -arm_diff, `Trajectory` = -path_ratio, `Spin Similarity` = spin_sim)]

CUES <- c("Arm Angle","Trajectory","Spin Similarity")
RES  <- c(res_shape = "1. Shape only",
          res_old   = "2. Shape + post-hoc location fit (what Fig 11 published)",
          res_new   = "3. Location & handedness inside the model")

out <- rbindlist(lapply(c("Breaking","Offspeed"), function(g)
  rbindlist(lapply(CUES, function(cu)
    rbindlist(lapply(names(RES), function(rv) {
      B <- d[group == g, .(n = .N, cx = mean(get(cu), na.rm = TRUE),
                           ry = 100*mean(get(rv), na.rm = TRUE)),
             by = .(pitcher, pitch_type)][n >= MINP & is.finite(cx) & is.finite(ry)]
      ct <- suppressWarnings(cor.test(B$cx, B$ry, method = "spearman", exact = FALSE))
      se <- sqrt(1.06/(nrow(B) - 3)); r <- unname(ct$estimate)
      data.table(group = g, cue = cu, resid = RES[[rv]], np = nrow(B), r = r,
                 lo = tanh(atanh(r) - 1.96*se), hi = tanh(atanh(r) + 1.96*se),
                 p = ct$p.value)
    }))))))
out[, `:=`(cue = factor(cue, levels = CUES), resid = factor(resid, levels = RES))]
fwrite(out, file.path(AST, "ext_cue_comparison_rebuilt.csv"))
print(out[, .(group, cue, stage = substr(resid,1,1), np, r = round(r,3),
              ci = sprintf("[%+.2f, %+.2f]", lo, hi), p = signif(p,2))], row.names = FALSE)

GREY <- "#8d99ae"; POS <- "#2a9d8f"
g <- ggplot(out, aes(cue, r, fill = group)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = .4) +
  geom_col(position = position_dodge(width = .72), width = .64, colour = "black",
           linewidth = .25) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = .72),
                width = .15, linewidth = .4, colour = "grey25") +
  geom_text(aes(y = fifelse(r >= 0, hi, lo), label = sprintf("%+.2f", r),
                vjust = fifelse(r >= 0, -1.5, 2.4)),
            position = position_dodge(width = .72), size = 3.3, fontface = "bold") +
  geom_text(aes(y = fifelse(r >= 0, hi, lo), vjust = fifelse(r >= 0, -0.35, 1.1),
                label = fifelse(p < .05, sprintf("p=%.2g", p), "n.s.")),
            position = position_dodge(width = .72), size = 2.6,
            colour = fifelse(out$p < .05, "grey20", "grey55")) +
  scale_fill_manual(values = c(Breaking = GREY, Offspeed = POS), name = NULL) +
  facet_wrap(~ resid, nrow = 1) +
  coord_cartesian(ylim = c(-0.34, 0.34)) +
  labs(title = "At the pitcher level, with location modeled properly, no look-alike cue predicts overperformance",
       subtitle = paste0("Between-pitcher Spearman correlation of each cue with whiff overperformance, one point per pitcher x pitch type, min 200 pitches (500 breaking, 244 offspeed).\n",
                         "The original Figure 11 ran this per pitch, but path_ratio and arm_diff are stored as pitcher-season aggregates and never vary within a pitcher - path_ratio varies\n",
                         "in 0 of 5,301 pitcher x type x season groups - so the per-pitch version inflated significance without adding information. Panels 1 and 2 show why the residual\n",
                         "definition mattered so much: two cues look significant on a shape-only residual and a DIFFERENT pair does on the published post-hoc version, with offspeed\n",
                         "trajectory even flipping sign between them. Panel 3 is the honest answer - across 18 tests, one bar clears p<.05, which is what chance alone produces."),
       x = "Look-alike cue vs primary fastball",
       y = "Between-pitcher correlation with overperformance") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 8.2), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(), legend.position = "top",
        strip.text = element_text(face = "bold", size = 9.5),
        panel.spacing.x = unit(14, "pt"))
ggsave(file.path(AST, "fig11_cue_comparison_rebuilt.png"), g, width = 13.5, height = 6.8, dpi = 150)

## ---- same thing on a significance axis ----------------------------------------
# Correlation size and statistical support tell different stories here, so this plots the
# evidence directly. -log10(p) is used so taller = stronger; p is unsigned, so the
# direction of each effect is printed on the bar.
out[, `:=`(lp = -log10(p), dir = fifelse(r >= 0, "more similar helps", "more DIFFERENT helps"))]
BONF <- -log10(.05/nrow(out))

g2 <- ggplot(out, aes(cue, lp, fill = group)) +
  geom_hline(yintercept = -log10(.05), linetype = "dashed", colour = "#c0392b", linewidth = .45) +
  geom_hline(yintercept = BONF, linetype = "dotted", colour = "#7d3c98", linewidth = .5) +
  geom_col(position = position_dodge(width = .72), width = .64, colour = "black",
           linewidth = .25) +
  geom_text(aes(label = sprintf("p=%.3g", p)), vjust = -1.55,
            position = position_dodge(width = .72), size = 2.9, fontface = "bold") +
  geom_text(aes(label = sprintf("r=%+.2f", r)), vjust = -0.35,
            position = position_dodge(width = .72), size = 2.5, colour = "grey35") +
  annotate("text", x = 0.5, y = -log10(.05) + .09, label = "p = 0.05", hjust = 0,
           size = 2.9, colour = "#c0392b") +
  annotate("text", x = 0.5, y = BONF + .09, label = "p = 0.05 after Bonferroni over all 18 tests",
           hjust = 0, size = 2.9, colour = "#7d3c98") +
  scale_fill_manual(values = c(Breaking = GREY, Offspeed = POS), name = NULL) +
  facet_wrap(~ resid, nrow = 1) +
  coord_cartesian(ylim = c(0, 5.2)) +
  labs(title = "Statistical support for each look-alike cue, and how completely it depends on the residual",
       subtitle = paste0("Between-pitcher tests, one point per pitcher x pitch type, min 200 pitches. Bars are -log10(p), so taller = stronger evidence; the correlation and its sign are printed\n",
                         "on each bar because a p-value has no direction. Read it left to right. The only two bars that clear Bonferroni are breaking trajectory in panels 1 and 2 - the panels\n",
                         "where location is mishandled - and both are NEGATIVE, claiming looser tunnels miss more bats. The significant bars also relocate between panels: shape-only flags\n",
                         "offspeed arm angle, the published residual flags offspeed trajectory instead, which moves from p=0.52 to p=0.0089 while its sign reverses. In panel 3 nothing clears\n",
                         "even the uncorrected line except arm angle at p=0.035, which flips sign out of sample. Evidence that moves whenever the nuisance variable is handled differently is\n",
                         "evidence about the nuisance variable."),
       x = "Look-alike cue vs primary fastball", y = expression(-log[10](p))) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 8.2), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(), legend.position = "top",
        strip.text = element_text(face = "bold", size = 9.5),
        panel.spacing.x = unit(14, "pt"))
ggsave(file.path(AST, "fig11b_cue_significance.png"), g2, width = 13.5, height = 6.8, dpi = 150)

cat("\n=== SIGNIFICANCE BY STAGE ===\n")
print(dcast(out, cue + group ~ resid, value.var = "p")[, lapply(.SD, function(z)
  if (is.numeric(z)) signif(z, 2) else z)], row.names = FALSE)
cat("\nwrote fig11_cue_comparison_rebuilt.png and fig11b_cue_significance.png\n")
