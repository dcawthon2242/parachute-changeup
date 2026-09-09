#!/usr/bin/env Rscript

# Velocity separation against whiff residual, inside the searched-optimum bin versus every other
# changeup season.
#
# The correct null here is zero, in both groups. Velocity separation is one of the features the
# residual model trains on, so anything the model has priced correctly leaves no correlation
# behind. A non-zero slope means the model is systematically wrong about separation for that
# group - which is the whole claim being tested.
#
# Residuals come from the arm-aware model, since the bin gates on arm slot.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
set.seed(3); options(width = 205)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")

R  <- readRDS(file.path(MDIR, "whiff_tjstuff.rds"))
AS <- readRDS(file.path(MDIR, "active_spin_long.rds"))
S <- R[, .(nsw = .N, axis = mean(axis_diff), arm = mean(arm_angle), velo_sep = -mean(speed_diff),
           w = 100*mean(r4)), by = .(pitcher, player_name, season)]
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, as_ch = active_spin)], by = c("pitcher","season"))
S <- merge(S, AS[pitch_type == "FF", .(pitcher, season, as_fb = active_spin)], by = c("pitcher","season"))
S <- S[nsw >= 60]
S[, name := trimws(paste(sub(".*,\\s*","",player_name), sub(",.*","",player_name)))]
S[, bin := as_fb >= .85 & as_ch >= .89 & axis <= 8 & arm >= 35]
cat(sprintf("bin %d seasons | rest %d seasons\n\n", sum(S$bin), sum(!S$bin)))

rep <- function(D, lab) { ct <- cor.test(D$velo_sep, D$w)
  f <- lm(w ~ velo_sep, D); s <- summary(f)$coefficients["velo_sep",]
  data.table(group = lab, n = nrow(D), r = ct$estimate, ci = sprintf("[%+.2f, %+.2f]",
             ct$conf.int[1], ct$conf.int[2]), slope = s[1], se = s[2], p = ct$p.value) }
cat("=== correlation and slope of whiff residual on velocity separation ===\n")
print(rbind(rep(S[bin == TRUE],  "searched-optimum bin"),
            rep(S[bin == FALSE], "all other changeups"),
            rep(S,               "everything pooled"))[, .(group, n, r = round(r,3), ci,
            slope_pp_per_mph = round(slope,3), p = round(p,4))], row.names = FALSE)

cat("\n=== mean velocity separation, for context ===\n")
cat(sprintf("  bin  %.1f mph (sd %.1f)   rest  %.1f mph (sd %.1f)\n",
    mean(S[bin==TRUE]$velo_sep), sd(S[bin==TRUE]$velo_sep),
    mean(S[bin==FALSE]$velo_sep), sd(S[bin==FALSE]$velo_sep)))

# Fisher z on the difference of the two correlations, plus the regression interaction, which
# answers the same question with the pooled error term.
z <- function(r) .5*log((1+r)/(1-r))
r1 <- cor(S[bin==TRUE]$velo_sep, S[bin==TRUE]$w); n1 <- sum(S$bin)
r2 <- cor(S[bin==FALSE]$velo_sep, S[bin==FALSE]$w); n2 <- sum(!S$bin)
zs <- (z(r1)-z(r2))/sqrt(1/(n1-3) + 1/(n2-3))
cat(sprintf("\nFisher z on the difference in correlation: z = %.2f, p = %.4f\n", zs, 2*pnorm(-abs(zs))))
for (wt in c(FALSE, TRUE)) { f <- if (wt) lm(w ~ velo_sep*bin, S, weights = nsw) else lm(w ~ velo_sep*bin, S)
  d <- summary(f)$coefficients["velo_sep:binTRUE",]
  cat(sprintf("interaction (%s): %+.2f pp per mph, se %.2f, p = %.4f\n",
      if (wt) "weighted by swings" else "unweighted", d[1], d[2], d[4])) }

# The bin holds twenty seasons and two of them sit far out on the x-axis, so leverage has to be
# checked before the slope means anything.
cat("\n=== leverage: is the bin slope carried by Cease and Skubal? ===\n")
B <- S[bin == TRUE]
for (drop in list(character(0), "Dylan Cease", "Tarik Skubal", c("Dylan Cease","Tarik Skubal"))) {
  D <- B[!name %in% drop]; ct <- cor.test(D$velo_sep, D$w)
  cat(sprintf("  drop %-28s n=%2d  r=%+.3f  slope %+.2f  p=%.4f\n",
      if (length(drop)) paste(drop, collapse=" + ") else "nothing", nrow(D), ct$estimate,
      coef(lm(w ~ velo_sep, D))[2], ct$p.value)) }
cat(sprintf("\n  Spearman rank correlation inside the bin (outlier-insensitive): %+.3f\n",
            cor(B$velo_sep, B$w, method = "spearman")))

## ---- figure ---------------------------------------------------------------------------------
S[, grp := fifelse(bin, sprintf("Searched-optimum bin (n = %d)", n1),
                        sprintf("All other changeups (n = %s)", format(n2, big.mark=",")))]
S[, grp := factor(grp, levels = c(unique(grp[bin]), unique(grp[!bin])))]
lab <- S[bin == TRUE][order(-abs(w))][1:6]
gg <- ggplot(S, aes(velo_sep, w, colour = grp)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_point(data = S[bin == FALSE], aes(size = nsw), alpha = .16) +
  geom_smooth(data = S[bin == FALSE], method = "lm", se = TRUE, colour = "grey35",
              fill = "grey70", linewidth = .7, formula = y ~ x) +
  geom_point(data = S[bin == TRUE], aes(size = nsw), alpha = .9) +
  geom_smooth(data = S[bin == TRUE], method = "lm", se = TRUE, linewidth = .9, formula = y ~ x) +
  geom_text_repel(data = lab, aes(label = sprintf("%s '%02d", name, season %% 100)),
                  size = 2.8, seed = 3, min.segment.length = 0, show.legend = FALSE) +
  scale_colour_manual(values = setNames(c("#1d7870","grey55"), levels(S$grp)), name = NULL) +
  scale_size_continuous(range = c(.6, 5), guide = "none") +
  labs(title = "The split looks decisive - r = +0.58 inside the bin against +0.00 outside - but two seasons are carrying it",
       subtitle = paste0("Each point is a pitcher-season, sized by changeup swings. Zero slope is the correct null in both groups: velocity separation is already one of the 24 features\n",
                         "the residual model trains on, so a flat line means the model has priced it correctly and a tilted line means it has not. The grey cloud delivers that null\n",
                         "almost exactly, at r = +0.002 over 1,066 seasons, which is a useful check that the residual is behaving rather than a result in itself. The twenty seasons\n",
                         "inside the bin tilt up at +1.05 points of whiff per mph, and the gap between the two correlations reaches p = .007 by Fisher z. The leverage check is what\n",
                         "should temper this: Cease and Skubal sit alone at the right edge of the x-axis, and dropping both leaves r = +0.19 at p = .46 on the remaining eighteen.\n",
                         "The rank correlation, which ignores how far out those two sit, is +0.24 against a Pearson of +0.58 - the size of that gap is itself the diagnosis. Both\n",
                         "means are 8.0 mph, so the bin is not simply a group that throws more separation; the claim is that separation pays off differently inside it. That claim\n",
                         "is consistent with everything else in this work and is not established by this figure."),
       x = "Velocity separation from the four-seamer (mph slower)",
       y = "Whiff above model (percentage points)",
       caption = "Source: Statcast 2020-2026 - four-seam anchor - residuals from a shape, location and arm-angle model - 60+ changeup swings per season") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11.5), plot.subtitle = element_text(size = 8.2),
        panel.grid.minor = element_blank(), legend.position = "top",
        plot.caption = element_text(size = 7.5, colour = "grey40"))
ggsave(file.path(AST, "fig38_velo_corr_searched.png"), gg, width = 11.5, height = 7.2, dpi = 150)
cat("\nwrote fig38_velo_corr_searched.png\n")
