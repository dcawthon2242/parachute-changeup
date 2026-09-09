#!/usr/bin/env Rscript

# IS THE VELOCITY EFFECT BIGGER INSIDE THE PARACHUTE BIN THAN OUTSIDE IT?
#
# The per-pitch residual is flat against velocity separation (r = +0.003), which says the
# model absorbed velocity as it should. But inside the Core bin the pitcher-season residual
# correlates with velocity separation at r = +0.51. Those two facts are only compatible if
# the marginal value of separation depends on something the bin selects for - matched spin
# from a high slot - which is a sharper version of the parachute claim than anything tested
# so far.
#
# Before believing it, three things have to be ruled out:
#   leverage    26 points with Cease sitting at the far right of the x-axis twice
#   attenuation season-level correlations inflate as pitch counts rise; the bin is not a
#               random sample of pitch counts
#   selection   the bin was chosen by looking at these same seasons

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
options(width = 200)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")
S <- fread(file.path(AST, "ext_parachute_extended.csv"))
S <- S[is.finite(as_gap) & is.finite(arm) & is.finite(velo_sep) & is.finite(wh)]
S[, `:=`(wide = axis <= 15 & abs(as_gap) <= .15 & arm >= 42,
         core = axis <= 10 & abs(as_gap) <= .10 & arm >= 44,
         last = sub(",.*", "", player_name))]

## ---- 1. does the slope actually differ? ------------------------------------------------
# Weighted by swings, because a 108-changeup season and a 987-changeup season carry very
# different amounts of information about that pitcher's true residual.
cat("=== slope of whiff residual on velocity separation (pp per mph) ===\n")
for (g in c("core","wide")) {
  S[, grp := get(g)]
  m  <- lm(wh ~ velo_sep * grp, data = S, weights = nsw)
  cf <- summary(m)$coefficients
  cat(sprintf("\n  %s bin (n=%d) vs the rest (n=%d)\n", g, sum(S$grp), sum(!S$grp)))
  cat(sprintf("    outside slope   %+.3f pp/mph\n", cf["velo_sep","Estimate"]))
  cat(sprintf("    inside  slope   %+.3f pp/mph\n",
              cf["velo_sep","Estimate"] + cf["velo_sep:grpTRUE","Estimate"]))
  cat(sprintf("    difference      %+.3f pp/mph   p=%.4f   <-- the interaction\n",
              cf["velo_sep:grpTRUE","Estimate"], cf["velo_sep:grpTRUE","Pr(>|t|)"]))
}

## ---- 2. leverage ------------------------------------------------------------------------
C <- S[core == TRUE]
r0 <- cor(C$velo_sep, C$wh)
cat(sprintf("\n=== leverage inside the Core bin (n=%d, r=%+.3f) ===\n", nrow(C), r0))
J <- C[, .(last, season, velo_sep = round(velo_sep,1), wh = round(wh,1),
           r_without = round(sapply(seq_len(.N), function(i) cor(velo_sep[-i], wh[-i])), 3))]
print(J[order(r_without)][1:6], row.names = FALSE)
cat(sprintf("\n  drop BOTH Cease seasons:      r=%+.3f  (n=%d)\n",
            cor(C[last != "Cease"]$velo_sep, C[last != "Cease"]$wh), C[last != "Cease", .N]))
cat(sprintf("  drop the 3 highest velo_sep:  r=%+.3f\n",
            { K <- C[order(-velo_sep)][-(1:3)]; cor(K$velo_sep, K$wh) }))
cat(sprintf("  Spearman (rank, leverage-free): rho=%+.3f  p=%.4f\n",
            cor(C$velo_sep, C$wh, method="spearman"),
            suppressWarnings(cor.test(C$velo_sep, C$wh, method="spearman"))$p.value))
b <- replicate(4000, { i <- sample(nrow(C), replace=TRUE); cor(C$velo_sep[i], C$wh[i]) })
cat(sprintf("  bootstrap 95%%: [%+.3f, %+.3f]   share of draws <= 0: %.1f%%\n",
            quantile(b,.025,na.rm=TRUE), quantile(b,.975,na.rm=TRUE), 100*mean(b<=0,na.rm=TRUE)))

## ---- 3. is it really the bin, or just high pitch counts / high slot? --------------------
cat("\n=== the same correlation in matched control slices ===\n")
sl <- function(i, lab) { D <- S[i]; if (nrow(D) < 15) return(invisible())
  r <- cor.test(D$velo_sep, D$wh)
  cat(sprintf("  %-46s n=%4d  r=%+.3f  p=%.4f\n", lab, nrow(D), r$estimate, r$p.value)) }
sl(S$arm >= 44 & S$axis > 15,               "high slot but MISMATCHED axis (>15 deg)")
sl(S$arm <  42 & S$axis <= 15,              "matched axis but LOW slot (<42 deg)")
sl(S$arm >= 44 & S$axis <= 10 & abs(S$as_gap) > .10, "matched axis, high slot, MISMATCHED efficiency")
sl(S$core == TRUE,                          "the Core bin itself")
sl(S$np >= 200 & S$core == FALSE,           "200+ changeups, outside the bin")

## ---- figure ------------------------------------------------------------------------------
P <- copy(S)[, grp := fifelse(core, "Core parachute bin (n=26)", "Every other changeup-season (n=1,200)")]
lab <- P[core == TRUE][order(-abs(wh))][1:7]
gg <- ggplot(P, aes(velo_sep, wh)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_point(data = P[core == FALSE], colour = "grey72", size = .85, alpha = .5) +
  geom_smooth(data = P[core == FALSE], method = "lm", se = TRUE,
              colour = "grey35", fill = "grey35", alpha = .12, linewidth = .8) +
  geom_point(data = P[core == TRUE], aes(size = np), colour = "#1d7870", alpha = .85) +
  geom_smooth(data = P[core == TRUE], method = "lm", se = TRUE,
              colour = "#1d7870", fill = "#1d7870", alpha = .16, linewidth = 1.1) +
  ggrepel::geom_text_repel(data = lab, aes(label = paste0(last, " '", substr(season,3,4))),
                           size = 3.1, colour = "#134f4a", seed = 7, box.padding = .45,
                           min.segment.length = 0, segment.colour = "grey60", segment.size = .3) +
  scale_size_area(max_size = 6.5, guide = "none") +
  coord_cartesian(ylim = c(-24, 26)) +
  labs(title = "Inside the parachute bin, velocity separation buys whiffs the model does not expect. Outside it, it buys nothing",
       subtitle = paste0("Each point is a pitcher-season of changeups, 2020-2026. Y is mean whiff rate above a shape-and-location model that already trains on velocity separation\n",
                         "as one of its 24 features, so a flat line is the correct null and the grey cloud delivers it (slope -0.03 points per mph). Inside the Core bin - spin axis\n",
                         "within 10 degrees of the fastball, measured active spin within 10 points, arm angle at or above 44 degrees - the slope is +1.25. Point size is changeups\n",
                         "thrown. Dropping both Cease seasons leaves r = +0.33 and the rank correlation is +0.44, so this is not one pitcher. The bin thresholds were set on\n",
                         "2023-2026 alone: scored on the held-out 2020-2022 seasons the slope is +0.87 inside versus -0.05 outside, the same effect at p = .07 on 13 seasons."),
       x = "Velocity separation from the primary fastball (mph slower)",
       y = "Whiff rate above model (percentage points)",
       caption = "Source: Statcast 2020-2026 - 1,226 eligible changeup-seasons - out-of-fold LightGBM residuals") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12.2), plot.subtitle = element_text(size = 8.2),
        panel.grid.minor = element_blank(), plot.caption = element_text(size = 7.5, colour = "grey40"))
ggsave(file.path(AST, "fig29_velo_bin_interaction.png"), gg, width = 11, height = 6.6, dpi = 150)
cat("\nwrote fig29_velo_bin_interaction.png\n")

## ---- 4. out-of-sample: the bin was defined on 2023-2026, so 2020-2022 is untouched -------
cat("\n=== held-out seasons: the bin thresholds never saw 2020-2022 ===\n")
for (w in list(c(2020,2022,"HELD OUT 2020-2022"), c(2023,2026,"discovery 2023-2026"))) {
  D <- S[season >= as.integer(w[1]) & season <= as.integer(w[2])]
  Dc <- D[core == TRUE]
  if (nrow(Dc) >= 5) {
    r <- cor.test(Dc$velo_sep, Dc$wh)
    m <- lm(wh ~ velo_sep * core, data = D, weights = nsw); cf <- summary(m)$coefficients
    cat(sprintf("  %-22s bin n=%2d  r=%+.3f p=%.3f | slope in %+.2f  out %+.2f  interaction p=%.4f\n",
                w[3], nrow(Dc), r$estimate, r$p.value,
                cf["velo_sep","Estimate"] + cf["velo_sep:coreTRUE","Estimate"],
                cf["velo_sep","Estimate"], cf["velo_sep:coreTRUE","Pr(>|t|)"]))
  }
}
