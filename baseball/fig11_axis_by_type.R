#!/usr/bin/env Rscript

# FIG 11, REBUILT - spin AXIS similarity to the fastball, split by pitch type.
#
# Replaces the old three-cue chart's spin bar. That bar used spin_sim, a Gaussian kernel
# on active-spin gap and axis gap tuned for offspeed; on breaking balls it saturates into
# 1e-19..1e-6 and carries no interpretable meaning (Kershaw's gyro slider, 46% active
# spin against 85% on his fastball, scored in its top decile). Here the cue is the raw
# spin-axis difference from the primary fastball, oriented as similarity:
#     axis_sim = -axis_diff   ->  positive r = an axis closer to the fastball's helps
#                                 negative r = an axis closer to a MIRROR helps
#
# Measured at the pitcher x pitch-type level, because that is where the cue is
# interpretable, against location-adjusted whiff overperformance. Error bars are a
# pitcher bootstrap; the diamond is the same correlation refit in 2026 alone.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")
set.seed(42)

TYPES <- c("SL","ST","CU","KC","CH","FS")
LAB <- c(SL="Slider", ST="Sweeper", CU="Curveball", KC="Knuckle-\ncurve",
         CH="Changeup", FS="Splitter")
MIN_N <- 200; MIN_26 <- 120

oof <- readRDS(file.path(MDIR, "oof_whiff_resid.rds"))
oof[, res_loc := NA_real_]
for (g in c("breaking","offspeed")) {
  i <- which(oof$grp == g); s <- oof[i]
  oof$res_loc[i] <- residuals(lm(wres ~ poly(plate_x,3)*poly(plate_z,3)+below_zone+VAA+HAA, data=s))
}
d <- readRDS(file.path(MDIR, "miss_grade_data_activespin.rds"))
key <- c("game_pk","at_bat_number","pitch_number","season")
setkeyv(d, key); setkeyv(oof, key)
b <- d[, .(game_pk, at_bat_number, pitch_number, season, axis_diff)][oof]
b <- b[pitch_type %in% TYPES & is.finite(res_loc) & is.finite(axis_diff)]
b[, axis_sim := -axis_diff]

agg <- function(D, minn) D[, .(n = .N, axis_sim = mean(axis_sim), axis_diff = mean(axis_diff),
                               resid_pp = 100*mean(res_loc)),
                           by = .(pitcher, pitch_type)][n >= minn]
ALL <- agg(b, MIN_N); T26 <- agg(b[season == 2026], MIN_26)
TR  <- agg(b[season <= 2025], MIN_N)
J <- merge(TR[, .(pitcher, pitch_type, cue_tr = axis_sim)],
           T26[, .(pitcher, pitch_type, resid26 = resid_pp)], by = c("pitcher","pitch_type"))

S <- rbindlist(lapply(TYPES, function(ty) {
  a <- ALL[pitch_type == ty]; a26 <- T26[pitch_type == ty]; j <- J[pitch_type == ty]
  rs <- replicate(2000, { i <- sample.int(nrow(a), replace = TRUE)
                          suppressWarnings(cor(a$axis_sim[i], a$resid_pp[i], method="spearman")) })
  ct <- cor.test(a$axis_sim, a$resid_pp, method="spearman", exact=FALSE)
  data.table(pitch_type = ty, np = nrow(a),
             r = unname(ct$estimate), p = ct$p.value,
             lo = unname(quantile(rs,.025)), hi = unname(quantile(rs,.975)),
             mean_axis = mean(a$axis_diff),
             r26 = if (nrow(a26) >= 10) cor(a26$axis_sim, a26$resid_pp, method="spearman") else NA_real_,
             np26 = nrow(a26),
             r_pred = if (nrow(j) >= 12) cor(j$cue_tr, j$resid26, method="spearman") else NA_real_,
             np_pred = nrow(j))
}))
S[, verdict := fifelse(lo > 0, "matches the fastball",
                fifelse(hi < 0, "mirrors the fastball", "nothing reliable"))]
S[, lab := factor(LAB[pitch_type], levels = unname(LAB[TYPES]))]
print(S[, .(pitch = LAB[pitch_type], np, mean_axis = round(mean_axis),
            r = round(r,3), ci = sprintf("[%+.2f, %+.2f]", lo, hi), p = signif(p,2),
            r_2026 = round(r26,3), np26, r_predictive = round(r_pred,3), np_pred, verdict)])
fwrite(S, file.path(AST, "ext_axis_similarity_by_pitchtype.csv"))

theme_set(theme_minimal(base_size = 12) +
  theme(plot.title=element_text(face="bold"), panel.grid.major.x=element_blank(),
        plot.subtitle=element_text(size=9.5)))
VC <- c("matches the fastball" = "#2c7fb8", "mirrors the fastball" = "#c0392b",
        "nothing reliable" = "#9aa5b1")
ann <- S[!is.na(r_pred) & lo > 0,
         .(lab, y = 0.70, txt = sprintf("holds out of sample:\n2023-25 cue vs 2026 miss, r=%+.2f", r_pred))]

g <- ggplot(S, aes(lab, r, fill = verdict)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = .4) +
  geom_col(width = .62, colour = "black", linewidth = .25) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .16, linewidth = .45, colour = "grey20") +
  geom_point(aes(y = r26), shape = 23, size = 2.6, fill = "white", colour = "grey15",
             na.rm = TRUE) +
  geom_text(aes(y = ifelse(r >= 0, hi, lo), label = sprintf("%+.2f", r),
                vjust = ifelse(r >= 0, -0.8, 1.8)), size = 3.6, fontface = "bold") +
  geom_text(data = ann, inherit.aes = FALSE, aes(lab, y, label = txt),
            size = 2.8, colour = "#1f5f8b", lineheight = .95) +
  geom_text(aes(y = -0.86, label = paste0("n=", np, " pitchers")), size = 2.9, colour = "grey40") +
  geom_text(aes(y = -0.95, label = paste0(round(mean_axis), "\u00b0 avg gap")), size = 2.9,
            colour = "grey40") +
  scale_fill_manual(values = VC, name = NULL) +
  coord_cartesian(ylim = c(-1.02, 0.82)) +
  labs(title = "Sliders and sweepers want the fastball's spin axis; curveballs show nothing reliable",
       subtitle = paste0("Pitcher x pitch-type Spearman correlation of spin-axis SIMILARITY to the primary fastball with location-adjusted whiff\n",
                         "overperformance (min 200 pitches, 2023H2-2026). Positive = an axis closer to the fastball's goes with beating expectation;\n",
                         "negative = closer to a 180\u00b0 mirror does. Bars are the pooled estimate, whiskers a 2,000-rep pitcher bootstrap, diamonds the\n",
                         "same fit in 2026 alone. Replaces the old spin_sim cue, which had no resolution on breaking balls."),
       x = NULL, y = "correlation with overperformance (pitcher level)") +
  theme(legend.position = "top")
ggsave(file.path(AST, "fig11_spin_axis_by_type.png"), g, width = 12, height = 6.6, dpi = 150)
cat("\nwrote fig11_spin_axis_by_type.png\n")
