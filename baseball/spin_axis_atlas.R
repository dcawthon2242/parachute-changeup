#!/usr/bin/env Rscript

# AN ATLAS OF SPIN-AXIS RELATIONSHIPS ACROSS THE ARSENAL
#
# Everything so far has measured a secondary pitch against the primary fastball. This
# looks at every pitch-type PAIR a pitcher throws, so the geometry is visible at once:
# which pairs spin alike, which are mirrors, and which sit in the gyro middle.
#
# Spin axis is circular, so all averaging is circular (atan2 of mean sin/cos) and every
# pairwise gap is folded to 0-180: 0 = identical spin direction, 180 = perfect mirror,
# ~90 = orthogonal, which in practice means one of the two is gyro-dominant.
#
# Handedness is mirrored (LHP axes reflected) so arm side means the same thing in both.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(grid) })
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")

TYPES <- c("FF","SI","FC","SL","ST","SV","CU","KC","CH","FS")
LAB <- c(FF="4-Seam", SI="Sinker", FC="Cutter", SL="Slider", ST="Sweeper", SV="Slurve",
         CU="Curveball", KC="Knuckle-curve", CH="Changeup", FS="Splitter")
MIN_PT <- 50     # pitches of a type, in a season, for that pitcher to count
MIN_PAIR <- 40   # pitcher-seasons a pair needs before it gets reported

d <- readRDS(file.path(MDIR, "miss_grade_data_activespin.rds"))
d <- d[pitch_type %in% TYPES & is.finite(spin_axis)]
# p_throws was dropped upstream; release_pos_x recovers it cleanly (positive = LHP).
d[, throws := fifelse(mean(release_pos_x, na.rm = TRUE) < 0, "R", "L"), by = pitcher]
# Mirror LHP so "arm side" is the same direction for everyone. Statcast measures the
# axis clockwise from 0, so reflecting a lefty is 360 - axis.
d[, ax_m := fifelse(throws == "R", spin_axis, (360 - spin_axis) %% 360)]
cat("handedness check - circular mean 4-seam axis before mirroring:\n")
print(d[pitch_type == "FF", .(mean_axis = round((atan2(mean(sin(spin_axis*pi/180)),
        mean(cos(spin_axis*pi/180)))*180/pi) %% 360)), by = throws])
d[, `:=`(sx = sin(ax_m*pi/180), cx = cos(ax_m*pi/180))]

pt <- d[, .(n = .N, sx = mean(sx), cx = mean(cx), active = mean(active_spin, na.rm = TRUE)),
        by = .(pitcher, season, pitch_type)][n >= MIN_PT]
pt[, axis := (atan2(sx, cx)*180/pi) %% 360]

## ---- every within-arsenal pair ------------------------------------------------
pr <- merge(pt[, .(pitcher, season, a = pitch_type, ax_a = axis, act_a = active, n_a = n)],
            pt[, .(pitcher, season, b = pitch_type, ax_b = axis, act_b = active, n_b = n)],
            by = c("pitcher","season"), allow.cartesian = TRUE)
pr <- pr[match(a, TYPES) < match(b, TYPES)]
pr[, gap := abs(ax_a - ax_b)][gap > 180, gap := 360 - gap]

M <- pr[, .(pairs = .N, med_gap = median(gap), q25 = quantile(gap,.25), q75 = quantile(gap,.75),
            pct_under_30 = 100*mean(gap < 30), pct_over_150 = 100*mean(gap > 150)),
        by = .(a, b)][pairs >= MIN_PAIR]
setorder(M, med_gap)

cat("=== MOST SPIN-ALIKE PITCH-TYPE PAIRS (median axis gap, deg) ===\n")
print(head(M[, .(pair = paste(LAB[a], "+", LAB[b]), pairs, med_gap = round(med_gap,1),
                 IQR = sprintf("%.0f-%.0f", q25, q75),
                 pct_within_30deg = round(pct_under_30))], 12))
cat("\n=== MOST MIRRORED PAIRS ===\n")
print(head(M[order(-med_gap), .(pair = paste(LAB[a], "+", LAB[b]), pairs,
                 med_gap = round(med_gap,1), IQR = sprintf("%.0f-%.0f", q25, q75),
                 pct_beyond_150deg = round(pct_over_150))], 12))
fwrite(M, file.path(AST, "ext_spin_axis_pair_atlas.csv"))

## ---- the matrix ---------------------------------------------------------------
FULL <- rbind(M[, .(a, b, med_gap, pairs)], M[, .(a = b, b = a, med_gap, pairs)])
ord <- TYPES[TYPES %in% unique(c(FULL$a, FULL$b))]
FULL[, `:=`(a = factor(LAB[a], levels = LAB[ord]), b = factor(LAB[b], levels = rev(LAB[ord])))]

theme_set(theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face="bold"), panel.grid = element_blank(),
        plot.subtitle = element_text(size = 9.5)))

g1 <- ggplot(FULL, aes(a, b, fill = med_gap)) +
  geom_tile(colour = "white", linewidth = .8) +
  geom_text(aes(label = sprintf("%.0f\u00b0", med_gap),
                colour = med_gap > 60 & med_gap < 140), size = 3.5, fontface = "bold") +
  scale_colour_manual(values = c("TRUE"="grey20","FALSE"="white"), guide = "none") +
  scale_fill_gradientn(colours = c("#1a5fb4","#7fb3e0","#f2f2f2","#e8a33d","#c0392b"),
                       values = scales::rescale(c(0, 45, 90, 140, 180)), limits = c(0,180),
                       name = "median axis\ngap (deg)") +
  coord_fixed() +
  labs(title = "How every pair of pitches in an arsenal relates by spin axis",
       subtitle = paste0("Median circular gap between the two pitches' spin axes, over pitcher-seasons where both were thrown 50+ times (2023H2-2026).\n",
                         "Blue = the two spin alike. Red = mirror images. Grey middle (~90\u00b0) = orthogonal, which in practice means one pitch is gyro-dominant.\n",
                         "Handedness mirrored so arm side reads the same for lefties and righties."),
       x = NULL, y = NULL)
ggsave(file.path(AST, "fig12_spin_axis_matrix.png"), g1, width = 11.5, height = 9.4, dpi = 150)

## ---- distributions vs the primary fastball ------------------------------------
fbp <- pt[pitch_type %in% c("FF","SI","FC")][, pr := match(pitch_type, c("FF","SI","FC"))][
  order(pitcher, season, pr)][, .SD[1], by = .(pitcher, season)][, .(pitcher, season, fb_axis = axis,
                                                                     fb_type = pitch_type)]
sec <- merge(pt[!pitch_type %in% c("FF")], fbp, by = c("pitcher","season"))
sec <- sec[pitch_type != fb_type]
sec[, gap := abs(axis - fb_axis)][gap > 180, gap := 360 - gap]
sec[, lab := LAB[pitch_type]]
sec <- sec[, if (.N >= 40) .SD, by = lab]

med <- sec[, .(m = median(gap), n = .N, sprd = IQR(gap)), by = lab][order(m)]
med[, ttl := sprintf("%s   \u2014 spread %.0f\u00b0", lab, sprd)]
sec <- merge(sec, med[, .(lab, m, ttl)], by = "lab")
sec[, ttl := factor(ttl, levels = med$ttl)]
med[, ttl := factor(ttl, levels = med$ttl)]

g2 <- ggplot(sec, aes(gap)) +
  geom_histogram(aes(y = after_stat(density), fill = ttl), bins = 36, colour = NA,
                 alpha = .85, show.legend = FALSE) +
  geom_vline(data = med, aes(xintercept = m), colour = "grey15", linewidth = .5) +
  geom_label(data = med, aes(x = fifelse(m > 90, 178, 2), y = Inf,
                             label = sprintf("median %.0f\u00b0\nn=%d", m, n),
                             hjust = fifelse(m > 90, 1, 0)),
             vjust = 1.05, size = 2.9, colour = "grey15", label.size = 0,
             fill = alpha("white", .75), lineheight = .95) +
  facet_wrap(~ ttl, ncol = 4, scales = "free_y") +
  scale_x_continuous(breaks = c(0, 90, 180), limits = c(-4, 184),
                     labels = c("0\u00b0\nidentical","90\u00b0\northogonal","180\u00b0\nmirror")) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Every secondary pitch's spin axis, measured against its own primary fastball",
       subtitle = paste0("One point per pitcher-season, both pitches thrown 50+ times, 2023H2-2026. Panels run from spin-alike to mirrored; \"spread\" is the interquartile width.\n",
                         "Almost every pitch type has its axis fixed by the label itself - a changeup is always ~20\u00b0 off the fastball, a curveball is always ~170\u00b0 off. The slider\n",
                         "is the lone exception, spanning the full range, and that is exactly why axis similarity separates good sliders from bad ones and nothing else."),
       x = "spin-axis gap from the primary fastball", y = NULL) +
  theme(axis.text.y = element_blank(), strip.text = element_text(face = "bold", size = 10),
        panel.spacing.x = unit(14, "pt"))
ggsave(file.path(AST, "fig12b_axis_gap_distributions.png"), g2, width = 13, height = 6.6, dpi = 150)

cat("\n=== SPREAD OF EACH SECONDARY VS ITS FASTBALL ===\n")
print(sec[, .(pitcher_seasons = .N, median_gap = round(median(gap)),
              iqr = sprintf("%.0f-%.0f", quantile(gap,.25), quantile(gap,.75)),
              spread_iqr_deg = round(IQR(gap))), by = lab][order(median_gap)])
cat("\nwrote fig12_spin_axis_matrix.png and fig12b_axis_gap_distributions.png\n")
