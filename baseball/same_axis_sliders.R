#!/usr/bin/env Rscript

# WHAT A "SAME-AXIS SLIDER" ACTUALLY LOOKS LIKE
#
# The atlas said sliders are the one secondary whose spin axis is not fixed by the label.
# This pulls the actual pitchers at both ends so the shape is concrete: same-axis sliders
# next to mirror-axis sliders, with movement, velocity and spin efficiency attached.
#
# Key to reading it: Statcast's spin_axis is the direction the ball is DEFLECTED, so a
# same-axis slider is not a pitch that moves like a fastball. It is a pitch whose small
# residual movement points the same way the fastball's large movement points. The break
# is drained out by gyro spin rather than redirected. Efficiency, not direction.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(grid) })
options(width = 200)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")

COLS <- c("pitch_type","player_name","pitcher","p_throws","release_speed","pfx_x","pfx_z",
          "spin_axis","release_spin_rate","description","game_year")
read_yr <- function(f) {
  d <- fread(f, select = COLS, showProgress = FALSE)
  d[pitch_type %in% c("FF","SI","FC","SL","ST") & is.finite(spin_axis) & is.finite(pfx_x)]
}
d <- rbindlist(lapply(c("data/statcast_2025/statcast_2025_all.csv",
                        "data/statcast_2026/statcast_2026_all.csv"), read_yr))

# Mirror lefties into a right-handed frame: axis reflects, horizontal movement flips.
# Raw pfx_x is catcher's view, so arm side is negative for a righty and positive for a
# lefty; both get sent to positive here.
d[, `:=`(ax_m  = fifelse(p_throws == "R", spin_axis, (360 - spin_axis) %% 360),
         hb    = fifelse(p_throws == "R", -1, 1) * pfx_x * 12,   # + = arm side
         ivb   = pfx_z * 12)]
d[, `:=`(swing = description %in% c("swinging_strike","swinging_strike_blocked","foul",
                                    "foul_tip","hit_into_play","foul_bunt","missed_bunt"),
         whiff = description %in% c("swinging_strike","swinging_strike_blocked","foul_tip",
                                    "missed_bunt"))]

agg <- d[, .(n = .N, velo = mean(release_speed), hb = mean(hb), ivb = mean(ivb),
             spin = mean(release_spin_rate, na.rm = TRUE),
             sx = mean(sin(ax_m*pi/180)), cx = mean(cos(ax_m*pi/180)),
             sw = sum(swing), wh = sum(whiff), name = player_name[1], hand = p_throws[1]),
         by = .(pitcher, game_year, pitch_type)]
agg[, `:=`(axis = (atan2(sx, cx)*180/pi) %% 360, whiff_rate = 100*wh/pmax(sw,1))]
agg[, mv := sqrt(hb^2 + ivb^2)]

fb <- agg[pitch_type %in% c("FF","SI") & n >= 200][order(pitcher, game_year, -n)][
  , .SD[1], by = .(pitcher, game_year)]
sl <- agg[pitch_type %in% c("SL","ST") & n >= 200]

m <- merge(sl, fb[, .(pitcher, game_year, fb_type = pitch_type, fb_axis = axis, fb_velo = velo,
                      fb_hb = hb, fb_ivb = ivb, fb_spin = spin, fb_n = n)],
           by = c("pitcher","game_year"))
m[, gap := abs(axis - fb_axis)][gap > 180, gap := 360 - gap]
m[, `:=`(velo_gap = fb_velo - velo, mv_ratio = mv / sqrt(fb_hb^2 + fb_ivb^2))]

# Measured active spin, so "gyro" is a fact rather than an inference.
as_l <- readRDS(file.path(MDIR, "active_spin_long.rds"))
setDT(as_l)
m <- merge(m, as_l[, .(pitcher, season, pitch_type, active_spin)],
           by.x = c("pitcher","game_year","pitch_type"),
           by.y = c("pitcher","season","pitch_type"), all.x = TRUE)

# Overperformance vs the whiff model, which already prices in raw stuff. Two versions,
# because they disagree and the disagreement is the story: wres is the plain residual,
# res_loc strips the systematic location effect out of it first. Fig 11b's axis result
# lives entirely in res_loc - same-axis sliders get thrown to different spots, and the
# plain residual charges them for that.
rs <- readRDS(file.path(MDIR, "oof_whiff_resid.rds"))
setDT(rs)
rs[, res_loc := NA_real_]
for (g in c("breaking","offspeed")) {
  i <- which(rs$grp == g)
  rs$res_loc[i] <- residuals(lm(wres ~ poly(plate_x,3)*poly(plate_z,3) + below_zone + VAA + HAA,
                                data = rs[i]))
}
m <- merge(m, rs[, .(wres = 100*mean(wres), res_loc = 100*mean(res_loc),
                     plate_z = mean(plate_z), below = 100*mean(below_zone)),
                by = .(pitcher, season, pitch_type)],
           by.x = c("pitcher","game_year","pitch_type"),
           by.y = c("pitcher","season","pitch_type"), all.x = TRUE)

fmt <- function(x) x[, .(pitcher = substr(name, 1, 20), yr = game_year, T = hand,
                         typ = pitch_type, n = n,
                         gap = sprintf("%.0f", gap),
                         SL_mv = sprintf("%+.0f/%+.0f", hb, ivb),
                         FB_mv = sprintf("%+.0f/%+.0f", fb_hb, fb_ivb),
                         vsFB = sprintf("%.0f%%", 100*mv_ratio),
                         dvelo = sprintf("%.1f", -velo_gap),
                         actsp = fifelse(is.na(active_spin), "--", sprintf("%.0f%%", 100*active_spin)),
                         whiff = sprintf("%.1f", whiff_rate),
                         raw_vs_exp = fifelse(is.na(wres), "--", sprintf("%+.1f", wres)),
                         loc_adj = fifelse(is.na(res_loc), "--", sprintf("%+.1f", res_loc)))]

setorder(m, gap)
cat("=== SLIDERS THAT SPIN LIKE THE FASTBALL (smallest axis gap, 200+ thrown) ===\n")
print(fmt(head(m, 20)), row.names = FALSE)

cat("\n=== FOR CONTRAST: SLIDERS THAT MIRROR THE FASTBALL (largest axis gap) ===\n")
print(fmt(head(m[order(-gap)], 12)), row.names = FALSE)

cat("\n=== WHAT CHANGES ACROSS THE AXIS-GAP RANGE (true sliders only) ===\n")
m[, band := cut(gap, c(-1, 60, 120, 181),
                labels = c("same axis (<60\u00b0)","gyro middle (60-120\u00b0)","mirror (>120\u00b0)"))]
print(m[pitch_type == "SL", .(pitcher_seasons = .N, mean_gap = round(mean(gap)),
            SL_movement_in = round(mean(mv),1),
            pct_of_FB_movement = round(100*mean(mv_ratio)),
            mean_HB = round(mean(hb),1), mean_IVB = round(mean(ivb),1),
            active_spin = round(100*mean(active_spin, na.rm = TRUE)),
            velo_gap = round(mean(velo_gap),1),
            mean_plate_z = round(mean(plate_z, na.rm = TRUE),2),
            pct_below_zone = round(mean(below, na.rm = TRUE)),
            whiff = round(mean(whiff_rate),1),
            raw_vs_exp = round(mean(wres, na.rm = TRUE),2),
            loc_adj_vs_exp = round(mean(res_loc, na.rm = TRUE),2)), by = band][order(mean_gap)],
      row.names = FALSE)
cat("\nSweepers are excluded above: 91% of them sit past 120 deg, so there is no same-axis\n",
    "sweeper to compare against.\n", sep = "")

fwrite(m, file.path(AST, "ext_same_axis_sliders.csv"))

## ---- picture: movement plot for six arsenals ----------------------------------
pick <- function(who, yr) m[grepl(who, name) & game_year == yr][1]
ex <- rbindlist(list(
  pick("Kershaw", 2025), pick("Cease", 2025), pick("Skubal", 2025),
  m[gap < 25 & n >= 400][order(-whiff_rate)][1],
  m[gap > 150 & n >= 400][order(-whiff_rate)][1],
  m[gap %between% c(80,100) & n >= 400][order(-whiff_rate)][1]
), fill = TRUE)
ex <- unique(ex[!is.na(pitcher)], by = c("pitcher","game_year"))
ex[, ttl := sprintf("%s %d  \u2014  %.0f\u00b0 gap", sub(",.*","",name), game_year, gap)]
ex[, ttl := factor(ttl, levels = ex[order(gap)]$ttl)]

pts <- rbind(ex[, .(ttl, gap, lab = "fastball", hb = fb_hb, ivb = fb_ivb, velo = fb_velo)],
             ex[, .(ttl, gap, lab = pitch_type,  hb, ivb, velo)])

g <- ggplot(pts, aes(hb, ivb)) +
  geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_segment(aes(x = 0, y = 0, xend = hb, yend = ivb, colour = lab == "fastball"),
               arrow = arrow(length = unit(7,"pt"), type = "closed"), linewidth = 1.1,
               show.legend = FALSE) +
  geom_text(aes(label = sprintf("%s %.0f", lab, velo), colour = lab == "fastball",
                hjust = fifelse(hb > 0, -0.10, 1.10)), size = 3.1, fontface = "bold",
            show.legend = FALSE) +
  scale_colour_manual(values = c("TRUE" = "#c0392b", "FALSE" = "#1a5fb4")) +
  facet_wrap(~ ttl, nrow = 2) +
  coord_fixed(xlim = c(-26, 26), ylim = c(-8, 24), clip = "off") +
  labs(title = "A same-axis slider is a fastball with the break drained out, not redirected",
       subtitle = paste0("Arrows are average movement in inches, lefties mirrored so + horizontal is always arm side. The direction an arrow points is what Statcast's\n",
                         "spin_axis measures, so a small gap means the two arrows point the same way - and that can only happen if the slider arrow is short. Pitchers at\n",
                         "a large gap separate the pitches by pointing the arrow the other way; pitchers at a small gap separate them by shortening it to almost nothing."),
       x = "horizontal movement (in, + = arm side)", y = "induced vertical break (in)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 9),
        strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())
ggsave(file.path(AST, "fig12c_same_axis_slider_examples.png"), g, width = 13, height = 6.4, dpi = 150)
cat("\nwrote fig12c_same_axis_slider_examples.png\n")
