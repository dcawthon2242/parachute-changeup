#!/usr/bin/env Rscript

# THE CLOCK FACE IS NOT THE BREAK
#
# spin_axis survives every test for being a real Hawk-Eye measurement of the ball's
# rotation rather than the movement vector relabelled (see spin_axis_provenance.R). What
# it reports is the 2D projection of that rotation onto the plane facing the hitter - the
# clock face, the thing a hitter actually reads out of the hand.
#
# Gyro spin points along the direction of travel. It is invisible on the clock face and it
# is the entire reason the ball does not move. So a pitch can present the fastball's clock
# and still break nothing like it. This draws both arrows to show the split.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(grid) })
options(width = 200)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")

COLS <- c("pitch_type","player_name","pitcher","p_throws","release_speed","pfx_x","pfx_z",
          "spin_axis","game_year")
d <- rbindlist(lapply(c("data/statcast_2025/statcast_2025_all.csv",
                        "data/statcast_2026/statcast_2026_all.csv"),
                      function(f) fread(f, select = COLS, showProgress = FALSE)))
d <- d[pitch_type %in% c("FF","SI","SL") & is.finite(spin_axis) & is.finite(pfx_x)]
d[, `:=`(hb = fifelse(p_throws=="R",-1,1)*pfx_x*12, ivb = pfx_z*12)]

a <- d[, .(n=.N, velo=mean(release_speed), hb=mean(hb), ivb=mean(ivb),
           sx=mean(sin(spin_axis*pi/180)), cx=mean(cos(spin_axis*pi/180)),
           name=player_name[1], hand=p_throws[1]),
       by=.(pitcher, game_year, pitch_type)][n >= 200]
a[, axis := (atan2(sx,cx)*180/pi) %% 360]
# The clock face as a direction in the same plane as the movement arrows. spin_axis 180 is
# pure backspin, i.e. straight up, so the implied direction is axis - 90 (then mirrored).
a[, th := (axis - 90)*pi/180]
a[, `:=`(ux = fifelse(hand=="R",-1,1)*cos(th), uz = sin(th))]

fb <- a[pitch_type %in% c("FF","SI")][order(pitcher,game_year,-n)][, .SD[1], by=.(pitcher,game_year)]
sl <- a[pitch_type == "SL"]
m  <- merge(sl, fb[, .(pitcher, game_year, fb_axis=axis, fb_hb=hb, fb_ivb=ivb,
                       fb_ux=ux, fb_uz=uz, fb_velo=velo)], by=c("pitcher","game_year"))
m[, gap := abs(axis-fb_axis)][gap > 180, gap := 360-gap]

pick <- function(who, yr) m[grepl(who, name) & game_year == yr][1]
ex <- rbindlist(list(pick("Ray, Robbie",2025), pick("Lee, Dylan",2025),
                     pick("Kershaw",2025), pick("Miller, Mason",2025)), fill = TRUE)
ex <- ex[!is.na(pitcher)]
ex[, who := factor(sprintf("%s  \u2014  %.0f\u00b0 axis gap", sub(",.*","",name), gap),
                   levels = sprintf("%s  \u2014  %.0f\u00b0 axis gap",
                                    sub(",.*","",name[order(gap)]), gap[order(gap)]))]

R <- 14  # clock arrows are direction only, so give them all the same length
rows <- rbind(
  ex[, .(who, row = "What the hitter reads out of the hand\n(measured spin axis, direction only)",
         lab = "fastball", x = R*fb_ux, y = R*fb_uz)],
  ex[, .(who, row = "What the hitter reads out of the hand\n(measured spin axis, direction only)",
         lab = "slider",   x = R*ux,    y = R*uz)],
  ex[, .(who, row = "What the ball actually does\n(movement, inches)",
         lab = "fastball", x = fb_hb,   y = fb_ivb)],
  ex[, .(who, row = "What the ball actually does\n(movement, inches)",
         lab = "slider",   x = hb,      y = ivb)])
rows[, row := factor(row, levels = unique(row))]

g <- ggplot(rows, aes(x, y)) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  # Fastball drawn wide and translucent underneath so a slider sitting on top of it still
  # reads as two arrows rather than one.
  geom_segment(data = rows[lab == "fastball"], aes(x = 0, y = 0, xend = x, yend = y,
               colour = lab), arrow = arrow(length = unit(11,"pt"), type = "closed"),
               linewidth = 3.2, alpha = .45) +
  geom_segment(data = rows[lab == "slider"], aes(x = 0, y = 0, xend = x, yend = y,
               colour = lab), arrow = arrow(length = unit(6,"pt"), type = "closed"),
               linewidth = 1.1) +
  scale_colour_manual(values = c(fastball = "#c0392b", slider = "#1a5fb4"), name = NULL) +
  facet_grid(row ~ who, switch = "y") +
  coord_fixed(xlim = c(-20, 20), ylim = c(-20, 20)) +
  labs(title = "A gyro slider can show the fastball's clock face and still break nothing like it",
       subtitle = paste0("Top row is Hawk-Eye's measured spin axis, drawn as a direction only - this is the orientation of the ball's rotation, the cue a hitter reads. Bottom row is\n",
                         "the same two pitches' actual movement. For Ray, Lee and Kershaw the top arrows sit nearly on top of each other while the bottom arrows do not: the slider\n",
                         "spins the way the fastball spins but most of that spin is gyro, pointed at the hitter, so it produces no break. Miller is the conventional case, where the\n",
                         "slider announces itself by spinning the opposite way. Lefties mirrored; + horizontal is arm side."),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face="bold"), plot.subtitle = element_text(size = 8.6),
        strip.text.x = element_text(face="bold"), strip.text.y.left = element_text(angle = 90, size = 8.5),
        panel.grid.minor = element_blank(), legend.position = "top",
        axis.text = element_text(size = 7.5))
ggsave(file.path(AST, "fig12d_clock_vs_break.png"), g, width = 12, height = 7.6, dpi = 150)

cat("=== THE FOUR EXAMPLES ===\n")
print(ex[order(gap), .(pitcher = sub(",.*","",name), yr = game_year, axis_gap = round(gap),
                       FB_axis = round(fb_axis), SL_axis = round(axis),
                       FB_move = sprintf("%+.0f/%+.0f", fb_hb, fb_ivb),
                       SL_move = sprintf("%+.0f/%+.0f", hb, ivb),
                       velo_gap = round(fb_velo - velo, 1))], row.names = FALSE)
cat("\nwrote fig12d_clock_vs_break.png\n")
