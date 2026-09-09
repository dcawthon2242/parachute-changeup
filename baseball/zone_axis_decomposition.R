#!/usr/bin/env Rscript

# DOES HORIZONTAL LOCATION BUY THE SWING AND VERTICAL LOCATION BUY THE MISS?
#
# Fig 3c-2 showed that the high-miss zone's HEIGHT band carries almost all of the miss
# while its SIDE band carries almost none -- but that the height band draws far fewer
# swings. That suggests the two axes do different jobs: side gets the hitter to commit,
# height gets him to miss. This tests it directly.
#
# For each breaking-ball type, three outcomes are regressed on binned plate location:
#   swing decision (all pitches), whiff | swing, miss distance | swing.
# Each outcome is fit on horizontal bins alone, vertical bins alone, and both, so the
# explained variance can be split into what only the horizontal axis knows, what only
# the vertical axis knows, and what they share. If the hypothesis holds, the swing model
# should be horizontal-heavy and the two miss models vertical-heavy.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

MDIR <- file.path("data","statcast_model")
ODIR <- file.path(MDIR, "tunnel_location")
AST  <- file.path(MDIR, "article_assets")

TYPES <- c("SL","CU","ST")
TYPE_LAB <- c(SL="Slider", CU="Curveball (CU+KC)", ST="Sweeper")

# Rulebook zone, league-average edges, used only to define a chase.
SZ_BOT <- 1.6; SZ_TOP <- 3.4; SZ_HALF <- 0.83

p <- readRDS(file.path(MDIR, "tunnel_pairs.rds"))[is_primary_setup == TRUE & brk_type %in% TYPES]
mir <- function(x, h) fifelse(h == "R", x, -x)
p[, `:=`(bb_x = mir(bb_plate_x, p_throws), bb_z = bb_plate_z)]
p <- p[is_swing == FALSE | is.finite(miss_distance)]
p <- p[is.finite(bb_x) & is.finite(bb_z) &
       bb_x >= -2 & bb_x <= 2 & bb_z >= -1 & bb_z <= 4.5]

BW <- 0.25
p[, `:=`(xbin = factor(floor(bb_x/BW)), zbin = factor(floor(bb_z/BW)))]
p[, in_rulebook := abs(bb_x) <= SZ_HALF & bb_z >= SZ_BOT & bb_z <= SZ_TOP]

################################################################################
## variance decomposition: what does each axis know?
################################################################################
# Linear-probability / linear models on binned location. Bins (not splines) so the fit
# is free to be non-monotonic -- swing rate peaks in the middle of the plate on both
# axes, so a linear term would report nothing.
r2 <- function(d, y, rhs) {
  f <- as.formula(paste(y, "~", rhs))
  s <- summary(lm(f, data = d))
  s$r.squared
}

decomp <- rbindlist(lapply(TYPES, function(ty) {
  d  <- p[brk_type == ty]
  ds <- d[is_swing == TRUE]
  dc <- d[in_rulebook == FALSE]          # chase = swing at a pitch out of the zone
  rbindlist(list(
    data.table(outcome = "Swing decision",   n = nrow(d),
               rx = r2(d,"is_swing","xbin"),  rz = r2(d,"is_swing","zbin"),
               rb = r2(d,"is_swing","xbin+zbin")),
    data.table(outcome = "Chase (pitch out of zone)", n = nrow(dc),
               rx = r2(dc,"is_swing","xbin"), rz = r2(dc,"is_swing","zbin"),
               rb = r2(dc,"is_swing","xbin+zbin")),
    data.table(outcome = "Whiff | swing",    n = nrow(ds),
               rx = r2(ds,"is_whiff","xbin"), rz = r2(ds,"is_whiff","zbin"),
               rb = r2(ds,"is_whiff","xbin+zbin")),
    data.table(outcome = "Miss distance | swing", n = nrow(ds),
               rx = r2(ds,"miss_distance","xbin"), rz = r2(ds,"miss_distance","zbin"),
               rb = r2(ds,"miss_distance","xbin+zbin"))
  ))[, brk_type := ty][]
}))
# Unique = what the other axis cannot recover; shared = the part both axes can see,
# which exists because where a pitch finishes horizontally is correlated with its height.
decomp[, `:=`(uniq_x = rb - rz, uniq_z = rb - rx, shared = rx + rz - rb)]
decomp[, vert_share := 100 * uniq_z / (uniq_x + uniq_z)]

cat("=== VARIANCE EXPLAINED BY PLATE LOCATION, SPLIT BY AXIS ===\n")
print(decomp[, .(brk_type, outcome, n,
                 R2_horiz = round(100*rx,2), R2_vert = round(100*rz,2),
                 R2_both  = round(100*rb,2),
                 uniq_horiz = round(100*uniq_x,2), uniq_vert = round(100*uniq_z,2),
                 shared = round(100*shared,2),
                 pct_of_unique_that_is_vertical = round(vert_share))])
fwrite(decomp, file.path(AST, "ext_zone_axis_decomposition.csv"))

################################################################################
## marginal profiles, for the picture
################################################################################
prof <- rbindlist(lapply(TYPES, function(ty) {
  d <- p[brk_type == ty]
  bx <- d[, .(n = .N, swing = mean(is_swing),
              whiff = mean(is_whiff[is_swing]), miss = mean(miss_distance[is_swing])),
          by = .(loc = (as.integer(as.character(xbin)) + .5) * BW)][, axis := "Horizontal (mirrored, + = arm side)"][]
  bz <- d[, .(n = .N, swing = mean(is_swing),
              whiff = mean(is_whiff[is_swing]), miss = mean(miss_distance[is_swing])),
          by = .(loc = (as.integer(as.character(zbin)) + .5) * BW)][, axis := "Vertical (height above plate)"][]
  rbind(bx, bz)[, brk_type := ty][]
}))
prof <- prof[n >= 150]
pl <- melt(prof, id.vars = c("brk_type","axis","loc","n"),
           measure.vars = c("swing","miss"), variable.name = "outcome")
pl[, outcome := factor(fifelse(outcome == "swing", "Swing rate", "Miss distance | swing"),
                       levels = c("Swing rate","Miss distance | swing"))]
pl[, value := fifelse(outcome == "Swing rate", 100*value, value)]
pl[, brk_lab := TYPE_LAB[brk_type]]

theme_set(theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face="bold"), strip.text = element_text(face="bold")))
COLS <- c("Curveball (CU+KC)" = "#c0392b", "Slider" = "#2c3e50", "Sweeper" = "#e08e0b")

g <- ggplot(pl, aes(loc, value, colour = brk_lab)) +
  geom_line(linewidth = 1) + geom_point(size = 1.1) +
  facet_grid(outcome ~ axis, scales = "free") +
  scale_colour_manual(values = COLS, name = NULL) +
  labs(title = "Both axes govern the swing decision; only height governs the miss",
       subtitle = paste0("2-strike breaking balls after the primary fastball, 2023H2-2026, in 3-inch location bins (bins with <150 pitches dropped).\n",
                         "Swing rate traces the same hump on either axis - hitters swing at what is over the plate and at what is belt-high. Miss distance is\n",
                         "the asymmetric one: it climbs steeply as the pitch gets lower, but stays near zero across the whole horizontal range until the far edge."),
       x = "plate location (ft)", y = NULL) +
  theme(legend.position = "top")
ggsave(file.path(AST, "fig3c3_axis_jobs.png"), g, width = 12, height = 6.8, dpi = 150)
ggsave(file.path(ODIR, "fig3c3_axis_jobs.png"), g, width = 12, height = 6.8, dpi = 150)

################################################################################
## the decomposition itself, as bars
################################################################################
db <- melt(decomp, id.vars = c("brk_type","outcome"),
           measure.vars = c("uniq_x","uniq_z","shared"), variable.name = "part")
db[, part := factor(c(uniq_x = "horizontal only", uniq_z = "vertical only",
                      shared = "shared by both")[as.character(part)],
                    levels = c("vertical only","shared by both","horizontal only"))]
db[, value := 100*value]
db[, brk_lab := TYPE_LAB[brk_type]]
db[, outcome := factor(outcome, levels = c("Swing decision","Chase (pitch out of zone)",
                                           "Whiff | swing","Miss distance | swing"))]

hl <- decomp[, .(brk_lab = TYPE_LAB[brk_type], outcome,
                 y = 100*rb + 2.5,
                 txt = sprintf("%.0f%% H", 100*uniq_x/(uniq_x + uniq_z)))]
hl[, outcome := factor(outcome, levels = levels(db$outcome))]

g2 <- ggplot(db, aes(outcome, value, fill = part)) +
  geom_col(width = .68) +
  geom_text(data = hl, inherit.aes = FALSE, aes(outcome, y, label = txt),
            size = 3.2, fontface = "bold", colour = "#2c3e50") +
  facet_wrap(~ brk_lab) +
  scale_fill_manual(values = c("vertical only" = "#c0392b", "shared by both" = "#b9c3cc",
                               "horizontal only" = "#2c3e50"), name = NULL) +
  labs(title = "Height decides the miss; the side of the plate does its work earlier, at the swing decision",
       subtitle = paste0("Share of each outcome's variance explained by 3-inch plate-location bins, split into what only one axis knows and what both do.\n",
                         "Labels give the horizontal axis's share of the two unique pieces. That share falls by half or more moving from the swing decision to\n",
                         "the miss for every pitch type - but height still outweighs side at both stages, except for the sweeper's swing decision.\n",
                         "2-strike breaking balls after the primary fastball, 2023H2-2026."),
       x = NULL, y = "variance explained by location (%)") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 20, hjust = 1),
        panel.grid.major.x = element_blank())
ggsave(file.path(AST, "fig3c4_axis_decomposition.png"), g2, width = 12, height = 6.4, dpi = 150)
ggsave(file.path(ODIR, "fig3c4_axis_decomposition.png"), g2, width = 12, height = 6.4, dpi = 150)

cat("\nwrote fig3c3_axis_jobs.png and fig3c4_axis_decomposition.png\n")
