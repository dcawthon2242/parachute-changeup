#!/usr/bin/env Rscript

# HOW LOCATION-SENSITIVE IS EACH PITCH TYPE'S TUNNEL?
#
# Folds the pitcher's ability to land the breaking ball in the high-miss zone into the
# pitch-characteristic tunneling picture, two ways:
#
#  (1) High-Miss% as a DRIVER alongside the fastball-relative shape gaps, all measured
#      at the pitcher x pitch-type level so the comparison is apples-to-apples.
#      High-Miss% is the share of a pitcher's breaking balls that land in the high-miss
#      zone -- the top 20% of the chase-adjusted miss surface. It is NOT a strike rate;
#      the zone sits mostly below the strike zone.
#
#  (2) the INTERACTION that actually defines location sensitivity: does a tight tunnel
#      pay off more when the pitch is located? Miss is split by tunnel tercile x in/out
#      of the high-miss zone. The gap between in-zone and out-of-zone, and how much that
#      gap widens for tight tunnels, is the location sensitivity of that pitch type.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

MDIR <- file.path("data","statcast_model")
ODIR <- file.path(MDIR, "tunnel_location")
AST  <- file.path(MDIR, "article_assets")

TYPES <- c("SL","CU","ST")
TYPE_LAB <- c(SL="Slider", CU="Curveball (CU+KC)", ST="Sweeper")

GX <- seq(-2.0, 2.0, by = 0.10); GZ <- seq(-0.6, 4.4, by = 0.10)
cell_of <- function(x, z) {
  ix <- round((x - GX[1])/0.1) + 1; iz <- round((z - GZ[1])/0.1) + 1
  out <- rep(NA_integer_, length(x))
  ok <- ix >= 1 & ix <= length(GX) & iz >= 1 & iz <= length(GZ)
  out[ok] <- (iz[ok]-1)*length(GX) + ix[ok]; out
}

p <- readRDS(file.path(MDIR, "tunnel_pairs.rds"))[is_primary_setup == TRUE & brk_type %in% TYPES]
mir <- function(x, h) fifelse(h == "R", x, -x)
p[, `:=`(bb_x = mir(bb_plate_x, p_throws), fb_x = mir(fb_plate_x, p_throws),
         bb_z = bb_plate_z, fb_z = fb_plate_z)]
p[, matchup := fifelse(p_throws == stand, "same", "opp")]
p <- p[is_swing == FALSE | is.finite(miss_distance)]
p[, miss_pp := fifelse(is_swing == TRUE, miss_distance, 0)]

zs <- readRDS(file.path(ODIR, "zone_surfaces.rds"))
p[, cell := cell_of(bb_x, bb_z)]
p <- merge(p, zs[, .(brk_type, matchup, cell = cell_of(gx, gz), inzone_pp)],
           by = c("brk_type","matchup","cell"), all.x = TRUE)
p[is.na(inzone_pp), inzone_pp := FALSE]

# Split hitting the zone into its two axes: getting the HEIGHT right and getting the SIDE
# right, each ignoring the other. Bands are the 5th-95th percentile of the zone's own
# cells, so "height" means the pitch finished in the zone's vertical range at any x.
bands <- zs[inzone_pp == TRUE, .(zlo = quantile(gz, .05), zhi = quantile(gz, .95),
                                 xlo = quantile(gx, .05), xhi = quantile(gx, .95)),
            by = .(brk_type, matchup)]
p <- merge(p, bands, by = c("brk_type","matchup"), all.x = TRUE)
p[, `:=`(inzone_z = bb_z >= zlo & bb_z <= zhi,
         inzone_x = bb_x >= xlo & bb_x <= xhi)]
cat("=== high-miss zone bands (mirrored ft) and how often each axis is hit ===\n")
print(merge(bands, p[, .(hit_both = round(100*mean(inzone_pp)), hit_z = round(100*mean(inzone_z)),
                         hit_x = round(100*mean(inzone_x))), by = .(brk_type, matchup)],
            by = c("brk_type","matchup"))[order(brk_type, matchup)])

################################################################################
## (1) High-Miss% as a driver of the tunnel, at pitcher x type level
################################################################################
md <- readRDS(file.path(MDIR, "miss_grade_data_activespin.rds"))
fbA <- md[grp == "fastball", .(fb_active = mean(active_spin, na.rm = TRUE)), by = .(pitcher, season)]
md <- merge(md, fbA, by = c("pitcher","season"), all.x = TRUE)
keep <- c("game_pk","at_bat_number","pitch_number","season","speed_diff","ax_diff","az_diff",
          "axis_diff","active_spin","fb_active","release_extension","release_pos_x","release_pos_z")
d <- merge(p, md[, ..keep], by = c("game_pk","at_bat_number","pitch_number","season"))

# IMPORTANT: path_ratio = early separation / plate separation, and the high-miss zone
# sits far from the fastball, so landing there mechanically shrinks the denominator
# (cor(path_ratio, plate_sep) = -0.54; cor(High-Miss%, plate_sep) = +0.45). Correlating
# location ability against path_ratio would therefore be largely an artifact. The clean
# tunneling quantity is the NUMERATOR alone -- how close the two pitches stayed through
# the hitter's decision window -- so that is the headline target here.
agg <- d[, .(
  n            = .N,
  early_sep    = mean(tunnel),
  path_ratio   = mean(path_ratio),
  plate_sep    = mean(plate_sep),
  `High-Miss%`            = mean(inzone_pp),
  `High-Miss% (height)`   = mean(inzone_z),
  `High-Miss% (side)`     = mean(inzone_x),
  `Velocity gap`           = mean(abs(speed_diff), na.rm=TRUE),
  `Movement gap`           = mean(sqrt(ax_diff^2 + az_diff^2), na.rm=TRUE),
  `Spin-axis gap`          = mean(axis_diff, na.rm=TRUE),
  `Active-spin gap`        = mean(abs(active_spin - fb_active), na.rm=TRUE)
), by = .(pitcher, player_name, brk_type)][n >= 40]

DRV <- c("High-Miss%","High-Miss% (height)","High-Miss% (side)",
         "Velocity gap","Movement gap","Spin-axis gap","Active-spin gap")
LOC <- DRV[1:3]
drv <- rbindlist(lapply(TYPES, function(ty) {
  a <- agg[brk_type == ty]
  rbindlist(lapply(DRV, function(v) {
    ok <- is.finite(a[[v]])
    data.table(brk_type = ty, driver = v, n = sum(ok),
               r = cor(a[[v]][ok], a$early_sep[ok]),
               r_pathratio = cor(a[[v]][ok], a$path_ratio[ok]))
  }))
}))
cat("=== (1) DRIVERS OF THE TUNNEL, pitcher x type level ===\n")
cat("  r = vs early-flight separation (clean; + = pitches split sooner)\n")
print(dcast(drv, driver ~ brk_type, value.var = "r")[order(-SL)])
cat("\n  for reference, the same drivers vs path_ratio (denominator-confounded):\n")
print(dcast(drv, driver ~ brk_type, value.var = "r_pathratio")[order(-SL)])
fwrite(drv, file.path(AST, "ext_tunnel_drivers_with_location.csv"))

# Early separation is overwhelmingly a velocity story (r ~ .9), so any location signal has
# to be checked against that: does it survive holding the velocity gap fixed?
pcor <- function(a, b, c) {
  rab <- cor(a,b); rac <- cor(a,c); rbc <- cor(b,c)
  (rab - rac*rbc) / sqrt((1-rac^2)*(1-rbc^2))
}
cat("\n  location drivers vs early separation, raw and net of the velocity gap:\n")
print(rbindlist(lapply(TYPES, function(ty) {
  a <- agg[brk_type == ty]
  rbindlist(lapply(LOC, function(v) {
    ct <- cor.test(a[[v]], a$early_sep)
    data.table(brk_type = ty, driver = v, n_pitchers = nrow(a),
               r = round(unname(ct$estimate), 3), p = signif(ct$p.value, 2),
               r_net_of_velo = round(pcor(a[[v]], a$early_sep, a$`Velocity gap`), 3))
  }))
})))

################################################################################
## (2) location sensitivity of the tunnel: tunnel tercile x zone
################################################################################
# Terciles of early-flight separation, NOT path_ratio: splitting on the ratio would sort
# pitches by how far they finished from the fastball, which is the very thing the zone
# flag measures, and would manufacture the interaction we are testing for.
p[, tun_terc := cut(tunnel, quantile(tunnel, c(0,1/3,2/3,1), na.rm=TRUE),
                    labels = c("tight tunnel","middle","loose tunnel"),
                    include.lowest = TRUE), by = brk_type]

cell2 <- p[!is.na(tun_terc), .(
  pitches   = .N,
  miss_sw   = mean(miss_distance[is_swing], na.rm = TRUE),   # E[miss | swing]
  miss_pp   = mean(miss_pp),                                 # chase-adjusted
  swing_pct = 100*mean(is_swing)),
  by = .(brk_type, tun_terc, inzone_pp)]
setorder(cell2, brk_type, tun_terc, -inzone_pp)

cat("\n=== (2) MISS BY TUNNEL TERCILE x HIGH-MISS ZONE ===\n")
print(cell2[, .(brk_type, tunnel = tun_terc, zone = fifelse(inzone_pp,"in zone","outside"),
                pitches, `E[miss|swing]` = round(miss_sw,2),
                `miss/pitch` = round(miss_pp,2), swing_pct = round(swing_pct,1))])

# Same cut, but with the zone broken into its two axes so we can see which one carries the
# miss: finishing in the zone's height band, its side band, both, or neither.
p[, zone4 := factor(fifelse(inzone_z & inzone_x, "both bands",
                     fifelse(inzone_z, "height band only",
                      fifelse(inzone_x, "side band only", "neither"))),
                    levels = c("both bands","height band only","side band only","neither"))]
cell4 <- p[!is.na(tun_terc), .(
  pitches = .N, miss_sw = mean(miss_distance[is_swing], na.rm = TRUE),
  miss_pp = mean(miss_pp), swing_pct = 100*mean(is_swing)),
  by = .(brk_type, tun_terc, zone4)]
setorder(cell4, brk_type, tun_terc, zone4)
cat("\n=== (2b) MISS BY TUNNEL TERCILE x ZONE AXIS ===\n")
print(cell4[, .(brk_type, tunnel = tun_terc, zone4, pitches,
                `E[miss|swing]` = round(miss_sw,2), `miss/pitch` = round(miss_pp,2),
                swing_pct = round(swing_pct,1))])
pooled <- cell4[, .(m = weighted.mean(miss_sw, pitches)), by = .(brk_type, zone4)]
pooled[, vs_neither := round(m - m[zone4 == "neither"], 2), by = brk_type]
cat("\n  axis value, pooled over terciles (inches of E[miss|swing] vs 'neither'):\n")
print(dcast(pooled, brk_type ~ zone4, value.var = "vs_neither"))

sens <- dcast(cell2, brk_type + tun_terc ~ inzone_pp, value.var = c("miss_sw","miss_pp"))
setnames(sens, c("miss_sw_FALSE","miss_sw_TRUE","miss_pp_FALSE","miss_pp_TRUE"),
         c("out_sw","in_sw","out_pp","in_pp"))
sens[, `:=`(sens_sw = in_sw - out_sw, sens_pp = in_pp - out_pp)]

cat("\n=== LOCATION SENSITIVITY = in-zone minus outside (inches) ===\n")
print(sens[, .(brk_type, tunnel = tun_terc,
               sensitivity_cond = round(sens_sw,2), sensitivity_perpitch = round(sens_pp,2))])

idx <- sens[tun_terc %in% c("tight tunnel","loose tunnel"),
  .(tight = sens_sw[tun_terc=="tight tunnel"], loose = sens_sw[tun_terc=="loose tunnel"]),
  by = brk_type]
idx[, amplification := round(tight/loose, 2)]
cat("\n=== HOW MUCH THE TUNNEL AMPLIFIES LOCATION (E[miss|swing]) ===\n")
print(idx[, .(brk_type, tight_tunnel_sens = round(tight,2), loose_tunnel_sens = round(loose,2),
              amplification)][order(-amplification)])
fwrite(sens, file.path(AST, "ext_tunnel_location_sensitivity.csv"))

################################################################################
## figures
################################################################################
theme_set(theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face="bold"), panel.grid.major.x = element_blank(),
        strip.text = element_text(face="bold")))
ACC <- "#c0392b"; BLU <- "#2c3e50"

c2 <- copy(cell2)[!is.na(tun_terc)]
c2[, zone := factor(fifelse(inzone_pp, "in the high-miss zone", "outside the zone"),
                    levels = c("in the high-miss zone","outside the zone"))]
c2[, brk_lab := TYPE_LAB[brk_type]]
lab <- idx[, .(brk_lab = TYPE_LAB[brk_type],
               txt = sprintf("location value, tight vs loose tunnel: %.2fx", amplification))]

g1 <- ggplot(c2, aes(tun_terc, miss_sw, fill = zone)) +
  geom_col(position = position_dodge(width = .72), width = .68) +
  geom_text(aes(label = sprintf("%.1f", miss_sw)),
            position = position_dodge(width = .72), vjust = -0.4, size = 3.1) +
  geom_text(data = lab, inherit.aes = FALSE, aes(x = 2, y = 8.5, label = txt),
            size = 3.5, fontface = "italic") +
  facet_wrap(~ brk_lab) +
  scale_fill_manual(values = c("in the high-miss zone" = ACC, "outside the zone" = BLU), name = NULL) +
  coord_cartesian(ylim = c(0, 9)) +
  labs(title = "How much a tunnel depends on location is a property of the pitch type",
       subtitle = paste0("2-strike breaking balls after the primary fastball, 2023H2-2026. The red/navy gap is what landing in the high-miss zone is worth.\n",
                         "Curveballs: that gap widens as the tunnel tightens (6.4 vs 4.8 in) - a tunneled curveball lives or dies on its spot. Sliders run the\n",
                         "other way and sweepers are flat, so a tight slider tunnel keeps missing bats even when the location is imperfect."),
       x = "trajectory tunnel (terciles of early-flight separation)", y = "E[miss | swing]  (inches)") +
  theme(legend.position = "top")
ggsave(file.path(AST, "fig3c_tunnel_location_sensitivity.png"), g1, width = 12, height = 6.4, dpi = 150)
ggsave(file.path(ODIR, "fig3c_tunnel_location_sensitivity.png"), g1, width = 12, height = 6.4, dpi = 150)

# Same chart with the zone opened up into its two axes.
c4 <- copy(cell4)[!is.na(tun_terc)]
c4[, brk_lab := TYPE_LAB[brk_type]]
AXC <- c("both bands" = ACC, "height band only" = "#e08e0b",
         "side band only" = "#7f8c9b", "neither" = BLU)
g1b <- ggplot(c4, aes(tun_terc, miss_sw, fill = zone4)) +
  geom_col(position = position_dodge(width = .82), width = .78) +
  geom_text(aes(label = sprintf("%.1f", miss_sw)),
            position = position_dodge(width = .82), vjust = -0.4, size = 2.7) +
  facet_wrap(~ brk_lab) +
  scale_fill_manual(values = AXC, name = NULL) +
  coord_cartesian(ylim = c(0, 9)) +
  labs(title = "Getting a breaking ball to the right height is most of what the high-miss zone is",
       subtitle = paste0("The high-miss zone, opened into its two axes: the pitch finished in the zone's height band, its side band, both, or neither.\n",
                         "Height carries nearly all of it - side alone is worth 0.2-0.6 in. Curveballs need only the height (+5.9 in on its own, vs +5.5 for\n",
                         "both); the sweeper is the one type that truly needs both (+4.3 vs +1.9). Caveat: height-only pitches miss sideways and draw half\n",
                         "as many swings, so per pitch both bands still wins. 2-strike breaking balls after the primary fastball, 2023H2-2026."),
       x = "trajectory tunnel (terciles of early-flight separation)", y = "E[miss | swing]  (inches)") +
  theme(legend.position = "top")
ggsave(file.path(AST, "fig3c2_zone_axis_split.png"), g1b, width = 12, height = 6.4, dpi = 150)
ggsave(file.path(ODIR, "fig3c2_zone_axis_split.png"), g1b, width = 12, height = 6.4, dpi = 150)

dv <- copy(drv); dv[, brk_lab := TYPE_LAB[brk_type]]
dv[, kind := fifelse(driver %in% LOC, "location", "shape / spin gap")]
ord <- dv[, .(m = mean(r)), by = driver][order(m)]$driver
dv[, driver := factor(driver, levels = ord)]
g2 <- ggplot(dv, aes(driver, r, fill = kind)) +
  geom_col(width = .68) + geom_hline(yintercept = 0, colour = "black", linewidth = .3) +
  geom_text(aes(label = sprintf("%+.2f", r), hjust = ifelse(r < 0, 1.15, -0.15)), size = 3.2) +
  coord_flip(ylim = c(-1, 1.28)) + facet_wrap(~ brk_lab) +
  scale_fill_manual(values = c(location = ACC, `shape / spin gap` = BLU), name = NULL) +
  labs(title = "What loosens a tunnel, with location ability folded in",
       subtitle = paste0("Pitcher x pitch-type level (min 40 pairs). + = the trait goes with pitches splitting SOONER, - = staying together longer.\n",
                         "High-Miss% is the share of breaking balls landing in the put-away spot - the low glove-side region where whiff rate, swinging-strike\n",
                         "rate and miss distance all peak, a foot-plus below the chase peak. It is split into its height and side bands here.\n",
                         "Measured against early-flight separation, not path_ratio: the spot sits far from the fastball, so it would shrink that ratio mechanically."),
       x = NULL, y = "correlation with early-flight separation (decision window)") +
  theme(legend.position = "top", panel.grid.major.y = element_blank())
for (f in c("fig3d_tunnel_drivers_with_location.png", "fig3d_tunnel_drivers_putaway.png")) {
  ggsave(file.path(AST, f),  g2, width = 12, height = 6.0, dpi = 150)
  ggsave(file.path(ODIR, f), g2, width = 12, height = 6.0, dpi = 150)
}

cat("\nwrote fig3c_tunnel_location_sensitivity.png, fig3c2_zone_axis_split.png and fig3d_tunnel_drivers_putaway.png\n")
