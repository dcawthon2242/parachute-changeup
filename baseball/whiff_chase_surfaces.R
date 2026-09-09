#!/usr/bin/env Rscript

# The high-miss heatmaps, refit on whiff and on chase.
#
# Same pairs, same kernel, same top-20% zone rule as fig_highmiss_zone.png, but with two
# other targets:
#   whiff | swing : P(swing and miss) among swings
#   whiff per pitch : P(swing) x P(whiff | swing), i.e. swinging-strike rate over all
#                     pitches -- the product of the two surfaces below, and the one a
#                     pitcher actually gets paid for
#   chase         : P(swing) among pitches outside the rulebook zone
# The point is where the zones sit relative to each other. The conditional whiff zone and
# the chase zone are the two ends of the tradeoff; the per-pitch whiff zone is where the
# arithmetic actually lands.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

MDIR <- file.path("data","statcast_model")
ODIR <- file.path(MDIR, "tunnel_location")
AST  <- file.path(MDIR, "article_assets")

p <- readRDS(file.path(MDIR, "tunnel_pairs.rds"))[is_primary_setup == TRUE]
mir <- function(x, h) fifelse(h == "R", x, -x)
p[, `:=`(bb_x = mir(bb_plate_x, p_throws), bb_z = bb_plate_z)]
p[, matchup := fifelse(p_throws == stand, "same", "opp")]
p <- p[is_swing == FALSE | is.finite(miss_distance)]

SZ_BOT <- 1.6; SZ_TOP <- 3.4; SZ_HALF <- 0.83
p[, in_rulebook := abs(bb_x) <= SZ_HALF & bb_z >= SZ_BOT & bb_z <= SZ_TOP]

GX <- seq(-2.0, 2.0, by = 0.10); GZ <- seq(-0.6, 4.4, by = 0.10)
grid <- as.data.table(expand.grid(gx = GX, gz = GZ))
H <- 0.35; K_LEAGUE <- 20; REL_SUP <- 0.12; CORE_SUP <- 0.25
ZONE_Q <- 0.80

kern_sums <- function(x, z, y, G = grid) {
  sw <- numeric(nrow(G)); sy <- numeric(nrow(G)); n <- length(x)
  for (i in seq(1, n, by = 4000)) {
    j <- i:min(i + 3999, n)
    W <- exp(-0.5 * (outer(G$gx, x[j], "-")^2 + outer(G$gz, z[j], "-")^2) / H^2)
    sw <- sw + rowSums(W); sy <- sy + as.numeric(W %*% y[j])
  }
  list(w = sw, y = sy)
}
fit_surface <- function(x, z, y, G = grid) {
  k <- kern_sums(x, z, y, G)
  loc <- k$y / pmax(k$w, 1e-9)
  fit <- (k$w * loc + K_LEAGUE * mean(y)) / (k$w + K_LEAGUE)
  ok  <- k$w >= max(12, min(40, 0.02*length(x))) & k$w >= REL_SUP * max(k$w)
  fit[!ok] <- NA_real_
  list(fit = fit, dens = k$w / max(k$w))
}
argmax_cell <- function(fit, dens, G = grid) {
  core <- !is.na(fit) & dens >= CORE_SUP
  if (!any(core)) core <- !is.na(fit)
  v <- fit; v[!core] <- NA_real_; i <- which.max(v)
  list(x = G$gx[i], z = G$gz[i], val = v[i])
}

TRAIN <- p[season <= 2025]          # same training window as the miss zone
TYPES <- c("SL","CU","ST")
TYPE_LAB <- c(SL="Slider", CU="Curveball (CU+KC)", ST="Sweeper")

################################################################################
## fit
################################################################################
TRAIN[, whiff_pp := as.numeric(is_swing == TRUE & is_whiff == TRUE)]

out <- list(); opt <- list()
for (ty in TYPES) for (mu in c("same","opp")) {
  dsw <- TRAIN[brk_type == ty & matchup == mu & is_swing == TRUE]
  dal <- TRAIN[brk_type == ty & matchup == mu]
  dch <- dal[in_rulebook == FALSE]
  if (nrow(dsw) < 400) next
  sw <- fit_surface(dsw$bb_x, dsw$bb_z, as.numeric(dsw$is_whiff))
  pp <- fit_surface(dal$bb_x, dal$bb_z, dal$whiff_pp)
  ch <- fit_surface(dch$bb_x, dch$bb_z, as.numeric(dch$is_swing))
  zw <- !is.na(sw$fit) & sw$fit >= quantile(sw$fit, ZONE_Q, na.rm = TRUE)
  zp <- !is.na(pp$fit) & pp$fit >= quantile(pp$fit, ZONE_Q, na.rm = TRUE)
  zc <- !is.na(ch$fit) & ch$fit >= quantile(ch$fit, ZONE_Q, na.rm = TRUE)
  aw <- argmax_cell(sw$fit, sw$dens); ap <- argmax_cell(pp$fit, pp$dens)
  ac <- argmax_cell(ch$fit, ch$dens)
  # Read the swinging-strike surface at the other two peaks: is the chase spot or the
  # whiff spot the better place to actually get a swinging strike?
  idx_of <- function(x, z) which(abs(grid$gx - x) < 1e-6 & abs(grid$gz - z) < 1e-6)[1]
  out[[paste(ty,mu)]] <- data.table(brk_type=ty, matchup=mu, gx=grid$gx, gz=grid$gz,
    whiff = sw$fit, whiffpp = pp$fit, chase = ch$fit,
    zone_whiff = zw, zone_whiffpp = zp, zone_chase = zc)
  opt[[length(opt)+1]] <- data.table(brk_type=ty, matchup=mu,
    n_swings = nrow(dsw), n_pitches = nrow(dal), n_outzone = nrow(dch),
    whiff_peak = sprintf("(%+.1f,%+.1f)", aw$x, aw$z), whiff_val = round(100*aw$val,1),
    swstr_peak = sprintf("(%+.1f,%+.1f)", ap$x, ap$z), swstr_val = round(100*ap$val,1),
    chase_peak = sprintf("(%+.1f,%+.1f)", ac$x, ac$z), chase_val = round(100*ac$val,1),
    whiff_peak_z = aw$z, swstr_peak_z = ap$z, chase_peak_z = ac$z,
    whiff_peak_x = aw$x, swstr_peak_x = ap$x, chase_peak_x = ac$x,
    swstr_at_chase_peak = round(100*pp$fit[idx_of(ac$x, ac$z)], 1),
    swstr_at_whiff_peak = round(100*pp$fit[idx_of(aw$x, aw$z)], 1))
}
S <- rbindlist(out); O <- rbindlist(opt)

cat("=== PEAK OF EACH SURFACE (mirrored ft, + x = glove side) ===\n")
print(O[, .(brk_type, matchup, n_pitches,
            chase_peak, `chase%` = chase_val,
            swstr_peak, `swstr%` = swstr_val,
            whiff_peak, `whiff|sw%` = whiff_val)])
cat("\nwhere the per-pitch whiff peak sits between the other two (ft):\n")
print(O[, .(brk_type, matchup,
            below_chase_ft = round(chase_peak_z - swstr_peak_z, 2),
            above_condwhiff_ft = round(swstr_peak_z - whiff_peak_z, 2),
            chase_to_whiff_ft  = round(chase_peak_z - whiff_peak_z, 2))])
cat("\nswinging strikes per pitch if you aim at the chase spot vs the whiff spot (%):\n")
print(O[, .(brk_type, matchup, at_chase_spot = swstr_at_chase_peak,
            at_whiff_spot = swstr_at_whiff_peak, best = swstr_val)])
cat("\nmean separation between the chase peak and the conditional whiff peak:\n")
print(O[, .(height_ft = round(mean(chase_peak_z - whiff_peak_z),2),
            side_ft   = round(mean(abs(chase_peak_x - whiff_peak_x)),2)), by = brk_type])

# The miss surfaces were already fit with this same kernel and training window by
# tunnel_command.R: `fit` is E[miss | swing], `fit_pp` is E[miss per pitch].
zs <- readRDS(file.path(ODIR, "zone_surfaces.rds"))
S <- merge(S, zs[, .(brk_type, matchup, gx, gz, dens, misscond = fit, misspp = fit_pp,
                     zone_misscond = inzone, zone_miss = inzone_pp)],
           by = c("brk_type","matchup","gx","gz"), all.x = TRUE)
S[is.na(zone_miss), zone_miss := FALSE]; S[is.na(zone_misscond), zone_misscond := FALSE]

mp <- S[, {
  core <- !is.na(misspp) & dens >= CORE_SUP
  v <- misspp; v[!core] <- NA_real_; i <- which.max(v)
  cc <- !is.na(misscond) & dens >= CORE_SUP
  w <- misscond; w[!cc] <- NA_real_; j <- which.max(w)
  .(misspp_peak = sprintf("(%+.1f,%+.1f)", gx[i], gz[i]), misspp_val = round(v[i],2),
    misspp_peak_x = gx[i], misspp_peak_z = gz[i],
    misscond_peak = sprintf("(%+.1f,%+.1f)", gx[j], gz[j]), misscond_val = round(w[j],2))
}, by = .(brk_type, matchup)]
O <- merge(O, mp, by = c("brk_type","matchup"))
cat("\n=== MISS DISTANCE, PER PITCH VS PER SWING ===\n")
print(O[, .(brk_type, matchup, misspp_peak, `miss/pitch_in` = misspp_val,
            misscond_peak, `miss|swing_in` = misscond_val,
            swstr_peak, `swstr%` = swstr_val)])
cat("\nhow far the per-pitch miss peak sits from the other peaks (ft, + = higher):\n")
print(O[, .(brk_type, matchup,
            chase_above_misspp = round(chase_peak_z - misspp_peak_z, 2),
            swstr_above_misspp = round(swstr_peak_z - misspp_peak_z, 2),
            misspp_above_condwhiff = round(misspp_peak_z - whiff_peak_z, 2))])
cat("\noverlap of the top-20% zones (% of the union that both cover):\n")
print(S[, .(whiff_vs_chase  = round(100*sum(zone_whiff & zone_chase)/sum(zone_whiff | zone_chase),1),
            swstr_vs_chase  = round(100*sum(zone_whiffpp & zone_chase)/sum(zone_whiffpp | zone_chase),1),
            swstr_vs_whiff  = round(100*sum(zone_whiffpp & zone_whiff)/sum(zone_whiffpp | zone_whiff),1),
            swstr_vs_miss   = round(100*sum(zone_whiffpp & zone_miss )/sum(zone_whiffpp | zone_miss ),1)),
        by = .(brk_type, matchup)])
fwrite(O, file.path(AST, "ext_whiff_chase_peaks.csv"))

################################################################################
## figures, in the same six-panel format as fig_highmiss_zone.png
################################################################################
theme_set(theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face="bold"), panel.grid = element_blank(),
        strip.text = element_text(face="bold", size=10)))
SZBOX <- data.table(xmin=-0.83, xmax=0.83, ymin=1.6, ymax=3.4)
source(file.path("baseball", "batter_silhouette.R"))
MLEV <- c("vs same-handed hitter", "vs opposite-handed hitter")
mlab_of <- function(m) factor(fifelse(m == "same", MLEV[1], MLEV[2]), levels = MLEV)
tlab_of <- function(t) factor(TYPE_LAB[t], levels = unname(TYPE_LAB[TYPES]))
XL <- c(-3.05, 3.05); YL <- c(-0.6, 6.2)
batter_layer <- function() c(
  batter_layers(list(list(side = -1, data = data.table(mlab = mlab_of("same"))),
                     list(side =  1, data = data.table(mlab = mlab_of("opp"))))),
  list(geom_hline(yintercept = 0, colour = "grey35", linewidth = .4)))
zone_grid <- function() list(
  facet_grid(mlab ~ tlab, switch = "y"),
  scale_x_continuous(breaks = c(-1,0,1)),
  coord_fixed(xlim = XL, ylim = YL, expand = FALSE),
  theme(strip.text.y.left = element_text(angle = 90, face = "bold")))

S[, `:=`(mlab = mlab_of(matchup), tlab = tlab_of(brk_type))]

heat <- function(val, zone, ttl, sub, legend, scale = 100) {
  ggplot(S, aes(gx, gz)) +
    geom_raster(aes(fill = scale * .data[[val]]), interpolate = TRUE) +
    geom_contour(aes(z = as.numeric(.data[[zone]])), breaks = 0.5,
                 colour = "white", linewidth = .6) +
    geom_rect(data=SZBOX, inherit.aes=FALSE, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
              fill=NA, colour="grey20", linetype="dashed", linewidth=.4) +
    batter_layer() + zone_grid() +
    scale_fill_viridis_c(option="inferno", na.value="white", name=legend) +
    labs(title = ttl, subtitle = sub, x = "mirrored plate_x (ft)", y = "plate_z (ft)")
}

g1 <- heat("whiff", "zone_whiff",
  "Where a 2-strike breaking ball gets swung through",
  "P(whiff | swing), fit on 2023H2-2025. White outline = top 20% of the surface. Mirrored: + x = glove side, so the hitter switches boxes between rows.",
  "whiff rate\ngiven a swing (%)")
ggsave(file.path(ODIR, "fig_whiff_zone.png"), g1, width = 13.5, height = 10.0, dpi = 150)
ggsave(file.path(AST,  "fig_whiff_zone.png"), g1, width = 13.5, height = 10.0, dpi = 150)

g1b <- heat("whiffpp", "zone_whiffpp",
  "Where a 2-strike breaking ball actually produces a swinging strike",
  "Whiffs per pitch, not per swing: P(swing) x P(whiff | swing), fit on 2023H2-2025. White outline = top 20%. Mirrored: + x = glove side.",
  "swinging strikes\nper pitch (%)")
ggsave(file.path(ODIR, "fig_swstr_zone.png"), g1b, width = 13.5, height = 10.0, dpi = 150)
ggsave(file.path(AST,  "fig_swstr_zone.png"), g1b, width = 13.5, height = 10.0, dpi = 150)

g1c <- heat("misspp", "zone_miss",
  "Miss distance per pitch, with takes counted as zero",
  "E[miss | swing] x P(swing), fit on 2023H2-2025 - the chase-adjusted companion to fig_highmiss_zone.png, which conditions on a swing.\nWhite outline = top 20%. Mirrored: + x = glove side.",
  "miss distance\nper pitch (in)", scale = 1)
ggsave(file.path(ODIR, "fig_misspp_zone.png"), g1c, width = 13.5, height = 10.0, dpi = 150)
ggsave(file.path(AST,  "fig_misspp_zone.png"), g1c, width = 13.5, height = 10.0, dpi = 150)

g2 <- heat("chase", "zone_chase",
  "Where a 2-strike breaking ball gets chased",
  "P(swing) among pitches outside the rulebook zone, fit on 2023H2-2025. White outline = top 20%. Mirrored: + x = glove side.",
  "chase rate (%)")
ggsave(file.path(ODIR, "fig_chase_zone.png"), g2, width = 13.5, height = 10.0, dpi = 150)
ggsave(file.path(AST,  "fig_chase_zone.png"), g2, width = 13.5, height = 10.0, dpi = 150)

# The comparison figure: all three zones as outlines on one set of panels.
ZLEV <- c("chase (swing rate)", "swinging strikes per pitch",
          "whiff given a swing", "miss distance per pitch (takes = 0)")
ZC <- c("#2c7fb8", "#1a9850", "#c0392b", "#f0a202"); names(ZC) <- ZLEV
g3 <- ggplot(S, aes(gx, gz)) +
  geom_raster(aes(fill = 100*whiffpp), interpolate = TRUE, alpha = .6) +
  geom_contour(aes(z = as.numeric(zone_chase),   colour = ZLEV[1]), breaks = .5, linewidth = .9) +
  geom_contour(aes(z = as.numeric(zone_whiffpp), colour = ZLEV[2]), breaks = .5, linewidth = 1.1) +
  geom_contour(aes(z = as.numeric(zone_whiff),   colour = ZLEV[3]), breaks = .5, linewidth = .9) +
  geom_contour(aes(z = as.numeric(zone_miss),    colour = ZLEV[4]), breaks = .5,
               linewidth = .8, linetype = "22") +
  geom_rect(data=SZBOX, inherit.aes=FALSE, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            fill=NA, colour="grey20", linetype="dashed", linewidth=.4) +
  batter_layer() + zone_grid() +
  scale_fill_gradient(low = "white", high = "grey60", na.value = "white", guide = "none") +
  scale_colour_manual(values = ZC, breaks = ZLEV, name = "top 20% zone for") +
  guides(colour = guide_legend(nrow = 1)) +
  labs(title = "The swinging strike lives where the whiff lives, not where the chase lives",
       subtitle = paste0("Top-20% regions of four surfaces on the same panels, fit on 2023H2-2025. Grey shading is the swinging-strike surface underneath.\n",
                         "The chase zone (blue) sits at the bottom of the strike zone and the conditional whiff zone (red) more than a foot below it, sharing none\n",
                         "of their area. Whiffs per pitch (green) is the outcome that counts, and it sides with the whiff zone: it overlaps that zone 35-61% but\n",
                         "the chase zone only 0-13%. The higher miss rate down there more than pays for the swings you give up by leaving the chase zone."),
       x = "mirrored plate_x (ft)", y = "plate_z (ft)") +
  theme(legend.position = "top")
ggsave(file.path(ODIR, "fig_zone_overlay.png"), g3, width = 13.5, height = 10.4, dpi = 150)
ggsave(file.path(AST,  "fig_zone_overlay.png"), g3, width = 13.5, height = 10.4, dpi = 150)

cat("\nwrote fig_whiff_zone.png, fig_swstr_zone.png, fig_misspp_zone.png, fig_chase_zone.png, fig_zone_overlay.png\n")
