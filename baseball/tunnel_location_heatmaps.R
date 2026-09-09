#!/usr/bin/env Rscript

# PHASE 2b - Verification heatmaps for the tunnel-location surfaces.
#
#   Surface A: where the setup FASTBALL goes  -> expected miss on the next breaking ball
#   Surface B: where the BREAKING BALL goes   -> expected miss, given a well-tunneled FB
#
# Colour = expected miss distance in inches PER PITCH THROWN (chase-adjusted:
# P(swing) * E[miss|swing]); white X = optimal cell. All plots are catcher's view,
# so plate_x > 0 is to the catcher's right. Dashed box = nominal strike zone.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

MDIR <- file.path("data","statcast_model")
ODIR <- file.path(MDIR, "tunnel_location")
AST  <- file.path(MDIR, "article_assets")

league <- readRDS(file.path(ODIR, "league_surfaces.rds"))
psurf  <- readRDS(file.path(ODIR, "pitcher_surfaces.rds"))
optima <- fread(file.path(ODIR, "ext_tunnel_optima.csv"))

TYPE_LAB <- c(SL = "Slider", CU = "Curveball (CU+KC)", ST = "Sweeper")
ZONE <- data.table(xmin = -0.83, xmax = 0.83, ymin = 1.6, ymax = 3.4)

theme_set(theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 10)))

base_map <- function(dt, optdt) {
  ggplot(dt, aes(gx, gz)) +
    geom_raster(aes(fill = val), interpolate = TRUE) +
    geom_rect(data = ZONE, inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = NA, colour = "grey20", linetype = "dashed", linewidth = .4) +
    geom_point(data = optdt, inherit.aes = FALSE, aes(ox, oz),
               shape = 4, size = 3, stroke = 1.4, colour = "white") +
    scale_fill_viridis_c(option = "inferno", na.value = "white",
                         name = "E[miss | swing]\n(inches)") +
    scale_x_continuous(breaks = c(-1, 0, 1)) +
    scale_y_continuous(breaks = 1:4) +
    geom_hline(yintercept = 0, colour = "grey35", linewidth = .3) +
    coord_fixed(xlim = c(-2.0, 2.0), ylim = c(-0.6, 4.4), expand = FALSE) +
    labs(x = "plate_x (ft, catcher's view)", y = "plate_z (ft)")
}

## ---------------------------------------------------------------- league ----
lg_opt <- optima[scope == "league"]
for (ty in names(TYPE_LAB)) {
  d <- league[brk_type == ty]
  if (!nrow(d)) next
  long <- rbindlist(list(
    d[, .(gx, gz, p_throws, stand, val = fitA, surface = "A: setup FASTBALL location")],
    d[, .(gx, gz, p_throws, stand, val = fitB, surface = "B: BREAKING BALL location (well-tunneled)")]))
  long[, panel := sprintf("%sHP vs %sHB", p_throws, stand)]

  o <- lg_opt[brk_type == ty]
  oo <- rbindlist(list(
    o[, .(panel = sprintf("%sHP vs %sHB", p_throws, stand),
          surface = "A: setup FASTBALL location", ox = fb_opt_x, oz = fb_opt_z)],
    o[, .(panel = sprintf("%sHP vs %sHB", p_throws, stand),
          surface = "B: BREAKING BALL location (well-tunneled)", ox = bb_opt_x, oz = bb_opt_z)]))

  p <- base_map(long, oo) +
    facet_grid(surface ~ panel) +
    labs(title = sprintf("Tunnel location map - %s (2-strike counts, after the primary fastball)", TYPE_LAB[ty]),
         subtitle = "League aggregate, 2023H2-2026. Colour = expected miss distance in inches given a swing. White X = optimal spot.")
  ggsave(file.path(ODIR, sprintf("fig_league_%s.png", ty)), p, width = 13, height = 7.2, dpi = 150)
  message("wrote fig_league_", ty, ".png")
}

## --------------------------------------------------- per-pitcher dashboard ----
# Rank 2026 arms by pooled sample so the examples are the best-supported ones.
rank26 <- unique(psurf[, .(pitcher, player_name, brk_type, p_throws, stand, n_A, n_B)])
rank26 <- rank26[, .(n_tot = sum(n_A), n_wt = sum(n_B)), by = .(pitcher, player_name, brk_type, p_throws)]
setorder(rank26, brk_type, -n_tot)
fwrite(rank26, file.path(ODIR, "pitcher_lookup_index.csv"))

draw_pitcher <- function(pid, ty) {
  d <- psurf[pitcher == pid & brk_type == ty]
  if (!nrow(d)) return(invisible(NULL))
  nm <- d$player_name[1]
  long <- rbindlist(list(
    d[, .(gx, gz, stand, n_A, n_B, val = fitA, surface = "A: setup FASTBALL location")],
    d[, .(gx, gz, stand, n_A, n_B, val = fitB, surface = "B: BREAKING BALL location (well-tunneled)")]))
  long[, panel := sprintf("vs %sHB  (n=%d)", stand, n_A)]

  o <- optima[scope == "pitcher" & pitcher == pid & brk_type == ty]
  oo <- rbindlist(list(
    o[, .(panel = sprintf("vs %sHB  (n=%d)", stand, n_A),
          surface = "A: setup FASTBALL location", ox = fb_opt_x, oz = fb_opt_z)],
    o[, .(panel = sprintf("vs %sHB  (n=%d)", stand, n_A),
          surface = "B: BREAKING BALL location (well-tunneled)", ox = bb_opt_x, oz = bb_opt_z)]))

  p <- base_map(long, oo) +
    facet_grid(surface ~ panel) +
    labs(title = sprintf("%s - %s tunnel location map", nm, TYPE_LAB[ty]),
         subtitle = "Pooled 2023H2-2026, EB-shrunk toward the league prior. White X = optimal spot.")
  fn <- file.path(ODIR, sprintf("fig_pitcher_%s_%d.png", ty, pid))
  ggsave(fn, p, width = 9, height = 7.2, dpi = 150)
  message("wrote ", basename(fn), "  [", nm, "]")
}

for (ty in names(TYPE_LAB)) {
  top <- rank26[brk_type == ty][1:2]
  for (i in seq_len(nrow(top))) draw_pitcher(top$pitcher[i], ty)
}

cat("\n=============== TOP 2026 ARMS BY POOLED SAMPLE ===============\n")
print(rank26[, head(.SD, 5), by = brk_type][, .(brk_type, player_name, p_throws, n_tot, n_wt)])
cat(sprintf("\nlookup index: %s\n", file.path(ODIR, "pitcher_lookup_index.csv")))
