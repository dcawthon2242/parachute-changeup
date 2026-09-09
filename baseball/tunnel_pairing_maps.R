#!/usr/bin/env Rscript

# Pairing-specific tunnel location maps: SI -> ST (sinker/sweeper) and FF -> CU.
#
# Restricted to pitchers with ELITE path ratios for that specific pairing (tightest
# tercile of pitcher-mean path_to_location_ratio), and compared against everyone else,
# to answer: where should an arm that actually tunnels this pair put the two pitches?
#
# Handedness is MIRRORED so both hands pool: plate_x is flipped for LHP, and the
# resulting axis reads + = pitcher's GLOVE side, - = arm side. (Verified against the
# per-hand fits: a RHP's putaway slider sits at positive plate_x, a LHP's at negative.)
# Batter side is expressed relative to the pitcher: same-handed vs opposite-handed.
#
# Target: E[miss | swing], in inches.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

MDIR <- file.path("data","statcast_model")
ODIR <- file.path(MDIR, "tunnel_location")
dir.create(ODIR, showWarnings = FALSE, recursive = TRUE)

pairs <- readRDS(file.path(MDIR, "tunnel_pairs.rds"))
pairs[, x_m_fb := fifelse(p_throws == "R", fb_plate_x, -fb_plate_x)]
pairs[, x_m_bb := fifelse(p_throws == "R", bb_plate_x, -bb_plate_x)]
pairs[, matchup := fifelse(p_throws == stand, "same-handed hitter", "opposite-handed hitter")]

PAIRINGS <- list(
  list(id = "SI_ST", setup = "SI", brk = "ST", lab = "Sinker -> Sweeper"),
  list(id = "FF_CU", setup = "FF", brk = "CU", lab = "Four-seam -> Curveball"))

MIN_PITCHER_N <- 25    # pairs needed for a pitcher to be rated on this pairing
H <- 0.35; K_LEAGUE <- 20; K_SHRINK <- 10; REL_SUP <- 0.12; CORE_SUP <- 0.25

GX <- seq(-2.0, 2.0, by = 0.10); GZ <- seq(-0.6, 4.4, by = 0.10)
grid <- as.data.table(expand.grid(gx = GX, gz = GZ))

kern <- function(x, z) exp(-0.5 * (outer(grid$gx, x, "-")^2 + outer(grid$gz, z, "-")^2) / H^2)

## E[miss|swing] surface over swings, masked to real support
fit_surface <- function(x, z, y, prior = NULL) {
  W <- kern(x, z); ns <- rowSums(W)
  loc <- as.numeric((W %*% y) / pmax(ns, 1e-9))
  fit <- if (is.null(prior)) (ns*loc + K_LEAGUE*mean(y)) / (ns + K_LEAGUE)
         else (ns*loc + K_SHRINK*prior) / (ns + K_SHRINK)
  ok <- ns >= max(12, min(40, 0.02*length(x))) & ns >= REL_SUP * max(ns)
  fit[!ok] <- NA_real_
  if (!is.null(prior)) fit[is.na(prior)] <- NA_real_
  list(fit = fit, dens = ns / max(ns))
}

argmax_cell <- function(fit, dens) {
  core <- !is.na(fit) & dens >= CORE_SUP
  if (!any(core)) core <- !is.na(fit)
  if (!any(core)) return(list(x=NA_real_, z=NA_real_, val=NA_real_, rng=NA_real_))
  v <- fit; v[!core] <- NA_real_; i <- which.max(v)
  list(x = grid$gx[i], z = grid$gz[i], val = v[i], rng = diff(range(v, na.rm=TRUE)))
}

ZONE <- data.table(xmin = -0.83, xmax = 0.83, ymin = 1.6, ymax = 3.4)
theme_set(theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face="bold"), panel.grid = element_blank(),
        strip.text = element_text(face="bold", size = 9.5)))

base_map <- function(dt, optdt) {
  ggplot(dt, aes(gx, gz)) +
    geom_raster(aes(fill = val), interpolate = TRUE) +
    geom_rect(data = ZONE, inherit.aes = FALSE,
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              fill = NA, colour = "grey20", linetype = "dashed", linewidth = .4) +
    geom_hline(yintercept = 0, colour = "grey35", linewidth = .3) +
    geom_point(data = optdt, inherit.aes = FALSE, aes(ox, oz),
               shape = 4, size = 3, stroke = 1.4, colour = "white") +
    scale_fill_viridis_c(option = "inferno", na.value = "white",
                         name = "E[miss | swing]\n(inches)") +
    scale_x_continuous(breaks = c(-1,0,1)) + scale_y_continuous(breaks = 1:4) +
    coord_fixed(xlim = c(-2.0, 2.0), ylim = c(-0.6, 4.4), expand = FALSE) +
    labs(x = "mirrored plate_x (ft):  + = glove side,  - = arm side", y = "plate_z (ft)")
}

all_opt <- list()

for (P in PAIRINGS) {
  d <- pairs[setup_type == P$setup & brk_type == P$brk]
  d <- d[is_swing == FALSE | is.finite(miss_distance)]

  ## ---- rate pitchers on this pairing's tunnel tightness ----
  pr <- d[, .(n = .N, n_sw = sum(is_swing), path = mean(path_ratio)), by = .(pitcher, player_name, p_throws)]
  pr <- pr[n >= MIN_PITCHER_N]
  cut_elite <- quantile(pr$path, 1/3)
  pr[, grp := fifelse(path <= cut_elite, "elite tunnel", "rest")]
  d <- merge(d, pr[, .(pitcher, grp, pitcher_path = path)], by = "pitcher")

  message(sprintf("\n%s: %d rated pitchers, elite cutoff path_ratio <= %.3f (%d elite)",
                  P$lab, nrow(pr), cut_elite, sum(pr$grp == "elite tunnel")))

  ## ---- group surfaces: elite vs rest, by matchup ----
  panels <- list(); popt <- list()
  for (g in c("elite tunnel","rest")) for (mu in unique(d$matchup)) {
    dd <- d[grp == g & matchup == mu]
    ss <- dd[is_swing == TRUE]
    if (nrow(ss) < 120) next
    sA <- fit_surface(ss$x_m_fb, ss$fb_plate_z, ss$miss_distance)
    sB <- fit_surface(ss$x_m_bb, ss$bb_plate_z, ss$miss_distance)
    pk <- sprintf("%s\n%s  (n=%d)", g, mu, nrow(ss))
    panels[[length(panels)+1]] <- rbindlist(list(
      data.table(gx=grid$gx, gz=grid$gz, val=sA$fit, pk=pk, surface="A: setup pitch location"),
      data.table(gx=grid$gx, gz=grid$gz, val=sB$fit, pk=pk, surface="B: breaking ball location")))
    aA <- argmax_cell(sA$fit, sA$dens); aB <- argmax_cell(sB$fit, sB$dens)
    popt[[length(popt)+1]] <- rbindlist(list(
      data.table(pk=pk, surface="A: setup pitch location", ox=aA$x, oz=aA$z),
      data.table(pk=pk, surface="B: breaking ball location", ox=aB$x, oz=aB$z)))
    all_opt[[length(all_opt)+1]] <- data.table(
      pairing=P$lab, scope=g, matchup=mu, player_name=NA_character_, pitcher=NA_integer_,
      n_swings=nrow(ss), pitcher_path=mean(dd$pitcher_path),
      setup_x=aA$x, setup_z=aA$z, setup_miss=aA$val,
      bb_x=aB$x, bb_z=aB$z, bb_miss=aB$val, bb_rng=aB$rng)
    if (g == "elite tunnel") assign(paste0("prior_", gsub("[^a-z]","",mu)), list(A=sA$fit, B=sB$fit))
  }
  pl <- rbindlist(panels); po <- rbindlist(popt)
  pl[, pk := factor(pk, levels = unique(po$pk))]; po[, pk := factor(pk, levels = levels(pl$pk))]

  p <- base_map(pl, po) + facet_grid(surface ~ pk) +
    labs(title = sprintf("%s tunnel location map - elite tunnelers vs everyone else", P$lab),
         subtitle = sprintf("2-strike counts, 2023H2-2026. Elite = tightest third of pitcher-mean path ratio (<= %.2f). White X = optimal spot.", cut_elite))
  ggsave(file.path(ODIR, sprintf("fig_pairing_%s.png", P$id)), p, width = 13.5, height = 7.4, dpi = 150)
  message("  wrote fig_pairing_", P$id, ".png")

  ## ---- individual elite arms, EB-shrunk toward the elite-group prior ----
  elite <- pr[grp == "elite tunnel"][order(-n_sw)]
  for (i in seq_len(min(3, nrow(elite)))) {
    pid <- elite$pitcher[i]; nm <- elite$player_name[i]
    dd  <- d[pitcher == pid & is_swing == TRUE]
    ipan <- list(); iopt <- list()
    for (mu in unique(dd$matchup)) {
      ss <- dd[matchup == mu]; if (nrow(ss) < 15) next
      ref <- d[grp == "elite tunnel" & matchup == mu & is_swing == TRUE]
      rA <- fit_surface(ref$x_m_fb, ref$fb_plate_z, ref$miss_distance)
      rB <- fit_surface(ref$x_m_bb, ref$bb_plate_z, ref$miss_distance)
      fA <- fit_surface(ss$x_m_fb, ss$fb_plate_z, ss$miss_distance, prior = rA$fit)
      fB <- fit_surface(ss$x_m_bb, ss$bb_plate_z, ss$miss_distance, prior = rB$fit)
      pk <- sprintf("%s  (n=%d)", mu, nrow(ss))
      ipan[[length(ipan)+1]] <- rbindlist(list(
        data.table(gx=grid$gx, gz=grid$gz, val=fA$fit, pk=pk, surface="A: setup pitch location"),
        data.table(gx=grid$gx, gz=grid$gz, val=fB$fit, pk=pk, surface="B: breaking ball location")))
      aA <- argmax_cell(fA$fit, rA$dens); aB <- argmax_cell(fB$fit, rB$dens)
      iopt[[length(iopt)+1]] <- rbindlist(list(
        data.table(pk=pk, surface="A: setup pitch location", ox=aA$x, oz=aA$z),
        data.table(pk=pk, surface="B: breaking ball location", ox=aB$x, oz=aB$z)))
      all_opt[[length(all_opt)+1]] <- data.table(
        pairing=P$lab, scope="elite pitcher", matchup=mu, player_name=nm, pitcher=pid,
        n_swings=nrow(ss), pitcher_path=elite$path[i],
        setup_x=aA$x, setup_z=aA$z, setup_miss=aA$val,
        bb_x=aB$x, bb_z=aB$z, bb_miss=aB$val, bb_rng=aB$rng)
    }
    if (!length(ipan)) next
    q <- base_map(rbindlist(ipan), rbindlist(iopt)) + facet_grid(surface ~ pk) +
      labs(title = sprintf("%s - %s", nm, P$lab),
           subtitle = sprintf("Elite tunneler (path ratio %.2f, league-pairing median %.2f). EB-shrunk toward the elite-group prior.",
                              elite$path[i], median(pr$path)))
    ggsave(file.path(ODIR, sprintf("fig_pairing_%s_%d.png", P$id, pid)), q,
           width = 9.5, height = 7.4, dpi = 150)
    message("  wrote fig_pairing_", P$id, "_", pid, ".png  [", nm, "]")
  }

  ## elite leaderboard for this pairing
  cat(sprintf("\n===== %s : ELITE-TUNNEL ARMS (tightest path ratio, min %d pairs) =====\n",
              P$lab, MIN_PITCHER_N))
  print(head(pr[grp == "elite tunnel"][order(path),
    .(player_name, hand = p_throws, pairs = n, swings = n_sw, path_ratio = round(path,3))], 12))
}

opt <- rbindlist(all_opt)
fwrite(opt, file.path(ODIR, "ext_pairing_optima.csv"))

cat("\n============ GROUP OPTIMA: ELITE vs REST ============\n")
print(opt[scope %in% c("elite tunnel","rest"),
  .(pairing, group = scope, matchup, n_swings,
    setup = sprintf("(%+.1f,%+.1f) %.1f", setup_x, setup_z, setup_miss),
    breaking = sprintf("(%+.1f,%+.1f) %.1f", bb_x, bb_z, bb_miss),
    bb_spread = round(bb_rng,2))][order(pairing, matchup, group)])
cat(sprintf("\nwritten: %s\n", file.path(ODIR,"ext_pairing_optima.csv")))
