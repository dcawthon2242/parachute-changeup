#!/usr/bin/env Rscript

# PHASE 2a - Fit the two miss-distance location surfaces.
#
#   Surface A: setup-FB plate location -> E[breaking-ball miss | swing]
#              (all qualifying primary-FB -> BB 2-strike pairs)
#   Surface B: breaking ball's OWN plate location -> E[miss | swing]
#              (restricted to WELL-TUNNELED pairs = tightest tercile of path_ratio)
#
# Two views:
#   league  : GAM te(x,z) per brk_type x pitcher-hand x batter-side  (the EB prior)
#   pitcher : Gaussian-kernel local mean shrunk toward the league prior
#             posterior = (n_eff*local + K*prior) / (n_eff + K)
#
# Surfaces are masked to the region with real data support so the argmax
# ("optimal location") can't land on an extrapolated grid corner.

suppressPackageStartupMessages({ library(data.table); library(mgcv) })

MDIR <- file.path("data","statcast_model")
ODIR <- file.path(MDIR, "tunnel_location")
dir.create(ODIR, showWarnings = FALSE, recursive = TRUE)

pairs <- readRDS(file.path(MDIR, "tunnel_pairs.rds"))
pairs <- pairs[is_primary_setup == TRUE]   # Phase 2 = primary fastball as the setup

# TARGET: E[miss | swing], consistent with the article's miss model. Note this is a
# swing-conditional quantity, so its high ground is the chase region -- swings at balls
# well out of the zone miss by a lot. The chase-adjusted expected miss per pitch thrown
# (miss_pp, = P(swing) * E[miss|swing]) is carried alongside as a diagnostic.
sw <- pairs[is_swing == FALSE | is.finite(miss_distance)]
sw[, miss_pp := fifelse(is_swing == TRUE, miss_distance, 0)]

## ---- grid over the plate region (catcher's view, feet) --------------------
# Grid runs well past the zone so the chase-adjusted optimum can turn over inside
# the domain rather than getting pinned to a boundary.
GX <- seq(-2.0, 2.0, by = 0.10)
GZ <- seq(-0.6, 4.4, by = 0.10)   # floor below 0 on purpose: ~10% of well-tunneled
                                  # curveballs are bounced, and the chase-adjusted
                                  # peak sits just above the dirt
grid <- as.data.table(expand.grid(gx = GX, gz = GZ))
NG   <- nrow(grid)

H        <- 0.35  # kernel bandwidth (ft)
K_LEAGUE <- 25    # league shrink toward the cell's global mean
K_SHRINK <- 10    # EB strength: effective pitches needed to move off the prior
MIN_SUP  <- 40    # absolute: mask cells with fewer than this many effective pitches
REL_SUP  <- 0.12  # relative: also mask cells below this fraction of peak support.
                  # Absolute support alone lets sparse distribution tails (e.g. a
                  # fastball at plate_z 0.8) survive and win the argmax on noise.

MIN_A <- 40       # min pooled pitches for a per-pitcher Surface A unit
MIN_B <- 15       # min pooled well-tunneled pitches for a per-pitcher Surface B unit

## kernel weight matrix: NG x n
kern <- function(x, z) {
  d2 <- outer(grid$gx, x, "-")^2 + outer(grid$gz, z, "-")^2
  exp(-0.5 * d2 / H^2)
}

## League surface: kernel local mean shrunk to the global mean, masked to real
## support. Non-parametric on purpose -- a GAM extrapolates past the data hull and
## invents impossible miss values in the grid corners.
## Primary surface = E[miss | swing], estimated over SWINGS only, so its support mask
## is driven by swing counts. Per-pitch miss and swing rate ride along as diagnostics.
fit_league <- function(x, z, y, swing) {
  sel <- which(swing)
  Ws  <- kern(x[sel], z[sel]); ns <- rowSums(Ws)
  cond <- as.numeric((Ws %*% y[sel]) / pmax(ns, 1e-9))
  cond <- (ns * cond + K_LEAGUE * mean(y[sel])) / (ns + K_LEAGUE)

  Wa  <- kern(x, z); na_ <- rowSums(Wa)
  pp  <- as.numeric((Wa %*% fifelse(swing, y, 0)) / pmax(na_, 1e-9))
  pp  <- (na_ * pp + K_LEAGUE * mean(fifelse(swing, y, 0))) / (na_ + K_LEAGUE)
  swr <- as.numeric((Wa %*% as.numeric(swing)) / pmax(na_, 1e-9))

  # Adaptive absolute floor: a fixed 40-effective-swing requirement wipes out the
  # smallest platoon cells entirely once the target is swing-conditional.
  min_abs <- max(12, min(MIN_SUP, 0.02 * length(sel)))
  ok <- ns >= min_abs & ns >= REL_SUP * max(ns)
  cond[!ok] <- NA_real_; pp[!ok] <- NA_real_; swr[!ok] <- NA_real_
  list(fit = cond, dens = ns / max(ns), swrate = swr, pp = pp, ok = ok)
}

## per-pitcher: kernel local mean shrunk toward the prior surface
fit_pitcher <- function(x, z, y, prior) {
  W     <- kern(x, z)
  n_eff <- rowSums(W)
  local <- as.numeric((W %*% y) / pmax(n_eff, 1e-9))
  post  <- (n_eff * local + K_SHRINK * prior) / (n_eff + K_SHRINK)
  post[is.na(prior)] <- NA_real_
  list(fit = post, dens = n_eff / max(n_eff), n_eff = n_eff)
}

CORE_SUP <- 0.25  # argmax only over well-sampled "core" cells, so the reported
                  # optimum is a real spot rather than a hot pixel on the rim

## argmax over core-supported cells; also reports how much the surface varies at all
## (spread = how much location matters, which is the sensitivity signal for Phase 3)
argmax_cell <- function(fit, dens) {
  core <- !is.na(fit) & dens >= CORE_SUP
  if (!any(core)) core <- !is.na(fit)
  if (!any(core)) return(list(x=NA_real_, z=NA_real_, val=NA_real_, rng=NA_real_, sd=NA_real_))
  v <- fit; v[!core] <- NA_real_
  i <- which.max(v)
  list(x = grid$gx[i], z = grid$gz[i], val = v[i],
       rng = diff(range(v, na.rm = TRUE)), sd = sd(v, na.rm = TRUE))
}

CELLS <- unique(sw[, .(brk_type, p_throws, stand)])[order(brk_type, p_throws, stand)]

league <- list(); optima <- list()

message("Fitting league surfaces (the EB prior)...")
for (i in seq_len(nrow(CELLS))) {
  ce <- CELLS[i]
  dA <- sw[brk_type == ce$brk_type & p_throws == ce$p_throws & stand == ce$stand]
  dB <- dA[well_tunneled == TRUE]
  if (nrow(dA) < 100 || nrow(dB) < 60) next

  sA <- fit_league(dA$fb_plate_x, dA$fb_plate_z, dA$miss_distance, dA$is_swing)
  sB <- fit_league(dB$bb_plate_x, dB$bb_plate_z, dB$miss_distance, dB$is_swing)

  league[[length(league)+1]] <- data.table(
    brk_type = ce$brk_type, p_throws = ce$p_throws, stand = ce$stand,
    gx = grid$gx, gz = grid$gz,
    fitA = sA$fit, densA = sA$dens, swA = sA$swrate, ppA = sA$pp,
    fitB = sB$fit, densB = sB$dens, swB = sB$swrate, ppB = sB$pp)

  aA <- argmax_cell(sA$fit, sA$dens); aB <- argmax_cell(sB$fit, sB$dens)
  # same argmax on the chase-adjusted per-pitch surface, for side-by-side comparison
  pA <- argmax_cell(sA$pp, sA$dens);  pB <- argmax_cell(sB$pp, sB$dens)
  optima[[length(optima)+1]] <- data.table(
    scope = "league", pitcher = NA_integer_, player_name = "LEAGUE",
    brk_type = ce$brk_type, p_throws = ce$p_throws, stand = ce$stand,
    n_A = sum(dA$is_swing), n_B = sum(dB$is_swing),
    fb_opt_x = aA$x, fb_opt_z = aA$z, fb_opt_miss = aA$val,
    fb_rng = aA$rng, fb_sd = aA$sd,
    bb_opt_x = aB$x, bb_opt_z = aB$z, bb_opt_miss = aB$val,
    bb_rng = aB$rng, bb_sd = aB$sd,
    fb_pp_x = pA$x, fb_pp_z = pA$z, fb_pp_miss = pA$val,
    bb_pp_x = pB$x, bb_pp_z = pB$z, bb_pp_miss = pB$val)
}
league <- rbindlist(league)
saveRDS(league, file.path(ODIR, "league_surfaces.rds"))
message(sprintf("  league cells fit: %d", uniqueN(league[, .(brk_type,p_throws,stand)])))

## ---- per-pitcher surfaces (pooled years, displayed for 2026 arms) ---------
act26 <- unique(sw[season == 2026, .(pitcher, brk_type)])
poolA <- sw[act26, on = .(pitcher, brk_type)]
qualA <- poolA[, .(N = sum(is_swing)), by = .(pitcher, brk_type)][N >= MIN_A]
poolA <- poolA[qualA[, .(pitcher, brk_type)], on = .(pitcher, brk_type)]

message(sprintf("Fitting per-pitcher surfaces for %d qualifying pitcher x type units...",
                nrow(qualA)))
pitch_surf <- list()
units <- unique(poolA[, .(pitcher, player_name, brk_type, p_throws)])
for (u in seq_len(nrow(units))) {
  uu <- units[u]
  for (st in c("L","R")) {
    pri <- league[brk_type == uu$brk_type & p_throws == uu$p_throws & stand == st]
    if (!nrow(pri)) next
    dA <- poolA[pitcher == uu$pitcher & brk_type == uu$brk_type & stand == st]
    dB <- dA[well_tunneled == TRUE]
    if (sum(dA$is_swing) < 8) next

    sA <- dA[is_swing == TRUE]; sB <- dB[is_swing == TRUE]
    fA <- fit_pitcher(sA$fb_plate_x, sA$fb_plate_z, sA$miss_distance, pri$fitA)
    fB <- if (nrow(sB) >= MIN_B)
            fit_pitcher(sB$bb_plate_x, sB$bb_plate_z, sB$miss_distance, pri$fitB)
          else list(fit = pri$fitB, dens = rep(NA_real_, NG))

    pitch_surf[[length(pitch_surf)+1]] <- data.table(
      pitcher = uu$pitcher, player_name = uu$player_name, brk_type = uu$brk_type,
      p_throws = uu$p_throws, stand = st, n_A = nrow(sA), n_B = nrow(sB),
      gx = grid$gx, gz = grid$gz,
      fitA = fA$fit, densA = fA$dens, fitB = fB$fit, densB = fB$dens)

    aA <- argmax_cell(fA$fit, pri$densA); aB <- argmax_cell(fB$fit, pri$densB)
    pA <- list(x=NA_real_, z=NA_real_, val=NA_real_); pB <- pA
    optima[[length(optima)+1]] <- data.table(
      scope = "pitcher", pitcher = uu$pitcher, player_name = uu$player_name,
      brk_type = uu$brk_type, p_throws = uu$p_throws, stand = st,
      n_A = nrow(sA), n_B = nrow(sB),
      fb_opt_x = aA$x, fb_opt_z = aA$z, fb_opt_miss = aA$val,
      fb_rng = aA$rng, fb_sd = aA$sd,
      bb_opt_x = aB$x, bb_opt_z = aB$z, bb_opt_miss = aB$val,
      bb_rng = aB$rng, bb_sd = aB$sd,
      fb_pp_x = pA$x, fb_pp_z = pA$z, fb_pp_miss = pA$val,
      bb_pp_x = pB$x, bb_pp_z = pB$z, bb_pp_miss = pB$val)
  }
  if (u %% 50 == 0) message(sprintf("   ...%d/%d units", u, nrow(units)))
}
pitch_surf <- rbindlist(pitch_surf)
saveRDS(pitch_surf, file.path(ODIR, "pitcher_surfaces.rds"))

optima <- rbindlist(optima)
fwrite(optima, file.path(ODIR, "ext_tunnel_optima.csv"))

cat("\n================ LEAGUE OPTIMA ================\n")
print(optima[scope == "league",
  .(brk_type, p_throws, stand, n_A, n_B,
    FB_opt = sprintf("(%.1f, %.1f)", fb_opt_x, fb_opt_z), FB_miss = round(fb_opt_miss,2),
    BB_opt = sprintf("(%.1f, %.1f)", bb_opt_x, bb_opt_z), BB_miss = round(bb_opt_miss,2))])

cat("\n======== HOW MUCH DOES LOCATION MATTER? (core-cell spread, inches) ========\n")
cat("  rng = best-minus-worst E[miss|swing] across well-sampled cells\n")
print(optima[scope == "league", .(
  FB_setup_rng = round(mean(fb_rng, na.rm=TRUE),3), FB_setup_sd = round(mean(fb_sd, na.rm=TRUE),3),
  BB_loc_rng   = round(mean(bb_rng, na.rm=TRUE),3), BB_loc_sd   = round(mean(bb_sd, na.rm=TRUE),3),
  ratio_BB_over_FB = round(mean(bb_rng, na.rm=TRUE)/mean(fb_rng, na.rm=TRUE),2)),
  by = brk_type][order(brk_type)])

cat("\n==== SAME OPTIMA UNDER BOTH OBJECTIVES (swing-conditional vs chase-adjusted) ====\n")
print(optima[scope == "league", .(brk_type, matchup = paste0(p_throws,"HP/",stand,"HB"),
  BB_cond = sprintf("(%+.1f,%+.1f) %.1f", bb_opt_x, bb_opt_z, bb_opt_miss),
  BB_perpitch = sprintf("(%+.1f,%+.1f) %.2f", bb_pp_x, bb_pp_z, bb_pp_miss))])
cat(sprintf("\nper-pitcher surface rows: %d  (%d pitcher x type x side panels)\n",
            nrow(pitch_surf), uniqueN(pitch_surf[, .(pitcher,brk_type,stand)])))
cat(sprintf("optima written: %s\n", file.path(ODIR,"ext_tunnel_optima.csv")))
