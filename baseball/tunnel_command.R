#!/usr/bin/env Rscript

# Does the optimal breaking-ball spot MOVE with where the setup fastball went, and
# how often does each pitcher actually get the breaking ball there?
#
#  1. conditional-target test : refit the BB miss surface inside each setup-FB zone and
#                               see whether the optimum shifts (formally: does letting
#                               the surface vary by FB zone beat one global surface on
#                               held-out 2026?)
#  2. optimal zone            : high-miss cells of the BB-location surface, fit on
#                               2023H2-2025 only so 2026 stays clean for validation
#  3. distribution            : where breaking balls actually go vs where they should
#  4. High-Miss%              : per-pitcher share of breaking balls landing in the
#                               high-miss zone, validated out-of-sample against 2026
#                               miss and whiff rate
#
# Handedness mirrored: + x = pitcher's glove side. Target E[miss | swing], inches.
# High-Miss% is a share of pitches landing in the high-miss zone, NOT a strike rate --
# that zone sits mostly below the strike zone.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

MDIR <- file.path("data","statcast_model")
ODIR <- file.path(MDIR, "tunnel_location")

p <- readRDS(file.path(MDIR, "tunnel_pairs.rds"))
p <- p[is_primary_setup == TRUE]
mir <- function(x, h) fifelse(h == "R", x, -x)
p[, `:=`(fb_x = mir(fb_plate_x, p_throws), bb_x = mir(bb_plate_x, p_throws))]
p[, `:=`(fb_z = fb_plate_z, bb_z = bb_plate_z)]
p[, `:=`(dx = bb_x - fb_x, dz = bb_z - fb_z)]
p[, matchup := fifelse(p_throws == stand, "same", "opp")]
p <- p[is_swing == FALSE | is.finite(miss_distance)]

GX <- seq(-2.0, 2.0, by = 0.10); GZ <- seq(-0.6, 4.4, by = 0.10)
grid <- as.data.table(expand.grid(gx = GX, gz = GZ)); NG <- nrow(grid)
H <- 0.35; K_LEAGUE <- 20; REL_SUP <- 0.12; CORE_SUP <- 0.25
ZONE_Q <- 0.80          # "high-miss zone" = top 20% of supported surface cells

## chunked kernel-weighted sums, so we never build a 2091 x 20000 matrix at once
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
  if (!any(core)) return(list(x=NA_real_, z=NA_real_, val=NA_real_))
  v <- fit; v[!core] <- NA_real_; i <- which.max(v)
  list(x = G$gx[i], z = G$gz[i], val = v[i])
}
cell_of <- function(x, z) {
  ix <- round((x - GX[1])/0.1) + 1; iz <- round((z - GZ[1])/0.1) + 1
  out <- rep(NA_integer_, length(x))
  ok <- ix >= 1 & ix <= length(GX) & iz >= 1 & iz <= length(GZ)
  out[ok] <- (iz[ok]-1)*length(GX) + ix[ok]
  out
}

TRAIN <- p[season <= 2025]; TEST <- p[season == 2026]
TYPES <- c("SL","CU","ST")
TYPE_LAB <- c(SL="Slider", CU="Curveball (CU+KC)", ST="Sweeper")

################################################################################
## 1. Does the optimal BB spot move with the setup fastball's location?
################################################################################
cat("========== 1. DOES THE TARGET MOVE WITH THE SETUP FASTBALL? ==========\n")
p[, fb_zone := paste0(fifelse(fb_z >= 2.9, "FB high", fifelse(fb_z >= 2.0, "FB mid", "FB low")),
                      " / ", fifelse(fb_x >= 0.35, "glove", fifelse(fb_x <= -0.35, "arm", "middle")))]
shift <- list()
for (ty in TYPES) {
  d <- p[brk_type == ty & is_swing == TRUE]
  g <- argmax_cell(fit_surface(d$bb_x, d$bb_z, d$miss_distance)$fit,
                   fit_surface(d$bb_x, d$bb_z, d$miss_distance)$dens)
  for (fz in unique(d$fb_zone)) {
    dd <- d[fb_zone == fz]; if (nrow(dd) < 400) next
    s <- fit_surface(dd$bb_x, dd$bb_z, dd$miss_distance)
    a <- argmax_cell(s$fit, s$dens)
    shift[[length(shift)+1]] <- data.table(brk_type = ty, fb_zone = fz, n = nrow(dd),
      opt_x = a$x, opt_z = a$z, opt_miss = a$val,
      global_x = g$x, global_z = g$z,
      dist_from_global = sqrt((a$x-g$x)^2 + (a$z-g$z)^2))
  }
}
shift <- rbindlist(shift)
print(shift[, .(brk_type, fb_zone, n, opt = sprintf("(%+.1f,%+.1f)", opt_x, opt_z),
                global = sprintf("(%+.1f,%+.1f)", global_x, global_z),
                shift_ft = round(dist_from_global,2))][order(brk_type, fb_zone)])
cat("\nspread of the optimal BB spot across setup-FB zones (ft):\n")
print(shift[, .(sd_opt_x = round(sd(opt_x),2), sd_opt_z = round(sd(opt_z),2),
                mean_shift = round(mean(dist_from_global),2)), by = brk_type])
fwrite(shift, file.path(ODIR, "target_shift_by_fb_zone.csv"))

################################################################################
## 2. Optimal zone, fit on TRAIN years only
################################################################################
cat("\n========== 2. HIGH-MISS ZONE (fit on 2023H2-2025) ==========\n")
# Two competing definitions of "the spot", so the choice can be settled by which one
# actually predicts future performance rather than by argument:
#   cond : top cells of E[miss | swing]      -- rewards balls nobody competitive swings at
#   pp   : top cells of E[miss per pitch]    -- chase-adjusted, P(swing) * E[miss|swing]
TRAIN[, miss_pp := fifelse(is_swing == TRUE, miss_distance, 0)]
p[, miss_pp := fifelse(is_swing == TRUE, miss_distance, 0)]

zones <- list(); surfs <- list()
for (ty in TYPES) for (mu in c("same","opp")) {
  trs <- TRAIN[brk_type == ty & matchup == mu & is_swing == TRUE]
  tra <- TRAIN[brk_type == ty & matchup == mu]
  if (nrow(trs) < 400) next
  sc <- fit_surface(trs$bb_x, trs$bb_z, trs$miss_distance)   # swing-conditional
  sp <- fit_surface(tra$bb_x, tra$bb_z, tra$miss_pp)         # chase-adjusted per pitch
  zc <- !is.na(sc$fit) & sc$fit >= quantile(sc$fit, ZONE_Q, na.rm = TRUE)
  zp <- !is.na(sp$fit) & sp$fit >= quantile(sp$fit, ZONE_Q, na.rm = TRUE)
  ac <- argmax_cell(sc$fit, sc$dens); ap <- argmax_cell(sp$fit, sp$dens)
  surfs[[paste(ty,mu)]] <- data.table(brk_type=ty, matchup=mu, gx=grid$gx, gz=grid$gz,
    fit=sc$fit, dens=sc$dens, inzone=zc, fit_pp=sp$fit, inzone_pp=zp)
  zones[[length(zones)+1]] <- data.table(brk_type=ty, matchup=mu, n_train=nrow(trs),
    cond_opt=sprintf("(%+.1f,%+.1f)", ac$x, ac$z), cond_cells=sum(zc),
    pp_opt=sprintf("(%+.1f,%+.1f)", ap$x, ap$z), pp_cells=sum(zp),
    overlap_pct=round(100*sum(zc & zp)/sum(zc|zp),1))
}
surfs <- rbindlist(surfs); zones <- rbindlist(zones)
saveRDS(surfs, file.path(ODIR, "zone_surfaces.rds"))   # shared with the sensitivity script
print(zones)

################################################################################
## 3. Where breaking balls actually go
################################################################################
cat("\n========== 3. HOW OFTEN DOES A BREAKING BALL LAND IN THE ZONE? ==========\n")
p[, cell := cell_of(bb_x, bb_z)]
zl <- surfs[, .(brk_type, matchup, cell = cell_of(gx, gz), inzone, inzone_pp)]
p <- merge(p, zl, by = c("brk_type","matchup","cell"), all.x = TRUE)
p[is.na(inzone), inzone := FALSE]; p[is.na(inzone_pp), inzone_pp := FALSE]

print(p[, .(pitches = .N,
            cond_zone_pct = round(100*mean(inzone),1),
            pp_zone_pct   = round(100*mean(inzone_pp),1),
            swing_in_cond = round(100*mean(is_swing[inzone]),1),
            swing_in_pp   = round(100*mean(is_swing[inzone_pp]),1),
            missPP_in_cond = round(mean(miss_pp[inzone]),2),
            missPP_in_pp   = round(mean(miss_pp[inzone_pp]),2),
            missPP_out     = round(mean(miss_pp[!inzone & !inzone_pp]),2)),
        by = brk_type][order(brk_type)])

cat("\nHigh-Miss% by where the setup fastball went (is the target harder off some FBs?):\n")
print(dcast(p[, .(pitches = .N, hit = round(100*mean(inzone_pp),1)),
              by = .(brk_type, fb_zone)], fb_zone ~ brk_type, value.var = "hit"))

################################################################################
## 4. Per-pitcher tunnel command + out-of-sample validation
################################################################################
cat("\n========== 4. TUNNEL COMMAND (per pitcher x breaking type) ==========\n")
MIN_TR <- 30; MIN_TE <- 20
cmd_tr <- p[season <= 2025, .(n_tr = .N, cmd_cond = mean(inzone), cmd_pp = mean(inzone_pp)),
            by = .(pitcher, player_name, brk_type)][n_tr >= MIN_TR]
out_te <- p[season == 2026, .(n_te = .N, n_sw = sum(is_swing),
                              missPP26 = mean(miss_pp),                       # per pitch
                              miss26   = mean(miss_distance[is_swing], na.rm=TRUE),
                              swstr26  = sum(is_whiff)/.N,                    # per pitch
                              whiff26  = sum(is_whiff)/pmax(sum(is_swing),1)),
            by = .(pitcher, brk_type)][n_te >= MIN_TE]
cmd_te <- p[season == 2026, .(cmd26_cond = mean(inzone), cmd26_pp = mean(inzone_pp)),
            by = .(pitcher, brk_type)]
v <- merge(merge(cmd_tr, out_te, by = c("pitcher","brk_type")), cmd_te, by = c("pitcher","brk_type"))

cat(sprintf("validation sample: %d pitcher x type units\n", nrow(v)))
cat("\nOUT-OF-SAMPLE r: 2023-25 command rate vs 2026 outcomes, BOTH zone definitions\n")
cr <- function(a, b) round(cor(a, b, use = "complete.obs"), 3)
print(rbindlist(list(
  data.table(zone = "cond  E[miss|swing]", n = nrow(v),
    vs_missPP26 = cr(v$cmd_cond, v$missPP26), vs_swstr26 = cr(v$cmd_cond, v$swstr26),
    vs_miss26 = cr(v$cmd_cond, v$miss26), vs_whiff26 = cr(v$cmd_cond, v$whiff26),
    reliability = cr(v$cmd_cond, v$cmd26_cond)),
  data.table(zone = "pp    E[miss/pitch]", n = nrow(v),
    vs_missPP26 = cr(v$cmd_pp, v$missPP26), vs_swstr26 = cr(v$cmd_pp, v$swstr26),
    vs_miss26 = cr(v$cmd_pp, v$miss26), vs_whiff26 = cr(v$cmd_pp, v$whiff26),
    reliability = cr(v$cmd_pp, v$cmd26_pp)))))

cat("\nby breaking type, chase-adjusted zone:\n")
print(v[, .(n = .N, vs_missPP26 = cr(cmd_pp, missPP26), vs_swstr26 = cr(cmd_pp, swstr26),
            reliability = cr(cmd_pp, cmd26_pp)), by = brk_type][order(brk_type)])

lead <- p[, .(n = .N, cmd = mean(inzone_pp), n_sw = sum(is_swing),
              missPP = mean(miss_pp), swstr = sum(is_whiff)/.N),
          by = .(pitcher, player_name, brk_type)][n >= 60]
setorder(lead, -cmd)
fwrite(lead, file.path(ODIR, "tunnel_command_leaderboard.csv"))
cat("\n===== TOP 15 TUNNEL COMMAND (chase-adjusted zone, min 60 breaking balls) =====\n")
print(head(lead[, .(player_name, brk_type, n, cmd_pct = round(100*cmd,1),
                    missPP = round(missPP,2), swstr_pct = round(100*swstr,1))], 15))
cat("\n===== BOTTOM 10 =====\n")
print(head(lead[order(cmd), .(player_name, brk_type, n, cmd_pct = round(100*cmd,1),
                    missPP = round(missPP,2), swstr_pct = round(100*swstr,1))], 10))

################################################################################
## 5. Within-pitcher vs between-pitcher: does hitting the spot MORE help YOU?
################################################################################
# The pitch-level effect is huge but the cross-pitcher correlation is ~0, which is the
# classic aggregation trap: arms who live in the chase zone differ in stuff and intent
# from arms who pitch in the zone. Demeaning by pitcher removes stuff and asks the
# development question directly -- when a given pitcher hits the spot more often than
# his own norm, does he miss more bats?
cat("\n========== 5. WITHIN- vs BETWEEN-PITCHER ==========\n")
panel <- p[, .(n = .N, zone = mean(inzone_pp), missPP = mean(miss_pp),
               swstr = sum(is_whiff)/.N),
           by = .(pitcher, player_name, brk_type, season)][n >= 25]
panel[, nseas := .N, by = .(pitcher, brk_type)]
pan <- panel[nseas >= 2]
pan[, `:=`(zone_c = zone - mean(zone), miss_c = missPP - mean(missPP),
           swstr_c = swstr - mean(swstr)), by = .(pitcher, brk_type)]
btw <- panel[, .(zone = mean(zone), missPP = mean(missPP), swstr = sum(swstr*n)/sum(n)),
             by = .(pitcher, brk_type)]

cat(sprintf("panel: %d pitcher x type x season rows, %d units with 2+ seasons\n",
            nrow(panel), uniqueN(pan[, .(pitcher, brk_type)])))
print(rbindlist(list(
  data.table(effect = "BETWEEN pitchers", n = nrow(btw),
             r_zone_vs_missPP = cr(btw$zone, btw$missPP),
             r_zone_vs_swstr  = cr(btw$zone, btw$swstr)),
  data.table(effect = "WITHIN pitcher (demeaned)", n = nrow(pan),
             r_zone_vs_missPP = cr(pan$zone_c, pan$miss_c),
             r_zone_vs_swstr  = cr(pan$zone_c, pan$swstr_c)))))

fe <- lm(missPP ~ zone + factor(paste(pitcher, brk_type)), data = pan)
cat(sprintf("\nfixed-effects slope (miss per pitch per +1.00 High-Miss%%): %.2f  (t = %.2f)\n",
            coef(fe)[["zone"]], summary(fe)$coefficients["zone","t value"]))
cat(sprintf("  -> +10 percentage points of High-Miss%% = %+.3f in of miss per pitch\n",
            0.10 * coef(fe)[["zone"]]))
cat("\nPITCH-LEVEL effect for reference (in zone vs out), miss per pitch:\n")
print(p[, .(in_zone = round(mean(miss_pp[inzone_pp]),2),
            out = round(mean(miss_pp[!inzone_pp]),2),
            ratio = round(mean(miss_pp[inzone_pp])/mean(miss_pp[!inzone_pp]),1)), by = brk_type])
fwrite(panel, file.path(ODIR, "tunnel_command_panel.csv"))

################################################################################
## figures
################################################################################
theme_set(theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face="bold"), panel.grid = element_blank(),
        strip.text = element_text(face="bold", size=10)))
SZ <- data.table(xmin=-0.83, xmax=0.83, ymin=1.6, ymax=3.4)

# Both six-panel figures use a pitch-type x handedness grid rather than a wrapped strip,
# so each pitch type's two platoon splits sit one above the other and can be read as a
# pair. A batter silhouette anchors which side of the plate is which.
source(file.path("baseball", "batter_silhouette.R"))
MLEV <- c("vs same-handed hitter", "vs opposite-handed hitter")
mlab_of <- function(m) factor(fifelse(m == "same", MLEV[1], MLEV[2]), levels = MLEV)
tlab_of <- function(t) factor(TYPE_LAB[t], levels = unname(TYPE_LAB[TYPES]))

XL <- c(-3.05, 3.05); YL <- c(-0.6, 6.2)
# one batter per panel, on the side the matchup puts him
batter_layer <- function() c(
  batter_layers(list(list(side = -1, data = data.table(mlab = mlab_of("same"))),
                     list(side =  1, data = data.table(mlab = mlab_of("opp"))))),
  list(geom_hline(yintercept = 0, colour = "grey35", linewidth = .4)))
zone_grid <- function() list(
  facet_grid(mlab ~ tlab, switch = "y"),
  scale_x_continuous(breaks = c(-1,0,1)),
  coord_fixed(xlim = XL, ylim = YL, expand = FALSE),
  theme(strip.text.y.left = element_text(angle = 90, face = "bold")))

sf <- copy(surfs)
sf[, `:=`(mlab = mlab_of(matchup), tlab = tlab_of(brk_type))]
g1 <- ggplot(sf, aes(gx, gz)) +
  geom_raster(aes(fill = fit), interpolate = TRUE) +
  geom_contour(aes(z = as.numeric(inzone)), breaks = 0.5, colour = "white", linewidth = .6) +
  geom_rect(data=SZ, inherit.aes=FALSE, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            fill=NA, colour="grey20", linetype="dashed", linewidth=.4) +
  batter_layer() + zone_grid() +
  scale_fill_viridis_c(option="inferno", na.value="white", name="E[miss | swing]\n(inches)") +
  labs(title = "The high-miss zone for a breaking ball after a fastball (2 strikes)",
       subtitle = "Fit on 2023H2-2025. White outline = top 20% of the surface. Mirrored: + x = glove side, so the hitter switches boxes between rows.",
       x = "mirrored plate_x (ft)", y = "plate_z (ft)")
ggsave(file.path(ODIR, "fig_highmiss_zone.png"), g1, width = 13.5, height = 10.0, dpi = 150)

sh <- copy(shift)
g2 <- ggplot(sh, aes(opt_x, opt_z)) +
  geom_rect(data=SZ, inherit.aes=FALSE, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            fill=NA, colour="grey60", linetype="dashed", linewidth=.4) +
  geom_hline(yintercept=0, colour="grey35", linewidth=.3) +
  # no batter here: this chart is about how tightly the dots cluster, and the wide limits
  # a silhouette needs shrink them to nothing
  geom_point(aes(colour = fb_zone), size = 3) +
  geom_point(aes(global_x, global_z), shape = 4, size = 5, stroke = 1.3, colour = "black") +
  facet_wrap(~ tlab_of(brk_type)) +
  coord_fixed(xlim = c(-1.6, 2.0), ylim = c(-0.2, 3.8)) +
  labs(title = "The optimal breaking-ball spot barely moves with the setup fastball",
       subtitle = "Coloured dot = optimum refit inside each setup-FB zone. Black X = one global optimum.",
       x = "mirrored plate_x (ft)", y = "plate_z (ft)", colour = "setup FB zone")
ggsave(file.path(ODIR, "fig_target_stability.png"), g2, width = 12, height = 5.8, dpi = 150)

g3 <- ggplot(v, aes(100*cmd_pp, missPP26)) +
  geom_point(aes(size = n_te), alpha = .5, colour = "#2c3e50") +
  geom_smooth(method = "lm", se = TRUE, colour = "#c0392b", linewidth = .8) +
  facet_wrap(~ brk_type, labeller = labeller(brk_type = TYPE_LAB), scales = "free_x") +
  scale_size_continuous(range = c(1,4), name = "2026 pitches") +
  labs(title = "High-Miss% does not carry over: hitting the spot in past years says nothing about future misses",
       subtitle = "2023-25 share of breaking balls landing in the chase-adjusted high-miss zone vs 2026 miss per pitch (226 pitcher x type units, out-of-sample r = -0.09)",
       x = "High-Miss%, 2023-25 (% of breaking balls landing in the high-miss zone)",
       y = "2026 miss distance per pitch (in)")
ggsave(file.path(ODIR, "fig_tunnel_command_validation.png"), g3, width = 12, height = 5, dpi = 150)

## ---- the offset framing, shown for completeness -----------------------------
# E[miss|swing] over the FB->BB separation vector itself, plus where separations
# actually land. The optimum here is a large downward gap, but the RMSE test
# (tunnel_delta_test.R) shows this is mostly a restatement of "the BB finished low".
OG <- as.data.table(expand.grid(gx = seq(-3, 3, 0.1), gz = seq(-4, 2, 0.1)))
osurf <- list(); oopt <- list()
for (ty in TYPES) {
  d <- p[brk_type == ty & is_swing == TRUE]
  s <- fit_surface(d$dx, d$dz, d$miss_distance, G = OG)
  a <- argmax_cell(s$fit, s$dens, G = OG)
  osurf[[ty]] <- data.table(brk_type = ty, gx = OG$gx, gz = OG$gz, fit = s$fit)
  oopt[[ty]]  <- data.table(brk_type = ty, ox = a$x, oz = a$z, val = a$val,
                            mx = mean(d$dx), mz = mean(d$dz))
}
osurf <- rbindlist(osurf); oopt <- rbindlist(oopt)

# Second row: take the best separation from the panel above, add it back to the fixed
# optimal BB spot, and chart the fastball location it implies onto the strike zone. The
# implied spot lands ABOVE the zone while the fastball surface's own best cell sits at
# the BOTTOM of it -- the two disagree, which is the tell that the separation optimum is
# just "put the breaking ball low" restated, and carries no fastball instruction.
fbsurf <- list(); fbmark <- list()
for (ty in TYPES) {
  d <- p[brk_type == ty & is_swing == TRUE]
  A <- fit_surface(d$fb_x, d$fb_z, d$miss_distance)
  B <- fit_surface(d$bb_x, d$bb_z, d$miss_distance)
  aA <- argmax_cell(A$fit, A$dens); aB <- argmax_cell(B$fit, B$dens)
  rngA <- diff(range(A$fit[!is.na(A$fit) & A$dens >= CORE_SUP]))
  rngB <- diff(range(B$fit[!is.na(B$fit) & B$dens >= CORE_SUP]))
  s <- oopt[brk_type == ty]
  fbsurf[[ty]] <- data.table(brk_type = ty, gx = grid$gx, gz = grid$gz, fit = A$fit,
                             lab = sprintf("E[miss] spans %.2f in across FB spots\n(%.2f in across BB spots)",
                                           rngA, rngB))
  fbmark[[ty]] <- data.table(brk_type = ty,
    x = c(aB$x - s$ox, aA$x, mean(d$fb_x)), z = c(aB$z - s$oz, aA$z, mean(d$fb_z)),
    lab = c("FB spot implied by the best separation", "best FB spot the surface itself finds",
            "average FB actually thrown"))
}
fbsurf <- rbindlist(fbsurf); fbmark <- rbindlist(fbmark)
MK <- c("FB spot implied by the best separation", "best FB spot the surface itself finds",
        "average FB actually thrown")
fbmark[, lab := factor(lab, levels = MK)]

FLIM <- c(0, 4)   # shared with the BB-location surface so the flatness is not a scale trick

g4 <- ggplot(osurf, aes(gx, gz)) +
  geom_raster(aes(fill = fit), interpolate = TRUE) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = .3) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = .3) +
  geom_point(data = oopt, aes(ox, oz), shape = 4, size = 3.5, stroke = 1.4, colour = "white") +
  geom_point(data = oopt, aes(mx, mz), shape = 21, size = 3, fill = "cyan", colour = "black") +
  scale_fill_viridis_c(option = "inferno", limits = FLIM, oob = scales::squish,
                       na.value = "white", guide = "none") +
  scale_x_continuous(breaks = c(-2,0,2)) +
  facet_wrap(~ tlab_of(brk_type)) +
  coord_fixed(xlim = c(-3,3), ylim = c(-4,2), expand = FALSE) +
  labs(title = "Why the setup fastball's location carries almost no information",
       subtitle = "Top: miss distance vs the FB-to-BB separation at the plate. White X = best separation, cyan dot = the league-average separation actually thrown.",
       x = "horizontal separation, BB - FB (ft)", y = "vertical separation, BB - FB (ft)")

g4b <- ggplot(fbsurf, aes(gx, gz)) +
  geom_raster(aes(fill = fit), interpolate = TRUE) +
  batter_layers(list(list(side = -1), list(side = 1))) +   # aggregated over matchup: both boxes
  geom_rect(data = SZ, inherit.aes = FALSE, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            fill = NA, colour = "white", linetype = "dashed", linewidth = .5) +
  geom_hline(yintercept = 0, colour = "grey35", linewidth = .4) +
  geom_point(data = fbmark, aes(x, z, shape = lab), colour = "white", size = 4, stroke = 1.5) +
  geom_text(data = unique(fbsurf[, .(brk_type, lab)]), inherit.aes = FALSE,
            aes(x = 0, y = 5.75, label = lab), size = 3.2, fontface = "bold", lineheight = .95) +
  scale_shape_manual(values = c(4, 2, 1), name = NULL, drop = FALSE) +
  scale_fill_viridis_c(option = "inferno", limits = FLIM, oob = scales::squish,
                       na.value = "white", name = "E[miss | swing] (inches)") +
  scale_x_continuous(breaks = c(-1,0,1)) +
  facet_wrap(~ tlab_of(brk_type)) +
  coord_fixed(xlim = XL, ylim = YL, expand = FALSE) +
  guides(fill = guide_colourbar(order = 1, barwidth = 12, barheight = .8),
         shape = guide_legend(order = 2, nrow = 3,
                              override.aes = list(colour = "black", size = 3.2))) +
  labs(subtitle = paste0("Bottom: the same target charted over where the SETUP FASTBALL was located, on the same colour scale. That surface is featureless.\n",
                         "The spot implied by the best separation (X) sits above the zone, nowhere near the spot the fastball surface itself prefers (triangle)."),
       x = "mirrored plate_x (ft)", y = "plate_z (ft)") +
  theme(legend.position = "bottom", legend.box = "horizontal")

png(file.path(ODIR, "fig_offset_surface.png"), width = 13.5, height = 12.4, units = "in", res = 150)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1,
                   heights = grid::unit(c(1, 1.34), "null"))))
print(g4,  vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(g4b, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
invisible(dev.off())

cat("\nsetup-FB location is near-irrelevant:\n")
print(merge(fbmark[, .(brk_type, lab, spot = sprintf("(%+.1f,%+.1f)", x, z))],
            unique(fbsurf[, .(brk_type, lab2 = gsub("\n", "  ", lab))]), by = "brk_type"))

## where breaking balls actually go, against the zone they should hit
dens <- p[brk_type %in% TYPES]
dens[, `:=`(mlab = mlab_of(matchup), tlab = tlab_of(brk_type))]
zoutline <- copy(surfs)[, `:=`(mlab = mlab_of(matchup), tlab = tlab_of(brk_type))]
hit <- p[, .(hit = 100*mean(inzone_pp)), by = .(brk_type, matchup)]
hit[, `:=`(mlab = mlab_of(matchup), tlab = tlab_of(brk_type),
           lab = sprintf("High-Miss%% = %.0f%%", hit))]

g5 <- ggplot(dens, aes(bb_x, bb_z)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", contour = TRUE, bins = 12) +
  geom_contour(data = zoutline, aes(gx, gz, z = as.numeric(inzone_pp)), breaks = 0.5,
               colour = "#e74c3c", linewidth = .8) +
  geom_rect(data = SZ, inherit.aes = FALSE, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            fill = NA, colour = "grey20", linetype = "dashed", linewidth = .4) +
  batter_layer() + zone_grid() +
  geom_label(data = hit, inherit.aes = FALSE, aes(x = 0, y = 5.85, label = lab),
             size = 3.4, fontface = "bold", label.size = 0, fill = "white", alpha = .8) +
  scale_fill_gradient(low = "#dfe6ec", high = "#34495e", name = "pitch density") +
  labs(title = "Where breaking balls actually go vs where they miss the most bats",
       subtitle = "Shading = actual 2-strike breaking-ball locations. Red outline = high-miss zone. Mirrored: + x = glove side.",
       x = "mirrored plate_x (ft)", y = "plate_z (ft)")
ggsave(file.path(ODIR, "fig_bb_distribution.png"), g5, width = 13.5, height = 10.0, dpi = 150)

cat("\nwrote fig_highmiss_zone.png, fig_target_stability.png, fig_tunnel_command_validation.png,\n")
cat("      fig_offset_surface.png, fig_bb_distribution.png\n")
print(oopt[, .(brk_type, best_sep = sprintf("(%+.1f,%+.1f)", ox, oz),
               actual_mean_sep = sprintf("(%+.1f,%+.1f)", mx, mz), best_miss = round(val,2))])
