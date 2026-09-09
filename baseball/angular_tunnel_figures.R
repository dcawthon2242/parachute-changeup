#!/usr/bin/env Rscript

# Figures explaining what the angular tunnel metric actually measures.
#
# The metric lives entirely in the hitter's visual field, so the figures are built in that
# frame rather than in field coordinates. Everything is recomputed from the stored
# trajectories so the plots show the real quantities, not a schematic.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
options(width = 200)
AST <- "data/statcast_model/article_assets"
dir.create(AST, showWarnings = FALSE, recursive = TRUE)
R_BALL <- (2.9/2)/12          # baseball radius, feet
DEG <- 180/pi
TAU <- 0.05

d <- readRDS("data/statcast_model/angular_tunnel_2026.rds")
d[, grp := fifelse(pitch_type == "SL", "Slider",
            fifelse(pitch_type == "ST", "Sweeper",
            fifelse(pitch_type %in% c("CU","KC"), "Curveball",
            fifelse(pitch_type == "CH", "Changeup",
            fifelse(pitch_type == "FS", "Splitter", NA_character_)))))]

pos <- function(r, v, a, t) r + v*t + 0.5*a*t^2

# Angular coordinates as seen from the eye. Forward points from the eye to the mound, so
# azimuth and elevation are the hitter's own left-right and up-down.
ang_frame <- function(row, nstep = 240) {
  E <- c(row$eye_x, row$eye_y, row$eye_z)
  fwd <- c(0, 1, 0)                      # toward the pitcher
  rgt <- c(1, 0, 0); upv <- c(0, 0, 1)
  Tmax <- max(row$t_react, row$p_t_react)
  tk <- seq(0, Tmax, length.out = nstep)
  one <- function(w) {
    p <- if (w == "current")
      cbind(pos(row$release_pos_x, row$vx0, row$ax, tk),
            pos(row$release_pos_y, row$vy0, row$ay, tk),
            pos(row$release_pos_z, row$vz0, row$az, tk))
    else
      cbind(pos(row$p_release_pos_x, row$p_vx0, row$p_ax, tk),
            pos(row$p_release_pos_y, row$p_vy0, row$p_ay, tk),
            pos(row$p_release_pos_z, row$p_vz0, row$p_az, tk))
    v <- sweep(p, 2, E)
    f <- as.numeric(v %*% fwd); rr <- as.numeric(v %*% rgt); uu <- as.numeric(v %*% upv)
    dist <- sqrt(rowSums(v^2))
    list(az = atan2(rr, f)*DEG, el = atan2(uu, f)*DEG,
         rho = asin(pmin(1, R_BALL/dist))*DEG, dist = dist)
  }
  A <- one("current"); P <- one("previous")
  w <- data.table(t = tk,
                  az_current = A$az, el_current = A$el, rho_current = A$rho, dist_current = A$dist,
                  az_previous = P$az, el_previous = P$el, rho_previous = P$rho, dist_previous = P$dist)
  # angular separation of the two centres, then edge-to-edge and size-difference cues
  w[, theta := {
    a1 <- az_current/DEG; e1 <- el_current/DEG
    a2 <- az_previous/DEG; e2 <- el_previous/DEG
    v1 <- cbind(sin(a1)*cos(e1), cos(a1)*cos(e1), sin(e1))
    v2 <- cbind(sin(a2)*cos(e2), cos(a2)*cos(e2), sin(e2))
    acos(pmin(1, pmax(-1, rowSums(v1*v2))))*DEG }]
  w[, `:=`(g_pos = theta - rho_current - rho_previous,
           g_diam = abs(2*rho_current - 2*rho_previous))]
  w[, frac := t/max(t)]
  w[]
}

## ============ pick a matched pair of examples from one pitcher ============
# Same pitcher, same two pitch types, one well tunneled and one not, so the contrast
# cannot be explained by who threw it.
cand <- d[grp == "Changeup" & consecutive == TRUE & is.finite(brk_any_005)]
who <- cand[, .N, by = .(pitcher, player_name)][N >= 60][order(-N)][1]
sub <- cand[pitcher == who$pitcher]
good <- sub[which.max(brk_any_005)]
bad  <- sub[which.min(brk_any_005)]
cat(sprintf("example pitcher: %s (%d changeup-after-fastball pairs)\n", who$player_name, who$N))
cat(sprintf("  tunneled pair    break fraction %.3f\n", good$brk_any_005))
cat(sprintf("  separated pair   break fraction %.3f\n", bad$brk_any_005))

G <- ang_frame(good)[, lab := sprintf("tunnelled\nbreak fraction %.2f", good$brk_any_005)]
B <- ang_frame(bad)[,  lab := sprintf("separates early\nbreak fraction %.2f", bad$brk_any_005)]
EX <- rbind(G, B)
EX[, lab := factor(lab, levels = c(unique(G$lab), unique(B$lab)))]

circ <- function(cx, cy, r, n = 120) {
  th <- seq(0, 2*pi, length.out = n)
  data.table(x = cx + r*cos(th), y = cy + r*sin(th))
}
COLS <- c("fastball (previous pitch)" = "#3B6FA0", "changeup (current pitch)" = "#C0392B")

## ---- FIG 0: how the angle is defined -----------------------------------
# Left: the two flight paths seen from above, with the sight lines that define the angle.
# Right: the same instant blown up to the scale the eye works at.
tsnap <- 0.60
gi <- which.min(abs(G$frac - tsnap))
traj <- rbindlist(lapply(c("current","previous"), function(w) {
  tk <- seq(0, max(good$t_react, good$p_t_react), length.out = 200)
  px <- if (w == "current") pos(good$release_pos_x, good$vx0, good$ax, tk)
        else                pos(good$p_release_pos_x, good$p_vx0, good$p_ax, tk)
  py <- if (w == "current") pos(good$release_pos_y, good$vy0, good$ay, tk)
        else                pos(good$p_release_pos_y, good$p_vy0, good$p_ay, tk)
  data.table(which = w, depth = py, lat = px)
}))
traj[, which := factor(which, levels = c("previous","current"), labels = names(COLS))]
ball <- rbindlist(lapply(c("current","previous"), function(w) {
  tt <- G$t[gi]
  data.table(which = w,
    depth = if (w=="current") pos(good$release_pos_y, good$vy0, good$ay, tt)
            else              pos(good$p_release_pos_y, good$p_vy0, good$p_ay, tt),
    lat   = if (w=="current") pos(good$release_pos_x, good$vx0, good$ax, tt)
            else              pos(good$p_release_pos_x, good$p_vx0, good$p_ax, tt))
}))
ball[, which := factor(which, levels = c("previous","current"), labels = names(COLS))]
sight <- ball[, .(which, x = c(good$eye_y), y = c(good$eye_x), xe = depth, ye = lat)]

fA <- ggplot() +
  geom_path(data = traj, aes(depth, lat, colour = which), linewidth = .9) +
  geom_segment(data = sight, aes(x = x, y = y, xend = xe, yend = ye),
               colour = "grey35", linewidth = .4, linetype = "22") +
  geom_point(data = ball, aes(depth, lat, colour = which), size = 3) +
  geom_point(aes(x = good$eye_y, y = good$eye_x), shape = 21, size = 3.2,
             fill = "black", colour = "black") +
  annotate("text", x = good$eye_y, y = good$eye_x + .12, label = "hitter's eye",
           hjust = .5, vjust = 0, size = 3.1) +
  annotate("text", x = 28, y = max(traj$lat) + .45,
           label = "the two sight lines differ by theta", hjust = .5, size = 3.1,
           colour = "grey25") +
  scale_colour_manual(values = COLS, name = NULL) +
  scale_x_reverse() +
  labs(title = "A. The angle, in field coordinates",
       subtitle = paste0("Seen from above, zoomed onto the flight paths. Trajectories come from the\n",
                         "Statcast nine-parameter fit; the eye point comes from OpenBiomechanics\n",
                         "stance data, scaled to this hitter's listed height. At this instant the two\n",
                         "sight lines differ by less than a fifth of a degree, which is why the angle\n",
                         "has to be measured rather than eyeballed."),
       x = "distance from home plate (ft)", y = "lateral position (ft)") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# Right panel: the same instant in angular units, centred on the midpoint of the two discs.
z <- G[gi]
mx <- mean(c(z$az_current, z$az_previous)); my <- mean(c(z$el_current, z$el_previous))
dd <- rbindlist(list(
  circ(z$az_previous - mx, z$el_previous - my, z$rho_previous)[, which := names(COLS)[1]],
  circ(z$az_current  - mx, z$el_current  - my, z$rho_current )[, which := names(COLS)[2]]))
ctr <- data.table(which = names(COLS),
                  x = c(z$az_previous - mx, z$az_current - mx),
                  y = c(z$el_previous - my, z$el_current - my))
LIM <- max(abs(c(dd$x, dd$y)))
fB <- ggplot() +
  geom_polygon(data = dd, aes(x, y, fill = which), alpha = .55) +
  annotate("segment", x = ctr$x[1], y = ctr$y[1], xend = ctr$x[2], yend = ctr$y[2],
           linewidth = .5, colour = "grey20") +
  geom_point(data = ctr, aes(x, y), size = 1.1, colour = "grey20") +
  annotate("text", x = 0, y = LIM*1.55,
           label = sprintf("theta = %.3f deg  (centre to centre)", z$theta), size = 3.1) +
  annotate("text", x = 0, y = -LIM*1.35,
           label = sprintf("edge-to-edge gap = theta - rho_A - rho_B = %+.3f deg", z$g_pos),
           size = 3.1, colour = "#1F4E79") +
  annotate("text", x = 0, y = -LIM*1.75,
           label = sprintf("angular size difference = %.3f deg", z$g_diam),
           size = 3.1, colour = "#8A6D00") +
  scale_fill_manual(values = COLS, guide = "none") +
  coord_equal(xlim = c(-2.1*LIM, 2.1*LIM), ylim = c(-2.1*LIM, 2.1*LIM)) +
  labs(title = "B. The same instant, in the visual field",
       subtitle = paste0("Each ball subtends a disc of angular radius rho = asin(1.45 in / distance).\n",
                         "A negative edge-to-edge gap means the two discs still overlap on the retina:\n",
                         "the hitter has no positional evidence that this is not the fastball."),
       x = "degrees", y = "degrees") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

ggsave(file.path(AST, "angtun_fig0_geometry.png"),
       fA + fB + plot_layout(widths = c(1.25, 1)), width = 12.5, height = 5.4, dpi = 150)

## ---- FIG 1: what the hitter sees, over the window ----------------------
# Centred on the midpoint of the pair so the separation is legible; the axis scale is
# shared across every panel, so disc growth and separation are directly comparable.
snaps <- c(0.10, 0.30, 0.50, 0.70, 0.90, 1.00)
disc <- rbindlist(lapply(levels(EX$lab), function(L) {
  zz <- EX[lab == L]
  rbindlist(lapply(snaps, function(s) {
    i <- which.min(abs(zz$frac - s))
    mx <- mean(c(zz$az_current[i], zz$az_previous[i]))
    my <- mean(c(zz$el_current[i], zz$el_previous[i]))
    rbindlist(lapply(c("previous","current"), function(w) {
      cc <- circ(zz[[paste0("az_", w)]][i] - mx, zz[[paste0("el_", w)]][i] - my,
                 zz[[paste0("rho_", w)]][i])
      cc[, `:=`(which = if (w == "previous") names(COLS)[1] else names(COLS)[2],
                snap = sprintf("%.0f%% of window", 100*s), lab = L, id = paste(w, s))]
    }))
  }))
}))
disc[, snap := factor(snap, levels = sprintf("%.0f%% of window", 100*snaps))]
disc[, lab := factor(lab, levels = levels(EX$lab))]
gap <- rbindlist(lapply(levels(EX$lab), function(L) {
  zz <- EX[lab == L]
  rbindlist(lapply(snaps, function(s) {
    i <- which.min(abs(zz$frac - s))
    data.table(lab = L, snap = sprintf("%.0f%% of window", 100*s),
               txt = sprintf("edge gap %+.2f\nsize diff %.2f", zz$g_pos[i], zz$g_diam[i]),
               sep = zz$g_pos[i] > TAU | zz$g_diam[i] > TAU)
  }))
}))
gap[, snap := factor(snap, levels = levels(disc$snap))]
gap[, lab := factor(lab, levels = levels(EX$lab))]

f1 <- ggplot(disc, aes(x, y, group = id, fill = which)) +
  geom_polygon(alpha = .6, colour = NA) +
  geom_text(data = gap, aes(x = 0, y = 0.78, label = txt, colour = sep),
            inherit.aes = FALSE, size = 2.7, fontface = "bold", vjust = 1, lineheight = .95) +
  facet_grid(lab ~ snap, switch = "y") +
  scale_fill_manual(values = COLS, name = NULL) +
  scale_colour_manual(values = c("FALSE" = "grey45", "TRUE" = "#B7410E"), guide = "none") +
  coord_equal(xlim = c(-.85, .85), ylim = c(-.85, 1.0)) +
  labs(title = "What the hitter's eye actually receives",
       subtitle = paste0(who$player_name,
         " - the two pitches as angular discs, centred on their midpoint, same scale in every panel.\n",
         "Discs grow because the ball is getting closer. Orange numbers mark moments when at least one cue ",
         "clears 0.05 deg; the break\nfraction is the first such moment. Separation is not monotone - the ",
         "bottom pair is already distinguishable at release, from release point alone."),
       x = "degrees of visual angle", y = "degrees of visual angle") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 8),
        strip.placement = "outside", panel.grid.minor = element_blank())
ggsave(file.path(AST, "angtun_fig1_visual_field.png"), f1, width = 11.5, height = 5.4, dpi = 150)

## ---- FIG 2: the two cues over the decision window ----------------------
long <- melt(EX[, .(frac, lab, `edge-to-edge angular gap` = g_pos,
                    `angular size difference (depth cue)` = g_diam)],
             id.vars = c("frac","lab"), variable.name = "cue", value.name = "deg")
brkpt <- EX[, {
  i <- which(g_pos > TAU | g_diam > TAU)
  .(frac = if (length(i)) frac[min(i)] else 1) }, by = lab]

f2 <- ggplot(long, aes(frac, deg, colour = cue)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = .3) +
  geom_hline(yintercept = TAU, linetype = "dashed", colour = "grey35") +
  annotate("text", x = .02, y = TAU, label = "acuity threshold tau = 0.05 deg",
           hjust = 0, vjust = -.6, size = 2.9, colour = "grey35") +
  geom_line(linewidth = .8) +
  geom_vline(data = brkpt, aes(xintercept = frac), colour = "#B7410E", linewidth = .6) +
  geom_text(data = brkpt, aes(x = frac, y = Inf, label = "tunnel ends"), inherit.aes = FALSE,
            hjust = -.08, vjust = 1.6, size = 2.9, colour = "#B7410E") +
  facet_wrap(~ lab) +
  scale_colour_manual(values = c("edge-to-edge angular gap" = "#1F4E79",
                                 "angular size difference (depth cue)" = "#C99700"), name = NULL) +
  labs(title = "The two cues, over the hitter's decision window",
       subtitle = paste0("Negative edge-to-edge gap means the discs physically overlap on the retina.\n",
                         "The tunnel ends the first moment EITHER cue clears the acuity threshold; ",
                         "the break fraction is that x-position."),
       x = "share of the decision window elapsed (release to commit point)",
       y = "degrees") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())
ggsave(file.path(AST, "angtun_fig2_cues.png"), f2, width = 10.5, height = 4.8, dpi = 150)

## ---- FIG 3: distribution of the break fraction by pitch type ----------
D3 <- d[!is.na(grp)]
ord <- D3[, .(m = mean(brk_any_005)), by = grp][order(m)]$grp
D3[, grp := factor(grp, levels = ord)]
f3 <- ggplot(D3, aes(brk_any_005, grp, fill = grp)) +
  geom_violin(scale = "width", colour = NA, alpha = .75) +
  stat_summary(fun = mean, geom = "point", size = 2, colour = "black") +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(title = "How long each pitch type stays indistinguishable from the fastball",
       subtitle = paste0("Break fraction at tau = 0.05 deg. 0 = separable the instant it leaves the hand; ",
                         "1 = still fused at the commit point.\nBlack dot is the mean. ",
                         "Sweepers and curveballs give themselves away earliest."),
       x = "break fraction", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
ggsave(file.path(AST, "angtun_fig3_distribution.png"), f3, width = 9.5, height = 4.6, dpi = 150)

## ---- FIG 4/5: the chase effect under the real specification ------------
# Reproduces the hardest specification from angular_tunnel_chase_robust.R exactly, so the
# figures cannot disagree with the tables.
cf <- list.files("data/statcast_2026/chunks", pattern = "csv$", full.names = TRUE)
ex <- unique(rbindlist(lapply(cf, function(f) fread(f, select = c(
  "game_pk","at_bat_number","pitch_number","sz_top","sz_bot","release_speed",
  "plate_x","plate_z"), showProgress = FALSE))), by = c("game_pk","at_bat_number","pitch_number"))
prev <- ex[, .(game_pk, at_bat_number, pitch_number = pitch_number + 1L,
               pp_x = plate_x, pp_z = plate_z)]
ex[, c("plate_x","plate_z") := NULL]
KEY <- c("game_pk","at_bat_number","pitch_number")
D4 <- merge(merge(d[!is.na(grp)], ex, by = KEY, all.x = TRUE), prev, by = KEY, all.x = TRUE)
D4[, `:=`(plate_sep = sqrt((plate_x - pp_x)^2 + (plate_z - pp_z)^2),
          zdist = sqrt(pmax(abs(plate_x) - 0.95, 0)^2 +
                       pmax(plate_z - sz_top, sz_bot - plate_z, 0)^2),
          chase = as.numeric(swing))]
D4[, pit_count := paste(pitcher, balls, strikes)]
D4 <- D4[zdist > 0 & is.finite(zdist)]
D4[, fam := fifelse(grp %in% c("Changeup","Splitter"), "Offspeed (changeup, splitter)",
             fifelse(grp %in% c("Slider","Sweeper"), "Horizontal (slider, sweeper)", "Curveball"))]

LOC  <- c("plate_x","plate_z","zdist","release_speed")
HARD <- c(LOC, "plate_sep", "pp_x", "pp_z")

# Partial out the fixed effect and the controls from both sides. By Frisch-Waugh the slope
# through the residuals equals the coefficient the regression reports.
resid_on <- function(dt, yv, ctl, fe) {
  y <- as.numeric(dt[[yv]]); g <- dt[[fe]]
  y <- y - ave(y, g)
  X <- vapply(ctl, function(v) { z <- as.numeric(dt[[v]]); z - ave(z, g) }, numeric(nrow(dt)))
  as.numeric(lm.fit(cbind(1, X), y)$residuals)
}
clust_slope <- function(yr, xr, cl) {
  X <- cbind(1, xr); XtXi <- solve(crossprod(X)); b <- XtXi %*% crossprod(X, yr)
  e <- as.vector(yr - X %*% b); meat <- matrix(0, 2, 2)
  for (ix in split(seq_along(yr), cl)) {
    u <- crossprod(X[ix, , drop = FALSE], e[ix]); meat <- meat + u %*% t(u) }
  V <- XtXi %*% meat %*% XtXi; nc <- length(unique(cl)); V <- V*(nc/(nc-1))
  c(b = b[2], se = sqrt(V[2,2]))
}

D4 <- D4[complete.cases(D4[, c("chase","brk_any_005", HARD, "pit_count","pitcher"), with = FALSE])]
D4[, `:=`(cr = resid_on(.SD, "chase", HARD, "pit_count"),
          br = resid_on(.SD, "brk_any_005", HARD, "pit_count")),
   by = fam, .SDcols = c("chase","brk_any_005", HARD, "pit_count")]
D4[, bin := cut(br, breaks = quantile(br, seq(0,1,length.out = 11), na.rm = TRUE),
                include.lowest = TRUE, labels = FALSE), by = fam]
B4 <- D4[!is.na(bin), .(x = mean(br), y = 100*mean(cr), n = .N,
                        se = 100*sd(cr)/sqrt(.N)), by = .(fam, bin)]
slope <- D4[, as.list(clust_slope(cr, br, pitcher)), by = fam]
slope[, txt := sprintf("slope %+.1f pp per full window (SE %.1f)", 100*b, 100*se)]
B4 <- merge(B4, slope[, .(fam, txt)], by = "fam")

f4 <- ggplot(B4, aes(x, y, colour = fam)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = .3) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = .3) +
  geom_errorbar(aes(ymin = y - 1.96*se, ymax = y + 1.96*se), width = 0, alpha = .45) +
  geom_point(size = 2.2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = .7, formula = y ~ x) +
  geom_text(aes(x = -Inf, y = Inf, label = txt), hjust = -.04, vjust = 1.6,
            size = 2.9, colour = "grey20", check_overlap = TRUE) +
  facet_wrap(~ fam) +
  scale_colour_brewer(palette = "Dark2", guide = "none") +
  labs(title = "Chase rate against tunnel length, under the full specification",
       subtitle = paste0("Out-of-zone pitches after a fastball, decile bins. Both axes are residuals after ",
                         "removing pitcher-by-count fixed effects,\nlocation, release speed, the plate ",
                         "separation between the two pitches, and the setup fastball's own location.\n",
                         "The fitted slope is the coefficient the regression reports; SEs are clustered by pitcher."),
       x = "break fraction (residual)", y = "chase rate (residual, pp)") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())
ggsave(file.path(AST, "angtun_fig4_chase_effect.png"), f4, width = 11.5, height = 4.8, dpi = 150)

## ---- FIG 5: the same coefficient as the controls get harder -----------
# This is the whole reason the horizontal result is reported as a null: the raw effect is
# real but it is plate separation wearing a costume.
SPECS <- list(
  list(nm = "1. pitcher only",                 ctl = character(0),                fe = "pitcher"),
  list(nm = "2. + location, velocity",         ctl = LOC,                         fe = "pitcher"),
  list(nm = "3. + plate separation",           ctl = c(LOC,"plate_sep"),          fe = "pitcher"),
  list(nm = "4. + setup pitch location",       ctl = c(LOC,"plate_sep","pp_x","pp_z"), fe = "pitcher"),
  list(nm = "5. + pitcher-by-count effects",   ctl = HARD,                        fe = "pit_count"))
lad <- rbindlist(lapply(SPECS, function(s) D4[, {
  yr <- resid_on(.SD, "chase", s$ctl, s$fe)
  xr <- resid_on(.SD, "brk_any_005", s$ctl, s$fe)
  z <- clust_slope(yr, xr, pitcher)
  .(spec = s$nm, b = 100*z[1], se = 100*z[2])
}, by = fam, .SDcols = c("chase","brk_any_005", HARD, "pitcher","pit_count")]))
lad[, spec := factor(spec, levels = rev(vapply(SPECS, `[[`, "", "nm")))]

f5 <- ggplot(lad, aes(b, spec, colour = fam)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = .4) +
  geom_errorbar(aes(xmin = b - 1.96*se, xmax = b + 1.96*se), orientation = "y",
                width = 0, linewidth = .6) +
  geom_point(size = 2.4) +
  facet_wrap(~ fam) +
  scale_colour_brewer(palette = "Dark2", guide = "none") +
  labs(title = "Why the horizontal breaking balls are reported as a null",
       subtitle = paste0("Chase effect of a full-window tunnel, in percentage points, as controls are added. ",
                         "95% intervals, clustered by pitcher.\nEvery family starts out looking huge. Adding ",
                         "location and release speed removes about four fifths of it everywhere, and for\n",
                         "sliders and sweepers it removes essentially all of it. The offspeed effect survives ",
                         "every control that follows."),
       x = "chase effect of a full-window tunnel (pp)", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())
ggsave(file.path(AST, "angtun_fig5_ladder.png"), f5, width = 12, height = 4.2, dpi = 150)

cat("\n=== fig 5: coefficient ladder (pp per full window) ===\n")
print(lad[order(fam, spec), .(fam, spec, b = round(b,2), se = round(se,2),
                              t = round(b/se,2))], row.names = FALSE)

cat("\nwrote:\n")
for (f in c("angtun_fig0_geometry.png","angtun_fig1_visual_field.png","angtun_fig2_cues.png",
            "angtun_fig3_distribution.png","angtun_fig4_chase_effect.png","angtun_fig5_ladder.png"))
  cat(" ", file.path(AST, f), "\n")
