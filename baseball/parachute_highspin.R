#!/usr/bin/env Rscript

# THE HIGH-ACTIVE-SPIN BIN: GATE OUT SEAM-SHIFTED CHANGEUPS INSTEAD OF MATCHING EFFICIENCY.
#
# The old bin required the changeup's active spin to sit within ten points of the four-seamer's.
# That is the wrong instrument for the idea. What we actually want to exclude is the changeup
# that gets its shape from a seam-shifted wake rather than from spin, because that pitch is
# doing something the fastball cannot imitate. A gap filter does not do that: it happily admits
# a 62-percent changeup paired with a 66-percent fastball, and it throws out a 99-percent
# changeup paired with a 98-percent fastball only if the arithmetic drifts.
#
# So the gate becomes a floor on BOTH pitches rather than a limit on the distance between them:
#
#     active spin >= 80 percent on the changeup AND on the four-seamer
#
# and the active-spin gap disappears from the definition entirely. With seam-shifted changeups
# already excluded, the axis gap can open up, which is the point - Alex Vesia runs a 99 percent
# changeup off a 98 percent four-seam but sits 24 to 32 degrees off it on the clock face, and
# the ten-degree rule was never going to find him.
#
# The swing floor is reported at both 60 and 40, because Vesia's changeup peaks at 48 swings in
# a season and the honest way to include that kind of pitcher is to say what the looser floor
# costs in precision rather than to quietly lower it.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
set.seed(1); options(width = 210)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")

L <- readRDS(file.path(MDIR, "parachute_ff_resid.rds"))
S <- L$CH[, .(np = .N, nsw = sum(is_swing), axis = mean(axis_diff), kill = -mean(az_diff),
              velo_sep = -mean(speed_diff), arm = mean(arm_angle, na.rm = TRUE),
              rv100 = 100*mean(rv_res)), by = .(pitcher, player_name, season)]
S <- merge(S, L$SW[, .(wh = 100*mean(wh_res)), by = .(pitcher, season)], by = c("pitcher","season"))
S <- merge(S, L$BP[, .(nbip = .N, gb = 100*mean(gb_res)), by = .(pitcher, season)],
           by = c("pitcher","season"))
AS <- readRDS(file.path(MDIR, "active_spin_long.rds"))
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, as_ch = active_spin)], by = c("pitcher","season"))
S <- merge(S, AS[pitch_type == "FF", .(pitcher, season, as_fb = active_spin)], by = c("pitcher","season"))
S <- S[nsw >= 40L & is.finite(axis) & is.finite(arm)]
S[, `:=`(last = sub(",.*", "", player_name), hs = as_ch >= .80 & as_fb >= .80)]

cat(sprintf("population at a 40-swing floor: %d pitcher-seasons (%d at 60)\n", nrow(S), sum(S$nsw >= 60)))
cat(sprintf("high active spin on both pitches: %d (%.0f%%)   ... and arm >= 44: %d\n\n",
            sum(S$hs), 100*mean(S$hs), sum(S$hs & S$arm >= 44)))

## ---- how far can the axis gap open? -----------------------------------------------------
# Every row is the same bin with one threshold moved. Effects are against everyone outside the
# bin in the same swing-floor population, so the comparison group widens as the bin does.
sweep <- function(D, arc = 44) rbindlist(lapply(seq(10, 45, 5), function(ax) {
  i <- D$hs & D$arm >= arc & D$axis <= ax
  if (sum(i) < 6) return(NULL)
  tt <- function(v) { t <- t.test(D[[v]][i], D[[v]][!i]); c(diff(rev(t$estimate)), t$p.value) }
  C <- D[i]; m <- lm(wh ~ velo_sep * i, D, weights = nsw); cf <- summary(m)$coefficients
  w <- tt("wh"); g <- tt("gb"); r <- tt("rv100")
  data.table(axis_max = ax, seasons = sum(i), pitchers = uniqueN(D$pitcher[i]),
             whiff = w[1], whiff_p = w[2], gb = g[1], gb_p = g[2], rv = r[1], rv_p = r[2],
             velo_r = cor(C$velo_sep, C$wh),
             slope_in = cf["velo_sep","Estimate"] + cf["velo_sep:iTRUE","Estimate"],
             inter_p = cf["velo_sep:iTRUE","Pr(>|t|)"]) }))

for (fl in c(60L, 40L)) {
  D <- S[nsw >= fl]
  cat(sprintf("=== axis-gap sweep, high-active-spin gate, arm >= 44, %d-swing floor ===\n", fl))
  W <- sweep(D)
  print(W[, .(axis_max, seasons, pitchers,
              whiff = sprintf("%+.2f (p=%.3f)", whiff, whiff_p),
              grounders = sprintf("%+.2f (p=%.3f)", gb, gb_p),
              rv100 = sprintf("%+.2f (p=%.3f)", rv, rv_p),
              velo_r = round(velo_r, 3),
              velo_slope = sprintf("%+.2f (p=%.4f)", slope_in, inter_p))], row.names = FALSE)
  cat("\n")
  if (fl == 60L) W60 <- W else W40 <- W
}

## ---- does dropping the efficiency filter actually change the membership? -----------------
cat("=== old gate vs new gate, at a 60-swing floor and a 15-degree axis limit ===\n")
D <- S[nsw >= 60]
D[, `:=`(old = abs(as_ch - as_fb) <= .15 & arm >= 42 & axis <= 15,
         new = hs & arm >= 44 & axis <= 15)]
cat(sprintf("  old (efficiency gap <= 15 pts, arm >= 42): %d seasons\n", sum(D$old)))
cat(sprintf("  new (both active spin >= 80%%, arm >= 44):  %d seasons\n", sum(D$new)))
cat(sprintf("  in both %d | only old %d | only new %d\n",
            sum(D$old & D$new), sum(D$old & !D$new), sum(!D$old & D$new)))
cat("\n  admitted ONLY by the new gate (the pitches the gap filter was throwing away):\n")
print(D[!old & new][order(-wh)][1:12, .(last, season, np, arm = round(arm,1), axis = round(axis,1),
      as_ch = round(as_ch,3), as_fb = round(as_fb,3), gap = round(as_ch-as_fb,3),
      velo_sep = round(velo_sep,1), wh = round(wh,1))], row.names = FALSE)
cat("\n  admitted ONLY by the old gate (kept despite a seam-shifted changeup):\n")
print(D[old & !new][order(as_ch)][1:8, .(last, season, np, arm = round(arm,1), axis = round(axis,1),
      as_ch = round(as_ch,3), as_fb = round(as_fb,3), wh = round(wh,1))], row.names = FALSE)

## ---- where does Vesia land? ---------------------------------------------------------------
cat("\n=== Alex Vesia under the new gate ===\n")
V <- S[grepl("Vesia", player_name)]
print(V[, .(season, np, nsw, arm = round(arm,1), axis = round(axis,1), as_ch = round(as_ch,3),
            as_fb = round(as_fb,3), velo_sep = round(velo_sep,1), wh = round(wh,1),
            gb = round(gb,1), rv100 = round(rv100,2),
            admits_at = sprintf("axis<=%d", 5*ceiling(axis/5)))], row.names = FALSE)

## ---- the roster at the widest defensible setting -------------------------------------------
AXC <- 30
R <- S[nsw >= 60 & hs & arm >= 44 & axis <= AXC][order(-wh)]
cat(sprintf("\n=== roster: both active spin >= 80%%, arm >= 44, axis gap <= %d, 60 swings (%d seasons) ===\n",
            AXC, nrow(R)))
print(R[, .(Pitcher = last, Season = season, CH = np, Arm = round(arm,1), Axis = round(axis,1),
            ActCH = round(as_ch,3), ActFF = round(as_fb,3), VeloSep = round(velo_sep,1),
            Whiff = round(wh,1), Grounders = round(gb,1), RV100 = round(rv100,2))][1:30],
      row.names = FALSE)
fwrite(S[nsw >= 40 & hs & arm >= 44 & axis <= AXC][order(-wh)],
       file.path(AST, "ext_parachute_highspin.csv"))

## ---- figure ---------------------------------------------------------------------------------
G <- rbind(cbind(W60, floor = "60-swing floor"), cbind(W40, floor = "40-swing floor"))
G[, floor := factor(floor, levels = c("60-swing floor","40-swing floor"))]
M <- melt(G, id.vars = c("axis_max","seasons","floor"),
          measure.vars = list(v = c("whiff","slope_in"), p = c("whiff_p","inter_p")),
          variable.name = "what")
M[, what := factor(what, labels = c("Whiff above model (pp)",
                                    "Velocity slope inside bin (pp per mph)"))]
M[, sig := fifelse(p < .05, "p < .05", "not significant")]

gg <- ggplot(M, aes(axis_max, v, colour = sig, group = floor)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_line(colour = "grey60", linewidth = .6) +
  geom_point(aes(size = seasons)) +
  facet_grid(what ~ floor, scales = "free_y", switch = "y") +
  scale_colour_manual(values = c("p < .05" = "#1d7870", "not significant" = "grey58"), name = NULL) +
  scale_size_area(max_size = 6, name = "seasons in bin") +
  scale_x_continuous(breaks = seq(10, 45, 5)) +
  labs(title = "Gating on high active spin instead of matched efficiency, then opening the axis gap",
       subtitle = paste0("The bin now requires 80 percent active spin on both the changeup and the four-seamer, which excludes seam-shifted changeups outright, and the\n",
                         "active-spin gap is gone from the definition. Each point moves the spin-axis limit one step wider. TOP: whiff above a shape-and-location model that\n",
                         "never sees the axis or the arm angle. BOTTOM: how much a mph of velocity separation is worth inside the bin, the one effect that has held up so far.\n",
                         "The left column keeps the 60-swing floor used everywhere else; the right drops it to 40, which is what it takes to admit relievers like Vesia."),
       x = "Spin-axis gap limit (degrees from the four-seamer)", y = NULL,
       caption = "Source: Statcast 2020-2026 - out-of-fold LightGBM residuals - four-seam anchor") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12.2), plot.subtitle = element_text(size = 8.2),
        strip.placement = "outside", strip.text = element_text(face = "bold", size = 9),
        panel.grid.minor = element_blank(), legend.position = "top",
        plot.caption = element_text(size = 7.5, colour = "grey40"))
ggsave(file.path(AST, "fig32_highspin_sweep.png"), gg, width = 11, height = 7, dpi = 150)
cat("\nwrote ext_parachute_highspin.csv, fig32_highspin_sweep.png\n")
