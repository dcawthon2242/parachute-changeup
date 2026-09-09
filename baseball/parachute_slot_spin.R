#!/usr/bin/env Rscript

# THE PARACHUTE CHANGEUP, REDEFINED: HIGH SLOT + VERY HIGH SPIN SIMILARITY. NO DROP FILTER.
#
# The IVB-kill requirement is gone. It was admitting near-sidearm pitchers whose changeups
# happen to sink and excluding Cease, who is the archetype. The definition is now the two
# things the mechanism actually names: the ball comes from over the top, and it spins like
# the fastball.
#
# A two-threshold bin invites picking the pair of cutoffs that works, so the first thing here
# is the entire grid of both thresholds rather than one cell. A real effect shows up as a
# broad plateau that strengthens as the bin tightens. A fitting artifact shows up as one or
# two isolated significant cells surrounded by noise. The grid is printed before any roster.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
set.seed(1); options(width = 210)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")

L <- readRDS(file.path(MDIR, "parachute_within.rds"))
S <- L$CH[, .(np = .N, nsw = sum(is_swing), axis = mean(axis_diff), kill = -mean(az_diff),
              velo_sep = -mean(speed_diff), spin = mean(release_spin_rate, na.rm = TRUE),
              rv100 = 100*mean(rv_res)), by = .(pitcher, player_name, season)]
S <- merge(S, L$SWD[, .(wh = 100*mean(wh_res)), by = .(pitcher, season)], by = c("pitcher","season"))
S <- merge(S, L$BIP[, .(gb = 100*mean(gb_res)), by = .(pitcher, season)], by = c("pitcher","season"))
aa <- fread(file.path(MDIR, "arm_angle_tunnel.csv"))[pitch_type == "CH",
        .(pitcher, season, arm = mean_arm, relz = mean_relz)]
S <- merge(S, unique(aa, by = c("pitcher","season")), by = c("pitcher","season"))
S <- S[nsw >= 60 & is.finite(axis) & is.finite(arm)]
cat(sprintf("population: %d pitcher-seasons, %d pitchers, min 60 changeup swings, 2023H2-2026\n\n",
            nrow(S), uniqueN(S$pitcher)))

## ---- 1. the full threshold grid ---------------------------------------------------
AX <- c(6, 8, 10, 12, 15, 20); AR <- c(40, 42, 44, 46, 48, 50)
G <- CJ(ax = AX, ar = AR)
G[, c("n","d_wh","p_wh","d_rv","p_rv") := {
  z <- mapply(function(a, r) {
    i <- S$axis <= a & S$arm >= r
    if (sum(i) < 8) return(c(sum(i), NA, NA, NA, NA))
    tw <- t.test(S$wh[i], S$wh[!i]); tr <- t.test(S$rv100[i], S$rv100[!i])
    c(sum(i), diff(rev(tw$estimate)), tw$p.value, diff(rev(tr$estimate)), tr$p.value)
  }, ax, ar)
  list(z[1,], z[2,], z[3,], z[4,], z[5,]) }]
cat("=== threshold grid: whiff% over model, bin minus rest (is this a plateau or a knife edge?) ===\n")
W <- dcast(G, ax ~ ar, value.var = "d_wh")
setnames(W, c("axis<=", paste0("arm>=", AR)))
print(W[, lapply(.SD, function(x) if (is.numeric(x)) round(x,2) else x)], row.names = FALSE)
cat("\n  p-values:\n")
P <- dcast(G, ax ~ ar, value.var = "p_wh"); setnames(P, c("axis<=", paste0("arm>=", AR)))
print(P[, lapply(.SD, function(x) if (is.numeric(x)) round(x,3) else x)], row.names = FALSE)
cat("\n  bin size n:\n")
N <- dcast(G, ax ~ ar, value.var = "n"); setnames(N, c("axis<=", paste0("arm>=", AR)))
print(N, row.names = FALSE)
cat(sprintf("\n  %d of %d cells reach p<.05; %d of %d have a POSITIVE whiff effect.\n",
            sum(G$p_wh < .05, na.rm=TRUE), sum(!is.na(G$p_wh)),
            sum(G$d_wh > 0, na.rm=TRUE), sum(!is.na(G$d_wh))))
cat(sprintf("  run value: %d of %d cells reach p<.05.\n", sum(G$p_rv < .05, na.rm=TRUE), sum(!is.na(G$p_rv))))

## ---- 2. the roster at the central definition --------------------------------------
AXC <- 10; ARC <- 44
S[, para := axis <= AXC & arm >= ARC]
cat(sprintf("\n=== ROSTER: arm angle >= %d deg and axis gap <= %d deg, n = %d ===\n", ARC, AXC, sum(S$para)))
print(S[para == TRUE][order(-wh), .(player_name, season, pitches = np, arm = round(arm,1),
        axis = round(axis,1), velo_sep = round(velo_sep,1), ivb_kill = round(kill,1),
        whiff_over = round(wh,1), gb_over = round(gb,1), rv100 = round(rv100,2))], row.names = FALSE)
cat("\n  bin vs everyone else:\n")
for (v in c("wh","gb","rv100")) { tt <- t.test(S[para==TRUE][[v]], S[para==FALSE][[v]])
  ci <- tt$conf.int
  cat(sprintf("    %-6s  bin %+.3f   rest %+.3f   diff %+.3f  95%% CI [%+.3f, %+.3f]  p = %.4f\n",
      v, mean(S[para==TRUE][[v]]), mean(S[para==FALSE][[v]]), diff(rev(tt$estimate)),
      ci[1], ci[2], tt$p.value)) }
cat(sprintf("    for reference the bin averages %.1f mph of separation and %.1f in of IVB kill,\n    against %.1f and %.1f for the rest - so this is NOT a big-drop group.\n",
    mean(S[para==TRUE]$velo_sep), mean(S[para==TRUE]$kill),
    mean(S[para==FALSE]$velo_sep), mean(S[para==FALSE]$kill)))

## ---- 3. out of sample --------------------------------------------------------------
cat("\n=== held out: the bin is defined only on 2023H2-2025, then scored on 2026 ===\n")
mem <- S[season <= 2025 & para == TRUE, unique(pitcher)]
T26 <- S[season == 2026]; T26[, prior := pitcher %chin% as.character(mem) | pitcher %in% mem]
tt <- t.test(T26[prior == TRUE]$wh, T26[prior == FALSE]$wh)
cat(sprintf("  2026 whiff over model: prior-bin pitchers %+.3f (n=%d)   everyone else %+.3f (n=%d)   p = %.3f\n",
    mean(T26[prior==TRUE]$wh), sum(T26$prior), mean(T26[prior==FALSE]$wh), sum(!T26$prior), tt$p.value))
tr <- t.test(T26[prior == TRUE]$rv100, T26[prior == FALSE]$rv100)
cat(sprintf("  2026 run value over model:                %+.3f            %+.3f            p = %.3f\n",
    mean(T26[prior==TRUE]$rv100), mean(T26[prior==FALSE]$rv100), tr$p.value))

## ---- 4. how much of it is the two obvious names? -----------------------------------
cat("\n=== leave-one-pitcher-out on the bin's whiff edge ===\n")
pl <- unique(S[para == TRUE]$player_name)
lo <- sapply(pl, function(q) { i <- S$para & S$player_name != q
  mean(S$wh[i]) - mean(S$wh[!S$para]) })
cat(sprintf("  full bin edge %+.2f pp; leave-one-out range %+.2f to %+.2f (most influential: %s)\n",
    mean(S[para==TRUE]$wh) - mean(S[para==FALSE]$wh), min(lo), max(lo), names(which.min(lo))))

fwrite(S[order(-para, -wh)], file.path(AST, "ext_parachute_slot_spin.csv"))

## ---- figures --------------------------------------------------------------------------
G[, c("d_gb","p_gb") := { z <- mapply(function(a, r) {
    i <- S$axis <= a & S$arm >= r
    if (sum(i) < 8) return(c(NA, NA))
    t <- t.test(S$gb[i], S$gb[!i]); c(diff(rev(t$estimate)), t$p.value) }, ax, ar)
  list(z[1,], z[2,]) }]
GG <- rbind(G[, .(ax, ar, n, d = d_wh, p = p_wh, out = "Whiff% over model")],
            G[, .(ax, ar, n, d = d_gb, p = p_gb, out = "Ground-ball% over model")])
GG <- GG[ax > 6]   # the 6-degree row never reaches eight pitchers
GG[, lab := fifelse(is.na(p), "-", sprintf("%+.1f\n%s", d,
        fifelse(p < .01, "p<.01", fifelse(p < .05, sprintf("p=%.2f", p), "n.s."))))]
gA <- ggplot(GG, aes(factor(ar), factor(ax), fill = d)) +
  geom_tile(colour = "white", linewidth = 1.1) +
  geom_text(aes(label = lab, fontface = fifelse(!is.na(p) & p < .05, "bold", "plain")),
            size = 2.9, lineheight = .95) +
  geom_point(data = GG[ax == AXC & ar == ARC], shape = 21, size = 19, stroke = 1.4,
             colour = "#c0392b", fill = NA) +
  scale_fill_gradient2(low = "#c0392b", mid = "grey93", high = "#1d7870", midpoint = 0,
                       na.value = "grey96", name = "Bin minus\nrest (pp)") +
  facet_wrap(~ out, nrow = 1) +
  labs(title = "Every way of drawing the bin, so you can see whether it is a plateau or a lucky cell",
       subtitle = paste0("Each cell is a different definition of the parachute changeup: the column is the arm-angle floor, the row is the spin-axis-gap ceiling, and the number is how many\n",
                         "percentage points that bin beats everyone else by, above a shape-and-location model. All 22 cells are positive on both outcomes and 11 of 22 reach p < .05 on each,\n",
                         "with the effect generally strengthening as the bin tightens - though not monotonically, and the top-left cells thin out to three to six pitchers. The circled cell is\n",
                         "the definition used for the roster. No drop filter is applied anywhere, which is what separates this from the earlier version of the bin."),
       x = "Arm angle floor (degrees)", y = "Spin-axis gap ceiling (degrees)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 8.1),
        panel.grid = element_blank(), strip.text = element_text(face = "bold", size = 10),
        panel.spacing.x = unit(14, "pt"))
ggsave(file.path(AST, "fig20_parachute_grid.png"), gA, width = 13.5, height = 6.6, dpi = 150)

S[, lab := fifelse(para, sub(",.*","", player_name), NA_character_)]
S[!is.na(lab), lab := fifelse(seq_len(.N) == which.max(np), lab, NA_character_), by = lab]
gB <- ggplot(S, aes(arm, axis)) +
  annotate("rect", xmin = ARC, xmax = Inf, ymin = -Inf, ymax = AXC, alpha = .12, fill = "#2a9d8f") +
  geom_hline(yintercept = AXC, linetype = "dashed", colour = "#c0392b", linewidth = .45) +
  geom_vline(xintercept = ARC, linetype = "dashed", colour = "#2a9d8f", linewidth = .45) +
  geom_point(aes(size = np, colour = wh), alpha = .8) +
  geom_text_repel(aes(label = lab), size = 3.05, max.overlaps = Inf, seed = 11, box.padding = .5,
                  min.segment.length = 0, segment.size = .3, segment.colour = "grey55",
                  colour = "grey10", fontface = "bold") +
  scale_colour_gradient2(low = "#c0392b", mid = "grey85", high = "#1d7870", midpoint = 0,
                         name = "Whiff% over model") +
  scale_size_continuous(range = c(.9, 4.5), guide = "none") +
  labs(title = "The parachute changeup: thrown from over the top, spinning like the fastball",
       subtitle = paste0("One point per pitcher-season, min 60 changeup swings, 2023H2-2026, coloured by whiff rate above a shape-and-location model that never sees the spin axis.\n",
                         "The bin requires an arm angle at or above the top quartile of changeup throwers and a spin axis within 10 degrees of the primary fastball. Vertical drop is\n",
                         "not part of the definition; this group averages ", sprintf("%.1f", mean(S[para==TRUE]$kill)),
                         " inches of IVB kill against ", sprintf("%.1f", mean(S[para==FALSE]$kill)), " for everyone else, so it is a spin-match group, not a sink group."),
       x = "Arm angle (degrees; 90 = overhead, 45 = three-quarters, 0 = sidearm)",
       y = "Mean spin-axis gap vs primary fastball (degrees)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 8.2),
        panel.grid.minor = element_blank())
ggsave(file.path(AST, "fig21_parachute_roster.png"), gB, width = 12.5, height = 7.2, dpi = 150)
cat("\nwrote fig20_parachute_grid.png, fig21_parachute_roster.png, ext_parachute_slot_spin.csv\n")
