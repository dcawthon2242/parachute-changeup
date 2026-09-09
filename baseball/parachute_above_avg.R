#!/usr/bin/env Rscript

# THE PARACHUTE BIN, REDEFINED AS FOUR "ABOVE AVERAGE" CUTS.
#
#   above-average arm slot
#   above-average spin-axis difference to the four-seamer
#   above-average four-seam active spin
#   above-average changeup active spin
#
# The second cut is ambiguous and the ambiguity reverses the bin, so both readings are built:
#
#   MATCHED  axis gap BELOW the mean - a better-than-average match, which is what "parachute"
#            has meant everywhere else in this project and what the tight bins screened for.
#   WIDE     axis gap ABOVE the mean - the literal reading of the words.
#
# One design point is not optional. The bin conditions on arm slot, and the residual model this
# is scored against must therefore contain arm angle. Every previous version of this analysis
# showed the high-slot cue clearing significance against models that lacked an arm-slot feature
# and collapsing the moment one was added, because the residual was simply carrying the part of
# whiff that arm angle explains. Scoring an arm-slot bin against an arm-blind model would
# manufacture the result. Both are reported so the gap between them stays visible.

suppressPackageStartupMessages({ library(data.table) })
options(width = 210)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")

R  <- readRDS(file.path(MDIR, "whiff_tjstuff.rds"))   # r3 = shape+location, r4 = + arm angle
AS <- readRDS(file.path(MDIR, "active_spin_long.rds"))

S <- R[, .(nsw = .N, axis = mean(axis_diff), arm = mean(arm_angle), velo_sep = -mean(speed_diff),
           no_arm = 100*mean(r3), with_arm = 100*mean(r4)), by = .(pitcher, player_name, season)]
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, as_ch = active_spin)], by = c("pitcher","season"))
S <- merge(S, AS[pitch_type == "FF", .(pitcher, season, as_fb = active_spin)], by = c("pitcher","season"))
S <- S[nsw >= 60]
S[, name := trimws(paste(sub(".*,\\s*","",player_name), sub(",.*","",player_name)))]

M <- S[, .(arm = mean(arm), axis = mean(axis), as_fb = mean(as_fb), as_ch = mean(as_ch))]
cat(sprintf("population: %d pitcher-seasons, %s changeup swings\n", nrow(S),
            format(sum(S$nsw), big.mark=",")))
cat(sprintf("league means -> arm slot %.1f deg | axis gap %.1f deg | FF active spin %.3f | CH active spin %.3f\n\n",
            M$arm, M$axis, M$as_fb, M$as_ch))

S[, base := arm > M$arm & as_fb > M$as_fb & as_ch > M$as_ch]
S[, `:=`(MATCHED = base & axis < M$axis, WIDE = base & axis > M$axis)]
cat(sprintf("three shared cuts alone: %d seasons\n", sum(S$base)))
cat(sprintf("  MATCHED (axis gap < %.1f): %d seasons\n", M$axis, sum(S$MATCHED)))
cat(sprintf("  WIDE    (axis gap > %.1f): %d seasons\n\n", M$axis, sum(S$WIDE)))

## ---- do either version overperform? -------------------------------------------------------
cat("=== whiff above model, bin vs everyone else (percentage points) ===\n")
out <- rbindlist(lapply(c("MATCHED","WIDE","base"), function(b) rbindlist(lapply(
  c(no_arm = "no_arm", with_arm = "with_arm"), function(v) {
    i <- S[[b]]; t <- t.test(S[[v]][i], S[[v]][!i])
    data.table(bin = b, model = if (v=="no_arm") "shape+location (arm blind)" else "shape+location+arm",
               n = sum(i), effect = diff(rev(t$estimate)),
               lo = t$conf.int[1], hi = t$conf.int[2], p = t$p.value) }))))
out[bin == "base", bin := "three cuts, axis free"]
print(out[, .(bin, model, n, effect = round(effect,2), lo = round(lo,2), hi = round(hi,2),
              p = round(p,4))], row.names = FALSE)

cat("\n=== the four cuts scored one at a time, against the arm-aware model ===\n")
print(rbindlist(lapply(list(c("arm slot > mean","arm > M$arm"), c("FF active spin > mean","as_fb > M$as_fb"),
      c("CH active spin > mean","as_ch > M$as_ch"), c("axis gap < mean (better match)","axis < M$axis")),
  function(x) { i <- S[, eval(parse(text = x[2]))]; t <- t.test(S$with_arm[i], S$with_arm[!i])
    data.table(cut = x[1], n = sum(i), effect = round(diff(rev(t$estimate)),2),
               p = round(t$p.value,4)) })), row.names = FALSE)

## ---- rosters ------------------------------------------------------------------------------
show <- function(b, lab) { cat(sprintf("\n=== %s: %d seasons, top 15 by whiff above the arm-aware model ===\n", lab, sum(S[[b]])))
  print(S[get(b) == TRUE][order(-with_arm)][1:min(15,.N), .(Pitcher = name, Season = season,
    Swings = nsw, Arm = round(arm,1), AxisGap = round(axis,1), ActFF = round(as_fb,2),
    ActCH = round(as_ch,2), VeloSep = round(velo_sep,1), WhiffAbove = round(with_arm,1),
    ArmBlind = round(no_arm,1))], row.names = FALSE) }
show("MATCHED", "MATCHED"); show("WIDE", "WIDE")

cat("\n=== velocity separation x bin (the one effect that has survived elsewhere) ===\n")
for (b in c("MATCHED","WIDE")) for (v in c("no_arm","with_arm")) {
  S[, bin := get(b)]; cf <- summary(lm(get(v) ~ velo_sep*bin, S, weights = nsw))$coefficients
  d <- cf["velo_sep:binTRUE",]
  cat(sprintf("  %-8s %-9s slope in %+.2f out %+.2f  difference %+.2f  p=%.4f  r_in=%+.3f\n",
      b, v, cf["velo_sep","Estimate"]+d[1], cf["velo_sep","Estimate"], d[1], d[4],
      cor(S[bin==TRUE]$velo_sep, S[bin==TRUE][[v]]))) }

## ---- how tight does the axis cut have to be before anything appears? -----------------------
#
# "Above average" admits 156 seasons, roughly one in seven of the population. Every effect this
# project has found lived in bins of about twenty. Sweeping the axis threshold down from the
# mean, with the other three cuts held at their means, shows whether the null above is the
# definition being too permissive or the effect not existing at any width.

cat("\n=== axis threshold swept down from the mean, other three cuts held at their means ===\n")
sw <- rbindlist(lapply(c(22.4, 20, 18, 16, 14, 12, 10, 8), function(ax) {
  S[, bin := base & axis < ax]
  if (sum(S$bin) < 6) return(NULL)
  t <- t.test(S$with_arm[S$bin], S$with_arm[!S$bin])
  cf <- summary(lm(with_arm ~ velo_sep*bin, S, weights = nsw))$coefficients
  d <- cf["velo_sep:binTRUE",]
  data.table(axis_max = ax, seasons = sum(S$bin), whiff = diff(rev(t$estimate)), p_whiff = t$p.value,
             velo_diff = d[1], p_velo = d[4],
             r_velo = cor(S[bin==TRUE]$velo_sep, S[bin==TRUE]$with_arm)) }))
print(sw[, .(axis_max, seasons, whiff = round(whiff,2), p_whiff = round(p_whiff,3),
             velo_slope_extra = round(velo_diff,2), p_velo = round(p_velo,4),
             r_velo = round(r_velo,3))], row.names = FALSE)

# The arm cut is the other loose one: the population mean slot is well below the 44-degree line
# the earlier bins used. Re-run the sweep on the stricter slot to separate the two.
cat("\n=== same sweep, but arm slot >= 44 instead of above the 37.9 mean ===\n")
S[, base44 := arm >= 44 & as_fb > M$as_fb & as_ch > M$as_ch]
sw2 <- rbindlist(lapply(c(22.4, 20, 18, 16, 14, 12, 10, 8), function(ax) {
  S[, bin := base44 & axis < ax]
  if (sum(S$bin) < 6) return(NULL)
  t <- t.test(S$with_arm[S$bin], S$with_arm[!S$bin])
  cf <- summary(lm(with_arm ~ velo_sep*bin, S, weights = nsw))$coefficients
  d <- cf["velo_sep:binTRUE",]
  data.table(axis_max = ax, seasons = sum(S$bin), whiff = diff(rev(t$estimate)), p_whiff = t$p.value,
             velo_diff = d[1], p_velo = d[4],
             r_velo = cor(S[bin==TRUE]$velo_sep, S[bin==TRUE]$with_arm)) }))
print(sw2[, .(axis_max, seasons, whiff = round(whiff,2), p_whiff = round(p_whiff,3),
              velo_slope_extra = round(velo_diff,2), p_velo = round(p_velo,4),
              r_velo = round(r_velo,3))], row.names = FALSE)

## ---- figure ---------------------------------------------------------------------------------
suppressPackageStartupMessages(library(ggplot2))
G <- rbind(cbind(sw, slot = "Arm slot above the 37.9 mean"),
           cbind(sw2, slot = "Arm slot >= 44"))
P <- rbind(G[, .(axis_max, seasons, slot, panel = "Whiff above model (pp)", v = whiff, p = p_whiff)],
           G[, .(axis_max, seasons, slot, panel = "Velocity slope inside the bin (r)", v = r_velo, p = p_velo)])
P[, panel := factor(panel, levels = c("Whiff above model (pp)", "Velocity slope inside the bin (r)"))]
gg <- ggplot(P, aes(axis_max, v, colour = slot)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_line(linewidth = .7) +
  geom_point(aes(size = seasons)) +
  geom_text(data = P[p < .05], aes(label = "p<.05"), vjust = -1.5, size = 2.7, show.legend = FALSE) +
  # The significance labels sit above the topmost points, so the panel needs headroom or they
  # are clipped at the strip.
  geom_blank(aes(y = v * 1.18)) +
  facet_wrap(~panel, scales = "free_y") +
  scale_x_reverse(breaks = c(22.4,20,18,16,14,12,10,8)) +
  scale_size_continuous(range = c(1.4, 5), name = "seasons in bin") +
  scale_colour_manual(values = c("Arm slot above the 37.9 mean" = "#1d7870",
                                 "Arm slot >= 44" = "#b4632a"), name = NULL) +
  labs(title = "\"Above average\" is far too loose a cut. The parachute effect is a tail phenomenon, not a median one",
       subtitle = paste0("Both panels hold three of the four cuts fixed - arm slot, four-seam active spin and changeup active spin all above their league means - and tighten only the\n",
                         "spin-axis gap, moving right to left. At the mean gap of 22.4 degrees the bin holds 156 pitcher-seasons and is dead flat: whiff above model is -0.04 points\n",
                         "and the velocity correlation inside the bin is +0.004. Both quantities rise monotonically as the gap closes and only separate from zero below about 12\n",
                         "degrees, roughly the 15th percentile of the axis distribution rather than the 50th. The gradient being smooth is the reassuring part, since a single\n",
                         "isolated spike would read as noise. The caution is sample size: the bins that clear p < .05 hold 8 to 18 seasons, these p-values are uncorrected for the\n",
                         "eight thresholds swept, and the thresholds were chosen after seeing the data. Tightening the slot to 44 degrees does not add anything the axis cut has\n",
                         "not already delivered - it mostly just shrinks the bin."),
       x = "Spin-axis gap ceiling (degrees, tightening to the right)", y = NULL,
       caption = "Source: Statcast 2020-2026 - four-seam anchor - residuals from a shape, location and arm-angle model - 60+ changeup swings per season") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11.5), plot.subtitle = element_text(size = 8.2),
        strip.text = element_text(face = "bold", size = 9.5), panel.grid.minor = element_blank(),
        legend.position = "top", plot.caption = element_text(size = 7.5, colour = "grey40"))
ggsave(file.path(AST, "fig36_above_avg_sweep.png"), gg, width = 11.5, height = 6.6, dpi = 150)

fwrite(S[MATCHED | WIDE][order(-with_arm)], file.path(AST, "ext_parachute_above_avg.csv"))
cat("\nwrote ext_parachute_above_avg.csv, fig36_above_avg_sweep.png\n")
