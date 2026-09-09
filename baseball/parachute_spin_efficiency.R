#!/usr/bin/env Rscript

# DOES THE PARACHUTE BIN ACCOUNT FOR SPIN EFFICIENCY? IT DOES NOT. THIS ASKS WHETHER IT SHOULD.
#
# The bin is arm angle >= 44 deg plus a spin-axis gap <= 10 deg from the primary fastball.
# spin_axis is the clock-face angle of rotation projected onto the plane facing the catcher.
# Gyro spin points along the direction of travel and therefore does not appear in that
# projection at all, so two changeups can match a fastball's axis to within a degree while one
# is 95 percent transverse and the other is half gyro. The axis criterion is silent about it.
#
# That matters for this specific claim, because "same spin, slower" is a statement about the
# whole spin vector, not its shadow. Four questions:
#
#   1. How much efficiency mismatch is the bin actually admitting?
#   2. Does the axis effect survive controlling for the active-spin gap?
#   3. Does efficiency match predict anything on its own?
#   4. Does requiring BOTH sharpen the bin or just shrink it?
#
# Active spin here is Baseball Savant's measured leaderboard value, not the value inferred
# from movement, which an earlier pass in this project found unreliable.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
set.seed(1); options(width = 205)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")
FASTBALLS <- c("FF","SI","FC")

L <- readRDS(file.path(MDIR, "parachute_within.rds"))
S <- L$CH[, .(np = .N, nsw = sum(is_swing), axis = mean(axis_diff), kill = -mean(az_diff),
              velo_sep = -mean(speed_diff), spin = mean(release_spin_rate, na.rm = TRUE),
              rv100 = 100*mean(rv_res)), by = .(pitcher, player_name, season)]
S <- merge(S, L$SWD[, .(wh = 100*mean(wh_res)), by = .(pitcher, season)], by = c("pitcher","season"))
S <- merge(S, L$BIP[, .(gb = 100*mean(gb_res)), by = .(pitcher, season)], by = c("pitcher","season"))
aa <- fread(file.path(MDIR, "arm_angle_tunnel.csv"))[pitch_type == "CH", .(pitcher, season, arm = mean_arm)]
S <- merge(S, unique(aa, by = c("pitcher","season")), by = c("pitcher","season"))

# Primary fastball per pitcher-season, same FF > SI > FC rule used to build the axis gap.
d <- readRDS(file.path(MDIR, "parachute_rv.rds"))
fbt <- d[pitch_type %in% FASTBALLS, .N, by = .(pitcher, season, pitch_type)][N >= 50]
fbt[, rk := match(pitch_type, FASTBALLS)]
fbt <- fbt[order(pitcher, season, rk)][, .SD[1], by = .(pitcher, season)][, .(pitcher, season, fb = pitch_type)]

AS <- readRDS(file.path(MDIR, "active_spin_long.rds"))
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, as_ch = active_spin)],
           by = c("pitcher","season"), all.x = TRUE)
S <- merge(S, fbt, by = c("pitcher","season"), all.x = TRUE)
S <- merge(S, AS[, .(pitcher, season, fb = pitch_type, as_fb = active_spin)],
           by = c("pitcher","season","fb"), all.x = TRUE)
S <- S[nsw >= 60 & is.finite(axis) & is.finite(arm)]
S[, `:=`(as_gap = as_ch - as_fb, eff_ok = abs(as_ch - as_fb) <= 0.10)]
S[, para := axis <= 10 & arm >= 44]
cat(sprintf("population %d pitcher-seasons; %d (%.0f%%) have measured active spin for both the\nchangeup and the primary fastball\n\n",
            nrow(S), sum(!is.na(S$as_gap)), 100*mean(!is.na(S$as_gap))))
A <- S[is.finite(as_gap)]

## ---- 1. what is the bin admitting? -----------------------------------------------
cat("=== 1. efficiency spread inside the bin ===\n")
B <- A[para == TRUE]
cat(sprintf("  bin      n=%2d  active spin: CH %.2f  FB %.2f  gap %+.3f  (SD %.3f, range %+.2f to %+.2f)\n",
            nrow(B), mean(B$as_ch), mean(B$as_fb), mean(B$as_gap), sd(B$as_gap),
            min(B$as_gap), max(B$as_gap)))
R <- A[para == FALSE]
cat(sprintf("  rest     n=%3d  active spin: CH %.2f  FB %.2f  gap %+.3f  (SD %.3f)\n",
            nrow(R), mean(R$as_ch), mean(R$as_fb), mean(R$as_gap), sd(R$as_gap)))
tt <- t.test(B$as_gap, R$as_gap)
cat(sprintf("  bin vs rest on efficiency gap: %+.3f, p = %.3f\n", diff(rev(tt$estimate)), tt$p.value))
cat(sprintf("  %d of %d bin members are efficiency-MISMATCHED by more than 10 points:\n",
            sum(!B$eff_ok), nrow(B)))
print(B[order(as_gap), .(player_name, season, arm = round(arm,1), axis = round(axis,1),
        as_ch = round(as_ch,2), as_fb = round(as_fb,2), as_gap = round(as_gap,2),
        matched = eff_ok, whiff_over = round(wh,1), gb_over = round(gb,1))], row.names = FALSE)
cr <- suppressWarnings(cor.test(A$axis, abs(A$as_gap), method = "spearman"))
cat(sprintf("\n  axis gap vs |efficiency gap| across all changeups: r = %+.3f (p = %.3g)\n",
            cr$estimate, cr$p.value))
cat("  -> the two criteria are close to independent, so the axis filter is NOT quietly\n     selecting efficiency-matched pitches.\n")

## ---- 2. does the axis effect survive controlling for efficiency? -----------------
cat("\n=== 2. axis x arm interaction, with the efficiency gap controlled ===\n")
for (m in list(c("wh ~ axis * arm", "no efficiency term"),
               c("wh ~ axis * arm + as_gap", "+ efficiency gap"),
               c("wh ~ axis * arm + as_gap + I(abs(as_gap))", "+ |efficiency gap|"),
               c("wh ~ axis * arm + as_gap + kill + velo_sep + spin", "+ efficiency, drop, velo, spin"),
               c("gb ~ axis * arm + as_gap + kill + velo_sep + spin", "same on ground balls"))) {
  co <- summary(lm(as.formula(m[1]), data = A))$coefficients
  cat(sprintf("  %-32s axis:arm  beta = %+.5f   p = %.4f\n", m[2], co["axis:arm","Estimate"],
              co["axis:arm","Pr(>|t|)"]))
}

## ---- 3. does efficiency match do anything by itself? -----------------------------
cat("\n=== 3. efficiency match on its own (all changeups, and within high slot) ===\n")
sp <- function(x,y){ z <- suppressWarnings(cor.test(x,y,method="spearman"))
                     sprintf("%+.3f (p=%.3f)", z$estimate, z$p.value) }
for (g in list(list(A, "all changeups"), list(A[arm >= 44], "high slot only"))) {
  q <- g[[1]]
  cat(sprintf("  %-16s n=%3d   |eff gap| vs whiff %s   vs grounders %s   vs RV %s\n",
      g[[2]], nrow(q), sp(abs(q$as_gap), q$wh), sp(abs(q$as_gap), q$gb), sp(abs(q$as_gap), q$rv100)))
}
co <- summary(lm(wh ~ I(abs(as_gap)) * arm + kill + velo_sep + spin, data = A))$coefficients
cat(sprintf("  |efficiency gap| x arm interaction on whiff: beta = %+.4f, p = %.3f  (compare the\n  axis version at p = %.4f - if efficiency were the real cue this would be the stronger one)\n",
    co["I(abs(as_gap)):arm","Estimate"], co["I(abs(as_gap)):arm","Pr(>|t|)"],
    summary(lm(wh ~ axis*arm + kill + velo_sep + spin, data = A))$coefficients["axis:arm","Pr(>|t|)"]))

## ---- 4. does requiring both sharpen or just shrink? ------------------------------
cat("\n=== 4. adding an efficiency requirement to the bin ===\n")
A[, grp := fifelse(para & eff_ok, "axis + efficiency",
           fifelse(para & !eff_ok, "axis only", "neither"))]
print(A[, .(n = .N, whiff_over = round(mean(wh),2), gb_over = round(mean(gb),2),
            rv100 = round(mean(rv100),2), eff_gap = round(mean(as_gap),3)),
        by = grp][order(-whiff_over)], row.names = FALSE)
for (g in c("axis + efficiency","axis only")) {
  q <- A[grp == g]; if (nrow(q) < 3) next
  for (v in c("wh","gb")) { t <- t.test(q[[v]], A[grp == "neither"][[v]])
    cat(sprintf("  %-18s vs neither, %-3s  %+.2f vs %+.2f   p = %.4f  (n=%d)\n",
        g, v, mean(q[[v]]), mean(A[grp=="neither"][[v]]), t$p.value, nrow(q))) } }

fwrite(A[order(-para, -wh)], file.path(AST, "ext_parachute_spin_efficiency.csv"))

## ---- figure -----------------------------------------------------------------------
A[, lab := fifelse(para, sub(",.*","", player_name), NA_character_)]
A[!is.na(lab), lab := fifelse(seq_len(.N) == which.max(np), lab, NA_character_), by = lab]
gg <- ggplot(A, aes(axis, as_gap)) +
  annotate("rect", xmin = -Inf, xmax = 10, ymin = -.10, ymax = .10, alpha = .12, fill = "#2a9d8f") +
  geom_hline(yintercept = 0, colour = "grey45", linewidth = .4) +
  geom_vline(xintercept = 10, linetype = "dashed", colour = "#c0392b", linewidth = .45) +
  geom_point(aes(size = np, colour = wh), alpha = .78) +
  geom_text_repel(aes(label = lab), size = 3, max.overlaps = Inf, seed = 5, box.padding = .5,
                  min.segment.length = 0, segment.size = .3, segment.colour = "grey55",
                  colour = "grey10", fontface = "bold") +
  scale_colour_gradient2(low = "#c0392b", mid = "grey85", high = "#1d7870", midpoint = 0,
                         name = "Whiff% over model") +
  scale_size_continuous(range = c(.9, 4.3), guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Spin similarity in these charts is axis only, and axis says nothing about efficiency",
       subtitle = paste0("One point per pitcher-season with measured active spin on both the changeup and the primary fastball, min 60 changeup swings. The x-axis is the clock-face\n",
                         "spin-axis gap that defines the bin; the y-axis is how much less transverse spin the changeup carries than the fastball. The two are nearly independent\n",
                         "(r = ", sprintf("%+.2f", cr$estimate), "), so passing the axis filter tells you almost nothing about whether the pitch actually spins like the fastball in three dimensions.\n",
                         "The shaded box is what a full same-spin definition would require. Labelled points are the current bin; several sit well outside the box."),
       x = "Spin-axis gap vs primary fastball (degrees)",
       y = "Active-spin gap vs primary fastball\n(negative = changeup is more gyro)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12.8), plot.subtitle = element_text(size = 8.1),
        panel.grid.minor = element_blank())
ggsave(file.path(AST, "fig22_axis_vs_efficiency.png"), gg, width = 12.5, height = 7, dpi = 150)
cat("\nwrote fig22_axis_vs_efficiency.png, ext_parachute_spin_efficiency.csv\n")
