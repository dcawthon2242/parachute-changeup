#!/usr/bin/env Rscript

# THE PARACHUTE BIN, ANCHORED ON A GENUINE FOUR-SEAMER.
#
# Everything up to now compared the changeup to the pitcher's "primary" fastball, defined as
# his four-seamer if he threw fifty of them and otherwise his sinker. That fallback is a
# problem for a bin defined on spin-axis match. A sinker's seam-shifted wake gives it movement
# its spin axis does not predict, so two pitches sharing a spin axis with a sinker need not
# look alike out of the hand, while two pitches sharing an axis with a four-seamer largely do.
#
# The fallback also biases selection. Sinker-anchored seasons are 9.4 percent of the eligible
# population but 15 percent of the Core bin, and their mean axis gap is 18.4 degrees against
# 22.3 for four-seam-anchored seasons - a sinker sits closer to a changeup on the clock face
# almost by construction, so the old bin was partly a sinkerballer detector.
#
# This refits everything against four-seamers only. Pitchers without a real four-seam are
# dropped rather than given a substitute reference. Both versions are reported side by side.

suppressPackageStartupMessages({ library(data.table); library(lightgbm); library(ggplot2) })
set.seed(1); options(width = 205)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")
MINSW <- 60L
RES <- file.path(MDIR, "parachute_ff_resid.rds")

FEAT <- c("release_speed","release_spin_rate","release_extension","release_pos_x",
          "release_pos_z","ax","az","speed_diff","ax_diff","az_diff","plate_x","plate_z",
          "plate_x_in","plate_x_arm","z_rel_bot","z_rel_top","VAA","HAA","HAA_in",
          "stand_R","throws_R","same_hand","balls","strikes")

oof <- function(D, y, obj, tag) {
  K <- 4; fold <- sample(rep(1:K, length.out = nrow(D))); p <- rep(NA_real_, nrow(D))
  for (f in 1:K) {
    tr <- D[fold != f]; n <- nrow(tr); vi <- sample(n, floor(.12*n))
    dtr <- lgb.Dataset(as.matrix(tr[-vi, ..FEAT]), label = y[fold != f][-vi])
    dva <- lgb.Dataset.create.valid(dtr, as.matrix(tr[vi, ..FEAT]), label = y[fold != f][vi])
    m <- lgb.train(params = list(objective = obj,
                   metric = if (obj == "binary") "binary_logloss" else "l2",
                   learning_rate = .06, num_leaves = 31, min_data_in_leaf = 300,
                   feature_fraction = .8, bagging_fraction = .8, bagging_freq = 1),
                   data = dtr, nrounds = 1500, valids = list(v = dva),
                   early_stopping_rounds = 50, verbose = -1)
    p[fold == f] <- predict(m, as.matrix(D[fold == f, ..FEAT]))
  }
  cat(sprintf("  %-9s n=%8s  R2=%.4f\n", tag, format(nrow(D), big.mark=","), 1 - var(y-p)/var(y)))
  y - p
}

if (!file.exists(RES) || nzchar(Sys.getenv("REFIT"))) {
  CH <- readRDS(file.path(MDIR, "parachute_ff.rds"))
  CH <- CH[is.finite(axis_diff) & stats::complete.cases(CH[, ..FEAT])]
  cat("=== out-of-fold models on four-seam-anchored data; axis gap and arm angle excluded ===\n")
  CH[, rv_res := oof(CH, CH$rv, "regression", "run value")]
  SW <- CH[is_swing == TRUE]; SW[, wh_res := oof(SW, SW$whiff, "binary", "whiff")]
  BP <- CH[is_bip == TRUE]; BP[, gb := as.integer(bb_type == "ground_ball")]
  BP[, gb_res := oof(BP, BP$gb, "binary", "grounder")]
  saveRDS(list(CH = CH[, .(pitcher, player_name, season, axis_diff, az_diff, speed_diff,
                           arm_angle, release_spin_rate, rv_res, is_swing, is_bip)],
               SW = SW[, .(pitcher, season, wh_res)],
               BP = BP[, .(pitcher, season, gb_res, launch_speed)]), RES)
}
L <- readRDS(RES)
S <- L$CH[, .(np = .N, nsw = sum(is_swing), axis = mean(axis_diff), kill = -mean(az_diff),
              velo_sep = -mean(speed_diff), arm = mean(arm_angle, na.rm = TRUE),
              spin = mean(release_spin_rate, na.rm = TRUE), rv100 = 100*mean(rv_res)),
          by = .(pitcher, player_name, season)]
S <- merge(S, L$SW[, .(wh = 100*mean(wh_res)), by = .(pitcher, season)], by = c("pitcher","season"))
S <- merge(S, L$BP[, .(gb = 100*mean(gb_res)), by = .(pitcher, season)], by = c("pitcher","season"))

AS <- readRDS(file.path(MDIR, "active_spin_long.rds"))
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, as_ch = active_spin)],
           by = c("pitcher","season"), all.x = TRUE)
S <- merge(S, AS[pitch_type == "FF", .(pitcher, season, as_fb = active_spin)],
           by = c("pitcher","season"), all.x = TRUE)
# Active spin is published to the nearest percent, so the gap is meaningful only to about that
# precision. Rounding it here stops boundary cases like 0.72 - 0.87 landing at -0.15000000000002
# and failing an <= 0.15 test in memory while passing it after a CSV round-trip.
S[, as_gap := round(as_ch - as_fb, 4)]
S <- S[nsw >= MINSW & is.finite(axis) & is.finite(arm) & is.finite(as_gap)]
S[, `:=`(core = axis <= 10 & abs(as_gap) <= .10 & arm >= 44,
         wide = axis <= 15 & abs(as_gap) <= .15 & arm >= 42,
         last = sub(",.*", "", player_name))]
fwrite(S, file.path(AST, "ext_parachute_ff.csv"))

## ---- what the anchor change did ---------------------------------------------------------
O <- fread(file.path(AST, "ext_parachute_extended.csv"))
O <- O[is.finite(as_gap) & is.finite(arm)]
O[, `:=`(core = axis <= 10 & abs(as_gap) <= .10 & arm >= 44,
         wide = axis <= 15 & abs(as_gap) <= .15 & arm >= 42)]
cat(sprintf("\n=== population ===\n  primary anchor %4d seasons  |  four-seam anchor %4d seasons  (%d lost)\n",
            nrow(O), nrow(S), nrow(O) - nrow(S)))
cat(sprintf("  mean axis gap  %.1f deg            |  %.1f deg\n", mean(O$axis), mean(S$axis)))

cat("\n=== bin sizes and effects, both anchors ===\n")
eff <- function(D, i, v) { t <- t.test(D[[v]][i], D[[v]][!i])
  sprintf("%+.2f (p=%.3f)", diff(rev(t$estimate)), t$p.value) }
cmp <- rbindlist(lapply(c("core","wide"), function(b) rbindlist(lapply(
  list(list("primary", O), list("four-seam", S)), function(z) {
    D <- z[[2]]; i <- D[[b]]
    data.table(bin = b, anchor = z[[1]], seasons = sum(i), pitchers = uniqueN(D$pitcher[i]),
               grounders = eff(D, i, "gb"), whiff = eff(D, i, "wh"), rv100 = eff(D, i, "rv100")) }))))
print(cmp, row.names = FALSE)

cat("\n=== who enters and leaves the Core bin ===\n")
ko <- O[core == TRUE, paste(pitcher, season)]; kn <- S[core == TRUE, paste(pitcher, season)]
cat(sprintf("  in both: %d   dropped by the four-seam rule: %d   newly admitted: %d\n",
            length(intersect(ko,kn)), length(setdiff(ko,kn)), length(setdiff(kn,ko))))
if (length(setdiff(ko,kn))) {
  cat("\n  DROPPED (no four-seam, or the gap widened once measured against one):\n")
  print(O[core == TRUE][paste(pitcher, season) %in% setdiff(ko,kn),
        .(player_name = sub(",.*","",player_name), season, fb_type, axis_old = round(axis,1),
          wh_old = round(wh,1), gb_old = round(gb,1))][order(-wh_old)], row.names = FALSE)
}
if (length(setdiff(kn,ko))) {
  cat("\n  NEWLY ADMITTED:\n")
  print(S[core == TRUE][paste(pitcher, season) %in% setdiff(kn,ko),
        .(player_name = last, season, axis_new = round(axis,1), velo_sep = round(velo_sep,1),
          wh_new = round(wh,1), gb_new = round(gb,1))][order(-wh_new)], row.names = FALSE)
}

## ---- does the velocity interaction survive? ---------------------------------------------
cat("\n=== velocity separation x bin interaction, both anchors ===\n")
for (z in list(list("primary", O), list("four-seam", S))) {
  D <- z[[2]]
  for (b in c("core","wide")) {
    D[, grp := get(b)]; C <- D[grp == TRUE]
    m <- lm(wh ~ velo_sep * grp, D, weights = nsw); cf <- summary(m)$coefficients
    r <- cor.test(C$velo_sep, C$wh)
    cat(sprintf("  %-10s %-5s n=%3d  r=%+.3f p=%.3f | slope in %+.2f out %+.2f  interaction p=%.4f\n",
                z[[1]], b, nrow(C), r$estimate, r$p.value,
                cf["velo_sep","Estimate"] + cf["velo_sep:grpTRUE","Estimate"],
                cf["velo_sep","Estimate"], cf["velo_sep:grpTRUE","Pr(>|t|)"]))
  }
}
cat("\n=== held-out 2020-2022 under the four-seam anchor ===\n")
for (w in list(c(2020,2022,"HELD OUT 2020-2022"), c(2023,2026,"discovery 2023-2026"))) {
  D <- S[season >= as.integer(w[1]) & season <= as.integer(w[2])]; C <- D[core == TRUE]
  if (nrow(C) >= 5) {
    m <- lm(wh ~ velo_sep * core, D, weights = nsw); cf <- summary(m)$coefficients
    cat(sprintf("  %-22s bin n=%2d  r=%+.3f | slope in %+.2f out %+.2f  interaction p=%.4f\n",
                w[3], nrow(C), cor(C$velo_sep, C$wh),
                cf["velo_sep","Estimate"] + cf["velo_sep:coreTRUE","Estimate"],
                cf["velo_sep","Estimate"], cf["velo_sep:coreTRUE","Pr(>|t|)"]))
  }
}

## ---- figure -------------------------------------------------------------------------------
P <- rbind(O[, .(velo_sep, wh, np, core, anchor = "Primary fastball (four-seam, else sinker)")],
           S[, .(velo_sep, wh, np, core, anchor = "Four-seam fastball only")])
P[, anchor := factor(anchor, levels = c("Primary fastball (four-seam, else sinker)",
                                        "Four-seam fastball only"))]
N <- P[core == TRUE, .N, by = anchor]
P <- merge(P, N, by = "anchor")
P[, facet := sprintf("%s\nCore bin: %d seasons", anchor, N)]
P[, facet := factor(facet, levels = unique(facet[order(anchor)]))]
gg <- ggplot(P, aes(velo_sep, wh)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_point(data = P[core == FALSE], colour = "grey75", size = .75, alpha = .45) +
  geom_smooth(data = P[core == FALSE], method = "lm", se = TRUE, colour = "grey35",
              fill = "grey35", alpha = .12, linewidth = .8) +
  geom_point(data = P[core == TRUE], aes(size = np), colour = "#1d7870", alpha = .85) +
  geom_smooth(data = P[core == TRUE], method = "lm", se = TRUE, colour = "#1d7870",
              fill = "#1d7870", alpha = .16, linewidth = 1.1) +
  facet_wrap(~facet) + scale_size_area(max_size = 5.5, guide = "none") +
  coord_cartesian(ylim = c(-24, 26)) +
  labs(title = "Requiring a real four-seamer kills the ground-ball effect and leaves the velocity effect standing",
       subtitle = paste0("Both panels rebuild the bin from scratch, refitting the residual model on the anchor shown. Dropping the sinker fallback removes 116 pitcher-seasons and\n",
                         "four Core-bin members - Blackburn '26, Perez '25 and both Brett Anderson seasons - and takes the ground-ball edge from +1.82 points (p = .26) to +0.75\n",
                         "(p = .65). Three of those four had strongly positive grounder residuals, which is what a sinkerballer looks like, so that effect was substantially a\n",
                         "sinker artifact and the earlier reading of it was wrong. The velocity slope inside the bin barely moves: +1.19 points per mph becomes +1.09, and the\n",
                         "interaction against the rest of the league holds at p = .0095. No pitcher-season was newly admitted by the change."),
       x = "Velocity separation from the anchor fastball (mph slower)",
       y = "Whiff rate above model (percentage points)",
       caption = "Source: Statcast 2020-2026 - out-of-fold LightGBM residuals - point size is changeups thrown") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12.2), plot.subtitle = element_text(size = 8.2),
        strip.text = element_text(face = "bold", size = 9.5), panel.grid.minor = element_blank(),
        plot.caption = element_text(size = 7.5, colour = "grey40"))
saveRDS(gg, file.path(MDIR, "fig31_gg.rds"))
ggsave(file.path(AST, "fig31_ff_anchor.png"), gg, width = 11, height = 5.8, dpi = 150)
cat("\nwrote ext_parachute_ff.csv, fig31_ff_anchor.png\n")
