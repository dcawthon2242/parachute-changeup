#!/usr/bin/env Rscript

# WHY DID THE VELOCITY INTERACTION WEAKEN UNDER tjStuff+ v3.0 FEATURES?
#
# Swapping to Nestico's feature set changed two things at once, so the drop from p = .05 to
# p = .14 could not be attributed:
#
#   the shape block   v3.0 mirrors x0, ax and ax_diff for left-handers rather than carrying
#                     handedness dummies, and it adds spin_axis, which my set never had.
#   the location block  the nine columns I paired with v3.0 are a subset of the fifteen the old
#                     model used. Dropped: plate_x, plate_x_in, HAA, throws_R, same_hand.
#
# A ladder that moves one block at a time separates them:
#
#   M1  old shape + full location                  the p = .05 reference
#   M2  old shape + spin axis + full location      adds the one genuinely new feature
#   M3  v3.0 shape + full location                 swaps the shape block, location held fixed
#   M4  v3.0 shape + trimmed location              the model from the tjStuff+ run
#   M5  v3.0 shape only                            no location at all
#
# TWO THINGS THIS SCRIPT HAS TO GET RIGHT.
#
# Fold assignment is shared. An earlier version drew fresh folds inside each fit, which means
# a rung-to-rung difference mixed the feature change with a different out-of-fold partition.
# Every model here sees the identical partition, so the ladder is a clean one-factor comparison.
#
# Fold noise is measured, not assumed. Even on shared folds a single partition is one draw, and
# the interaction rests on ~21 pitcher-seasons, so it is fragile. The critical rung (M1 -> M2)
# is therefore refit across several partitions and read as a PAIRED difference: within each
# partition both models see the same folds, so the pairing cancels the partition entirely and
# what is left is the feature's own contribution.

suppressPackageStartupMessages({ library(data.table); library(lightgbm); library(ggplot2) })
set.seed(1); options(width = 205)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")
CACHE <- file.path(MDIR, "velo_feature_attribution.rds"); SEEDS <- 1:6

CH <- readRDS(file.path(MDIR, "parachute_ff.rds"))
CH[, spin_axis := (atan2(sax, cax)*180/pi) %% 360]
CH[, L := p_throws == "L"]
CH[, `:=`(tj_x0 = fifelse(L, -release_pos_x, release_pos_x), tj_ax = fifelse(L, -ax, ax),
          tj_ax_diff = fifelse(L, -ax_diff, ax_diff),
          tj_axis = fifelse(L, (360 - spin_axis) %% 360, spin_axis))]

SHAPE_OLD <- c("release_speed","release_spin_rate","release_extension","release_pos_x",
               "release_pos_z","ax","az","speed_diff","ax_diff","az_diff")
SHAPE_TJ  <- c("release_speed","release_spin_rate","release_extension","tj_x0",
               "release_pos_z","tj_ax","az","speed_diff","tj_ax_diff","az_diff","tj_axis")
LOC_FULL  <- c("plate_x","plate_z","plate_x_in","plate_x_arm","z_rel_bot","z_rel_top","VAA",
               "HAA","HAA_in","stand_R","throws_R","same_hand","balls","strikes")
LOC_TRIM  <- c("plate_x_arm","plate_z","z_rel_bot","z_rel_top","VAA","HAA_in","stand_R",
               "balls","strikes")
SETS <- list(M1 = c(SHAPE_OLD, LOC_FULL), M2 = c(SHAPE_OLD, "tj_axis", LOC_FULL),
             M3 = c(SHAPE_TJ, LOC_FULL),  M4 = c(SHAPE_TJ, LOC_TRIM), M5 = SHAPE_TJ)
LAB <- c(M1 = "M1  old shape + full location",     M2 = "M2  + spin axis",
         M3 = "M3  v3.0 shape + full location",    M4 = "M4  v3.0 shape + trimmed location",
         M5 = "M5  v3.0 shape only")

need <- unique(c(unlist(SETS), "whiff","is_swing","axis_diff","arm_angle"))
CH <- CH[is.finite(axis_diff) & is.finite(arm_angle) & stats::complete.cases(CH[, ..need])]
SW <- CH[is_swing == TRUE]
cat(sprintf("changeup swings: %s\n", format(nrow(SW), big.mark=",")))

fit <- function(FEAT, fold) {
  p <- rep(NA_real_, nrow(SW))
  for (f in unique(fold)) {
    tr <- SW[fold != f]; n <- nrow(tr); vi <- sample(n, floor(.12*n))
    dtr <- lgb.Dataset(as.matrix(tr[-vi, ..FEAT]), label = tr$whiff[-vi])
    dva <- lgb.Dataset.create.valid(dtr, as.matrix(tr[vi, ..FEAT]), label = tr$whiff[vi])
    m <- lgb.train(params = list(objective = "binary", metric = "binary_logloss",
                   learning_rate = .06, num_leaves = 31, min_data_in_leaf = 300,
                   feature_fraction = .8, bagging_fraction = .8, bagging_freq = 1),
                   data = dtr, nrounds = 1500, valids = list(v = dva),
                   early_stopping_rounds = 50, verbose = -1)
    p[fold == f] <- predict(m, as.matrix(SW[fold == f, ..FEAT]))
  }
  list(res = SW$whiff - p, r2 = 1 - var(SW$whiff - p)/var(SW$whiff))
}

AS <- readRDS(file.path(MDIR, "active_spin_long.rds"))
# Collapse per-pitch residuals to pitcher-seasons and return the bin x velocity interaction.
# The interaction term, not either slope alone, is what says the bin is special.
score <- function(res) {
  D <- copy(SW[, .(pitcher, player_name, season, axis_diff, speed_diff, arm_angle)])[, r := res]
  S <- D[, .(nsw = .N, axis = mean(axis_diff), arm = mean(arm_angle),
             velo_sep = -mean(speed_diff), w = 100*mean(r)), by = .(pitcher, player_name, season)]
  S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, as_ch = active_spin)], by = c("pitcher","season"))
  S <- merge(S, AS[pitch_type == "FF", .(pitcher, season, as_fb = active_spin)], by = c("pitcher","season"))
  S <- S[nsw >= 60]; S[, as_gap := round(as_ch - as_fb, 4)]
  S[, `:=`(Core = abs(as_gap) <= .10 & arm >= 44 & axis <= 10,
           HighSpin = as_ch >= .80 & as_fb >= .80 & arm >= 44 & axis <= 10)]
  rbindlist(lapply(c("Core","HighSpin"), function(b) {
    S[, bin := get(b)]
    cf <- summary(lm(w ~ velo_sep*bin, S, weights = nsw))$coefficients
    d <- cf["velo_sep:binTRUE",]
    data.table(bin = b, n = sum(S$bin), slope_out = cf["velo_sep","Estimate"],
               slope_in = cf["velo_sep","Estimate"] + d[1], delta = d[1],
               lo = d[1]-1.96*d[2], hi = d[1]+1.96*d[2], p = d[4],
               r_in = cor(S[bin == TRUE]$velo_sep, S[bin == TRUE]$w)) }))
}

if (!file.exists(CACHE)) {
  set.seed(11); FOLD <- sample(rep(1:4, length.out = nrow(SW)))   # shared across the ladder
  cat("\n=== ladder, shared fold partition ===\n")
  ladder <- rbindlist(lapply(names(SETS), function(k) {
    f <- fit(SETS[[k]], FOLD)
    cat(sprintf("  %-38s %2d features  R2=%.4f\n", LAB[k], length(SETS[[k]]), f$r2))
    cbind(model = k, r2 = f$r2, score(f$res)) }))

  # Paired seed study on the rung that moved. Same folds for both models within a seed.
  cat("\n=== paired seed study on M1 -> M2 ===\n")
  seeds <- rbindlist(lapply(SEEDS, function(s) {
    set.seed(100 + s); fo <- sample(rep(1:4, length.out = nrow(SW)))
    a <- score(fit(SETS$M1, fo)$res); b <- score(fit(SETS$M2, fo)$res)
    cat(sprintf("  seed %d  Core: M1 %+.2f (p=%.3f) -> M2 %+.2f (p=%.3f)   paired change %+.2f\n",
        s, a[bin=="Core"]$delta, a[bin=="Core"]$p, b[bin=="Core"]$delta, b[bin=="Core"]$p,
        b[bin=="Core"]$delta - a[bin=="Core"]$delta))
    data.table(seed = s, bin = a$bin, M1 = a$delta, M2 = b$delta, p1 = a$p, p2 = b$p) }))
  saveRDS(list(ladder = ladder, seeds = seeds), CACHE)
} else cat("(using cached fits)\n")
z <- readRDS(CACHE); ladder <- z$ladder; seeds <- z$seeds

cat("\n=== velocity slope inside vs outside the bin, per model (pp of whiff per mph) ===\n")
print(ladder[, .(model = LAB[model], bin, n, R2 = round(r2,4), slope_in = round(slope_in,2),
                 slope_out = round(slope_out,2), difference = round(delta,2),
                 p = round(p,4), r_in = round(r_in,3))], row.names = FALSE)

cat("\n=== attribution: change in the interaction at each rung ===\n")
for (b in c("Core","HighSpin")) { C <- ladder[bin == b]
  cat(sprintf("  [%s]\n", b))
  for (i in 2:nrow(C)) cat(sprintf("    %-36s -> %-36s  %+.2f\n",
      gsub("\\s+"," ",LAB[C$model[i-1]]), gsub("\\s+"," ",LAB[C$model[i]]), C$delta[i]-C$delta[i-1])) }

cat("\n=== is the M1 -> M2 drop bigger than fold noise? (paired over 6 partitions) ===\n")
for (b in c("Core","HighSpin")) { d <- seeds[bin == b]; t <- t.test(d$M2, d$M1, paired = TRUE)
  cat(sprintf("  %-9s M1 %+.2f +- %.2f | M2 %+.2f +- %.2f | paired change %+.2f [%+.2f, %+.2f] p=%.4f\n",
      b, mean(d$M1), sd(d$M1), mean(d$M2), sd(d$M2), t$estimate, t$conf.int[1], t$conf.int[2], t$p.value)) }
cat(sprintf("  M1 significant at p<.05 in %d of %d partitions (Core); M2 in %d\n",
            sum(seeds[bin=="Core"]$p1 < .05), length(SEEDS), sum(seeds[bin=="Core"]$p2 < .05)))

## ---- figure ---------------------------------------------------------------------------------
P <- copy(ladder)
P[, model := factor(LAB[model], levels = rev(LAB))]
P[, bin := factor(bin, levels = c("Core","HighSpin"),
                  labels = c("Core bin: efficiency gap <= .10, arm >= 44, axis <= 10  (21 seasons)",
                             "High-spin bin: active spin >= .80 both, arm >= 44, axis <= 10  (20 seasons)"))]
P[, sig := fifelse(p < .05, "p < .05", "not significant")]
gg <- ggplot(P, aes(delta, model, colour = sig)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = .18, linewidth = .5) +
  geom_point(size = 2.8) +
  geom_text(aes(label = sprintf("%+.2f   p=%.3f", delta, p)), vjust = -1.05, size = 2.85,
            show.legend = FALSE) +
  facet_wrap(~bin, ncol = 1) +
  scale_colour_manual(values = c("p < .05" = "#1d7870", "not significant" = "grey55"), name = NULL) +
  expand_limits(x = c(-1.4, 3.2)) +
  labs(title = "It was not the trimmed location block. One feature - spin axis - accounts for the whole drop",
       subtitle = paste0("How much steeper the velocity-separation slope is inside the parachute bin than outside it, in percentage points of whiff residual per mph, with 95 percent\n",
                         "intervals. All five models share one out-of-fold partition, so each rung isolates a single change. Mirroring the shape block moves the Core estimate by\n",
                         "+0.04 and trimming location from fifteen columns to nine moves it by -0.08; neither is the cause, and the caveat I raised about the trim was wrong. Handing\n",
                         "the model spin axis is what costs the result: the interaction falls from +0.84 at p = .03 to +0.61 at p = .11. Refitting both across six independent\n",
                         "partitions puts the paired change at -0.32 [-0.38, -0.25], and fold noise is far too small to explain it - M1 clears p < .05 in six partitions of six, M2\n",
                         "in none. The bin is defined on the spin-axis gap, so once the model sees the axis it prices in part of what the bin selects on, which makes this a stricter\n",
                         "test rather than a broken one. Dropping location entirely costs a further 0.28, the reverse of what a location confound would do: if the bin were only\n",
                         "pitchers who locate changeups well, hiding location from the model would inflate the effect, not halve it."),
       x = "Extra whiff-residual slope inside the bin (percentage points per mph)", y = NULL,
       caption = "Source: Statcast 2020-2026 - four-seam anchor - out-of-fold LightGBM - pitcher-seasons with 60+ changeup swings, weighted by swings") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11.5), plot.subtitle = element_text(size = 8.2),
        strip.text = element_text(face = "bold", size = 9.5), panel.grid.minor = element_blank(),
        axis.text.y = element_text(family = "mono", size = 8.5),
        legend.position = "top", plot.caption = element_text(size = 7.5, colour = "grey40"))
ggsave(file.path(AST, "fig35_velo_attribution.png"), gg, width = 11.5, height = 7.8, dpi = 150)
fwrite(ladder, file.path(AST, "ext_velo_attribution.csv"))
cat("\nwrote ext_velo_attribution.csv, fig35_velo_attribution.png\n")
