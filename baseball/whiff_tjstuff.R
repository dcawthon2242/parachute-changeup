#!/usr/bin/env Rscript

# THE WHIFF MODEL, REBUILT ON tjStuff+ v3.0 FEATURES.
#
# Nestico's v3.0 feature set is eleven columns and nothing else:
#
#   start_speed  spin_rate  extension  ax  az  x0  z0  spin_axis  speed_diff  ax_diff  az_diff
#
# with two design choices that matter here:
#
#   mirroring   every x-dimension quantity is flipped for left-handers so both hands train on
#               one scale. That covers x0, ax and ax_diff, and it covers spin_axis too, which
#               reflects as 360 - axis. v3.0 also drops the explicit handedness dummies my
#               earlier feature set carried, because mirroring makes them redundant.
#   no location the model never sees plate_x, plate_z or the approach angles. Nestico is
#               explicit that this is why tjStuff+ is a poor descriptive metric. For us it
#               means the residual now contains command, so "above tjStuff" is a different and
#               broader quantity than "above shape and location" - it is closer to the way the
#               phrase "outperforming his stuff" is normally used.
#
# v3.0 also applies a RobustScaler. That is a per-feature affine transform and LightGBM splits
# on order statistics, so it cannot change a tree model's fit; it is omitted deliberately
# rather than by oversight.
#
# Nestico targets run value. The target here stays whiff-on-swing, because that is the question
# this project is asking; only the feature set is his.
#
# Four models are fit so the feature choice can be attributed rather than assumed.

suppressPackageStartupMessages({ library(data.table); library(lightgbm); library(ggplot2) })
set.seed(1); options(width = 205)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")
CACHE <- file.path(MDIR, "whiff_tjstuff.rds")

CH <- readRDS(file.path(MDIR, "parachute_ff.rds"))

# The build stores the spin axis as its sine and cosine; atan2 recovers the original degrees
# exactly, so there is no need to re-scrape for it.
CH[, spin_axis := (atan2(sax, cax)*180/pi) %% 360]
CH[, L := p_throws == "L"]
CH[, `:=`(tj_x0      = fifelse(L, -release_pos_x, release_pos_x),
          tj_ax      = fifelse(L, -ax, ax),
          tj_ax_diff = fifelse(L, -ax_diff, ax_diff),
          tj_axis    = fifelse(L, (360 - spin_axis) %% 360, spin_axis))]

TJ   <- c("release_speed","release_spin_rate","release_extension","tj_ax","az","tj_x0",
          "release_pos_z","tj_axis","speed_diff","tj_ax_diff","az_diff")
LOC  <- c("plate_x_arm","plate_z","z_rel_bot","z_rel_top","VAA","HAA_in","stand_R","balls","strikes")
SETS <- list("tjStuff+ v3.0"                = TJ,
             "tjStuff+ v3.0 + arm angle"    = c(TJ, "arm_angle"),
             "tjStuff+ v3.0 + location"     = c(TJ, LOC),
             "tjStuff+ v3.0 + arm + location" = c(TJ, "arm_angle", LOC))

need <- unique(c(unlist(SETS), "whiff","is_swing","axis_diff","arm_angle","speed_diff"))
CH <- CH[is.finite(axis_diff) & is.finite(arm_angle) & stats::complete.cases(CH[, ..need])]
SW <- CH[is_swing == TRUE]
cat(sprintf("changeup swings: %s   pitchers: %d\n\n", format(nrow(SW), big.mark=","), uniqueN(SW$pitcher)))

if (!file.exists(CACHE)) {
  fit <- function(FEAT, tag) {
    K <- 4; fold <- sample(rep(1:K, length.out = nrow(SW))); p <- rep(NA_real_, nrow(SW))
    for (f in 1:K) {
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
    cat(sprintf("  %-32s %2d features  R2=%.4f\n", tag, length(FEAT), 1 - var(SW$whiff-p)/var(SW$whiff)))
    SW$whiff - p
  }
  cat("=== out-of-fold whiff models ===\n")
  for (i in seq_along(SETS)) SW[[paste0("r", i)]] <- fit(SETS[[i]], names(SETS)[i])
  saveRDS(SW[, .(pitcher, player_name, season, axis_diff, speed_diff, arm_angle,
                 r1, r2, r3, r4)], CACHE)
} else cat("(using cached fits)\n")
R <- readRDS(CACHE)

# One importance run on the full sample, purely to report which of the eleven the model leans on.
FEAT <- TJ
mi <- lgb.train(params = list(objective = "binary", learning_rate = .06, num_leaves = 31,
                min_data_in_leaf = 300, feature_fraction = .8, bagging_fraction = .8,
                bagging_freq = 1), verbose = -1, nrounds = 400,
                data = lgb.Dataset(as.matrix(SW[, ..FEAT]), label = SW$whiff))
imp <- as.data.table(lgb.importance(mi))[order(-Gain)]
cat("\n=== which tjStuff+ features the whiff model leans on ===\n")
print(imp[, .(feature = Feature, gain = round(100*Gain,1), split_share = round(100*Frequency,1))],
      row.names = FALSE)

## ---- season table, one residual column per feature set ----------------------------------
AS <- readRDS(file.path(MDIR, "active_spin_long.rds"))
S <- R[, .(nsw = .N, axis = mean(axis_diff), arm = mean(arm_angle), velo_sep = -mean(speed_diff),
           w1 = 100*mean(r1), w2 = 100*mean(r2), w3 = 100*mean(r3), w4 = 100*mean(r4)),
       by = .(pitcher, player_name, season)]
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, as_ch = active_spin)], by = c("pitcher","season"))
S <- merge(S, AS[pitch_type == "FF", .(pitcher, season, as_fb = active_spin)], by = c("pitcher","season"))
S <- S[nsw >= 60]
S[, `:=`(hs = as_ch >= .80 & as_fb >= .80, hi = arm >= 44, last = sub(",.*","",player_name))]
cat(sprintf("\npitcher-seasons at a 60-swing floor: %d\n", nrow(S)))

NM <- names(SETS); VS <- paste0("w", 1:4)
cat("\n=== the cues, scored against each feature set (pp of whiff above model) ===\n")
res <- rbindlist(lapply(seq_along(VS), function(i) {
  v <- VS[i]
  ct <- function(cond, lab) { t <- t.test(S[[v]][cond], S[[v]][!cond])
    data.table(model = NM[i], cue = lab, est = diff(rev(t$estimate)), p = t$p.value) }
  cr <- function(x, lab, D = S) { r <- cor.test(D[[x]], D[[v]])
    data.table(model = NM[i], cue = lab, est = r$estimate, p = r$p.value) }
  rbind(ct(S$hi, "High slot >= 44 deg (pp)"),
        ct(S$hs, "High active spin both >= .80 (pp)"),
        ct(S$hs & S$hi & S$axis <= 15, "Full gate, axis <= 15 (pp)"),
        cr("axis", "Axis gap vs whiff, inside gate (r)", S[hs & hi]),
        cr("arm", "Arm angle as continuous cue (r)")) }))
W <- dcast(res, cue ~ model, value.var = c("est","p"))
setcolorder(W, c("cue", as.vector(rbind(paste0("est_", NM), paste0("p_", NM)))))
print(W[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], row.names = FALSE)

cat("\n=== velocity separation x bin interaction (axis <= 10), each feature set ===\n")
S[, bin := hs & hi & axis <= 10]
for (i in seq_along(VS)) {
  v <- VS[i]; C <- S[bin == TRUE]
  m <- lm(get(v) ~ velo_sep*bin, S, weights = nsw); cf <- summary(m)$coefficients
  cat(sprintf("  %-32s n=%2d  r=%+.3f | slope in %+.2f out %+.2f  interaction p=%.4f\n",
      NM[i], nrow(C), cor(C$velo_sep, C[[v]]),
      cf["velo_sep","Estimate"] + cf["velo_sep:binTRUE","Estimate"],
      cf["velo_sep","Estimate"], cf["velo_sep:binTRUE","Pr(>|t|)"]))
}

cat("\n=== biggest overperformers against tjStuff+ v3.0 ===\n")
print(S[order(-w1)][1:15, .(Pitcher = last, Season = season, Swings = nsw, Arm = round(arm,1),
      Axis = round(axis,1), VeloSep = round(velo_sep,1), ActCH = round(as_ch,3),
      tjStuff = round(w1,1), plus_arm = round(w2,1), plus_loc = round(w3,1))], row.names = FALSE)
fwrite(S, file.path(AST, "ext_whiff_tjstuff.csv"))

## ---- figure -------------------------------------------------------------------------------
R2 <- c(.0142, .0158, .1738, .1759)
P <- res[grepl("pp\\)", cue)]
P[, cue := factor(sub(" \\(pp\\)","",cue),
                  levels = rev(c("High slot >= 44 deg","High active spin both >= .80",
                                 "Full gate, axis <= 15")))]
P[, model := factor(sprintf("%s   (R2 = %.3f)", model, R2[match(model, NM)]),
                    levels = sprintf("%s   (R2 = %.3f)", NM, R2))]
P[, sig := fifelse(p < .05, "p < .05", "not significant")]
gg <- ggplot(P, aes(est, cue, colour = sig)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_point(size = 3) +
  # Labels sit on the outboard side of each point so the two negative estimates do not run off
  # the left edge of their panels.
  geom_text(aes(label = sprintf("%+.2f  p=%.3f", est, p), hjust = fifelse(est < 0, 1.15, -0.15)),
            size = 2.9, show.legend = FALSE) +
  facet_wrap(~model, ncol = 2) +
  scale_colour_manual(values = c("p < .05" = "#1d7870", "not significant" = "grey55"), name = NULL) +
  expand_limits(x = c(-2.6, 3.4)) +
  labs(title = "Scored against tjStuff+ v3.0, two of the three cues look real. Both are artifacts of what v3.0 leaves out",
       subtitle = paste0("Whiff above model, in percentage points, for each bin definition. The top-left panel uses Nestico's eleven v3.0 features exactly - velocity, spin rate, extension,\n",
                         "both accelerations, both release coordinates, spin axis, and the velocity and acceleration gaps to the four-seamer - with every x-dimension quantity mirrored\n",
                         "for left-handers. Against it the high slot is worth +1.91 points and high active spin costs 1.28, both comfortably significant. The other panels add back what\n",
                         "v3.0 omits, and each one removes a different result. Arm angle kills the slot cue, because tjStuff+ has no arm-slot feature and release height alone does not\n",
                         "stand in for one. Location kills the active-spin cue, so that penalty was never about spin: it is that these pitchers command the changeup worse. The spin-axis\n",
                         "gap the parachute bin is built on does nothing in any of the four. Note the R-squared column - v3.0 is a run-value model and explains 1.4 percent of whiff\n",
                         "variance on its own against 17.4 with location, so it is a demanding baseline for shape and a very weak one for whiffs."),
       x = "Whiff above model (percentage points)", y = NULL,
       caption = "Source: Statcast 2020-2026 - four-seam anchor - out-of-fold LightGBM - target is whiff on swing, features are tjStuff+ v3.0") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 8.2),
        strip.text = element_text(face = "bold", size = 9.5), panel.grid.minor = element_blank(),
        legend.position = "top", plot.caption = element_text(size = 7.5, colour = "grey40"))
ggsave(file.path(AST, "fig34_tjstuff_whiff.png"), gg, width = 11.5, height = 7.4, dpi = 150)
cat("\nwrote ext_whiff_tjstuff.csv, fig34_tjstuff_whiff.png\n")
