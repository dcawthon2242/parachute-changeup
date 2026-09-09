#!/usr/bin/env Rscript

# ARM ANGLE RECONSTRUCTION, SECOND ATTEMPT: PER PITCH RATHER THAN PER PITCHER-SEASON.
#
# The first attempt fitted season means and got R2 = 0.61, which is useless - it recovered
# only half the truly high-slot changeup seasons. That failure was self-inflicted: arm angle
# is an arctangent of the release point, and the mean of an arctangent is not the arctangent
# of the means, so averaging first threw away the relationship. Statcast carries a per-pitch
# arm_angle column, so this fits the geometry where it actually lives.
#
# If Savant computes arm angle from release point and a height-derived shoulder, a per-pitch
# fit should be nearly deterministic. Anything much below R2 = 0.95 means the calculation uses
# body-tracking information that release point cannot stand in for, and back-casting to
# pre-2023 seasons would be guesswork.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
set.seed(1); options(width = 200)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")

COLS <- c("pitcher","p_throws","arm_angle","release_pos_x","release_pos_z","release_pos_y",
          "release_extension","pitch_type")
d <- rbindlist(lapply(2023:2026, function(y) {
  k <- fread(sprintf("data/statcast_%d/statcast_%d_all.csv", y, y), select = COLS,
             showProgress = FALSE); k[, season := y]; k }))
H <- fread("data/pitcher_heights.csv")
d <- merge(d[is.finite(arm_angle) & is.finite(release_pos_x) & is.finite(release_pos_z)],
           H[, .(pitcher, ht_in)], by = "pitcher")
d <- d[is.finite(ht_in)]
d[, `:=`(ht_ft = ht_in/12, ax = fifelse(p_throws == "R", -release_pos_x, release_pos_x))]
cat(sprintf("%s pitches with per-pitch arm angle, release point and height\n",
            format(nrow(d), big.mark = ",")))

TR <- d[season <= 2025][sample(.N, min(.N, 600000))]; TE <- d[season == 2026]

## ---- physical form, shoulder constants fitted -------------------------------------
f <- function(k, D) atan2(D$release_pos_z - k[1]*D$ht_ft, D$ax - k[2]*D$ht_ft) * 180/pi
op <- optim(c(0.70, 0.10), function(k) sum((f(k, TR) - TR$arm_angle)^2), method = "BFGS")
cat(sprintf("fitted shoulder: vertical %.4f x height, lateral %.4f x height\n", op$par[1], op$par[2]))

met <- function(p, t, lab) {
  cat(sprintf("  %-30s R2 = %.4f   RMSE = %.2f deg   median |err| = %.2f deg\n", lab,
              1 - sum((p-t)^2)/sum((t-mean(t))^2), sqrt(mean((p-t)^2)), median(abs(p-t))))
  invisible(1 - sum((p-t)^2)/sum((t-mean(t))^2))
}
cat("\n=== per-pitch accuracy, held out on 2026 ===\n")
TE[, phys := f(op$par, TE)]
met(TE$phys, TE$arm_angle, "physical formula")

## ---- flexible model, in case the shoulder is not a pure height fraction -----------
suppressPackageStartupMessages(library(lightgbm))
FEAT <- c("release_pos_x","release_pos_z","release_pos_y","release_extension","ht_ft","ax")
TR[, thr := as.integer(p_throws == "R")]; TE[, thr := as.integer(p_throws == "R")]
FEAT <- c(FEAT, "thr")
TRc <- TR[stats::complete.cases(TR[, ..FEAT])]; TEc <- TE[stats::complete.cases(TE[, ..FEAT])]
vi <- sample(nrow(TRc), floor(.1*nrow(TRc)))
dtr <- lgb.Dataset(as.matrix(TRc[-vi, ..FEAT]), label = TRc$arm_angle[-vi])
dva <- lgb.Dataset.create.valid(dtr, as.matrix(TRc[vi, ..FEAT]), label = TRc$arm_angle[vi])
m <- lgb.train(params = list(objective="regression", metric="l2", learning_rate=.05,
               num_leaves=63, min_data_in_leaf=200, feature_fraction=.9, bagging_fraction=.8,
               bagging_freq=1), data = dtr, nrounds = 3000, valids = list(v=dva),
               early_stopping_rounds = 60, verbose = -1)
TEc[, gbm := predict(m, as.matrix(TEc[, ..FEAT]))]
met(TEc$gbm, TEc$arm_angle, "gradient boosting on geometry")

## ---- the test that matters: pitcher-season classification -------------------------
cat("\n=== aggregate to pitcher-season changeups and check the 44-degree line ===\n")
S <- TEc[pitch_type == "CH", .(n = .N, arm = mean(arm_angle), phys = mean(phys),
                               gbm = mean(gbm)), by = pitcher][n >= 60]
for (v in c("phys","gbm")) {
  tp <- sum(S$arm >= 44 & S[[v]] >= 44); fp <- sum(S$arm < 44 & S[[v]] >= 44)
  fn <- sum(S$arm >= 44 & S[[v]] < 44); tn <- sum(S$arm < 44 & S[[v]] < 44)
  cat(sprintf("  %-5s agreement %.1f%%   high-slot recovered %d/%d   false positives %d   season-mean R2 %.3f\n",
      v, 100*(tp+tn)/nrow(S), tp, tp+fn, fp,
      1 - sum((S[[v]]-S$arm)^2)/sum((S$arm-mean(S$arm))^2)))
}
cat(sprintf("  (n = %d changeup pitcher-seasons in 2026 with >= 60 pitches)\n", nrow(S)))

saveRDS(list(par = op$par, gbm = m, feat = FEAT), file.path(MDIR, "arm_angle_model2.rds"))
g <- ggplot(S, aes(arm, gbm)) +
  annotate("rect", xmin = 44, xmax = Inf, ymin = 44, ymax = Inf, alpha = .1, fill = "#2a9d8f") +
  geom_abline(slope = 1, intercept = 0, colour = "#c0392b", linetype = "dashed", linewidth = .5) +
  geom_hline(yintercept = 44, colour = "grey55", linewidth = .3) +
  geom_vline(xintercept = 44, colour = "grey55", linewidth = .3) +
  geom_point(alpha = .5, size = 1.6, colour = "#1d3557") +
  labs(title = "Arm angle reconstructed per pitch from release point and listed height",
       subtitle = sprintf(paste0("Fitted on 2023-2025 pitches, held out on 2026, then averaged to the pitcher-season level shown here. Season-mean R-squared %.3f.\n",
                          "The shaded corner is agreement on high-slot membership, which is the only property the parachute bin needs."),
                          1 - sum((S$gbm-S$arm)^2)/sum((S$arm-mean(S$arm))^2)),
       x = "Published Statcast arm angle (degrees)", y = "Reconstructed arm angle (degrees)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12.5), plot.subtitle = element_text(size = 8.3),
        panel.grid.minor = element_blank())
ggsave(file.path(AST, "fig25_arm_angle_reconstruction.png"), g, width = 8.5, height = 7, dpi = 150)
cat("\nwrote fig25_arm_angle_reconstruction.png and arm_angle_model2.rds\n")
