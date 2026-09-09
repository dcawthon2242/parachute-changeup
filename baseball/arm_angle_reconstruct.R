#!/usr/bin/env Rscript

# CAN ARM ANGLE BE RECONSTRUCTED FROM RELEASE POINT + LISTED HEIGHT?
#
# Statcast publishes arm angle only from 2023, which is the single thing stopping the
# parachute analysis from using more seasons. Savant computes it as the angle from the
# pitcher's shoulder to the release point, and estimates the shoulder from listed height.
# That is reproducible: fit the shoulder constants on 2023-2026, where both the release
# point and the true published angle exist, then apply them backwards.
#
# The fit is judged on the thing that actually matters here, which is not RMSE in degrees
# but whether the reconstruction puts the SAME pitchers above the 44-degree line that
# defines the high-slot bin. A model with a respectable R-squared that reshuffles bin
# membership is useless for this purpose.
#
# Validation is out-of-sample by season: fit on 2023-2025, test on 2026.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
set.seed(1); options(width = 200)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")

A <- fread(file.path(MDIR, "arm_angle_tunnel.csv"))[
  , .(pitcher, season, pitch_type, arm = mean_arm, relx = mean_relx, relz = mean_relz)]
A <- A[is.finite(arm) & is.finite(relx) & is.finite(relz)]
H <- fread("data/pitcher_heights.csv")
# p_throws decides which side of the rubber the arm is on, so |relx| needs the handedness.
d <- readRDS(file.path(MDIR, "parachute_rv.rds"))
TH <- unique(d[, .(pitcher, p_throws)], by = "pitcher")
A <- merge(merge(A, H[, .(pitcher, ht_in)], by = "pitcher"), TH, by = "pitcher")
A <- A[is.finite(ht_in)]
A[, `:=`(ht_ft = ht_in/12, ax = fifelse(p_throws == "R", -relx, relx))]
cat(sprintf("%d pitcher-season-pitchtype rows with published arm angle, release point and height\n",
            nrow(A)))

TR <- A[season <= 2025]; TE <- A[season == 2026]

## ---- the physical form, with the shoulder offsets fitted -------------------------
# arm angle = angle from shoulder to release point. Shoulder sits at some fraction of
# height vertically and some lateral offset toward the throwing side.
f <- function(k, D) atan2(D$relz - k[1]*D$ht_ft, D$ax - k[2]*D$ht_ft) * 180/pi
sse <- function(k) sum((f(k, TR) - TR$arm)^2)
op <- optim(c(0.70, 0.10), sse, method = "BFGS")
cat(sprintf("\nfitted shoulder: vertical %.4f x height, lateral %.4f x height\n", op$par[1], op$par[2]))

rep_metrics <- function(pred, truth, lab) {
  r2 <- 1 - sum((pred-truth)^2)/sum((truth-mean(truth))^2)
  cat(sprintf("  %-28s R2 = %.4f   RMSE = %.2f deg   median |err| = %.2f deg\n",
              lab, r2, sqrt(mean((pred-truth)^2)), median(abs(pred-truth))))
  invisible(r2)
}
cat("\n=== held out on 2026 ===\n")
TE[, phys := f(op$par, TE)]
rep_metrics(TE$phys, TE$arm, "physical formula")

# A small linear correction on top, in case the shoulder model is systematically off for
# very tall or very short pitchers.
lmf <- lm(arm ~ phys + ht_ft + I(phys*ht_ft), data = cbind(TR, phys = f(op$par, TR)))
TE[, phys_adj := predict(lmf, TE)]
rep_metrics(TE$phys_adj, TE$arm, "+ linear height correction")

## ---- the test that matters: does the 44-degree line survive? ---------------------
cat("\n=== does the reconstruction preserve high-slot membership? (CH only, the bin's population) ===\n")
CH <- TE[pitch_type == "CH"]
for (v in c("phys","phys_adj")) {
  tp <- sum(CH$arm >= 44 & CH[[v]] >= 44); fp <- sum(CH$arm < 44 & CH[[v]] >= 44)
  fn <- sum(CH$arm >= 44 & CH[[v]] < 44); tn <- sum(CH$arm < 44 & CH[[v]] < 44)
  cat(sprintf("  %-10s agreement %.1f%%   truly high-slot recovered %d/%d   false positives %d\n",
      v, 100*(tp+tn)/nrow(CH), tp, tp+fn, fp))
}
cat(sprintf("  (n = %d changeup pitcher-seasons in 2026)\n", nrow(CH)))

## ---- where does it break? --------------------------------------------------------
TE[, err := phys_adj - arm]
cat("\n=== error by true arm angle decile ===\n")
TE[, dec := cut(arm, quantile(arm, seq(0,1,.2)), include.lowest = TRUE)]
print(TE[, .(n = .N, mean_true = round(mean(arm),1), bias = round(mean(err),2),
             sd_err = round(sd(err),2)), by = dec][order(mean_true)], row.names = FALSE)

fwrite(A, file.path(MDIR, "arm_angle_reconstruct_input.csv"))
saveRDS(list(par = op$par, lmf = lmf), file.path(MDIR, "arm_angle_model.rds"))

g <- ggplot(TE, aes(arm, phys_adj)) +
  annotate("rect", xmin = 44, xmax = Inf, ymin = 44, ymax = Inf, alpha = .1, fill = "#2a9d8f") +
  geom_abline(slope = 1, intercept = 0, colour = "#c0392b", linetype = "dashed", linewidth = .5) +
  geom_hline(yintercept = 44, colour = "grey55", linewidth = .3) +
  geom_vline(xintercept = 44, colour = "grey55", linewidth = .3) +
  geom_point(alpha = .3, size = 1.2, colour = "#1d3557") +
  labs(title = "Arm angle reconstructed from release point and listed height, held out on 2026",
       subtitle = sprintf(paste0("Shoulder position fitted on 2023-2025 as a fraction of listed height, then applied to unseen 2026 data. R-squared %.3f, RMSE %.2f degrees.\n",
                          "The shaded corner is agreement on high-slot membership, which is what the parachute bin depends on. Points in the off-diagonal corners are\n",
                          "pitchers the reconstruction would misclassify."),
                          1 - sum((TE$phys_adj-TE$arm)^2)/sum((TE$arm-mean(TE$arm))^2),
                          sqrt(mean((TE$phys_adj-TE$arm)^2))),
       x = "Published Statcast arm angle (degrees)", y = "Reconstructed arm angle (degrees)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12.5), plot.subtitle = element_text(size = 8.3),
        panel.grid.minor = element_blank())
ggsave(file.path(AST, "fig25_arm_angle_reconstruction.png"), g, width = 8.5, height = 7, dpi = 150)
cat("\nwrote fig25_arm_angle_reconstruction.png and arm_angle_model.rds\n")
