#!/usr/bin/env Rscript

# WHICH OF THE FOUR CONDITIONS IS ACTUALLY DOING THE WORK?
#
# The bin test says the axis condition carries the chase effect and the other three mostly cost
# sample. That is worth believing only after the obvious confound is ruled out. A tight spin-axis
# gap between changeup and four-seamer is not independent of arm slot - an over-the-top pitcher
# releases both pitches from a similar orientation, so slot and axis gap are entangled by
# construction. The earlier whiff version of this bin turned out to be an arm-angle artifact for
# exactly that reason, and there is no reason to assume chase behaves differently.
#
# Three checks: each condition alone, the axis split held inside arm bands, and all cues entered
# jointly. The within-pitcher test is the strictest - it asks whether the same pitcher chases above
# model in his own seasons where the axis gap happens to be tighter, which no stable pitcher trait
# can explain.

suppressPackageStartupMessages(library(data.table)); options(width = 200)
MDIR <- "data/statcast_model"
F  <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(F); F <- F[pitch_type == "CH"]
AS <- readRDS(file.path(MDIR, "active_spin_long.rds")); setDT(AS)
R  <- readRDS(file.path(MDIR, "mlb_chase_resid.rds")); setDT(R)

S <- F[, .(axis = mean(axis_diff, na.rm = TRUE), arm = mean(arm_angle, na.rm = TRUE),
           velo_sep = -mean(speed_diff, na.rm = TRUE)), by = .(pitcher, season)]
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, as_ch = active_spin)], by = c("pitcher","season"))
S <- merge(S, AS[pitch_type == "FF", .(pitcher, season, as_fb = active_spin)], by = c("pitcher","season"))
S <- merge(S, R[, .(c_all = 100*mean(r_all), n = .N), by = .(pitcher, season)], by = c("pitcher","season"))
S <- S[n >= 40 & is.finite(arm)]
cat(sprintf("n = %d pitcher-seasons\n", nrow(S)))

cat("\n=== each condition on its own ===\n")
one <- function(l, i) { x <- S$c_all[i]; y <- S$c_all[!i]; t <- t.test(x, y)
  data.table(cut = l, n = length(x), bin = round(mean(x),2),
             diff = round(mean(x)-mean(y),2), p = round(t$p.value,4)) }
print(rbind(one("arm >= 44 only",     S$arm >= 44),
            one("axis <= 10 only",    S$axis <= 10),
            one("eff .85 both only",  S$as_ch >= .85 & S$as_fb >= .85)), row.names = FALSE)

cat("\n=== axis split held inside arm bands ===\n")
print(S[, .(n = .N, chase_above = round(mean(c_all),2), se = round(sd(c_all)/sqrt(.N),2)),
        by = .(arm_band = cut(arm, c(0,35,44,90), c("<35","35-44","44+")),
               axis = fifelse(axis <= 10, "axis<=10", "axis>10"))][order(arm_band, axis)],
      row.names = FALSE)

cat("\n=== all cues entered jointly, standardized ===\n")
S[, `:=`(z_arm = scale(arm)[,1], z_axis = scale(-axis)[,1], z_ch = scale(as_ch)[,1],
         z_fb = scale(as_fb)[,1], z_velo = scale(velo_sep)[,1])]
print(round(summary(lm(c_all ~ z_arm + z_axis + z_ch + z_fb + z_velo, S))$coefficients, 4))
cat(sprintf("\ncorrelation between axis gap and arm slot: %+.3f\n", cor(S$axis, S$arm)))

cat("\n=== within-pitcher: his own tighter-axis seasons ===\n")
M <- S[, if (.N >= 3) .(dax = axis - mean(axis), dc = c_all - mean(c_all)), by = pitcher]
ct <- cor.test(-M$dax, M$dc)
cat(sprintf("  %d season-deviations from %d pitchers: r = %+.3f (p = %.3f)\n",
            nrow(M), uniqueN(M$pitcher), ct$estimate, ct$p.value))
