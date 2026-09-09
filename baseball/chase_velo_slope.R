#!/usr/bin/env Rscript

# DOES VELOCITY SEPARATION STILL BUY CHASE ONCE YOU ARE INSIDE THE BIN?
#
# The chase model already has speed_diff as a feature, so on the full population the residual should
# carry no velocity slope at all - if it did, the model would simply be underfitting velocity. The
# interesting question is conditional: inside the matched-spin, high-slot subset, does extra velo
# separation pay off MORE than the model expects? That is what a parachute claim actually predicts.
# A changeup that looks identical out of the hand is supposed to convert its speed gap into swing
# decisions more efficiently than a changeup the hitter can pick up early.
#
# The test that matters is therefore not the within-bin correlation on its own - with a dozen
# seasons that number is nearly unreadable - but the INTERACTION against everyone else. A steep
# slope inside the bin means nothing if the slope outside is just as steep.
#
# Both leagues are run side by side because they fail differently: D1's bin is small but its arms
# are nearly all distinct, while MLB's is small AND repeats pitchers.

suppressPackageStartupMessages(library(data.table)); options(width = 200)
MDIR <- "data/statcast_model"
GATE <- 40L

## ---- D1 ----------------------------------------------------------------------------------------
ARM <- as.data.table(readRDS(file.path(MDIR, "ncaa_armangle.rds")))
PR  <- as.data.table(readRDS(file.path(MDIR, "ncaa_spineff_pairs.rds")))
PR  <- merge(PR, ARM[, .(PitcherId, season, arm_hat)], by = c("PitcherId","season"))
setnames(PR, "axis_gap", "axis", skip_absent = TRUE)
CR  <- readRDS(file.path(MDIR, "ncaa_chase_resid.rds")); setDT(CR)
D <- merge(PR, CR[!is.na(PitcherId), .(y = 100*mean(r_all, na.rm = TRUE), n = .N),
                  by = .(PitcherId, season)], by = c("PitcherId","season"))[n >= GATE]
D[, `:=`(id = PitcherId, vs = velo_sep,
         b_eff  = eff_ch >= .85 & eff_ff >= .85,
         b_arm  = eff_ch >= .85 & eff_ff >= .85 & arm_hat >= 44,
         b_full = eff_ch >= .85 & eff_ff >= .85 & arm_hat >= 44 & axis <= 10)]

## ---- MLB ---------------------------------------------------------------------------------------
F  <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(F); F <- F[pitch_type == "CH"]
AS <- readRDS(file.path(MDIR, "active_spin_long.rds")); setDT(AS)
R  <- readRDS(file.path(MDIR, "mlb_chase_resid.rds")); setDT(R)
S <- F[, .(axis = mean(axis_diff, na.rm = TRUE), arm_hat = mean(arm_angle, na.rm = TRUE),
           vs = -mean(speed_diff, na.rm = TRUE)), by = .(pitcher, season)]
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, eff_ch = active_spin)], by = c("pitcher","season"))
S <- merge(S, AS[pitch_type == "FF", .(pitcher, season, eff_ff = active_spin)], by = c("pitcher","season"))
S <- merge(S, R[, .(y = 100*mean(r_all), n = .N), by = .(pitcher, season)], by = c("pitcher","season"))
S <- S[n >= GATE & is.finite(arm_hat)]
S[, `:=`(id = pitcher,
         b_eff  = eff_ch >= .85 & eff_ff >= .85,
         b_arm  = eff_ch >= .85 & eff_ff >= .85 & arm_hat >= 44,
         b_full = eff_ch >= .85 & eff_ff >= .85 & arm_hat >= 44 & axis <= 10)]

## ---- slopes -------------------------------------------------------------------------------------
slope <- function(X, lab, col) {
  I <- X[[col]]; A <- X[I]; B <- X[!I]
  if (nrow(A) < 4) return(data.table(cut = lab, n = nrow(A), r = NA_real_, slope = NA_real_,
                                     p_in = NA_real_, slope_out = NA_real_, p_interaction = NA_real_))
  ct <- cor.test(A$vs, A$y)
  fi <- summary(lm(y ~ vs, A))$coefficients
  # The baseline row compares the population to nothing, so there is no outside group to contrast.
  has_out <- nrow(B) >= 4
  fo <- if (has_out) summary(lm(y ~ vs, B))$coefficients[, "Estimate"]["vs"] else NA_real_
  ix <- if (has_out) { m <- summary(lm(y ~ vs*I, X))$coefficients; m[nrow(m), "Pr(>|t|)"] } else NA_real_
  data.table(cut = lab, n = nrow(A), r = round(ct$estimate,3), slope = round(fi["vs","Estimate"],3),
             p_in = round(ct$p.value,4), slope_out = round(fo,3), p_interaction = round(ix,4))
}
for (lg in c("D1","MLB")) {
  X <- if (lg == "D1") D else S
  cat(sprintf("\n=== %s: velocity separation vs chase above model (gate %d, n = %d) ===\n",
              lg, GATE, nrow(X)))
  cat(sprintf("  velo separation: median %.1f mph, sd %.1f\n", median(X$vs), sd(X$vs)))
  X[, ok := TRUE]
  print(rbind(slope(X, "everyone (baseline)", "ok"),
              slope(X, "eff .85 both",        "b_eff"),
              slope(X, "eff + arm 44",        "b_arm"),
              slope(X, "full bin + axis<=10", "b_full")), row.names = FALSE)
}

## ---- the bin's own seasons, listed ---------------------------------------------------------------
# With a dozen points a correlation is one observation away from changing sign, so the raw pairs are
# worth seeing rather than trusting the summary statistic.
cat("\n=== D1 full-bin seasons: velo separation against chase above model ===\n")
print(D[b_full == TRUE][order(-vs), .(Pitcher = name, Season = season, OOZ_CH = n,
        VeloSep = round(vs,1), ArmSlot = round(arm_hat,1), AxisGap = round(axis,1),
        Chase_above = round(y,2))], row.names = FALSE)
