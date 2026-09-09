#!/usr/bin/env Rscript

# MLB AT ITS OWN TOP-THIRD ARM SLOT, AXIS SWEPT.
#
# D1 was run at 40 degrees, which is close to its top third once the proxy's compression is taken
# into account. The equivalent cut in MLB's measured distribution is 42.6. Matching on percentile
# rather than degrees is the only way to ask the same question of both leagues, since a fixed
# threshold selects 17 percent of D1 and 27 percent of MLB.
#
# MLB is the more informative side of this comparison by some distance. Every condition is measured
# rather than estimated, the pool is comparable in size, and the axis is Hawk-Eye's reading of the
# ball's rotation rather than a quantity reconstructed from the break - which in D1 turned out to be
# entangled with efficiency at r = -.37 and so was partly duplicating a condition already applied.

suppressPackageStartupMessages(library(data.table)); options(width = 205)
MDIR <- "data/statcast_model"; GATE <- 40L; EFF <- .85

F  <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(F); F <- F[pitch_type == "CH"]
AS <- readRDS(file.path(MDIR, "active_spin_long.rds")); setDT(AS)
R  <- readRDS(file.path(MDIR, "mlb_chase_resid.rds")); setDT(R)
S <- F[, .(axis = mean(axis_diff, na.rm = TRUE), arm = mean(arm_angle, na.rm = TRUE),
           vs = -mean(speed_diff, na.rm = TRUE)), by = .(pitcher, season)]
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, ec = active_spin)], by = c("pitcher","season"))
S <- merge(S, AS[pitch_type == "FF", .(pitcher, season, ef = active_spin)], by = c("pitcher","season"))
S <- merge(S, R[, .(y = 100*mean(r_all), n = .N), by = .(pitcher, season)], by = c("pitcher","season"))
S <- S[n >= GATE & is.finite(arm)]
THR <- as.numeric(quantile(S$arm, 2/3))
cat(sprintf("MLB top-third arm slot threshold: %.1f degrees (pool of %d pitcher-seasons)\n", THR, nrow(S)))

P <- S[ec >= EFF & ef >= EFF & arm >= THR]
cat(sprintf("base pool after efficiency and slot: %d seasons from %d arms, axis gap median %.1f\n",
            nrow(P), uniqueN(P$pitcher), median(P$axis)))

cat("\n=== axis gap swept inside eff >= .85 both + top-third slot ===\n")
print(rbindlist(lapply(c(8, 10, 12, 15, 20, 25, 30, 999), function(a) {
  i <- S[, ec >= EFF & ef >= EFF & arm >= THR & axis <= a]
  x <- S$y[i]; y <- S$y[!i]; if (length(x) < 3) return(NULL)
  t <- t.test(x, y)
  A <- data.table(id = S$pitcher[i], v = x)[, .(m = mean(v)), by = id]
  data.table(axis_max = if (a == 999) NA_integer_ else a, n = length(x), arms = nrow(A),
             bin = round(mean(x),2), diff = round(mean(x)-mean(y),2), p = round(t$p.value,4),
             se = round(sd(x)/sqrt(length(x)),2),
             se_by_arm = round(sd(A$m)/sqrt(nrow(A)),2)) })), row.names = FALSE)

cat("\n=== disjoint axis bands inside that same pool ===\n")
print(P[, .(n = .N, arms = uniqueN(pitcher), chase_above = round(mean(y),2),
            se = round(sd(y)/sqrt(.N),2)),
        by = .(band = cut(axis, c(-1,10,20,30,45,360), c("0-10","10-20","20-30","30-45","45+")))
       ][order(band)], row.names = FALSE)
ct <- cor.test(-P$axis, P$y)
cat(sprintf("  axis match vs chase above model inside the pool: r = %+.3f (p = %.3f, n = %d)\n",
            ct$estimate, ct$p.value, nrow(P)))

# The comparison that matters against D1: same percentile slot, same axis cut, side by side.
cat("\n=== the headline cell, both leagues at their own top third ===\n")
i <- S[, ec >= EFF & ef >= EFF & arm >= THR & axis <= 10]
t <- t.test(S$y[i], S$y[!i])
cat(sprintf("  MLB  n = %2d from %2d arms | %+.2f vs %+.2f | diff %+.2f, p = %.3f | 95%% CI %+.2f to %+.2f\n",
            sum(i), uniqueN(S$pitcher[i]), mean(S$y[i]), mean(S$y[!i]),
            mean(S$y[i])-mean(S$y[!i]), t$p.value, t$conf.int[1], t$conf.int[2]))
cat("  D1   n = 24 from 21 arms | +1.50 vs +0.04 | diff +1.46, p = 0.201 | (from the percentile sweep)\n")

# With slot and efficiency held, does anything else in the bin predict overperformance?
cat("\n=== inside the top-third pool, what does predict chase above model? ===\n")
P[, `:=`(z_axis = scale(-axis)[,1], z_arm = scale(arm)[,1], z_velo = scale(vs)[,1],
         z_ch = scale(ec)[,1], z_ff = scale(ef)[,1])]
print(round(summary(lm(y ~ z_axis + z_arm + z_velo + z_ch + z_ff, P))$coefficients, 4))
