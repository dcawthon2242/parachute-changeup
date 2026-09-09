#!/usr/bin/env Rscript

# LOOSENING THE ARM SLOT, WITH THE AXIS AND EFFICIENCY CONDITIONS HELD.
#
# The axis condition is the expensive one and it is not negotiable - it is the whole hypothesis. Arm
# slot is the negotiable one: 44 degrees was chosen because it is roughly the league median, not
# because anything identified it. Sweeping it shows whether the effect is genuinely a property of
# over-the-top arms or whether 44 was simply where a small sample happened to look good.
#
# The diagnostic to watch is not the p-value, which will improve automatically as the sample grows.
# It is whether the effect SIZE decays smoothly as lower slots are admitted. A real slot dependence
# decays; a spurious one jumps around, and a threshold that was pure sample-size luck will show the
# effect collapsing the moment the bin grows.

suppressPackageStartupMessages(library(data.table)); options(width = 200)
MDIR <- "data/statcast_model"; GATE <- 40L; EFF <- .85; AXIS <- 10

ARM <- as.data.table(readRDS(file.path(MDIR, "ncaa_armangle.rds")))
PR  <- as.data.table(readRDS(file.path(MDIR, "ncaa_spineff_pairs.rds")))
PR  <- merge(PR, ARM[, .(PitcherId, season, arm_hat)], by = c("PitcherId","season"))
setnames(PR, "axis_gap", "axis", skip_absent = TRUE)
CR  <- readRDS(file.path(MDIR, "ncaa_chase_resid.rds")); setDT(CR)
D <- merge(PR, CR[!is.na(PitcherId), .(y = 100*mean(r_all, na.rm = TRUE), n = .N),
                  by = .(PitcherId, season)], by = c("PitcherId","season"))[n >= GATE]
D[, `:=`(id = PitcherId, arm = arm_hat)]

F  <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(F); F <- F[pitch_type == "CH"]
AS <- readRDS(file.path(MDIR, "active_spin_long.rds")); setDT(AS)
R  <- readRDS(file.path(MDIR, "mlb_chase_resid.rds")); setDT(R)
S <- F[, .(axis = mean(axis_diff, na.rm = TRUE), arm = mean(arm_angle, na.rm = TRUE)),
       by = .(pitcher, season)]
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, eff_ch = active_spin)], by = c("pitcher","season"))
S <- merge(S, AS[pitch_type == "FF", .(pitcher, season, eff_ff = active_spin)], by = c("pitcher","season"))
S <- merge(S, R[, .(y = 100*mean(r_all), n = .N), by = .(pitcher, season)], by = c("pitcher","season"))
S <- S[n >= GATE & is.finite(arm)][, id := pitcher]

sweep <- function(X, lg) {
  cat(sprintf("\n=== %s: eff >= %.2f both, axis gap <= %d, arm slot swept (gate %d, pool %d) ===\n",
              lg, EFF, AXIS, GATE, nrow(X)))
  cat(sprintf("  arm slot in this pool: median %.1f, quartiles %.1f / %.1f\n",
              median(X$arm), quantile(X$arm,.25), quantile(X$arm,.75)))
  print(rbindlist(lapply(c(0, 20, 25, 30, 35, 38, 41, 44, 47, 50), function(a) {
    i <- X[, eff_ch >= EFF & eff_ff >= EFF & axis <= AXIS & arm >= a]
    x <- X$y[i]; y <- X$y[!i]
    if (length(x) < 3) return(data.table(arm_min = a, n = length(x), arms = uniqueN(X$id[i]),
                                         bin = NA_real_, diff = NA_real_, p = NA_real_, se = NA_real_))
    t <- t.test(x, y)
    # Clustering by arm matters here because loosening the threshold pulls in repeat seasons from
    # pitchers already present, which inflates n without adding independent evidence.
    A <- data.table(id = X$id[i], v = x)[, .(m = mean(v)), by = id]
    data.table(arm_min = a, n = length(x), arms = nrow(A), bin = round(mean(x),2),
               diff = round(mean(x)-mean(y),2), p = round(t$p.value,4),
               se = round(sd(x)/sqrt(length(x)),2),
               se_by_arm = if (nrow(A) > 1) round(sd(A$m)/sqrt(nrow(A)),2) else NA_real_) })),
    row.names = FALSE)
}
sweep(D, "D1"); sweep(S, "MLB")

# The axis condition on its own, with no slot requirement at all, is the cleanest statement of the
# hypothesis and the largest sample it can have.
cat("\n=== axis gap alone, no efficiency and no slot condition ===\n")
for (lg in c("D1","MLB")) { X <- if (lg == "D1") D else S
  print(rbindlist(lapply(c(8, 10, 12, 15, 20), function(a) {
    i <- X$axis <= a; x <- X$y[i]; y <- X$y[!i]
    if (length(x) < 3) return(NULL); t <- t.test(x, y)
    data.table(league = lg, axis_max = a, n = length(x), arms = uniqueN(X$id[i]),
               bin = round(mean(x),2), diff = round(mean(x)-mean(y),2),
               p = round(t$p.value,4)) })), row.names = FALSE) }
