#!/usr/bin/env Rscript

# HOW MANY DISTINCT ARMS, NOT PITCHER-SEASONS, SIT IN EACH BIN?
#
# Every significance test so far has treated pitcher-seasons as the unit, which quietly assumes they
# are independent. They are not: the same pitcher recurring across four seasons contributes four
# rows that share his mechanics, his repertoire and most of his catchers. The distinct-arm count is
# the honest measure of how much independent evidence exists, and where it is much smaller than the
# season count, the reported standard errors are optimistic.

suppressPackageStartupMessages(library(data.table)); options(width = 190)
MDIR <- "data/statcast_model"

## ---- MLB -------------------------------------------------------------------------------------
F  <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(F); F <- F[pitch_type == "CH"]
AS <- readRDS(file.path(MDIR, "active_spin_long.rds")); setDT(AS)
R  <- readRDS(file.path(MDIR, "mlb_chase_resid.rds")); setDT(R)
S <- F[, .(axis = mean(axis_diff, na.rm = TRUE), arm = mean(arm_angle, na.rm = TRUE)),
       by = .(pitcher, season)]
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, as_ch = active_spin)], by = c("pitcher","season"))
S <- merge(S, AS[pitch_type == "FF", .(pitcher, season, as_fb = active_spin)], by = c("pitcher","season"))
S <- merge(S, R[, .(c_all = 100*mean(r_all), n = .N), by = .(pitcher, season)], by = c("pitcher","season"))
S <- S[n >= 40 & is.finite(arm)]

mlb <- function(lab, i) data.table(league = "MLB", cut = lab, seasons = sum(i),
                                   arms = uniqueN(S$pitcher[i]),
                                   seasons_per_arm = round(sum(i)/uniqueN(S$pitcher[i]), 2))
OUT <- rbind(
  mlb("qualifying pool",         rep(TRUE, nrow(S))),
  mlb("arm >= 44",               S$arm >= 44),
  mlb("axis <= 10",              S$axis <= 10),
  mlb("eff .85 + arm 44",        S$as_ch >= .85 & S$as_fb >= .85 & S$arm >= 44),
  mlb("full bin (.85+arm+axis)", S$as_ch >= .85 & S$as_fb >= .85 & S$arm >= 44 & S$axis <= 10),
  mlb("full bin at .90",         S$as_ch >= .90 & S$as_fb >= .90 & S$arm >= 44 & S$axis <= 10))

## ---- D1 --------------------------------------------------------------------------------------
# The same two tables the chase script joins - ncaa_parachute90.rds is a narrower, already-screened
# set and would understate the bin.
ARM <- as.data.table(readRDS(file.path(MDIR, "ncaa_armangle.rds")))
PR  <- as.data.table(readRDS(file.path(MDIR, "ncaa_spineff_pairs.rds")))
PR  <- merge(PR, ARM[, .(PitcherId, season, arm_hat)], by = c("PitcherId","season"))
setnames(PR, "axis_gap", "axis", skip_absent = TRUE)
CR <- readRDS(file.path(MDIR, "ncaa_chase_resid.rds")); setDT(CR)
D <- merge(PR, CR[!is.na(PitcherId), .(c_all = 100*mean(r_all, na.rm = TRUE), n = .N),
                  by = .(PitcherId, season)], by = c("PitcherId","season"))[n >= 40]

d1 <- function(lab, i) data.table(league = "D1", cut = lab, seasons = sum(i),
                                  arms = uniqueN(D$PitcherId[i]),
                                  seasons_per_arm = round(sum(i)/max(uniqueN(D$PitcherId[i]),1), 2))
OUT <- rbind(OUT,
  d1("qualifying pool",         rep(TRUE, nrow(D))),
  d1("arm >= 44",               D$arm_hat >= 44),
  d1("axis <= 10",              D$axis <= 10),
  d1("eff .85 + arm 44",        D$eff_ch >= .85 & D$eff_ff >= .85 & D$arm_hat >= 44),
  d1("full bin (.85+arm+axis)", D$eff_ch >= .85 & D$eff_ff >= .85 & D$arm_hat >= 44 & D$axis <= 10),
  d1("full bin at .90",         D$eff_ch >= .90 & D$eff_ff >= .90 & D$arm_hat >= 44 & D$axis <= 10))
print(OUT, row.names = FALSE)

## ---- what the repetition costs -----------------------------------------------------------------
# If the seasons inside a bin are near-copies of one another, the effective sample is closer to the
# arm count. Clustering the standard error by pitcher shows how much of the reported precision is
# real. A ratio well above 1 means the pitcher-season test was overstating its confidence.
cat("\n=== standard error, clustered by pitcher, on the full .85 bin ===\n")
for (lg in c("MLB","D1")) {
  X <- if (lg == "MLB") S[as_ch >= .85 & as_fb >= .85 & arm >= 44 & axis <= 10,
                          .(id = pitcher, y = c_all)]
       else D[eff_ch >= .85 & eff_ff >= .85 & arm_hat >= 44 & axis <= 10,
              .(id = PitcherId, y = c_all)]
  if (!nrow(X)) next
  naive <- sd(X$y)/sqrt(nrow(X))
  A <- X[, .(m = mean(y), k = .N), by = id]
  clus <- if (nrow(A) > 1) sd(A$m)/sqrt(nrow(A)) else NA_real_
  cat(sprintf("  %-3s  mean %+.2f  |  %d seasons from %d arms  |  naive SE %.2f, by-arm SE %.2f (x%.2f)\n",
              lg, mean(X$y), nrow(X), nrow(A), naive, clus, clus/naive))
}
