#!/usr/bin/env Rscript

# WHY THE SLOT INTERACTION DIED: ONE CHANGE AT A TIME.
#
# The axis x arm-slot interaction is beta = -0.0062, p = .002 (p = .005 arm-clustered) on the
# sample parachute_slot_audit.R used, and beta = -0.15 z-units, p = .48 on the frozen pool. Seven
# things differ between those two tables. Attributing the collapse to "the four-seam filter" or
# "the velocity adjustment" without checking is guessing, and two of the seven are differences
# nobody chose on purpose - the outcome model and the arm-angle source both changed underneath
# this result while it was being used as the project's one surviving finding.
#
#   arm source   arm_angle_tunnel.csv (2023+, measured) vs parachute_ff (2020+, backcast pre-2023)
#   axis source  parachute_within vs parachute_ff
#   outcome      parachute_within SWD wh_res vs mlb_whiff_locaware r_all - different whiff models
#   swing gate   60 vs 40
#   population   all arms vs four-seam primary only
#   residual     raw vs season-level velocity-adjusted
#   span         2023-2026 vs 2020-2026
#
# Two passes. Cumulative walks from the original to the frozen spec applying one change at a time.
# Leave-one-out starts at the frozen spec and reverts one change at a time, which catches a change
# that only matters in the presence of the others.

suppressPackageStartupMessages({ library(data.table) })
options(width = 205); set.seed(23)
MDIR <- "data/statcast_model"

crob <- function(f, D, key = "id") {
  D <- D[complete.cases(D[, c(all.vars(f), key), with = FALSE])]
  m <- lm(f, data = D); u <- residuals(m); X <- model.matrix(m); cl <- D[[key]]
  nc <- uniqueN(cl); b <- solve(crossprod(X))
  V <- b %*% crossprod(rowsum(X * u, cl)) %*% b * (nc/(nc-1)) * ((nrow(X)-1)/(nrow(X)-ncol(X)))
  k <- "axis:arm"; est <- coef(m)[k]; se <- sqrt(diag(V))[k]
  list(est = est, p = 2*pt(-abs(est/se), nc-1), n = nrow(X), arms = nc,
       zbeta = est * sd(D$axis) * sd(D$arm))
}

## ---- master table, every ingredient side by side ----------------------------------------------
F <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(F); FC <- F[pitch_type == "CH"]
M <- FC[, .(axis_ff = mean(axis_diff, na.rm = TRUE), arm_ff = mean(arm_angle, na.rm = TRUE),
            vs = -mean(speed_diff, na.rm = TRUE)), by = .(pitcher, season)]

RW <- readRDS(file.path(MDIR, "mlb_whiff_locaware.rds")); setDT(RW)
M <- merge(M, RW[, .(y_lock = 100*mean(r_all), nsw_lock = .N), by = .(pitcher, season)],
           by = c("pitcher","season"), all.x = TRUE)

LW <- readRDS(file.path(MDIR, "parachute_within.rds"))
CW <- as.data.table(LW$CH)[, .(axis_w = mean(axis_diff, na.rm = TRUE)), by = .(pitcher, season)]
SW <- as.data.table(LW$SWD)[, .(y_orig = 100*mean(wh_res), nsw_orig = .N), by = .(pitcher, season)]
M <- merge(M, CW, by = c("pitcher","season"), all.x = TRUE)
M <- merge(M, SW, by = c("pitcher","season"), all.x = TRUE)

aa <- unique(fread(file.path(MDIR, "arm_angle_tunnel.csv"))[pitch_type == "CH",
        .(pitcher, season, arm_tun = mean_arm)], by = c("pitcher","season"))
M <- merge(M, aa, by = c("pitcher","season"), all.x = TRUE)

AR <- rbindlist(lapply(2020:2026, function(y) {
  x <- fread(sprintf("data/savant/arsenal_%d.csv", y), showProgress = FALSE)
  setnames(x, names(x)[2], "pitcher"); x[, season := y]
  x[, .(pitcher, season, ff_use = suppressWarnings(as.numeric(n_ff)),
        si_use = suppressWarnings(as.numeric(n_si)))] }))
AR[is.na(ff_use), ff_use := 0][is.na(si_use), si_use := 0]
M <- merge(M, AR, by = c("pitcher","season"), all.x = TRUE)
M[, id := as.character(pitcher)]

## ---- one config -> one row ----------------------------------------------------------------------
run <- function(arm = "tun", axis = "w", out = "orig", gate = 60L,
                ffprim = FALSE, veloadj = FALSE, span = "2023") {
  D <- copy(M)
  D[, arm := if (arm == "tun") arm_tun else arm_ff]
  D[, axis := if (axis == "w") axis_w else axis_ff]
  D[, `:=`(y = if (out == "orig") y_orig else y_lock,
           nsw = if (out == "orig") nsw_orig else nsw_lock)]
  if (span == "2023") D <- D[season >= 2023]
  D <- D[is.finite(arm) & is.finite(axis) & is.finite(y) & is.finite(nsw) & nsw >= gate]
  if (ffprim) D <- D[is.finite(ff_use) & is.finite(si_use) & ff_use >= si_use]
  if (veloadj) D[, y := resid(lm(y ~ vs))]
  crob(y ~ axis * arm, D)
}

show <- function(tag, r) cat(sprintf("  %-46s beta %+.5f | z-beta %+.3f | p %.4f | %4d seasons, %3d arms\n",
                                     tag, r$est, r$zbeta, r$p, r$n, r$arms))

BASE <- list(arm = "tun", axis = "w",  out = "orig", gate = 60L, ffprim = FALSE, veloadj = FALSE, span = "2023")
FROZ <- list(arm = "ff",  axis = "ff", out = "lock", gate = 40L, ffprim = TRUE,  veloadj = TRUE,  span = "all")
STEPS <- c(arm = "arm angle from parachute_ff (backcast-capable)",
           axis = "axis gap from parachute_ff",
           out = "outcome = mlb_whiff_locaware model",
           gate = "swing gate 60 -> 40",
           ffprim = "four-seam-primary population filter",
           veloadj = "season velocity-adjusted residual",
           span = "span 2023-2026 -> 2020-2026")

cat("=== the two endpoints ===\n")
show("ORIGINAL (parachute_slot_audit)", do.call(run, BASE))
show("FROZEN (locked_spec conventions)", do.call(run, FROZ))

cat("\n=== cumulative: original -> frozen, one change at a time ===\n")
cfg <- BASE
for (s in names(STEPS)) { cfg[[s]] <- FROZ[[s]]; show(paste0("+ ", STEPS[[s]]), do.call(run, cfg)) }

cat("\n=== leave-one-out: frozen, with a single change reverted ===\n")
for (s in names(STEPS)) {
  cfg <- FROZ; cfg[[s]] <- BASE[[s]]
  show(paste0("revert ", STEPS[[s]]), do.call(run, cfg))
}

## ---- the two that are not analyst choices -------------------------------------------------------
cat("\n=== how much do the two silently-swapped inputs actually disagree? ===\n")
C <- M[season >= 2023 & is.finite(arm_tun) & is.finite(arm_ff)]
cat(sprintf("  arm angle:  r = %.4f over %d pitcher-seasons, mean abs gap %.2f deg\n",
            cor(C$arm_tun, C$arm_ff), nrow(C), mean(abs(C$arm_tun - C$arm_ff))))
C2 <- M[season >= 2023 & is.finite(y_orig) & is.finite(y_lock)]
cat(sprintf("  whiff resid: r = %.4f over %d pitcher-seasons (sd %.2f vs %.2f)\n",
            cor(C2$y_orig, C2$y_lock), nrow(C2), sd(C2$y_orig), sd(C2$y_lock)))
B <- M[season < 2023 & is.finite(arm_ff)]
cat(sprintf("  pre-2023 seasons carrying a backcast arm angle: %d of %d frozen-era MLB seasons\n",
            nrow(B), M[is.finite(arm_ff) & is.finite(y_lock), .N]))
