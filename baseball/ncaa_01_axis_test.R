#!/usr/bin/env Rscript

# IS TRACKMAN'S SpinAxis AN INDEPENDENT MEASUREMENT, OR IS IT DERIVED FROM THE BREAK?
#
# This decides whether the parachute analysis ports at all. The MLB work rests on Hawk-Eye
# imaging the seams, so a pitch can rotate on one axis and move along another; the discrepancy
# between them IS the gyro component. If TrackMan instead computes its axis from the observed
# movement, a "matched spin axis" bin is a matched-movement bin, and movement is already in the
# pitch-grade model - so the bin would select on something the model has fully priced.
#
# METHOD. Comparing two angle conventions directly is a trap: an arbitrary zero point or a
# handedness reflection produces a large constant offset that says nothing about whether the
# two carry the same information. So no convention is assumed. For each dataset the best
# rigid alignment is solved for - both a rotation and an optional reflection - and what gets
# reported is the CIRCULAR SPREAD that remains afterwards. That residual is the question:
# near zero means the axis is a relabelling of the break, and large means it is carrying
# something the break does not.

suppressPackageStartupMessages({ library(data.table) })
options(width = 200)
DL <- path.expand("~/Downloads")

SCHEMA <- list(
  "2023" = list(f = "pbp23tm.csv", pfxx = "pfx_x", pfxz = "pfx_z"),
  "2024" = list(f = "D1TM24.csv",  pfxx = "pfxx",  pfxz = "pfxz"),
  "2025" = list(f = "D1TM25.csv",  pfxx = "pfxx",  pfxz = "pfxz"))

# Circular spread of a - b after the best rotation, checking both orientations. Returns the
# residual sd in degrees (0 = perfectly determined, ~104 = uniform/no relationship).
align <- function(a, b) {
  best <- NULL
  for (s in c(1, -1)) {
    d <- (s*a - b) * pi/180
    Rbar <- Mod(mean(exp(1i*d)))                       # resultant length: 1 = identical, 0 = uniform
    sdc <- sqrt(-2*log(max(Rbar, 1e-12))) * 180/pi
    if (is.null(best) || sdc < best$sd) best <- list(sd = sdc, R = Rbar, flip = s,
      off = (Arg(mean(exp(1i*d)))*180/pi) %% 360)
  }
  best
}
rpt <- function(a, b, lab, n) { z <- align(a, b)
  cat(sprintf("  %-34s n=%9s  residual spread %6.1f deg   R=%.3f  %s\n", lab,
      format(n, big.mark=","), z$sd, z$R, if (z$flip < 0) "(reflected)" else "")) }

cat("=== NCAA TrackMan: reported SpinAxis vs the axis implied by the break ===\n")
POOL <- rbindlist(lapply(names(SCHEMA), function(yr) {
  s <- SCHEMA[[yr]]
  D <- fread(file.path(DL, s$f), showProgress = FALSE,
             select = c("Level","AutoPitchType","PitcherThrows","SpinAxis","SpinRate",
                        "RelSpeed", s$pfxx, s$pfxz))
  setnames(D, c(s$pfxx, s$pfxz), c("pfx_x","pfx_z"))
  D <- D[Level == "D1" & is.finite(SpinAxis) & is.finite(pfx_x) & is.finite(pfx_z) &
         (abs(pfx_x) + abs(pfx_z)) > .05]
  D[, mv := (atan2(pfx_x, pfx_z)*180/pi) %% 360]
  rpt(D$SpinAxis, D$mv, sprintf("%s, all D1 pitches", yr), nrow(D))
  cbind(season = yr, D[, .(AutoPitchType, PitcherThrows, SpinAxis, mv, SpinRate, RelSpeed)]) }))

cat("\n  by pitch type (pooled). A derived column is uniformly tight; a measured one loosens\n")
cat("  on the pitch types that carry the most gyro spin:\n")
for (pt in POOL[!is.na(AutoPitchType), .N, by = AutoPitchType][order(-N)]$AutoPitchType) {
  d <- POOL[AutoPitchType == pt]; rpt(d$SpinAxis, d$mv, paste("   ", pt), nrow(d)) }

cat("\n=== MLB Hawk-Eye, identical treatment, as the calibration ===\n")
M <- readRDS("data/statcast_model/parachute_ff.rds")
M[, spin_axis := (atan2(sax, cax)*180/pi) %% 360]
# Total az includes gravity; removing it leaves the Magnus component, which is what the break
# direction encodes.
M <- M[is.finite(spin_axis) & is.finite(ax) & is.finite(az)]
M[, mv := (atan2(ax, az + 32.174)*180/pi) %% 360]
rpt(M$spin_axis, M$mv, "MLB changeups (Hawk-Eye)", nrow(M))

# The decisive MLB check: if the axis-to-break discrepancy is real physics rather than noise,
# it must track measured active spin. A derived axis cannot possibly do this.
M[, d := ((spin_axis - mv + 180) %% 360) - 180]
AS <- readRDS("data/statcast_model/active_spin_long.rds")
MM <- merge(M[, .(gap = mean(abs(d)), n = .N), by = .(pitcher, season)],
            AS[pitch_type == "CH", .(pitcher, season, active_spin)], by = c("pitcher","season"))[n >= 50]
ct <- cor.test(MM$gap, MM$active_spin)
cat(sprintf("\n  MLB: axis-to-break discrepancy vs MEASURED active spin: r = %+.3f (p = %.2g, n = %d)\n",
            ct$estimate, ct$p.value, nrow(MM)))
cat("  A larger discrepancy going with lower active spin is the signature of a genuinely\n")
cat("  independent axis measurement. This is the property TrackMan has to reproduce.\n")

cat("\n=== can spin efficiency be recovered from TrackMan at all? ===\n")
cat(sprintf("columns matching eff/active/gyro in the NCAA files: %s\n",
    paste(setdiff(grep("eff|Eff|ctive|yro", names(fread(file.path(DL,"D1TM25.csv"), nrows=1)),
                       value=TRUE), character(0)), collapse=", ")))
