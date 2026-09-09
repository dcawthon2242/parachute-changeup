#!/usr/bin/env Rscript

# RECONNAISSANCE ON THE NCAA TRACKMAN FILES.
#
# One question decides whether the parachute analysis is portable at all. Everything in the MLB
# work rests on spin axis being an INDEPENDENT measurement of how the ball rotates - Hawk-Eye
# images the seams, so a pitch can spin on one axis and move on another, and the difference
# between those two is exactly what spin efficiency measures.
#
# TrackMan is Doppler radar. It measures spin RATE directly, but the axis it reports may be
# inferred from the observed break, in which case "spin axis" carries no information that the
# movement columns do not already carry, and a matched-axis bin is really a matched-movement
# bin wearing different clothes. That is a much weaker claim and it has to be established
# before anything is built on top of it.
#
# The test: reconstruct the axis implied by pfx_x and pfx_z and compare it to the reported
# SpinAxis. Hawk-Eye MLB data shows a loose relationship - that looseness IS the gyro component.
# If TrackMan shows a near-deterministic one, the column is derived.

suppressPackageStartupMessages({ library(data.table) })
options(width = 200)
DL <- "~/Downloads"
FILES <- c("2023" = "pbp23tm.csv", "2024" = "D1TM24.csv", "2025" = "D1TM25.csv")

for (yr in names(FILES)) {
  f <- file.path(DL, FILES[yr])
  cat(sprintf("\n================ %s : %s ================\n", yr, FILES[yr]))
  # nrows keeps this cheap; the goal is structure, not statistics.
  D <- fread(f, nrows = 400000, showProgress = FALSE)
  cat(sprintf("sampled %s rows of the file\n", format(nrow(D), big.mark=",")))
  cat("\nLevel values: "); print(D[, .N, by = Level][order(-N)])
  cat("\nPitchCall values:\n"); print(D[, .N, by = PitchCall][order(-N)])
  cat("\nAutoPitchType:\n"); print(D[, .N, by = AutoPitchType][order(-N)][1:12])
  cat(sprintf("\nany spin-efficiency-like column? %s\n",
      paste(grep("eff|Eff|active|Active|gyro|Gyro", names(D), value = TRUE), collapse = ", ")))
  cat(sprintf("SpinAxis non-missing: %.1f%%   SpinRate: %.1f%%   Extension: %.1f%%\n",
      100*mean(!is.na(D$SpinAxis)), 100*mean(!is.na(D$SpinRate)), 100*mean(!is.na(D$Extension))))
}

## ---- the decisive test, on 2025 D1 only ------------------------------------------------------
cat("\n\n================ IS SpinAxis DERIVED FROM MOVEMENT? ================\n")
D <- fread(file.path(DL, FILES["2025"]), showProgress = FALSE, nrows = 1500000,
           select = c("Level","AutoPitchType","PitcherThrows","SpinAxis","SpinRate",
                      "pfx_x","pfx_z","InducedVertBreak","HorzBreak","RelSpeed"))
D <- D[Level == "D1" & is.finite(SpinAxis) & is.finite(pfx_x) & is.finite(pfx_z)]
cat(sprintf("D1 pitches with axis and movement: %s\n", format(nrow(D), big.mark=",")))

# Statcast convention: 180 degrees is pure backspin. The movement-implied axis is the direction
# the Magnus force points, rotated to that same convention.
D[, mv_axis := (atan2(-pfx_x, -pfx_z)*180/pi) %% 360]
D[, d := ((SpinAxis - mv_axis + 180) %% 360) - 180]
cat(sprintf("\nSpinAxis vs movement-implied axis:  median gap %+.2f deg   IQR %.2f to %.2f   sd %.2f\n",
            median(D$d), quantile(D$d,.25), quantile(D$d,.75), sd(D$d)))
cat(sprintf("share within 5 degrees: %.1f%%   within 10: %.1f%%   within 20: %.1f%%\n",
            100*mean(abs(D$d) <= 5), 100*mean(abs(D$d) <= 10), 100*mean(abs(D$d) <= 20)))
cat("\nby pitch type (a real measurement should disagree most on gyro-heavy pitches):\n")
print(D[, .(n = .N, median_gap = round(median(d),1), sd_gap = round(sd(d),1),
            within10 = round(100*mean(abs(d) <= 10),1)), by = AutoPitchType][order(-n)][1:10],
      row.names = FALSE)

cat("\n--- the same test on MLB Hawk-Eye, for calibration ---\n")
M <- readRDS("data/statcast_model/parachute_ff.rds")
M[, spin_axis := (atan2(sax, cax)*180/pi) %% 360]
M <- M[is.finite(spin_axis) & is.finite(ax) & is.finite(az)]
# Statcast pfx is not stored here, but ax/az carry the same direction information.
M[, mv := (atan2(-ax, -(az + 32.174))*180/pi) %% 360]
M[, d := ((spin_axis - mv + 180) %% 360) - 180]
cat(sprintf("MLB changeups: median gap %+.2f, sd %.2f, within 10 deg %.1f%%\n",
            median(M$d), sd(M$d), 100*mean(abs(M$d) <= 10)))
