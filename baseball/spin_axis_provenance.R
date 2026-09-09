#!/usr/bin/env Rscript

# WHERE DOES spin_axis COME FROM?
#
# The first audit found the reported axis deviates from the movement-implied angle by a
# median of 14 deg. Two explanations, with opposite consequences:
#
#   A. spin_axis is an independent Hawk-Eye measurement of the ball's actual rotation.
#      Then the deviation is real signal and "spin similarity" is a genuine spin cue.
#   B. spin_axis IS derived from movement and the deviation is my mapping being sloppy
#      (constant offset) plus angle noise on low-movement pitches.
#
# B predicts: after removing the best constant offset, deviation shrinks toward zero and
# whatever remains scales as 1/movement - big only where the movement vector is tiny.
# A predicts: deviation survives on high-movement pitches too, and tracks spin efficiency,
# because that is exactly where an observed axis and an inferred one part ways.

suppressPackageStartupMessages({ library(data.table) })
options(width = 200)

COLS <- c("pitch_type","pitcher","p_throws","release_speed","pfx_x","pfx_z","ax","az",
          "spin_axis","release_spin_rate")
d <- fread("data/statcast_2025/statcast_2025_all.csv", select = COLS, showProgress = FALSE)
d <- d[is.finite(spin_axis) & is.finite(pfx_x) & is.finite(ax) &
       pitch_type %in% c("FF","SI","FC","SL","ST","CU","KC","CH","FS")]

d[, `:=`(mag_pfx = 12*sqrt(pfx_x^2 + pfx_z^2),                      # inches of movement
         imp = (atan2(az + 32.174, ax)*180/pi + 90) %% 360)]        # movement-implied axis
sgn <- function(a, b) { e <- (a - b + 180) %% 360 - 180; e }
d[, raw_dev := sgn(spin_axis, imp)]

# Best constant offset per pitch type, circularly. If B is right this is where the
# deviation goes.
off <- d[, .(offset = (atan2(mean(sin(raw_dev*pi/180)), mean(cos(raw_dev*pi/180)))*180/pi)),
         by = pitch_type]
d <- merge(d, off, by = "pitch_type")
d[, dev := abs(sgn(raw_dev, offset))]

cat("=== BEST CONSTANT OFFSET PER PITCH TYPE, AND WHAT IS LEFT AFTER REMOVING IT ===\n")
print(d[, .(n = .N, offset = round(offset[1],1), median_abs_dev = round(median(dev),2),
            p90 = round(quantile(dev,.90),1), pct_within_3deg = round(100*mean(dev<3),1),
            median_movement_in = round(median(mag_pfx),1)),
        by = pitch_type][order(median_abs_dev)], row.names = FALSE)
cat("\nA single global offset would suffice if spin_axis were computed from movement.\n",
    "Spread of the per-type offsets: ", round(diff(range(off$offset)),1), " deg\n", sep="")

## ---- does the leftover deviation scale like angle noise (1/movement)? ---------
cat("\n=== DEVIATION BY MOVEMENT MAGNITUDE (the noise test) ===\n")
cat("If spin_axis were movement-derived, deviation would vanish wherever movement is large.\n\n")
d[, mbin := cut(mag_pfx, c(0,4,8,12,16,20,100),
                labels = c("0-4 in","4-8","8-12","12-16","16-20","20+ in"))]
print(dcast(d[pitch_type %in% c("FF","SI","SL","ST","CU","CH")],
            pitch_type ~ mbin, value.var = "dev",
            fun.aggregate = function(z) round(median(z),1)), row.names = FALSE)

cat("\n=== THE DECISIVE CELL: 4-SEAMERS WITH 15+ INCHES OF MOVEMENT ===\n")
z <- d[pitch_type == "FF" & mag_pfx >= 15]
cat("  n =", nrow(z), " median deviation =", round(median(z$dev),2), "deg,",
    round(100*mean(z$dev < 3),1), "% within 3 deg\n")
cat("  A movement-derived field would read 0.00 here. The movement vector is 15+ inches\n",
    "  long, so there is no angle noise to hide behind.\n", sep="")

## ---- does the deviation track spin efficiency? --------------------------------
as_l <- readRDS("data/statcast_model/active_spin_long.rds"); setDT(as_l)
pa <- d[, .(n = .N, dev = median(dev), mag = median(mag_pfx)), by = .(pitcher, pitch_type)][n >= 200]
pa <- merge(pa, as_l[season == 2025, .(pitcher, pitch_type, active_spin)],
            by = c("pitcher","pitch_type"))
cat("\n=== DOES THE DEVIATION TRACK SPIN EFFICIENCY? ===\n")
cat("Observed and inferred axes agree when spin is efficient and part ways when it is gyro.\n\n")
for (ty in c("FF","SI","SL","ST","CU","CH")) {
  x <- pa[pitch_type == ty]
  if (nrow(x) < 25) next
  ct <- suppressWarnings(cor.test(x$active_spin, x$dev, method = "spearman", exact = FALSE))
  cat(sprintf("  %-3s n=%3d  median active spin %3.0f%%  r(active spin, deviation) = %+.3f (p=%.2g)\n",
              ty, nrow(x), 100*median(x$active_spin, na.rm=TRUE), ct$estimate, ct$p.value))
}
cat("\nNegative r = the less efficient the spin, the further the reported axis sits from\n",
    "the direction the ball actually breaks. That gap cannot exist in a derived field.\n", sep="")
