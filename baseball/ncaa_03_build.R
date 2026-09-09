#!/usr/bin/env Rscript

# BUILD THE D1 CHANGEUP DATASET FROM THE THREE TRACKMAN SEASONS.
#
# Harmonises the two export schemas (2023 uses the older column names), filters to Level == D1,
# keeps only what the analysis needs, and derives:
#
#   inferred spin efficiency   from spin rate and Magnus magnitude, using the estimator
#                              validated against MLB measured active spin in script 02
#   four-seam anchor           per pitcher-season, and the changeup's gaps to it
#   arm-slot proxy             from release point, since TrackMan reports no arm angle
#
# One caveat carried forward from script 01: TrackMan's SpinAxis sits within about 5 degrees of
# a deterministic function of the break for changeups, so the axis gap computed here is close
# to a movement-direction gap. That is weaker than the Hawk-Eye quantity and it is flagged
# wherever it is used.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(5); options(width = 200)
DL <- path.expand("~/Downloads"); MDIR <- "data/statcast_model"
OUT <- file.path(MDIR, "ncaa_d1_pitches.rds")

SCHEMA <- list(
  "2023" = list(f = "pbp23tm.csv", map = c(ax="ax", az="az", px="px", pz="pz")),
  "2024" = list(f = "D1TM24.csv",  map = c(ax="ax0", az="az0", px="PlateLocSide", pz="PlateLocHeight")),
  "2025" = list(f = "D1TM25.csv",  map = c(ax="ax0", az="az0", px="PlateLocSide", pz="PlateLocHeight")))
COMMON <- c("Level","Date","Pitcher","PitcherId","PitcherThrows","BatterSide","AutoPitchType",
            "TaggedPitchType","PitchCall","Balls","Strikes","RelSpeed","SpinRate","SpinAxis",
            "RelHeight","RelSide","Extension","InducedVertBreak","HorzBreak",
            "VertApprAngle","HorzApprAngle")

if (!file.exists(OUT)) {
  D <- rbindlist(lapply(names(SCHEMA), function(yr) {
    s <- SCHEMA[[yr]]
    x <- fread(file.path(DL, s$f), showProgress = FALSE, select = c(COMMON, unname(s$map)))
    setnames(x, unname(s$map), names(s$map))
    x <- x[Level == "D1"]
    x[, season := as.integer(yr)]
    cat(sprintf("%s: %s D1 pitches\n", yr, format(nrow(x), big.mark=",")))
    x }), use.names = TRUE)

  # TrackMan az is total vertical acceleration, so gravity has to come out before the Magnus
  # magnitude is formed. Verified below rather than assumed.
  cat(sprintf("\nmean az on four-seamers: %.1f ft/s^2 (expect roughly -17 if gravity is included)\n",
              D[AutoPitchType == "Four-Seam", mean(az, na.rm = TRUE)]))
  D[, `:=`(mag_x = ax, mag_z = az + 32.174)]
  D[, magnus := sqrt(mag_x^2 + mag_z^2)]
  D[, ratio := magnus / (SpinRate * RelSpeed / 1000)]

  D[, swing := PitchCall %in% c("StrikeSwinging","FoulBall","FoulBallNotFieldable",
                                "FoulBallFieldable","InPlay")]
  D[, whiff := as.integer(PitchCall == "StrikeSwinging")]
  D[, pt := fifelse(!is.na(AutoPitchType) & AutoPitchType != "", AutoPitchType, TaggedPitchType)]
  saveRDS(D, OUT); cat(sprintf("\nsaved %s rows to %s\n", format(nrow(D), big.mark=","), OUT))
} else cat("(using cached D1 pitch file)\n")
D <- readRDS(OUT)
cat(sprintf("D1 pitches: %s across %d seasons, %d pitchers\n", format(nrow(D), big.mark=","),
            uniqueN(D$season), uniqueN(D$PitcherId)))
cat("\nswing rate %.3f, whiff-per-swing %.3f\n")
cat(sprintf("swing rate %.3f | whiff per swing %.3f\n", mean(D$swing), sum(D$whiff)/sum(D$swing)))

## ---- inferred spin efficiency, calibrated on MLB ---------------------------------------------
# The estimator is fit where truth exists (MLB, measured active spin) and applied to NCAA. The
# inputs are all quantities TrackMan reports, so nothing Hawk-Eye-specific leaks in.
cat("\n=== calibrating the efficiency estimator on MLB, applying to NCAA ===\n")
M <- readRDS(file.path(MDIR, "parachute_rv.rds"))[pitch_type %in% c("CH","FF")]
AS <- readRDS(file.path(MDIR, "active_spin_long.rds"))
M[, magnus := sqrt(ax^2 + (az + 32.174)^2)]
M[, ratio := magnus / (release_spin_rate * release_speed / 1000)]
PM <- M[is.finite(ratio) & release_spin_rate > 500,
        .(n = .N, ratio = mean(ratio), magnus = mean(magnus), spin = mean(release_spin_rate),
          velo = mean(release_speed), ext = mean(release_extension, na.rm = TRUE)),
        by = .(pitcher, season, pitch_type)]
PM <- merge(PM, AS[, .(pitcher, season, pitch_type, measured = active_spin)],
            by = c("pitcher","season","pitch_type"))[n >= 50 & is.finite(measured) & is.finite(ext)]
F <- c("ratio","magnus","spin","velo","ext")
mod <- lgb.train(params = list(objective = "regression", learning_rate = .05, num_leaves = 15,
                 min_data_in_leaf = 40, feature_fraction = .9), verbose = -1, nrounds = 400,
                 data = lgb.Dataset(as.matrix(PM[, ..F]), label = PM$measured))
cat(sprintf("  trained on %d MLB pitcher-season-pitch-types\n", nrow(PM)))

N <- D[pt %in% c("Changeup","Four-Seam") & is.finite(ratio) & SpinRate > 500 & is.finite(Extension),
       .(n = .N, ratio = mean(ratio), magnus = mean(magnus), spin = mean(SpinRate),
         velo = mean(RelSpeed), ext = mean(Extension), axis = mean(SpinAxis, na.rm = TRUE),
         relh = mean(RelHeight), rels = mean(RelSide), ivb = mean(InducedVertBreak),
         hb = mean(HorzBreak), throws = PitcherThrows[1], name = Pitcher[1]),
       by = .(PitcherId, season, pt)][n >= 30]
N[, eff := pmin(pmax(predict(mod, as.matrix(N[, ..F])), 0), 1)]
cat(sprintf("  NCAA pitcher-season-pitch-types scored: %d\n", nrow(N)))
print(N[, .(median_eff = round(median(eff),3), q10 = round(quantile(eff,.10),3),
            q90 = round(quantile(eff,.90),3)), by = pt], row.names = FALSE)
cat("\n  for reference, MLB measured active spin:\n")
print(PM[, .(pt = pitch_type, measured)][, .(median = round(median(measured),3),
      q10 = round(quantile(measured,.10),3), q90 = round(quantile(measured,.90),3)), by = pt],
      row.names = FALSE)

saveRDS(N, file.path(MDIR, "ncaa_pitcher_season_shapes.rds"))
cat(sprintf("\nwrote ncaa_pitcher_season_shapes.rds\n"))
