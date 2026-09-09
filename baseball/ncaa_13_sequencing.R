#!/usr/bin/env Rscript

# REBUILD THE D1 CHANGEUP RESIDUALS WITH SEQUENCING ATTACHED.
#
# The existing residual table has no pitch identifier, so there is no way to ask what preceded any
# given changeup. That question is the whole of M2, so the pipeline is re-run from the raw exports
# with the plate-appearance structure carried through.
#
# Two schema details. The 2023 export writes GameId and the later two write GameID, the same
# capitalisation drift already handled elsewhere in this project. And all three carry PitchUID,
# which is a genuine unique key - preferable to relying on row order surviving a re-read, since a
# silent misalignment there would corrupt every downstream sequencing result without any visible
# symptom.
#
# The whiff model is refit rather than reused because the cached residuals cannot be joined back.
# Only the location-inclusive specification is fit, since that is the one the locked spec consumes.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(5); options(width = 200)
DL <- path.expand("~/Downloads"); MDIR <- "data/statcast_model"
RAW <- file.path(MDIR, "ncaa_d1_seq_pitches.rds")
OUT <- file.path(MDIR, "ncaa_whiff_resid_seq.rds")

SCHEMA <- list(
  "2023" = list(f = "pbp23tm.csv", map = c(ax="ax", az="az", px="px", pz="pz", game_id="GameId")),
  "2024" = list(f = "D1TM24.csv",  map = c(ax="ax0", az="az0", px="PlateLocSide",
                                           pz="PlateLocHeight", game_id="GameID")),
  "2025" = list(f = "D1TM25.csv",  map = c(ax="ax0", az="az0", px="PlateLocSide",
                                           pz="PlateLocHeight", game_id="GameID")))
COMMON <- c("Level","Date","Pitcher","PitcherId","PitcherThrows","BatterSide","AutoPitchType",
            "TaggedPitchType","PitchCall","Balls","Strikes","RelSpeed","SpinRate","SpinAxis",
            "RelHeight","RelSide","Extension","InducedVertBreak","HorzBreak",
            "VertApprAngle","HorzApprAngle","PitchUID","PitchNo","PAofInning","PitchofPA",
            "Inning","BatterId")

if (!file.exists(RAW)) {
  D <- rbindlist(lapply(names(SCHEMA), function(yr) {
    s <- SCHEMA[[yr]]
    x <- fread(file.path(DL, s$f), showProgress = FALSE, select = c(COMMON, unname(s$map)))
    setnames(x, unname(s$map), names(s$map))
    x <- x[Level == "D1"]; x[, season := as.integer(yr)]
    cat(sprintf("%s: %s D1 pitches\n", yr, format(nrow(x), big.mark = ","))); x }),
    use.names = TRUE)
  D[, `:=`(mag_x = ax, mag_z = az + 32.174)]
  D[, magnus := sqrt(mag_x^2 + mag_z^2)]
  D[, ratio := magnus / (SpinRate * RelSpeed / 1000)]
  D[, swing := PitchCall %in% c("StrikeSwinging","FoulBall","FoulBallNotFieldable",
                                "FoulBallFieldable","InPlay")]
  D[, whiff := as.integer(PitchCall == "StrikeSwinging")]
  D[, pt := fifelse(!is.na(AutoPitchType) & AutoPitchType != "", AutoPitchType, TaggedPitchType)]
  saveRDS(D, RAW)
} else cat("(using cached sequenced pitch file)\n")
D <- readRDS(RAW); setDT(D)
D[, PitcherId := as.character(PitcherId)]

cat(sprintf("\nD1 pitches with sequencing: %s | unique PitchUID: %s\n",
            format(nrow(D), big.mark = ","), format(uniqueN(D$PitchUID), big.mark = ",")))
cat(sprintf("row count matches the original build (4,171,926): %s\n",
            if (nrow(D) == 4171926) "yes" else sprintf("NO - %d", nrow(D))))

## ---- the previous pitch in the same plate appearance -------------------------------------------
# A plate appearance is identified by game, inning, its index within the inning, and the batter.
# BatterId is included because PAofInning restarts each half-inning, so without it the top and
# bottom of an inning would be merged into one sequence.
setorder(D, game_id, Inning, PAofInning, PitchofPA, PitchNo)
D[, prev_pt := shift(pt), by = .(game_id, Inning, PAofInning, BatterId)]
D[, prev_velo := shift(RelSpeed), by = .(game_id, Inning, PAofInning, BatterId)]
cat(sprintf("pitches with a known predecessor in the same PA: %s (%.1f%%)\n",
            format(sum(!is.na(D$prev_pt)), big.mark = ","), 100*mean(!is.na(D$prev_pt))))
print(D[pt == "Changeup" & !is.na(prev_pt), .N, by = prev_pt][order(-N)][1:6], row.names = FALSE)

## ---- rebuild the changeup residual, carrying the key -------------------------------------------
if (!file.exists(OUT)) {
  D[, L := PitcherThrows == "Left"]
  FF <- D[pt == "Four-Seam" & is.finite(RelSpeed),
          .(nff = .N, ff_speed = mean(RelSpeed), ff_ax = mean(ax, na.rm = TRUE),
            ff_az = mean(az, na.rm = TRUE), ff_axis = mean(SpinAxis, na.rm = TRUE)),
          by = .(PitcherId, season)][nff >= 30]
  C <- merge(D[pt == "Changeup"], FF, by = c("PitcherId","season"))
  C <- C[is.finite(RelSpeed) & is.finite(ax) & is.finite(az) & is.finite(SpinAxis) &
         is.finite(px) & is.finite(pz) & is.finite(SpinRate) & is.finite(Extension) &
         is.finite(VertApprAngle) & is.finite(HorzApprAngle)]
  C[, `:=`(speed_diff = RelSpeed - ff_speed, ax_diff = ax - ff_ax, az_diff = az - ff_az,
           axis_diff = pmin(abs(SpinAxis - ff_axis), 360 - abs(SpinAxis - ff_axis)))]
  C[, `:=`(tj_x0 = fifelse(L, -RelSide, RelSide), tj_ax = fifelse(L, -ax, ax),
           tj_ax_diff = fifelse(L, -ax_diff, ax_diff), tj_px = fifelse(L, -px, px),
           tj_haa = fifelse(L, -HorzApprAngle, HorzApprAngle),
           tj_axis = fifelse(L, (360 - SpinAxis) %% 360, SpinAxis),
           same_hand = as.integer((PitcherThrows == "Left") == (BatterSide == "Left")))]
  FEAT <- c("RelSpeed","SpinRate","Extension","tj_ax","az","tj_x0","RelHeight","tj_axis",
            "speed_diff","tj_ax_diff","az_diff",
            "tj_px","pz","VertApprAngle","tj_haa","same_hand","Balls","Strikes")
  SW <- C[swing == TRUE]
  cat(sprintf("\nchangeup swings for the model: %s   whiff rate %.3f\n",
              format(nrow(SW), big.mark = ","), mean(SW$whiff)))
  K <- 4; fold <- sample(rep(1:K, length.out = nrow(SW))); p <- rep(NA_real_, nrow(SW))
  for (f in 1:K) {
    tr <- SW[fold != f]; n <- nrow(tr); vi <- sample(n, floor(.12*n))
    dtr <- lgb.Dataset(as.matrix(tr[-vi, ..FEAT]), label = tr$whiff[-vi])
    dva <- lgb.Dataset.create.valid(dtr, as.matrix(tr[vi, ..FEAT]), label = tr$whiff[vi])
    m <- lgb.train(params = list(objective = "binary", metric = "binary_logloss",
                   learning_rate = .06, num_leaves = 31, min_data_in_leaf = 300,
                   feature_fraction = .8, bagging_fraction = .8, bagging_freq = 1),
                   data = dtr, nrounds = 1500, valids = list(v = dva),
                   early_stopping_rounds = 50, verbose = -1)
    p[fold == f] <- predict(m, as.matrix(SW[fold == f, ..FEAT]))
  }
  SW[, r_all := whiff - p]
  cat(sprintf("  out-of-fold R2 = %.4f\n", 1 - var(SW$r_all)/var(SW$whiff)))
  saveRDS(SW[, .(PitchUID, PitcherId, Pitcher, season, game_id, Inning, PAofInning, PitchofPA,
                 BatterId, BatterSide, Balls, Strikes, prev_pt, prev_velo, RelSpeed,
                 axis_diff, speed_diff, whiff, r_all)], OUT)
} else cat("\n(using cached sequenced residuals)\n")
R <- readRDS(OUT); setDT(R)
cat(sprintf("\nsequenced changeup swings: %s | with a known predecessor: %s (%.1f%%)\n",
            format(nrow(R), big.mark = ","), format(sum(!is.na(R$prev_pt)), big.mark = ","),
            100*mean(!is.na(R$prev_pt))))

## ---- does the refit reproduce the cached residual? ---------------------------------------------
# The model is refit here, so the season-level numbers must be checked against the ones the locked
# spec was built on. A large divergence would mean the two pipelines are not the same analysis.
OLD <- readRDS(file.path(MDIR, "ncaa_whiff_resid.rds")); setDT(OLD)
OLD[, PitcherId := as.character(PitcherId)]
A <- OLD[, .(old = 100*mean(r_all), n_old = .N), by = .(PitcherId, season)][n_old >= 40]
B <- R[, .(new = 100*mean(r_all), n_new = .N), by = .(PitcherId, season)][n_new >= 40]
M <- merge(A, B, by = c("PitcherId","season"))
cat(sprintf("\nagreement with the cached pipeline: %d shared pitcher-seasons, r = %.4f, mean abs diff %.3f pts\n",
            nrow(M), cor(M$old, M$new), mean(abs(M$old - M$new))))
cat(sprintf("wrote %s\n", OUT))
