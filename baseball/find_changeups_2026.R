#!/usr/bin/env Rscript

# Find 2026 pitchers whose changeup (CH) mirrors their four-seam fastball (FF)
# movement at much lower velocity.
#
# Criteria (FF is the four-seam baseline, per pitcher):
#   - velo diff  (FF - CH) > 10 mph
#   - |HB diff|            < 3 in
#   - |IVB diff|           < 8 in
#
# HB  = pfx_x * 12 (inches), IVB = pfx_z * 12 (inches) -- standard Statcast
# induced-break convention. FF and CH share arm side per pitcher, so absolute
# inch differences are directly comparable.

suppressPackageStartupMessages({
  library(data.table)
})

# ---- Config ---------------------------------------------------------------

in_csv  <- file.path("data", "statcast_2026", "statcast_2026_all.csv")
out_csv <- file.path("data", "statcast_2026", "changeup_matches_2026.csv")

MIN_PITCHES   <- 25     # min FF and min CH thrown to be included
REGULAR_ONLY  <- TRUE   # keep only regular-season pitches (game_type == "R")

VELO_DIFF_MIN <- 10     # (ff_velo - ch_velo) must exceed this
HB_DIFF_MAX   <- 3      # abs(ff_hb - ch_hb) must be below this
IVB_DIFF_MAX  <- 12     # abs(ff_ivb - ch_ivb) must be below this (allows more drop)

# ---- Load -----------------------------------------------------------------

if (!file.exists(in_csv)) {
  stop("Raw data not found at ", in_csv, ". Run scrape_statcast_2026.R first.")
}

dt <- data.table::fread(in_csv, showProgress = FALSE)

if (REGULAR_ONLY && "game_type" %in% names(dt)) {
  dt <- dt[game_type == "R"]
}

# Keep only the pitch types and columns we need, with valid movement data.
dt <- dt[pitch_type %in% c("FF", "CH") &
           !is.na(release_speed) & !is.na(pfx_x) & !is.na(pfx_z)]

dt[, hb  := pfx_x * 12]
dt[, ivb := pfx_z * 12]

# ---- Per-pitcher, per-pitch-type summary ----------------------------------

summ <- dt[, .(
  n    = .N,
  velo = mean(release_speed),
  hb   = mean(hb),
  ivb  = mean(ivb),
  spin = mean(release_spin_rate, na.rm = TRUE)
), by = .(pitcher, player_name, p_throws, pitch_type)]

ff <- summ[pitch_type == "FF",
           .(pitcher, player_name, p_throws,
             ff_n = n, ff_velo = velo, ff_hb = hb, ff_ivb = ivb, ff_spin = spin)]
ch <- summ[pitch_type == "CH",
           .(pitcher, ch_n = n, ch_velo = velo, ch_hb = hb, ch_ivb = ivb, ch_spin = spin)]

both <- merge(ff, ch, by = "pitcher")
both <- both[ff_n >= MIN_PITCHES & ch_n >= MIN_PITCHES]

both[, velo_diff := ff_velo - ch_velo]
both[, hb_diff   := abs(ff_hb - ch_hb)]
both[, ivb_diff  := abs(ff_ivb - ch_ivb)]

# ---- Apply the criteria ---------------------------------------------------

matches <- both[velo_diff > VELO_DIFF_MIN &
                  hb_diff  < HB_DIFF_MAX &
                  ivb_diff < IVB_DIFF_MAX]

setorder(matches, -velo_diff)

# Round for readability.
num_cols  <- c("ff_velo", "ff_hb", "ff_ivb", "ch_velo", "ch_hb", "ch_ivb",
               "velo_diff", "hb_diff", "ivb_diff")
matches[, (num_cols) := lapply(.SD, round, 2), .SDcols = num_cols]
spin_cols <- c("ff_spin", "ch_spin")
matches[, (spin_cols) := lapply(.SD, round, 0), .SDcols = spin_cols]

col_order <- c("player_name", "p_throws", "ff_n", "ch_n",
               "ff_velo", "ch_velo", "velo_diff",
               "ff_hb", "ch_hb", "hb_diff",
               "ff_ivb", "ch_ivb", "ivb_diff",
               "ff_spin", "ch_spin", "pitcher")
setcolorder(matches, intersect(col_order, names(matches)))

# ---- Output ---------------------------------------------------------------

data.table::fwrite(matches, out_csv)

cat(sprintf("Pitchers evaluated (>=%d FF and >=%d CH): %d\n",
            MIN_PITCHES, MIN_PITCHES, nrow(both)))
cat(sprintf("Matches (velo>%d, |HB|<%d, |IVB|<%d): %d\n",
            VELO_DIFF_MIN, HB_DIFF_MAX, IVB_DIFF_MAX, nrow(matches)))
cat(sprintf("Wrote %s\n\n", out_csv))

print(matches[, .(player_name, p_throws, ch_n,
                  velo_diff, hb_diff, ivb_diff, ch_spin, ff_spin)])
