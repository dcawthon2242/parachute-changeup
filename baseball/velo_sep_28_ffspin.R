#!/usr/bin/env Rscript

# Four-seam spin rate per pitcher-season, pulled from the raw Statcast files.
#
# The changeup tables carry the changeup's own spin but no fastball anchor for it, the way they do
# for velocity and movement. Spin-rate difference therefore has to be built from source. Cached so
# the search script does not re-read four million rows every run.

suppressPackageStartupMessages({ library(data.table) })
MDIR <- "data/statcast_model"; OUT <- file.path(MDIR, "ff_spin_seasons.rds")
COLS <- c("pitch_type", "pitcher", "game_year", "release_spin_rate", "release_speed")

R <- rbindlist(lapply(2020:2026, function(y) {
  f <- sprintf("data/statcast_%d/statcast_%d_all.csv", y, y)
  if (!file.exists(f)) { cat(sprintf("  %d missing\n", y)); return(NULL) }
  d <- fread(f, select = COLS, showProgress = FALSE)
  d <- d[pitch_type %in% c("FF","CH") & is.finite(release_spin_rate) & is.finite(release_speed)]
  cat(sprintf("  %d: %s rows kept\n", y, format(nrow(d), big.mark = ",")))
  d
}))
setnames(R, "game_year", "season")
S <- dcast(R[, .(spin = mean(release_spin_rate), velo = mean(release_speed), n = .N),
             by = .(pitcher, season, pitch_type)],
           pitcher + season ~ pitch_type, value.var = c("spin","velo","n"))
S <- S[is.finite(spin_FF) & is.finite(spin_CH) & n_FF >= 100 & n_CH >= 40]
S[, `:=`(id = as.character(pitcher), spindiff = spin_CH - spin_FF,
         spinratio = spin_CH/spin_FF)]
# Bauer units: spin per mph, the standard way to ask whether a pitch spins a lot for its speed.
S[, `:=`(bu_ch = spin_CH/velo_CH, bu_ff = spin_FF/velo_FF)]
S[, budiff := bu_ch - bu_ff]
saveRDS(S, OUT)
cat(sprintf("\n%d pitcher-seasons with both pitches. spin difference (CH minus FF):\n", nrow(S)))
print(round(quantile(S$spindiff, c(0, .1, .25, .5, .75, .9, 1)), 0))
cat("\nBauer-unit difference (spin per mph, CH minus FF):\n")
print(round(quantile(S$budiff, c(0, .1, .25, .5, .75, .9, 1)), 2))
cat(sprintf("\nwrote %s\n", OUT))
