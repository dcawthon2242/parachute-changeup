#!/usr/bin/env Rscript

# Compare the run value of the "fastball-mirroring" changeups (the pitchers
# flagged by find_changeups_2026.R) to the league-average changeup in 2026.
#
# Run value convention (matches Baseball Savant's pitch arsenal leaderboard):
#   RV      = -sum(delta_run_exp)          # runs saved from the pitcher's view
#   RV/100  = RV / pitches * 100           # positive = good for the pitcher
# delta_run_exp is the change in run expectancy on each pitch (batter's view),
# so negating it gives the pitcher's run value.

suppressPackageStartupMessages({
  library(data.table)
})

in_csv       <- file.path("data", "statcast_2026", "statcast_2026_all.csv")
matches_csv  <- file.path("data", "statcast_2026", "changeup_matches_2026.csv")
out_csv      <- file.path("data", "statcast_2026", "changeup_run_value_2026.csv")

REGULAR_ONLY <- TRUE

if (!file.exists(in_csv))      stop("Missing ", in_csv, "; run scrape_statcast_2026.R first.")
if (!file.exists(matches_csv)) stop("Missing ", matches_csv, "; run find_changeups_2026.R first.")

dt <- data.table::fread(
  in_csv,
  select = c("pitch_type", "game_type", "pitcher", "player_name", "delta_run_exp"),
  showProgress = FALSE
)

if (REGULAR_ONLY) dt <- dt[game_type == "R"]

# All changeups with a valid run-value figure.
ch <- dt[pitch_type == "CH" & !is.na(delta_run_exp)]

rv_summary <- function(d) {
  n  <- nrow(d)
  rv <- -sum(d$delta_run_exp)
  data.table(
    pitches   = n,
    run_value = round(rv, 1),
    rv_per_100 = round(rv / n * 100, 2)
  )
}

# ---- League-average changeup ----------------------------------------------

league <- rv_summary(ch)

# ---- The flagged pitchers' changeups --------------------------------------

matches <- data.table::fread(matches_csv, showProgress = FALSE)
match_ids <- matches$pitcher

group <- rv_summary(ch[pitcher %in% match_ids])

per_pitcher <- ch[pitcher %in% match_ids,
  {
    n <- .N; rv <- -sum(delta_run_exp)
    .(ch_pitches = n, run_value = round(rv, 1), rv_per_100 = round(rv / n * 100, 2))
  },
  by = .(pitcher, player_name)]
per_pitcher <- merge(per_pitcher,
                     matches[, .(pitcher, p_throws, velo_diff, hb_diff, ivb_diff)],
                     by = "pitcher", all.x = TRUE)
setorder(per_pitcher, -rv_per_100)

# ---- Output ---------------------------------------------------------------

cat("=== 2026 Changeup Run Value (pitcher perspective; + = runs saved) ===\n\n")
cat(sprintf("League-average CH:   %d pitches | RV %+.1f | RV/100 %+.2f\n",
            league$pitches, league$run_value, league$rv_per_100))
cat(sprintf("Flagged 4 CH pooled: %d pitches | RV %+.1f | RV/100 %+.2f\n\n",
            group$pitches, group$run_value, group$rv_per_100))
cat(sprintf("Difference vs league (RV/100): %+.2f\n\n", group$rv_per_100 - league$rv_per_100))

cat("Per pitcher:\n")
print(per_pitcher[, .(player_name, p_throws, ch_pitches, run_value, rv_per_100,
                      velo_diff, hb_diff, ivb_diff)])

out <- rbindlist(list(
  data.table(player_name = "LEAGUE AVERAGE CH", p_throws = NA, ch_pitches = league$pitches,
             run_value = league$run_value, rv_per_100 = league$rv_per_100,
             velo_diff = NA, hb_diff = NA, ivb_diff = NA, pitcher = NA),
  data.table(player_name = "FLAGGED 4 POOLED", p_throws = NA, ch_pitches = group$pitches,
             run_value = group$run_value, rv_per_100 = group$rv_per_100,
             velo_diff = NA, hb_diff = NA, ivb_diff = NA, pitcher = NA),
  per_pitcher[, .(player_name, p_throws, ch_pitches, run_value, rv_per_100,
                  velo_diff, hb_diff, ivb_diff, pitcher)]
), use.names = TRUE, fill = TRUE)
data.table::fwrite(out, out_csv)
cat(sprintf("\nWrote %s\n", out_csv))
