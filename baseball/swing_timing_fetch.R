#!/usr/bin/env Rscript

# Pull Savant's Swing Timing + Miss Distance leaderboard, pitcher view, split by
# pitch type and batter side, one file per season.
#
# Why the leaderboard and not the per-pitch feed. Savant defines the nine timing
# categories on the BARREL-to-ball distance at closest approach:
#   tied-up / centered / flail   centered  = barrel within +/-4 in horizontally
#   late / on-time / early       on time   = barrel within +/-7 MILLISECONDS
#   over / lined-up / under      lined up  = barrel within +/-2 in vertically
# None of those three components is published per pitch. The per-pitch feed has
# only the scalar miss_distance (whiffs only) plus the ball's intercept position
# relative to the BATTER'S BODY, which is a different quantity roughly 30 in out
# front and 40 in to the side -- not a barrel offset. So the categories cannot be
# reconstructed pitch by pitch, and the pre-aggregated leaderboard is the only
# Savant-matched source. Confirmed against the published conditional means: the
# +/-4, +/-7, +/-2 thresholds imply tail means of -6.0/+7.9, +/-10.6 and
# -3.4/+2.9, which is what the leaderboard reports.
#
# The undocumented query shape, recovered from the page's own links, is
#   ?type=pitcher&season[]=YYYY&split[]=api_pitch_type_group03&split[]=bat_side
# Three traps, all of which fail quietly rather than erroring.
#
# `year=` is ignored. The parameter is `season[]`, and passing `year` returns the
# CURRENT season for every request, so a four-season pull looks like it worked and
# is really the same file four times.
#
# `min` defaults to `q` (qualified), which drops about half the rows. `min=25` is
# the LOOSER setting and is used here to get the superset. It does not filter on
# n_swings, which still runs down to a single swing, so a sample floor belongs
# downstream rather than being trusted from the feed.
#
# `bat_side_formatted` is NOT the batter's side. In the pitcher view it carries the
# PITCHER's throwing hand (100% agreement with p_throws; its 74/26 R/L marginal is
# the league split of right-handed pitchers). The raw `bat_side` is the real split
# key. Joining on the formatted column pairs every row with the wrong platoon side
# while still matching, so both columns are preserved here and consumers must pick
# `bat_side` deliberately.
#
# Output: data/swing_timing/swing_timing_pitcher_<season>.csv

suppressPackageStartupMessages({ library(data.table) })

SEASONS <- 2023:2026   # swing timing begins in the second half of 2023
OUT_DIR <- file.path("data", "swing_timing")
BASE <- paste0("https://baseballsavant.mlb.com/leaderboard/bat-tracking/",
               "swing-timing-miss-distance")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

for (yr in SEASONS) {
  url <- sprintf(paste0("%s?type=pitcher&season%%5B%%5D=%d",
                        "&split%%5B%%5D=api_pitch_type_group03",
                        "&split%%5B%%5D=bat_side&min=25&csv=true"), BASE, yr)
  dest <- file.path(OUT_DIR, sprintf("swing_timing_pitcher_%d.csv", yr))
  # Savant serves the CSV only to a browser-ish agent.
  code <- system2("curl", c("-s", "-H", shQuote("User-Agent: Mozilla/5.0"),
                            shQuote(url), "-o", shQuote(dest)))
  if (code != 0 || !file.exists(dest)) { cat(sprintf("  %d FAILED\n", yr)); next }

  d <- fread(dest, showProgress = FALSE)
  setnames(d, gsub("^\ufeff", "", names(d)))
  # A silently-changed schema would poison everything downstream, so fail loudly.
  need <- c("id", "name", "year", "api_pitch_type", "bat_side_formatted", "n_swings",
            "tied_up_percent", "centered_percent", "flailed_percent",
            "early_percent", "on_time_percent", "late_percent",
            "over_percent", "lined_up_percent", "under_percent",
            "miss_distance", "whiff_rate", "flawed_percent", "perfect_percent")
  miss <- setdiff(need, names(d))
  if (length(miss)) stop(sprintf("%d: leaderboard is missing %s", yr,
                                 paste(miss, collapse = ", ")))
  if (!all(d$year == yr)) stop(sprintf("%d: got seasons %s -- season[] was ignored",
                                       yr, paste(unique(d$year), collapse = "/")))
  fwrite(d, dest)
  cat(sprintf("  %d  rows %5d  swings %9s  -> %s\n", yr, nrow(d),
              format(sum(d$n_swings), big.mark = ","), basename(dest)))
}
cat("done\n")
