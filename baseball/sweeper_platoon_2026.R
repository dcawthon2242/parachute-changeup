#!/usr/bin/env Rscript

# Which 2026 pitchers throw their sweeper (pitch_type "ST") MORE to
# opposite-handed hitters than to same-handed hitters?
#
# "More" is measured as usage rate -- sweepers as a share of that pitcher's
# total pitches against each batter side -- because raw counts are driven by
# how many opposite- vs same-handed batters a pitcher happens to face
# (relievers in particular see very lopsided platoon splits).
#
# Reads the cached Savant chunk files directly rather than the combined
# statcast_2026_all.csv so the script works while a backfill is still running.

suppressPackageStartupMessages({
  library(data.table)
})

chunks_dir <- file.path("data", "statcast_2026", "chunks")
out_csv    <- file.path("data", "statcast_2026", "sweeper_platoon_2026.csv")

keep_cols <- c("pitch_type", "pitch_name", "game_date", "game_type",
               "player_name", "pitcher", "stand", "p_throws")

chunk_files <- list.files(chunks_dir, pattern = "^chunk_.*\\.csv$", full.names = TRUE)
if (length(chunk_files) == 0) stop("No chunk files found in ", chunks_dir)

message(sprintf("Reading %d chunk files...", length(chunk_files)))
pitches <- rbindlist(
  lapply(chunk_files, function(f) fread(f, select = keep_cols, showProgress = FALSE)),
  use.names = TRUE, fill = TRUE
)

# Regular season only, and drop rows with no classified pitch or missing handedness.
pitches <- pitches[game_type == "R" &
                     !is.na(pitch_type) & pitch_type != "" &
                     stand %in% c("L", "R") & p_throws %in% c("L", "R")]

message(sprintf("Regular-season pitches: %s (%s to %s)",
                format(nrow(pitches), big.mark = ","),
                min(pitches$game_date), max(pitches$game_date)))

pitches[, side := ifelse(stand == p_throws, "same", "opp")]
pitches[, is_st := pitch_type == "ST"]

# Per-pitcher totals by batter side.
agg <- pitches[, .(pitches = .N, sweepers = sum(is_st)), by = .(pitcher, player_name, p_throws, side)]

wide <- dcast(agg, pitcher + player_name + p_throws ~ side,
              value.var = c("pitches", "sweepers"), fill = 0)

setnames(wide,
         c("pitches_same", "pitches_opp", "sweepers_same", "sweepers_opp"),
         c("pit_same", "pit_opp", "sw_same", "sw_opp"))

wide[, `:=`(
  sw_total  = sw_same + sw_opp,
  pit_total = pit_same + pit_opp,
  rate_same = 100 * sw_same / pit_same,
  rate_opp  = 100 * sw_opp  / pit_opp
)]
wide[, `:=`(
  usage_total = 100 * sw_total / pit_total,
  diff        = rate_opp - rate_same
)]

# Require a real sweeper and a real sample against both sides.
qualified <- wide[sw_total >= 50 & pit_same >= 100 & pit_opp >= 100]
message(sprintf("Qualified sweeper-throwers (>=50 ST, >=100 pitches vs each side): %d",
                nrow(qualified)))

reverse <- qualified[diff > 0][order(-diff)]

fwrite(qualified[order(-diff)], out_csv)

fmt <- function(dt) {
  dt[, .(
    Pitcher = player_name,
    T = p_throws,
    ST = sw_total,
    `Use%` = sprintf("%.1f", usage_total),
    `vsOpp%` = sprintf("%.1f", rate_opp),
    `vsSame%` = sprintf("%.1f", rate_same),
    Diff = sprintf("%+.1f", diff),
    `ST opp/same` = sprintf("%d/%d", sw_opp, sw_same)
  )]
}

cat("\n==== Sweeper usage HIGHER vs opposite-handed hitters (2026, qualified) ====\n")
print(fmt(reverse), nrows = 200)

cat(sprintf("\n%d of %d qualified sweeper-throwers (%.0f%%) use it more vs opposite hand.\n",
            nrow(reverse), nrow(qualified), 100 * nrow(reverse) / nrow(qualified)))

cat("\n==== For contrast: biggest same-handed skews ====\n")
print(head(fmt(qualified[order(diff)]), 10))

cat(sprintf("\nWrote %s\n", out_csv))
