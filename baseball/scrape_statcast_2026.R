#!/usr/bin/env Rscript

# Scrape the full 2026 MLB Statcast pitch-level dataset from Baseball Savant.
#
# We hit Savant's statcast_search CSV endpoint directly (the same source the
# baseballr package uses) rather than baseballr::statcast_search(), because
# baseballr 1.6.0 hardcodes a column-name vector that no longer matches
# Savant's current 119-column schema ("Can't assign 92 names to a 119-column
# data.table"). Reading the CSV directly lets the header define the columns.
#
# Savant caps a single response at ~25k rows and a busy MLB day is ~4-5k
# pitches, so we request the season in small (~4-day) windows and stitch the
# chunks together.

# ---- Packages -------------------------------------------------------------

required_pkgs <- c("dplyr", "data.table", "curl")
missing_pkgs <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]
if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(curl)
})

# Savant blocks/rate-limits requests without a browser-like User-Agent, so we
# always send one.
USER_AGENT <- paste0(
  "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) ",
  "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/126.0 Safari/537.36"
)

# ---- Config ---------------------------------------------------------------

season_start <- as.Date("2026-03-18")  # skip deep-spring windows (no served data)
season_end   <- Sys.Date()             # today (2026-07-31 at time of writing)
window_days  <- 4                      # keeps each response under Savant's 25k cap
sleep_secs   <- 8                      # Savant throttles bursts; keep spaced
max_retries  <- 6
empty_backoff_base <- 40               # seconds; multiplied by attempt on throttle
cap_warn     <- 24000                  # warn if a window approaches the 25k cap

# Windows ending on/after this date are expected to contain regular-season
# pitches, so a 0-row response there almost certainly means Savant throttled
# us (it returns HTTP 200 with an empty body) and we should retry with backoff.
# Windows before this are spring training / off-season where empty is legit.
expect_data_from <- as.Date("2026-03-26")

out_dir <- file.path("data", "statcast_2026")
out_csv <- file.path(out_dir, "statcast_2026_all.csv")
progress_file <- file.path("data", "statcast_2026_progress.log")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Append a line to the progress file (opens/closes each call, so it flushes
# immediately and can be tailed while the scrape runs).
log_progress <- function(...) {
  cat(sprintf(...), "\n", file = progress_file, append = TRUE)
}
cat("", file = progress_file)  # truncate

# ---- Helper: build the Savant CSV URL for a date window -------------------
# hfGT=R|PO|S includes regular season, postseason, and spring training so we
# capture every Statcast-tracked pitch. The analysis script can filter by
# game_type later.

savant_url <- function(start_date, end_date) {
  paste0(
    "https://baseballsavant.mlb.com/statcast_search/csv?",
    "all=true",
    "&hfGT=R%7CPO%7CS%7C",
    "&hfSea=2026%7C",
    "&player_type=pitcher",
    "&game_date_gt=", start_date,
    "&game_date_lt=", end_date,
    "&min_pitches=0&min_results=0",
    "&group_by=name&sort_col=pitches",
    "&player_event_sort=api_p_release_speed&sort_order=desc",
    "&min_pas=0&type=details&"
  )
}

# ---- Helper: fetch one window with retries --------------------------------

is_empty_response <- function(res) {
  nrow(res) == 0 ||
    (nrow(res) == 1 && !("pitch_type" %in% names(res))) ||
    (nrow(res) == 1 && "pitch_type" %in% names(res) && all(is.na(res[["pitch_type"]])))
}

# expect_data: if TRUE, treat a 0-row response as a likely throttle and retry.
fetch_window <- function(start_date, end_date, expect_data) {
  url <- savant_url(start_date, end_date)
  for (attempt in seq_len(max_retries)) {
    res <- tryCatch({
      h <- curl::new_handle()
      curl::handle_setheaders(h,
        "User-Agent" = USER_AGENT,
        "Accept"     = "text/csv,*/*")
      tmp <- tempfile(fileext = ".csv")
      req <- curl::curl_fetch_disk(url, tmp, handle = h)
      if (req$status_code != 200) {
        stop(sprintf("HTTP %d", req$status_code))
      }
      dt <- data.table::fread(tmp, showProgress = FALSE)
      unlink(tmp)
      dt
    }, error = function(e) {
      message(sprintf("  attempt %d failed: %s", attempt, conditionMessage(e)))
      NULL
    })

    if (!is.null(res) && is.data.frame(res)) {
      if (is_empty_response(res)) {
        if (expect_data && attempt < max_retries) {
          backoff <- empty_backoff_base * attempt
          message(sprintf("  empty response (likely throttled), retrying in %ds", backoff))
          Sys.sleep(backoff)
          next
        }
        return(res[0])  # genuinely empty (off-season) or out of retries
      }
      # Return all rows as-is. (Do NOT filter on game_date: fread parses it as a
      # Date, so comparing it to "" yields all-NA and would drop every row.)
      return(res)
    }
    Sys.sleep(sleep_secs * attempt)
  }
  message(sprintf("  giving up on %s -> %s", start_date, end_date))
  NULL
}

# ---- Loop over the season in windows --------------------------------------
# Each window is saved to its own chunk file so a killed/interrupted run can be
# resumed simply by re-running (completed chunks are skipped). Empty in-season
# windows are NOT cached (so they get retried on the next run); genuinely empty
# off-season windows are marked done with a sentinel file.

chunks_dir <- file.path(out_dir, "chunks")
dir.create(chunks_dir, recursive = TRUE, showWarnings = FALSE)

windows <- seq(season_start, season_end, by = window_days)
failed  <- character(0)

for (i in seq_along(windows)) {
  w_start <- windows[i]
  w_end   <- min(w_start + (window_days - 1), season_end)
  expect  <- w_end >= expect_data_from
  chunk_path  <- file.path(chunks_dir, sprintf("chunk_%s.csv", w_start))
  empty_path  <- file.path(chunks_dir, sprintf("chunk_%s.empty", w_start))

  # Resume: skip windows already fetched (non-empty chunk or off-season empty).
  if (file.exists(chunk_path)) {
    message(sprintf("[%d/%d] %s -> %s : cached, skipping", i, length(windows), w_start, w_end))
    log_progress("[%d/%d] %s -> %s : cached", i, length(windows), w_start, w_end)
    next
  }
  if (!expect && file.exists(empty_path)) {
    message(sprintf("[%d/%d] %s -> %s : cached empty, skipping", i, length(windows), w_start, w_end))
    next
  }

  message(sprintf("[%d/%d] %s -> %s", i, length(windows), w_start, w_end))

  df <- fetch_window(as.character(w_start), as.character(w_end), expect_data = expect)
  if (is.null(df)) {
    failed <- c(failed, sprintf("%s -> %s", w_start, w_end))
    message("  FAILED (all retries)")
    log_progress("[%d/%d] %s -> %s : FAILED", i, length(windows), w_start, w_end)
  } else if (nrow(df) > 0) {
    data.table::fwrite(df, chunk_path)
    message(sprintf("  got %d rows -> %s", nrow(df), basename(chunk_path)))
    log_progress("[%d/%d] %s -> %s : %d rows", i, length(windows), w_start, w_end, nrow(df))
    if (nrow(df) >= cap_warn) {
      message(sprintf("  WARNING: %d rows is near Savant's 25k cap; consider a smaller window_days", nrow(df)))
      log_progress("  WARNING near cap for %s -> %s (%d rows)", w_start, w_end, nrow(df))
    }
  } else if (expect) {
    # In-season window that stayed empty through all retries -> suspicious.
    # Do NOT cache, so a re-run will try again.
    failed <- c(failed, sprintf("%s -> %s (stayed empty)", w_start, w_end))
    message("  0 rows (in-season, flagged)")
    log_progress("[%d/%d] %s -> %s : 0 rows (flagged)", i, length(windows), w_start, w_end)
  } else {
    file.create(empty_path)  # off-season empty: cache the sentinel
    message("  0 rows")
    log_progress("[%d/%d] %s -> %s : 0 rows", i, length(windows), w_start, w_end)
  }
  Sys.sleep(sleep_secs)
}

if (length(failed) > 0) {
  message("\nWARNING: the following windows failed and are missing (re-run to retry):")
  for (f in failed) message("  ", f)
}

# ---- Combine all chunk files ----------------------------------------------

chunk_files <- list.files(chunks_dir, pattern = "^chunk_.*\\.csv$", full.names = TRUE)
if (length(chunk_files) == 0) {
  stop("No data chunks found. Check network access to baseballsavant.mlb.com.")
}
message(sprintf("Combining %d chunk files...", length(chunk_files)))
all_pitches <- data.table::rbindlist(
  lapply(chunk_files, function(f) data.table::fread(f, showProgress = FALSE)),
  use.names = TRUE, fill = TRUE
)

# De-dup in case overlapping windows returned the same pitch. game_pk +
# at_bat_number + pitch_number uniquely identifies a pitch.
dedup_keys <- intersect(c("game_pk", "at_bat_number", "pitch_number"),
                        names(all_pitches))
if (length(dedup_keys) == 3) {
  before <- nrow(all_pitches)
  all_pitches <- unique(all_pitches, by = dedup_keys)
  message(sprintf("De-duplicated %d -> %d rows", before, nrow(all_pitches)))
}

# ---- Save -----------------------------------------------------------------

data.table::fwrite(all_pitches, out_csv)
message(sprintf("Wrote %d pitches (%d columns) to %s",
                nrow(all_pitches), ncol(all_pitches), out_csv))
