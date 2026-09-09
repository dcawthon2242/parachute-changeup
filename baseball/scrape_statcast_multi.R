#!/usr/bin/env Rscript

# Scrape full-season MLB Statcast pitch-level data for one or more seasons from
# Baseball Savant. This is a parameterized generalization of scrape_statcast_2026.R
# used to pull the 2023-2025 training seasons for the miss-distance pitch grade model.
#
# Usage:
#   Rscript baseball/scrape_statcast_multi.R 2023 2024 2025
# (defaults to 2023 2024 2025 if no args given)
#
# Each season lands in data/statcast_<YEAR>/ with per-window chunk files, so an
# interrupted run resumes by re-running (completed chunks are skipped). The
# type=details&all=true response already carries miss_distance + bat tracking +
# full kinematics, so no query change is needed beyond the season.

suppressPackageStartupMessages({
  library(data.table)
  library(curl)
})

USER_AGENT <- paste0(
  "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) ",
  "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/126.0 Safari/537.36"
)

# ---- Config ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
seasons <- if (length(args) > 0) as.integer(args) else c(2023L, 2024L, 2025L)

window_days  <- 4
sleep_secs   <- 5
max_retries  <- 4
empty_backoff_base <- 20
cap_warn     <- 24000

# Regular season roughly late March -> end of September. We scan a slightly wider
# window; off-season windows return empty and are cached as sentinels. A 0-row
# response is only treated as a throttle (worth retrying) INSIDE [expect_from,
# expect_to]; trailing post-season windows return empty legitimately and are
# cached without burning time on backoff retries.
season_bounds <- function(yr) {
  # 2020 was a 60-game season starting July 23. Without this the April-July windows are
  # treated as throttled rather than legitimately empty and each burns four retries with
  # escalating backoff, which costs about an hour of wall clock for nothing.
  if (yr == 2020L)
    return(list(start = as.Date("2020-07-20"), end = as.Date("2020-09-30"),
                expect_from = as.Date("2020-07-24"), expect_to = as.Date("2020-09-27")))
  list(start = as.Date(sprintf("%d-03-15", yr)),
       end   = as.Date(sprintf("%d-10-02", yr)),
       expect_from = as.Date(sprintf("%d-04-01", yr)),
       expect_to   = as.Date(sprintf("%d-09-29", yr)))
}

savant_url <- function(start_date, end_date, yr) {
  paste0(
    "https://baseballsavant.mlb.com/statcast_search/csv?",
    "all=true",
    "&hfGT=R%7CPO%7CS%7C",
    "&hfSea=", yr, "%7C",
    "&player_type=pitcher",
    "&game_date_gt=", start_date,
    "&game_date_lt=", end_date,
    "&min_pitches=0&min_results=0",
    "&group_by=name&sort_col=pitches",
    "&player_event_sort=api_p_release_speed&sort_order=desc",
    "&min_pas=0&type=details&"
  )
}

is_empty_response <- function(res) {
  nrow(res) == 0 ||
    (nrow(res) == 1 && !("pitch_type" %in% names(res))) ||
    (nrow(res) == 1 && "pitch_type" %in% names(res) && all(is.na(res[["pitch_type"]])))
}

fetch_window <- function(url, expect_data) {
  for (attempt in seq_len(max_retries)) {
    res <- tryCatch({
      h <- curl::new_handle()
      curl::handle_setheaders(h,
        "User-Agent" = USER_AGENT,
        "Accept"     = "text/csv,*/*")
      tmp <- tempfile(fileext = ".csv")
      req <- curl::curl_fetch_disk(url, tmp, handle = h)
      if (req$status_code != 200) stop(sprintf("HTTP %d", req$status_code))
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
        return(res[0])
      }
      return(res)
    }
    Sys.sleep(sleep_secs * attempt)
  }
  NULL
}

# ---- Scrape one season ----------------------------------------------------
scrape_season <- function(yr) {
  b <- season_bounds(yr)
  out_dir    <- file.path("data", sprintf("statcast_%d", yr))
  out_csv    <- file.path(out_dir, sprintf("statcast_%d_all.csv", yr))
  chunks_dir <- file.path(out_dir, "chunks")
  progress_file <- file.path("data", sprintf("statcast_%d_progress.log", yr))
  dir.create(chunks_dir, recursive = TRUE, showWarnings = FALSE)
  cat("", file = progress_file)
  log_progress <- function(...) cat(sprintf(...), "\n", file = progress_file, append = TRUE)

  windows <- seq(b$start, b$end, by = window_days)
  failed  <- character(0)
  message(sprintf("\n===== SEASON %d : %d windows =====", yr, length(windows)))

  for (i in seq_along(windows)) {
    w_start <- windows[i]
    w_end   <- min(w_start + (window_days - 1), b$end)
    expect  <- w_end >= b$expect_from && w_start <= b$expect_to
    chunk_path <- file.path(chunks_dir, sprintf("chunk_%s.csv", w_start))
    empty_path <- file.path(chunks_dir, sprintf("chunk_%s.empty", w_start))

    if (file.exists(chunk_path)) {
      log_progress("[%d/%d] %s -> %s : cached", i, length(windows), w_start, w_end)
      next
    }
    if (!expect && file.exists(empty_path)) next

    message(sprintf("[%d/%d] %s -> %s", i, length(windows), w_start, w_end))
    df <- fetch_window(savant_url(as.character(w_start), as.character(w_end), yr), expect)

    if (is.null(df)) {
      failed <- c(failed, sprintf("%s -> %s", w_start, w_end))
      log_progress("[%d/%d] %s -> %s : FAILED", i, length(windows), w_start, w_end)
    } else if (nrow(df) > 0) {
      data.table::fwrite(df, chunk_path)
      message(sprintf("  got %d rows", nrow(df)))
      log_progress("[%d/%d] %s -> %s : %d rows", i, length(windows), w_start, w_end, nrow(df))
      if (nrow(df) >= cap_warn)
        log_progress("  WARNING near cap for %s -> %s (%d rows)", w_start, w_end, nrow(df))
    } else if (expect) {
      failed <- c(failed, sprintf("%s -> %s (stayed empty)", w_start, w_end))
      log_progress("[%d/%d] %s -> %s : 0 rows (flagged)", i, length(windows), w_start, w_end)
    } else {
      file.create(empty_path)
    }
    Sys.sleep(sleep_secs)
  }

  chunk_files <- list.files(chunks_dir, pattern = "^chunk_.*\\.csv$", full.names = TRUE)
  if (length(chunk_files) == 0) {
    message(sprintf("SEASON %d: no chunks fetched (network?).", yr))
    return(invisible(NULL))
  }
  message(sprintf("SEASON %d: combining %d chunks...", yr, length(chunk_files)))
  all_pitches <- data.table::rbindlist(
    lapply(chunk_files, function(f) data.table::fread(f, showProgress = FALSE)),
    use.names = TRUE, fill = TRUE
  )
  dedup_keys <- intersect(c("game_pk", "at_bat_number", "pitch_number"), names(all_pitches))
  if (length(dedup_keys) == 3) all_pitches <- unique(all_pitches, by = dedup_keys)
  data.table::fwrite(all_pitches, out_csv)
  message(sprintf("SEASON %d: wrote %d pitches (%d cols) -> %s",
                  yr, nrow(all_pitches), ncol(all_pitches), out_csv))
  if (length(failed) > 0)
    message(sprintf("SEASON %d: %d windows failed (re-run to retry).", yr, length(failed)))
}

for (yr in seasons) scrape_season(yr)
message("\nAll requested seasons processed.")
