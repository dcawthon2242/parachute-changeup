#!/usr/bin/env Rscript

# Follow-up to sweeper_platoon_2026.R: is the sweeper actually EFFECTIVE when
# thrown to opposite-handed hitters, and specifically for the pitchers who
# lean on it against the opposite hand?
#
# Run value is taken from delta_pitcher_run_exp, which is signed so that
# positive = good for the pitcher (verified: strikeouts +0.22, HR -1.42).
# RV/100 is therefore runs saved per 100 pitches.

suppressPackageStartupMessages(library(data.table))

chunks_dir <- file.path("data", "statcast_2026", "chunks")

keep_cols <- c("pitch_type", "game_type", "player_name", "pitcher",
               "stand", "p_throws", "description", "zone", "woba_value",
               "woba_denom", "estimated_woba_using_speedangle",
               "delta_pitcher_run_exp", "launch_speed")

chunk_files <- list.files(chunks_dir, pattern = "^chunk_.*\\.csv$", full.names = TRUE)
p <- rbindlist(lapply(chunk_files, function(f)
  fread(f, select = keep_cols, showProgress = FALSE)), use.names = TRUE, fill = TRUE)

p <- p[game_type == "R" & stand %in% c("L", "R") & p_throws %in% c("L", "R")]
p[, side := ifelse(stand == p_throws, "same", "opp")]

swing_desc <- c("swinging_strike", "swinging_strike_blocked", "foul", "foul_tip",
                "hit_into_play", "hit_into_play_score", "hit_into_play_no_out",
                "foul_bunt", "missed_bunt", "bunt_foul_tip")
whiff_desc <- c("swinging_strike", "swinging_strike_blocked", "foul_tip",
                "missed_bunt", "bunt_foul_tip")

p[, `:=`(
  swing    = description %in% swing_desc,
  whiff    = description %in% whiff_desc,
  in_zone  = zone >= 1 & zone <= 9,
  xwoba_num = ifelse(is.na(estimated_woba_using_speedangle),
                     woba_value, estimated_woba_using_speedangle)
)]

# Metrics for an arbitrary grouping.
score <- function(dt, by) {
  dt[, .(
    pitches   = .N,
    usage     = NA_real_,
    whiff_pct = 100 * sum(whiff) / sum(swing),
    chase_pct = 100 * sum(swing & !in_zone) / sum(!in_zone),
    zone_pct  = 100 * sum(in_zone) / .N,
    xwoba     = sum(xwoba_num[woba_denom > 0], na.rm = TRUE) /
                sum(woba_denom[woba_denom > 0], na.rm = TRUE),
    rv100     = 100 * sum(delta_pitcher_run_exp, na.rm = TRUE) / .N
  ), by = by]
}

st <- p[pitch_type == "ST"]

cat("=== 1. League-wide 2026: sweeper by batter side ===\n")
lg <- score(st, "side")[order(side)]
print(lg[, .(side, pitches, whiff_pct = round(whiff_pct, 1),
             chase_pct = round(chase_pct, 1), xwoba = round(xwoba, 3),
             rv100 = round(rv100, 2))])

cat("\n=== 2. All pitch types, for reference (RV/100 by side) ===\n")
ref <- score(p[pitch_type %in% c("FF", "SL", "ST", "CU", "CH", "SI", "FC", "FS")],
             c("pitch_type", "side"))
ref <- dcast(ref, pitch_type ~ side, value.var = c("rv100", "xwoba", "whiff_pct"))
print(ref[, .(pitch_type,
              rv_opp = round(rv100_opp, 2), rv_same = round(rv100_same, 2),
              xw_opp = round(xwoba_opp, 3), xw_same = round(xwoba_same, 3),
              wh_opp = round(whiff_pct_opp, 1), wh_same = round(whiff_pct_same, 1))][order(-rv_opp)])

# --- The reverse-platoon group from the usage analysis ----------------------
plat <- fread(file.path("data", "statcast_2026", "sweeper_platoon_2026.csv"))
rev_ids <- plat[diff > 0, pitcher]
sig_ids <- plat[player_name %in% c("Bazardo, Eduard", "Senga, Kodai",
                                   "Wilcox, Cole", "Soto, Gregory"), pitcher]

cat("\n=== 3. Reverse-platoon sweeper users: their ST vs opposite hand ===\n")
grp <- score(st[pitcher %in% rev_ids], c("pitcher", "player_name", "side"))
w <- dcast(grp, pitcher + player_name ~ side,
           value.var = c("pitches", "whiff_pct", "xwoba", "rv100"))
w <- merge(w, plat[, .(pitcher, diff, sw_total)], by = "pitcher")
out <- w[order(-diff), .(
  Pitcher = player_name,
  `ST vs opp` = pitches_opp,
  `whiff%` = round(whiff_pct_opp, 1),
  xwOBA = round(xwoba_opp, 3),
  `RV/100` = round(rv100_opp, 2),
  `RV/100 same` = round(rv100_same, 2),
  `xwOBA same` = round(xwoba_same, 3),
  sig = ifelse(pitcher %in% sig_ids, "*", "")
)]
print(out)

cat("\n=== 4. Group aggregate: reverse-platoon users' ST vs opp, ")
cat("compared to everyone else's ST vs opp ===\n")
st[, grp := ifelse(pitcher %in% rev_ids, "reverse-platoon users", "all other pitchers")]
agg <- score(st[side == "opp"], "grp")
print(agg[, .(grp, pitches, whiff_pct = round(whiff_pct, 1),
              chase_pct = round(chase_pct, 1), xwoba = round(xwoba, 3),
              rv100 = round(rv100, 2))])

cat("\n=== 5. The four significant cases, both sides side-by-side ===\n")
sig <- score(st[pitcher %in% sig_ids], c("player_name", "side"))
print(sig[order(player_name, side), .(player_name, side, pitches,
          whiff_pct = round(whiff_pct, 1), chase_pct = round(chase_pct, 1),
          xwoba = round(xwoba, 3), rv100 = round(rv100, 2))])
