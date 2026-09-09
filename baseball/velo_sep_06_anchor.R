#!/usr/bin/env Rscript

# A PER-OUTING FASTBALL ANCHOR.
#
# Every separation number in this project so far is measured against a season-constant four-seam
# velocity, which is why sd(speed_diff) turned out to be sd(release_speed) exactly. That conflates
# two different things:
#
#   a pitcher who cannot repeat the changeup within an outing        - a command problem
#   a pitcher whose whole arsenal drifts down from April to September - a workload story
#
# The first should hurt; the second is not really inconsistency at all, because if the fastball
# drops with the changeup the hitter still sees the same gap. Anchoring to the four-seam velocity
# the pitcher actually had THAT DAY separates them.
#
# Built from the raw season CSVs because parachute_ff holds changeups only and has no fastball rows
# and no game key. Only nine columns are read per season.
#
# Two anchors, because the choice is not obvious:
#   game    mean four-seam velocity in that pitcher's outing. Stable, needs enough fastballs.
#   recent  mean of the pitcher's nearest four-seamers earlier in the same outing, so it is causal
#           and reacts within a start. Falls back to the game mean before enough have been thrown.
#
# MIN_FF guards the anchor: a game mean off two fastballs is noise dressed as precision.

suppressPackageStartupMessages({ library(data.table) })
options(width = 205)
MDIR <- "data/statcast_model"; OUT <- file.path(MDIR, "velo_anchor_game.rds")
MIN_FF <- 5L; NREC <- 10L

if (!file.exists(OUT)) {
  COLS <- c("pitch_type","game_date","release_speed","pitcher","p_throws",
            "inning","game_pk","at_bat_number","pitch_number")
  A <- rbindlist(lapply(2020:2026, function(y) {
    f <- sprintf("data/statcast_%d/statcast_%d_all.csv", y, y)
    x <- fread(f, select = COLS, showProgress = FALSE)
    x <- x[pitch_type %chin% c("FF","CH") & is.finite(release_speed)]
    x[, season := y]
    cat(sprintf("  %d: %s four-seams and changeups\n", y, format(nrow(x), big.mark = ",")))
    x }), use.names = TRUE)

  setorder(A, pitcher, season, game_pk, at_bat_number, pitch_number)
  A[, seq := seq_len(.N), by = .(pitcher, game_pk)]

  # Game anchor: mean four-seam velocity in the outing, only where the pitcher threw enough.
  G <- A[pitch_type == "FF", .(fb_game = mean(release_speed), n_ff = .N),
         by = .(pitcher, season, game_pk)][n_ff >= MIN_FF]

  # Recent anchor: expanding mean of four-seamers thrown EARLIER in the same outing. Causal by
  # construction - a pitch is never anchored to fastballs the pitcher had not thrown yet.
  FF <- A[pitch_type == "FF", .(pitcher, game_pk, seq, release_speed)]
  setorder(FF, pitcher, game_pk, seq)
  FF[, `:=`(cs = cumsum(release_speed), ci = seq_len(.N)), by = .(pitcher, game_pk)]
  FF[, fb_recent := fifelse(ci >= NREC,
        (cs - shift(cs, NREC, fill = 0)) / NREC, cs / ci), by = .(pitcher, game_pk)]
  FF[, k := ci]
  saveRDS(list(game = G, ff = FF[, .(pitcher, game_pk, seq, fb_recent, k)],
               ch = A[pitch_type == "CH", .(pitcher, season, game_pk, at_bat_number, pitch_number,
                                            seq, ch_velo = release_speed)]), OUT)
  cat(sprintf("wrote %s\n", OUT))
} else cat("(using cached anchor build)\n")

L <- readRDS(OUT); G <- L$game; FF <- L$ff; CH <- L$ch
setDT(G); setDT(FF); setDT(CH)

# For each changeup, the most recent earlier four-seam anchor in the same outing.
setkey(FF, pitcher, game_pk, seq)
CH[, rid := .I]
M <- FF[CH, on = .(pitcher, game_pk, seq), roll = TRUE,
        .(rid = i.rid, fb_recent, k)]
CH <- merge(CH, M, by = "rid", all.x = TRUE)
CH <- merge(CH, G[, .(pitcher, game_pk, fb_game, n_ff)], by = c("pitcher","game_pk"), all.x = TRUE)
CH[, `:=`(sep_game = fb_game - ch_velo, sep_recent = fb_recent - ch_velo)]

cat(sprintf("\nchangeups: %s | with a game anchor: %s (%.1f%%) | with an earlier-fastball anchor: %s (%.1f%%)\n",
            format(nrow(CH), big.mark = ","),
            format(sum(is.finite(CH$sep_game)), big.mark = ","),
            100*mean(is.finite(CH$sep_game)),
            format(sum(is.finite(CH$sep_recent)), big.mark = ","),
            100*mean(is.finite(CH$sep_recent))))

## ---- how different is the outing anchor from the season constant? --------------------------------
CH[, fb_season := mean(fb_game, na.rm = TRUE), by = .(pitcher, season)]
CH[, sep_season := fb_season - ch_velo]
V <- CH[is.finite(sep_game) & is.finite(sep_season),
        .(np = .N, games = uniqueN(game_pk),
          sd_season = sd(sep_season), sd_game = sd(sep_game),
          sd_fb_across_games = sd(fb_game), mean_sep = mean(sep_game)),
        by = .(pitcher, season)][np >= 40 & games >= 5]
cat(sprintf("\n%d pitcher-seasons with 40+ changeups and 5+ anchored outings\n", nrow(V)))
cat("\n=== how much of 'inconsistent separation' was really day-to-day fastball drift? ===\n")
cat(sprintf("  sd of separation, season anchor : %.3f mph (median)\n", median(V$sd_season)))
cat(sprintf("  sd of separation, outing anchor : %.3f mph (median)\n", median(V$sd_game)))
cat(sprintf("  the pitcher's own four-seam moves %.3f mph sd across his outings (median)\n",
            median(V$sd_fb_across_games)))
cat(sprintf("  correlation between the two sd measures: %+.3f\n", cor(V$sd_season, V$sd_game)))
cat(sprintf("  outing anchoring changes the sd by %+.1f%% on average, and reorders arms:\n",
            100*mean(V$sd_game/V$sd_season - 1)))
cat(sprintf("  spearman rank correlation of the two = %+.3f\n",
            cor(V$sd_season, V$sd_game, method = "spearman")))

saveRDS(CH[, .(pitcher, season, game_pk, at_bat_number, pitch_number,
               ch_velo, fb_game, fb_recent, sep_game, sep_recent, sep_season, n_ff)],
        file.path(MDIR, "velo_anchor_pitches.rds"))
saveRDS(V, file.path(MDIR, "velo_anchor_seasons.rds"))
cat("\nwrote velo_anchor_pitches.rds and velo_anchor_seasons.rds\n")
