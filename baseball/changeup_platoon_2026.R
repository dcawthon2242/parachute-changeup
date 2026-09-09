#!/usr/bin/env Rscript

# Mirror of sweeper_platoon_2026.R, for the changeup.
#
# The changeup is the one pitch type that is genuinely BETTER against
# opposite-handed hitters, so nearly everyone throws it that way. Question here
# is the inverse: are there pitchers who use it MORE against same-handed hitters?
#
# "More" is usage rate -- changeups as a share of that pitcher's total pitches
# against each batter side -- because raw counts just track how many lefties vs
# righties a pitcher happens to face.

suppressPackageStartupMessages(library(data.table))
options(width = 210)

chunks_dir <- file.path("data", "statcast_2026", "chunks")
out_csv    <- file.path("data", "statcast_2026", "changeup_platoon_2026.csv")

cols <- c("pitch_type", "game_type", "player_name", "pitcher", "stand", "p_throws",
          "description", "zone", "woba_denom", "woba_value",
          "estimated_woba_using_speedangle", "delta_pitcher_run_exp", "release_speed")

cf <- list.files(chunks_dir, pattern = "^chunk_.*\\.csv$", full.names = TRUE)
p <- rbindlist(lapply(cf, function(f) fread(f, select = cols, showProgress = FALSE)),
               use.names = TRUE, fill = TRUE)
p <- p[game_type == "R" & !is.na(pitch_type) & pitch_type != "" &
         stand %in% c("L", "R") & p_throws %in% c("L", "R")]
p[, side := ifelse(stand == p_throws, "same", "opp")]

wd <- c("swinging_strike", "swinging_strike_blocked", "foul_tip", "missed_bunt", "bunt_foul_tip")
sd <- c(wd, "foul", "hit_into_play", "hit_into_play_score", "hit_into_play_no_out", "foul_bunt")
p[, `:=`(swing = description %in% sd, whiff = description %in% wd,
         inz = zone >= 1 & zone <= 9,
         xw = ifelse(is.na(estimated_woba_using_speedangle),
                     woba_value, estimated_woba_using_speedangle))]

TARGET <- "CH"

# ---- usage ---------------------------------------------------------------
agg <- p[, .(pitches = .N, ch = sum(pitch_type == TARGET)), by = .(pitcher, player_name, p_throws, side)]
w <- dcast(agg, pitcher + player_name + p_throws ~ side, value.var = c("pitches", "ch"), fill = 0)
setnames(w, c("pitches_same","pitches_opp","ch_same","ch_opp"),
         c("pit_same","pit_opp","ch_same","ch_opp"))
w[, `:=`(ch_total = ch_same + ch_opp, pit_total = pit_same + pit_opp,
         rate_same = 100*ch_same/pit_same, rate_opp = 100*ch_opp/pit_opp)]
w[, `:=`(usage_total = 100*ch_total/pit_total, diff_same = rate_same - rate_opp)]

q <- w[ch_total >= 50 & pit_same >= 100 & pit_opp >= 100]
cat(sprintf("Regular-season pitches: %s\n", format(nrow(p), big.mark = ",")))
cat(sprintf("Qualified changeup-throwers (>=50 CH, >=100 pitches vs each side): %d\n", nrow(q)))

rev <- q[diff_same > 0][order(-diff_same)]
rev[, pval := mapply(function(a, b, c, e) fisher.test(matrix(c(a, b-a, c, e-c), 2))$p.value,
                     ch_same, pit_same, ch_opp, pit_opp)]
cat(sprintf("%d of %d (%.0f%%) use the changeup MORE against same-handed hitters.\n\n",
            nrow(rev), nrow(q), 100*nrow(rev)/nrow(q)))

# ---- quality of those changeups, by side ---------------------------------
sc <- function(dt, by) dt[, .(
  n = .N, swings = sum(swing),
  whiff = round(100*sum(whiff)/sum(swing), 1),
  chase = round(100*sum(swing & !inz)/sum(!inz), 1),
  xwoba = round(sum(xw[woba_denom > 0 & !is.na(woba_denom)], na.rm = TRUE) /
                sum(woba_denom[woba_denom > 0 & !is.na(woba_denom)], na.rm = TRUE), 3),
  rv100 = round(100*sum(delta_pitcher_run_exp, na.rm = TRUE)/.N, 2)), by = by]

ch <- p[pitch_type == TARGET]

cat("=== Changeup usage HIGHER vs SAME-handed hitters (2026, qualified) ===\n")
print(rev[, .(Pitcher = player_name, T = p_throws, CH = ch_total,
              `Use%` = sprintf("%.1f", usage_total),
              `vsSame%` = sprintf("%.1f", rate_same),
              `vsOpp%` = sprintf("%.1f", rate_opp),
              Diff = sprintf("%+.1f", diff_same),
              `CH same/opp` = sprintf("%d/%d", ch_same, ch_opp),
              p = sprintf("%.4f", pval),
              sig = ifelse(pval < 0.05, "*", ""))], nrows = 100)

cat("\n=== League baseline: changeup by batter side ===\n")
print(sc(ch, "side")[order(side)])

cat("\n=== Those pitchers' changeups: same-hand vs opposite-hand results ===\n")
g <- sc(ch[pitcher %in% rev$pitcher], c("pitcher", "player_name", "side"))
gw <- dcast(g, pitcher + player_name ~ side,
            value.var = c("n", "swings", "whiff", "xwoba", "rv100"))
gw <- merge(gw, rev[, .(pitcher, diff_same, pval)], by = "pitcher")
print(gw[order(-diff_same), .(Pitcher = player_name,
        `CH same` = n_same, `whiff same` = whiff_same, `xwOBA same` = xwoba_same,
        `RV/100 same` = rv100_same, `CH opp` = n_opp, `whiff opp` = whiff_opp,
        `xwOBA opp` = xwoba_opp, `RV/100 opp` = rv100_opp,
        sig = ifelse(pval < 0.05, "*", ""))], nrows = 100)

cat("\n=== Group aggregate: same-hand-leaning group's CH vs same hand, ")
cat("compared to everyone else's CH vs same hand ===\n")
ch[, grp := ifelse(pitcher %in% rev$pitcher, "same-hand leaning", "all other pitchers")]
print(sc(ch[side == "same"], "grp"))

fwrite(q[order(-diff_same)], out_csv)
cat(sprintf("\nWrote %s\n", out_csv))
