#!/usr/bin/env Rscript

# Ten fastball-to-breaking-ball pairs per breaking ball type, ranked by the angular
# tunnel metric. Knuckle-curves are folded into curveballs as requested.
#
# Columns: fb is whichever fastball the pitcher most often sets the pitch up with. brk is
# the mean break fraction at tau=0.05 (share of the decision window the pair survives
# before the hitter can separate them; higher is a longer tunnel). pct is the percentile
# within that pitch type. overlap is the share of the window where the two discs actually
# intersect on the retina. chase and whiff are that pitcher's own rates on these pairs,
# and old_rk is where the previous weighted-distance metric ranked him, included because
# the two metrics disagree often enough that it is worth seeing.
#
# Rates on 25-100 pitches are noisy and are shown as description, not evidence.

suppressPackageStartupMessages(library(data.table))
options(width = 200)
MIN_PAIRS <- 25

d <- readRDS("data/statcast_model/angular_tunnel_2026.rds")
cf <- list.files("data/statcast_2026/chunks", pattern = "csv$", full.names = TRUE)
ex <- unique(rbindlist(lapply(cf, function(f) fread(f, select = c(
  "game_pk","at_bat_number","pitch_number","sz_top","sz_bot"), showProgress = FALSE))),
  by = c("game_pk","at_bat_number","pitch_number"))
d <- merge(d, ex, by = c("game_pk","at_bat_number","pitch_number"), all.x = TRUE)
d[, zdist := sqrt(pmax(abs(plate_x) - 0.95, 0)^2 + pmax(plate_z - sz_top, sz_bot - plate_z, 0)^2)]
d[, out_zone := zdist > 0 & is.finite(zdist)]
d[, grp := fifelse(pitch_type == "SL", "Slider",
            fifelse(pitch_type == "ST", "Sweeper",
            fifelse(pitch_type %in% c("CU","KC"), "Curveball", NA_character_)))]
d <- d[!is.na(grp)]

for (g in c("Slider","Sweeper","Curveball")) {
  s <- d[grp == g]
  a <- s[, .(pairs = .N,
             fb = names(sort(table(p_pitch_type), decreasing = TRUE))[1],
             brk = mean(brk_any_005), overlap = mean(overlap_frac),
             chase = 100*mean(swing[out_zone]), n_oz = sum(out_zone),
             whiff = 100*mean(whiff[swing]), n_sw = sum(swing),
             old = mean(tunnel_old)),
         by = .(pitcher, player_name, p_throws)][pairs >= MIN_PAIRS]
  setorder(a, old); a[, old_rk := .I]
  setorder(a, -brk); a[, rk := .I]
  a[, pct := round(100*(1 - (rk - 1)/.N))]
  cat(sprintf("\n=== %s: top 10 fastball-to-%s tunnels, 2026 (min %d pairs, %d qualified) ===\n",
              toupper(g), tolower(g), MIN_PAIRS, nrow(a)))
  cat(sprintf("    league mean break fraction %.3f | mean chase %.1f%% | mean whiff %.1f%%\n",
              mean(a$brk), 100*mean(s$swing[s$out_zone]), 100*mean(s$whiff[s$swing])))
  print(a[1:10, .(rk, pitcher = player_name, T = p_throws, fb, pairs,
                  brk = round(brk,3), pct, overlap = round(overlap,3),
                  chase = round(chase,1), n_oz, whiff = round(whiff,1), n_sw,
                  old_rk)], row.names = FALSE)
  cat(sprintf("    rank agreement with the old metric: Spearman %.3f\n",
              cor(a$rk, a$old_rk, method = "spearman")))
}
