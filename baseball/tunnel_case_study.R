#!/usr/bin/env Rscript

# Fastball -> sweeper tunneling case study for one pitcher, against
# OPPOSITE-handed hitters. Generalizes baseball/wilcox_tunnel_2026.R:
#
#   Rscript baseball/tunnel_case_study.R "Soto, Gregory"
#   Rscript baseball/tunnel_case_study.R "Wilcox, Cole"     # reproduces the original
#
# Three questions, in order:
#   1. path-to-location ratio on FB -> ST sequences vs the opposite hand, and
#      where that sits among same-handed sweeper-throwers league-wide
#   2. miss distance on sweepers that FOLLOW a fastball vs sweepers that do not
#   3. chase rate on sweepers that follow a fastball vs sweepers that do not
#
# Tunnel definition is copied exactly from driver_gain_all_types.R so the number
# is comparable to the rest of the project: weighted 3D separation between the
# two trajectories, integrated from release to the commit point (plate time minus
# a 150ms reaction window), divided by how far apart the two pitches finish.
# LOWER path ratio = the pitches stayed together longer relative to how far they
# ended up apart = tighter tunnel. Because the denominator is "how far apart they
# finish", a pitcher whose two pitches end up close together is penalised on the
# ratio even with a tight early path, so the numerator is reported separately --
# the same caveat tunnel_location_sensitivity.R raises.

suppressPackageStartupMessages(library(data.table))
options(width = 205)

args <- commandArgs(trailingOnly = TRUE)
TARGET_PITCHER <- if (length(args) >= 1) args[1] else "Soto, Gregory"
MIN_PAIRS <- 20   # league leaderboard qualifier

FASTBALLS <- c("FF", "SI", "FC"); yf <- 17/12
WX <- 1.2; WY <- 0.4; WZ <- 1.4; REACT <- 0.150; NSTEP <- 40

cols <- c("game_pk","at_bat_number","pitch_number","pitch_type","game_type","pitcher",
          "player_name","p_throws","stand","balls","strikes","plate_x","plate_z",
          "release_pos_x","release_pos_y","release_pos_z","vx0","vy0","vz0","ax","ay","az",
          "description","zone","miss_distance","release_speed","delta_pitcher_run_exp",
          "woba_denom","woba_value","estimated_woba_using_speedangle")

cf <- list.files(file.path("data","statcast_2026","chunks"), pattern = "^chunk_.*csv$", full.names = TRUE)
k <- rbindlist(lapply(cf, function(f) fread(f, select = cols, showProgress = FALSE)))
k <- k[game_type == "R" & !is.na(vx0) & !is.na(release_pos_y) & pitch_type != ""]
k <- unique(k, by = c("game_pk","at_bat_number","pitch_number"))
setorder(k, game_pk, at_bat_number, pitch_number)

prim <- k[pitch_type %in% FASTBALLS, .N, by = .(pitcher, pitch_type)][N >= 50]
prim[, rk := match(pitch_type, FASTBALLS)]
prim <- prim[order(pitcher, rk)][, .SD[1], by = pitcher][, .(pitcher, primary_fb = pitch_type)]

k[, t_plate := (-vy0 - sqrt(vy0^2 - 2*ay*(release_pos_y - yf)))/ay]
k[, t_react := pmax(t_plate - REACT, 0.05)]
lagc <- c("release_pos_x","release_pos_y","release_pos_z","vx0","vy0","vz0","ax","ay","az",
          "t_react","pitch_type","pitch_number","plate_x","plate_z")
for (cc in lagc) k[, (paste0("p_", cc)) := shift(get(cc)), by = .(game_pk, at_bat_number)]
k <- merge(k, prim, by = "pitcher", all.x = TRUE)

k[, consecutive := !is.na(p_pitch_number) & (pitch_number - p_pitch_number == 1L)]
k[, after_fb := consecutive & p_pitch_type %in% FASTBALLS]

pos <- function(r, v, a, t) r + v*t + 0.5*a*t^2
k[, `:=`(pair_path = NA_real_, tunnel = NA_real_, plate_sep = NA_real_)]
idx <- which(k$consecutive)
Tmax <- pmax(k$t_react[idx], k$p_t_react[idx]); acc <- numeric(length(idx))
for (s in 1:NSTEP) {
  tk <- (s - 0.5)/NSTEP * Tmax
  dx <- pos(k$release_pos_x[idx], k$vx0[idx], k$ax[idx], tk) - pos(k$p_release_pos_x[idx], k$p_vx0[idx], k$p_ax[idx], tk)
  dy <- pos(k$release_pos_y[idx], k$vy0[idx], k$ay[idx], tk) - pos(k$p_release_pos_y[idx], k$p_vy0[idx], k$p_ay[idx], tk)
  dz <- pos(k$release_pos_z[idx], k$vz0[idx], k$az[idx], tk) - pos(k$p_release_pos_z[idx], k$p_vz0[idx], k$p_az[idx], tk)
  acc <- acc + sqrt(WX*dx^2 + WY*dy^2 + WZ*dz^2) * (Tmax/NSTEP)
}
k[idx, tunnel := acc]
k[idx, plate_sep := pmax(sqrt((plate_x - p_plate_x)^2 + (plate_z - p_plate_z)^2), 0.1)]
k[idx, pair_path := tunnel / plate_sep]

wd <- c("swinging_strike","swinging_strike_blocked","foul_tip","missed_bunt","bunt_foul_tip")
sd <- c(wd, "foul","hit_into_play","hit_into_play_score","hit_into_play_no_out","foul_bunt")
k[, `:=`(swing = description %in% sd, whiff = description %in% wd, inz = zone >= 1 & zone <= 9,
         xw = ifelse(is.na(estimated_woba_using_speedangle), woba_value, estimated_woba_using_speedangle))]

W <- k[player_name == TARGET_PITCHER]
if (!nrow(W)) stop("no 2026 pitches found for ", TARGET_PITCHER)
HAND <- W$p_throws[1]
OPP  <- if (HAND == "R") "L" else "R"
OPPL <- if (OPP == "L") "LHH" else "RHH"
SAML <- if (OPP == "L") "RHH" else "LHH"
cat(sprintf("%s -- %sHP, %d pitches in 2026, primary fastball = %s\n",
            TARGET_PITCHER, HAND, nrow(W), W$primary_fb[1]))
cat(sprintf("Opposite-handed side = %s\n\n", OPPL))

## ---------------------------------------------------------------------------
cat(sprintf("=== 1. PATH-TO-LOCATION RATIO, %s FB -> ST ===\n", TARGET_PITCHER))
cat("    (lower = tighter tunnel; tunnel = early separation, plate_sep = how far apart they finish)\n\n")
print(W[pitch_type == "ST" & after_fb == TRUE & is.finite(pair_path),
        .(n = .N, path_ratio = round(mean(pair_path), 3), tunnel = round(mean(tunnel), 3),
          plate_sep_ft = round(mean(plate_sep), 2)),
        by = .(side = ifelse(stand == p_throws, sprintf("vs %s (same)", SAML),
                                               sprintf("vs %s (opp)", OPPL)))][order(side)])

cat(sprintf("\n  by setup pitch, vs %s only:\n", OPPL))
print(W[pitch_type == "ST" & stand == OPP & consecutive & is.finite(pair_path),
        .(n = .N, path_ratio = round(mean(pair_path), 3), tunnel = round(mean(tunnel), 3),
          plate_sep_ft = round(mean(plate_sep), 2)), by = .(setup = p_pitch_type)][order(-n)])

lg <- k[pitch_type == "ST" & after_fb == TRUE & p_throws == HAND & stand == OPP & is.finite(pair_path)]
agg <- lg[, .(n = .N, path_ratio = mean(pair_path), tunnel = mean(tunnel),
              plate_sep = mean(plate_sep)), by = .(pitcher, player_name)][n >= MIN_PAIRS]
setorder(agg, path_ratio); agg[, rank := .I]
agg[order(tunnel), rank_tunnel := .I]
agg[order(-plate_sep), rank_sep := .I]

cat(sprintf("\n  League context (%sHP, ST after a FB, vs %s, min %d pairs): %d qualifying pitchers\n",
            HAND, OPPL, MIN_PAIRS, nrow(agg)))
wr <- agg[player_name == TARGET_PITCHER]
if (nrow(wr)) {
  cat(sprintf("  path ratio     : league mean %.3f, median %.3f | %s %.3f -> rank %d/%d (%.0fth pctile)\n",
      mean(agg$path_ratio), median(agg$path_ratio), TARGET_PITCHER, wr$path_ratio,
      wr$rank, nrow(agg), 100*(1 - (wr$rank - 1)/nrow(agg))))
  cat(sprintf("  early sep only : league mean %.3f, median %.3f | %s %.3f -> rank %d/%d (%.0fth pctile)\n",
      mean(agg$tunnel), median(agg$tunnel), TARGET_PITCHER, wr$tunnel,
      wr$rank_tunnel, nrow(agg), 100*(1 - (wr$rank_tunnel - 1)/nrow(agg))))
  cat(sprintf("  plate sep (ft) : league mean %.2f, median %.2f | %s %.2f -> rank %d/%d (1 = finishes farthest apart)\n",
      mean(agg$plate_sep), median(agg$plate_sep), TARGET_PITCHER, wr$plate_sep, wr$rank_sep, nrow(agg)))
} else cat(sprintf("  (%s did not reach %d qualifying pairs)\n", TARGET_PITCHER, MIN_PAIRS))

cat("\n  tightest 10 by path ratio:\n")
print(head(agg[, .(rank, player_name, n, path_ratio = round(path_ratio, 3),
                   tunnel = round(tunnel, 3), plate_sep = round(plate_sep, 2))], 10))
if (nrow(wr)) {
  cat("\n  neighbourhood by path ratio:\n")
  print(agg[abs(rank - wr$rank) <= 3, .(rank, player_name, n, path_ratio = round(path_ratio, 3),
                                        tunnel = round(tunnel, 3), plate_sep = round(plate_sep, 2))])
}
cat("\n  tightest 10 by EARLY SEPARATION alone (the unconfounded tunnel):\n")
print(head(agg[order(tunnel), .(rank_tunnel, player_name, n, tunnel = round(tunnel, 3),
                                plate_sep = round(plate_sep, 2), path_ratio = round(path_ratio, 3))], 10))

## ---------------------------------------------------------------------------
seqlab <- function(dt) ifelse(dt$after_fb == TRUE, "after fastball",
                       ifelse(dt$consecutive == TRUE, "after non-fastball", "first pitch of AB"))
st_o <- W[pitch_type == "ST" & stand == OPP]
st_o[, seq3 := seqlab(st_o)]

cat(sprintf("\n\n=== 2. MISS DISTANCE on %s sweepers vs %s: after a FB vs not ===\n", TARGET_PITCHER, OPPL))
print(st_o[, .(pitches = .N, swings = sum(swing), whiffs = sum(whiff),
               `whiff%` = round(100*sum(whiff)/sum(swing), 1),
               n_miss = sum(swing & is.finite(miss_distance)),
               mean_miss_in = round(mean(miss_distance[swing & is.finite(miss_distance)], na.rm = TRUE), 2)),
           by = .(seq = ifelse(after_fb == TRUE, "after fastball", "not after fastball"))][order(-pitches)])
cat("\n  same, split three ways:\n")
print(st_o[, .(pitches = .N, swings = sum(swing), whiffs = sum(whiff),
               n_miss = sum(swing & is.finite(miss_distance)),
               mean_miss_in = round(mean(miss_distance[swing & is.finite(miss_distance)], na.rm = TRUE), 2)),
           by = seq3][order(-pitches)])

md_fb <- st_o[seq3 == "after fastball" & swing & is.finite(miss_distance), miss_distance]
md_nf <- st_o[seq3 != "after fastball" & swing & is.finite(miss_distance), miss_distance]
if (length(md_fb) >= 3 && length(md_nf) >= 3)
  cat(sprintf("\n  after-FB vs not: %.2f vs %.2f in, t-test p = %.4f (n = %d vs %d)\n",
              mean(md_fb), mean(md_nf), t.test(md_fb, md_nf)$p.value, length(md_fb), length(md_nf)))

lst <- k[pitch_type == "ST" & p_throws == HAND & stand == OPP]
lst[, seq3 := seqlab(lst)]
cat(sprintf("\n  league baseline, %sHP sweepers vs %s:\n", HAND, OPPL))
print(lst[, .(pitches = .N, swings = sum(swing), `whiff%` = round(100*sum(whiff)/sum(swing), 1),
              n_miss = sum(swing & is.finite(miss_distance)),
              mean_miss_in = round(mean(miss_distance[swing & is.finite(miss_distance)], na.rm = TRUE), 2)),
          by = seq3][order(-pitches)])

## ---------------------------------------------------------------------------
cat(sprintf("\n\n=== 3. CHASE RATE on %s sweepers vs %s: after a FB vs not ===\n", TARGET_PITCHER, OPPL))
cat("    (chase = swing at a pitch outside the Savant zone / pitches outside the zone)\n")
print(st_o[, .(pitches = .N, out_of_zone = sum(!inz), chases = sum(swing & !inz),
               `chase%` = round(100*sum(swing & !inz)/sum(!inz), 1),
               `zone%` = round(100*sum(inz)/.N, 1), `swing%` = round(100*sum(swing)/.N, 1)),
           by = .(seq = ifelse(after_fb == TRUE, "after fastball", "not after fastball"))][order(-pitches)])
ct <- st_o[, .(oz = sum(!inz), ch = sum(swing & !inz)),
           by = .(seq = ifelse(after_fb == TRUE, "after fastball", "not after fastball"))]
if (nrow(ct) == 2 && all(ct$oz > 0)) {
  m <- matrix(c(ct$ch[1], ct$oz[1]-ct$ch[1], ct$ch[2], ct$oz[2]-ct$ch[2]), 2)
  cat(sprintf("\n  Fisher exact p on the chase difference = %.4f\n", fisher.test(m)$p.value))
  bt <- binom.test(ct$ch[ct$seq == "after fastball"], ct$oz[ct$seq == "after fastball"])
  cat(sprintf("  after-FB chase = %d/%d, 95%% CI = %.0f%% to %.0f%%\n",
              ct$ch[ct$seq == "after fastball"], ct$oz[ct$seq == "after fastball"],
              100*bt$conf.int[1], 100*bt$conf.int[2]))
}
cat("\n  same, split three ways:\n")
print(st_o[, .(pitches = .N, out_of_zone = sum(!inz), chases = sum(swing & !inz),
               `chase%` = round(100*sum(swing & !inz)/sum(!inz), 1),
               `zone%` = round(100*sum(inz)/.N, 1)), by = seq3][order(-pitches)])
cat(sprintf("\n  league baseline, %sHP sweepers vs %s:\n", HAND, OPPL))
print(lst[, .(pitches = .N, out_of_zone = sum(!inz),
              `chase%` = round(100*sum(swing & !inz)/sum(!inz), 1),
              `zone%` = round(100*sum(inz)/.N, 1)), by = seq3][order(-pitches)])

## ---------------------------------------------------------------------------
sc <- function(dt, by) dt[, .(
  pit = .N, swings = sum(swing), whiff = round(100*sum(whiff)/sum(swing), 1),
  chase = round(100*sum(swing & !inz)/sum(!inz), 1), zone = round(100*sum(inz)/.N, 1),
  xwoba = round(sum(xw[woba_denom > 0 & !is.na(woba_denom)], na.rm = TRUE) /
                sum(woba_denom[woba_denom > 0 & !is.na(woba_denom)], na.rm = TRUE), 3),
  rv100 = round(100*sum(delta_pitcher_run_exp, na.rm = TRUE)/.N, 2),
  velo = round(mean(release_speed, na.rm = TRUE), 1)), by = by]

cat(sprintf("\n\n=== 4. ARSENAL CONTEXT: %s vs %s ===\n", TARGET_PITCHER, OPPL))
a <- sc(W[stand == OPP], "pitch_type"); a[, usage := round(100*pit/sum(pit), 1)]
print(a[order(-pit), .(pitch_type, pit, usage, velo, swings, whiff, chase, zone, xwoba, rv100)])
