#!/usr/bin/env Rscript

# Use Statcast's new "flail rate" (swing where the ball is OUTSIDE the barrel
# horizontally -- hitter reaching/fooled) as the target, and test whether the
# matched-look + velo-kill deception profile predicts it.
# Flail data: Savant swing-timing/miss-distance leaderboard (min 25 swings).

suppressPackageStartupMessages({ library(data.table) })
dir <- file.path("data", "statcast_2026")

fl <- fread(file.path(dir, "savant_swing_timing_2026.csv"), showProgress = FALSE)
fl <- fl[, .(pitcher = id, s_type = api_pitch_type, name, flail = flailed_percent,
             flail_in = avg_x_flail, miss_dist = miss_distance, whiff_rt = whiff_rate,
             over = over_percent, under = under_percent, early = early_percent,
             late = late_percent, st_n = n_swings)]

pp <- fread(file.path(dir, "pitch_pair_all_types_2026.csv"), showProgress = FALSE)
m  <- merge(pp, fl, by = c("pitcher","s_type"))

wm <- function(x, w) round(weighted.mean(x, w, na.rm = TRUE), 3)

# ---- league flail baselines by pitch type ---------------------------------
cat("=== League flail% baseline by pitch type (min 25 swings) ===\n")
base <- fl[, .(n_pitchers=.N, swings=sum(st_n), flail=wm(flail, st_n),
               miss_dist=wm(miss_dist, st_n)), by=s_type][order(-flail)]
print(base)

# ---- Profile changeups (matched look + big velo kill) vs other CH ----------
chL <- fl[s_type=="CH"]
prof_ids <- pp[s_type=="CH" & big_kill==TRUE & whiff_resid>0, unique(pitcher)]
chL[, group := fifelse(pitcher %in% prof_ids, "profile CH (matched+kill, overperf)", "other CH")]
cat("\n=== Flail: profile changeups vs the rest of MLB changeups ===\n")
print(chL[, .(n=.N, swings=sum(st_n), flail=wm(flail, st_n),
              flail_inches=wm(flail_in, st_n), miss_dist=wm(miss_dist, st_n),
              whiff=wm(whiff_rt, st_n)), by=group])

# ---- Does the mechanism predict flail? (changeups) ------------------------
ch <- m[s_type=="CH" & st_n>=25]
cat(sprintf("\n=== Changeup flail correlations (n=%d) ===\n", nrow(ch)))
cat(sprintf("  flail ~ velo_diff : r = %+.2f\n", cor(ch$flail, ch$velo_diff, use="complete.obs")))
cat(sprintf("  flail ~ axis_diff : r = %+.2f\n", cor(ch$flail, ch$axis_diff, use="complete.obs")))
cat(sprintf("  flail ~ eff_diff  : r = %+.2f\n", cor(ch$flail, ch$eff_diff,  use="complete.obs")))
cat(sprintf("  flail ~ whiff_resid(over stuff): r = %+.2f\n", cor(ch$flail, ch$whiff_resid, use="complete.obs")))
fit <- lm(flail ~ velo_diff + axis_diff + eff_diff, data=ch, weights=st_n)
cat("\n  flail ~ velo_diff + axis_diff + eff_diff (weighted):\n")
print(round(summary(fit)$coefficients, 4))

# ---- Flail as target across ALL matched-look secondaries -------------------
mm <- m[big_kill==TRUE & st_n>=25]
cat(sprintf("\n=== All matched-look + big-kill secondaries: flail vs velo kill (n=%d) ===\n", nrow(mm)))
cat(sprintf("  flail ~ velo_diff : r = %+.2f | flail ~ whiff_resid : r = %+.2f\n",
    cor(mm$flail, mm$velo_diff, use="complete.obs"),
    cor(mm$flail, mm$whiff_resid, use="complete.obs")))

# ---- leaderboard: profile changeups ranked by flail ------------------------
cat("\n=== Profile changeups ranked by FLAIL rate ===\n")
top <- ch[big_kill==TRUE & whiff_resid>0][order(-flail)]
print(top[, .(name, FB=a_type, sw=st_n, Dvelo=round(velo_diff,1), axis=round(axis_diff,1),
  flail=round(flail,2), flail_in=round(flail_in,1), miss=round(miss_dist,1),
  whiff=round(whiff_rt,2), wResid=round(whiff_resid,1))], nrows=40)

fwrite(m, file.path(dir, "deception_with_flail_2026.csv"))
cat(sprintf("\nWrote %s\n", file.path(dir, "deception_with_flail_2026.csv")))
