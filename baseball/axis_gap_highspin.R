#!/usr/bin/env Rscript

# What is the typical spin-axis gap between a changeup and the four-seamer it plays off, when
# both pitches are genuinely spin-efficient (measured active spin >= 80 percent on each)?
#
# The efficiency floor matters for this question specifically. Active spin below ~80 percent
# means a meaningful share of the rotation is gyro, and a reported clock-face axis describes
# less and less of what the ball is actually doing as that share grows - so the gap between two
# low-efficiency pitches is a noisier quantity than the same number between two clean ones.
# Restricting to high-efficiency pairs is the closest thing to an apples-to-apples axis gap.

suppressPackageStartupMessages({ library(data.table) })
options(width = 200)
MDIR <- "data/statcast_model"

CH <- readRDS(file.path(MDIR, "parachute_ff.rds"))
AS <- readRDS(file.path(MDIR, "active_spin_long.rds"))
CH <- CH[is.finite(axis_diff)]

# axis_diff is already the unsigned separation to the four-seam anchor, wrapped to 0-180.
P <- merge(CH[, .(pitcher, player_name, season, p_throws, axis_diff)],
           AS[pitch_type == "CH", .(pitcher, season, as_ch = active_spin)], by = c("pitcher","season"))
P <- merge(P, AS[pitch_type == "FF", .(pitcher, season, as_fb = active_spin)], by = c("pitcher","season"))
P[, hs := as_ch >= .80 & as_fb >= .80]

S <- P[, .(n = .N, axis = mean(axis_diff), as_ch = as_ch[1], as_fb = as_fb[1], hs = hs[1],
           throws = p_throws[1]), by = .(pitcher, player_name, season)][n >= 50]

f <- function(D, lab) data.table(group = lab, seasons = nrow(D),
  pitches = sum(D$n), mean = mean(D$axis), median = median(D$axis), sd = sd(D$axis),
  p10 = quantile(D$axis,.10), p90 = quantile(D$axis,.90),
  pct_under_10 = 100*mean(D$axis <= 10), pct_under_20 = 100*mean(D$axis <= 20))

cat("=== spin-axis gap, changeup vs four-seamer (pitcher-seasons, 50+ changeups) ===\n")
print(rbind(f(S[hs == TRUE],  "both active spin >= .80"),
            f(S[hs == FALSE], "at least one below .80"),
            f(S,              "all seasons"))[, lapply(.SD, function(x)
              if (is.numeric(x)) round(x,1) else x)], row.names = FALSE)

cat("\n=== same cut, weighting each season by changeups thrown ===\n")
cat(sprintf("  both >= .80 : %.1f deg over %s changeups\n",
            weighted.mean(S[hs==TRUE]$axis, S[hs==TRUE]$n), format(sum(S[hs==TRUE]$n), big.mark=",")))
cat(sprintf("  otherwise   : %.1f deg over %s changeups\n",
            weighted.mean(S[hs==FALSE]$axis, S[hs==FALSE]$n), format(sum(S[hs==FALSE]$n), big.mark=",")))

t <- t.test(S[hs==TRUE]$axis, S[hs==FALSE]$axis)
cat(sprintf("\ndifference high-spin minus rest: %+.1f deg [%+.1f, %+.1f]  p=%.2g\n",
            diff(rev(t$estimate)), t$conf.int[1], t$conf.int[2], t$p.value))

cat("\n=== by handedness, high-spin pairs only ===\n")
print(S[hs == TRUE, .(seasons = .N, mean = round(mean(axis),1), median = round(median(axis),1)),
        by = .(throws)][order(throws)], row.names = FALSE)

cat("\n=== decile of the gap within high-spin pairs ===\n")
H <- S[hs == TRUE][order(axis)]
H[, decile := ceiling(10*seq_len(.N)/.N)]
print(H[, .(seasons = .N, from = round(min(axis),1), to = round(max(axis),1)), by = decile],
      row.names = FALSE)

cat("\n=== tightest 12 high-spin pairs ===\n")
print(H[1:12, .(Pitcher = trimws(paste(sub(".*,\\s*","",player_name), sub(",.*","",player_name))),
                Season = season, CH = n, AxisGap = round(axis,1),
                ActCH = round(as_ch,2), ActFF = round(as_fb,2))], row.names = FALSE)
