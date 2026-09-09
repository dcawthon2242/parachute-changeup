#!/usr/bin/env Rscript

# CHANGEUP AGAINST THE SAME PITCHER'S OWN SINKER, IN THE SAME COUNTS.
#
# The pooled comparison in the last run is confounded twice over: sinkers and changeups are thrown by
# different pitchers in different counts, and the changeup population in the study is quality-selected
# by the swing floor. Both go away if the comparison is made inside a pitcher-season, with the count
# distribution held fixed.
#
# Every pitcher-season with 200+ changeups and 200+ sinkers. Run value per 100 is computed for each
# pitch type using the SAME count weights - that season's own combined distribution - so a difference
# cannot come from one pitch being saved for two-strike counts.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 200); MDIR <- "data/statcast_model"

A <- readRDS(file.path(MDIR, "parachute_rv.rds")); setDT(A)
A <- A[pitch_type %in% c("CH","SI") & is.finite(rv)]
A[, `:=`(id = as.character(pitcher), cnt = paste0(balls, "-", strikes))]
A[, gb := as.integer(bb_type == "ground_ball")]
A[, hard := as.integer(launch_speed >= 95)]

N <- A[, .(n = .N), by = .(id, season, pitch_type)]
N <- dcast(N, id + season ~ pitch_type, value.var = "n", fill = 0)
KEEP <- N[CH >= 200 & SI >= 200]
A <- merge(A, KEEP[, .(id, season)], by = c("id","season"))
cat(sprintf("%d pitcher-seasons throw 200+ of both, %d arms, %d pitches.\n",
            nrow(KEEP), uniqueN(KEEP$id), nrow(A)))

# count weights: each season's own combined distribution, applied to both pitch types
W <- A[, .(w = .N), by = .(id, season, cnt)]
CM <- A[, .(m = mean(rv), n = .N), by = .(id, season, cnt, pitch_type)]
CM <- merge(CM, W, by = c("id","season","cnt"))
CM <- CM[n >= 3]                                    # a cell needs a few pitches to mean anything
RVP <- CM[, .(rv100 = 100*sum(m*w)/sum(w), cells = .N), by = .(id, season, pitch_type)]
RVP <- dcast(RVP[cells >= 6], id + season ~ pitch_type, value.var = "rv100")
setnames(RVP, c("CH","SI"), c("ch_rv","si_rv"))

Q <- A[!is.na(bb_type) & bb_type != "", .(bip = .N, gb = 100*mean(gb),
        hard = 100*mean(hard, na.rm = TRUE), ev = mean(launch_speed, na.rm = TRUE),
        xw = 1000*mean(estimated_woba_using_speedangle, na.rm = TRUE)),
       by = .(id, season, pitch_type)]
Q <- dcast(Q[bip >= 30], id + season ~ pitch_type,
           value.var = c("bip","gb","hard","ev","xw"))
J <- merge(RVP, Q, by = c("id","season"))
J <- merge(J, unique(A[, .(id, season, nm = player_name)]), by = c("id","season"))
J <- merge(J, KEEP[, .(id, season, nch = CH, nsi = SI)], by = c("id","season"))
J <- J[complete.cases(J[, .(ch_rv, si_rv, gb_CH, gb_SI, xw_CH, xw_SI)])]
cat(sprintf("%d of them clear the 30-ball-in-play floor on both pitches.\n\n", nrow(J)))

# separation, from the study file, to split the seasons the same way as before
P <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(P)
SEP <- P[is.finite(speed_diff), .(sep = mean(-speed_diff), slot = mean(arm_angle, na.rm = TRUE)),
         by = .(id = as.character(pitcher), season)]
J <- merge(J, SEP, by = c("id","season"))
J[, t3 := cut(sep, quantile(sep, 0:3/3), include.lowest = TRUE,
              labels = c("T1 least sep","T2","T3 most sep"))]

pt <- function(a, b) { t <- t.test(a, b, paired = TRUE); sprintf("%+.1f%s", mean(a-b),
  ifelse(t$p.value < .01, "**", ifelse(t$p.value < .05, "*", ""))) }

cat("=== CHANGEUP MINUS THE SAME PITCHER'S SINKER, COUNT-MATCHED ===\n\n")
T <- J[, .(seasons = .N, sep = mean(sep),
           ch_rv = mean(ch_rv), si_rv = mean(si_rv), d_rv = mean(ch_rv - si_rv),
           ch_gb = mean(gb_CH), si_gb = mean(gb_SI),
           ch_hard = mean(hard_CH), si_hard = mean(hard_SI),
           ch_xw = mean(xw_CH), si_xw = mean(xw_SI)), by = t3][order(t3)]
print(T[, .(t3, seasons, sep = round(sep,1),
            ch_rv100 = round(ch_rv,2), si_rv100 = round(si_rv,2), diff_rv = round(d_rv,2),
            ch_gb = round(ch_gb,1), si_gb = round(si_gb,1),
            ch_hard = round(ch_hard,1), si_hard = round(si_hard,1),
            ch_xwcon = round(ch_xw,0), si_xwcon = round(si_xw,0))], row.names = FALSE)

cat("\n  paired differences with significance (changeup minus own sinker):\n\n")
cat(sprintf("  %-16s %8s %10s %10s %10s %10s\n", "group", "seasons", "run value", "GB%", "hard%", "xwOBAcon"))
for (g in c(levels(J$t3), "all")) {
  X <- if (g == "all") J else J[t3 == g]
  cat(sprintf("  %-16s %8d %10s %10s %10s %10s\n", g, nrow(X),
              pt(X$ch_rv, X$si_rv), pt(X$gb_CH, X$gb_SI), pt(X$hard_CH, X$hard_SI),
              pt(X$xw_CH, X$xw_SI)))
}

cat("\n=== IS THE LOW-SEPARATION CHANGEUP A BETTER SINKER OR A WORSE CHANGEUP? ===\n\n")
cat("  each pitcher-season's changeup, ranked against the LEAGUE sinker baseline and against\n")
cat("  the league changeup baseline. share of seasons that beat each.\n\n")
lgSI <- J[, .(gb = weighted.mean(gb_SI, bip_SI), hard = weighted.mean(hard_SI, bip_SI),
              xw = weighted.mean(xw_SI, bip_SI))]
lgCH <- J[, .(gb = weighted.mean(gb_CH, bip_CH), hard = weighted.mean(hard_CH, bip_CH),
              xw = weighted.mean(xw_CH, bip_CH))]
cat(sprintf("  league sinker: GB %.1f  hard %.1f  xwOBAcon %.0f\n", lgSI$gb, lgSI$hard, lgSI$xw))
cat(sprintf("  league changeup: GB %.1f  hard %.1f  xwOBAcon %.0f\n\n", lgCH$gb, lgCH$hard, lgCH$xw))
B <- J[, .(seasons = .N,
           beats_SI_xw = 100*mean(xw_CH < lgSI$xw), beats_CH_xw = 100*mean(xw_CH < lgCH$xw),
           beats_SI_gb = 100*mean(gb_CH > lgSI$gb), beats_SI_hard = 100*mean(hard_CH < lgSI$hard)),
       by = t3][order(t3)]
print(B[, .(t3, seasons, `% better xwOBAcon than a sinker` = round(beats_SI_xw),
            `% better xwOBAcon than a changeup` = round(beats_CH_xw),
            `% more grounders than a sinker` = round(beats_SI_gb),
            `% softer than a sinker` = round(beats_SI_hard))], row.names = FALSE)

cat("\n  the ten lowest-separation seasons with both pitches:\n\n")
print(J[order(sep)][1:10, .(pitcher = nm, season, sep = round(sep,1), slot = round(slot,0),
        ch = nch, si = nsi, ch_rv100 = round(ch_rv,2), si_rv100 = round(si_rv,2),
        ch_gb = round(gb_CH,1), si_gb = round(gb_SI,1), ch_xwcon = round(xw_CH,0),
        si_xwcon = round(xw_SI,0))], row.names = FALSE)
saveRDS(J, file.path(MDIR, "ch_vs_si_paired.rds"))
