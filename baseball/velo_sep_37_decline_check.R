#!/usr/bin/env Rscript

# IS THE DOWNGRADE RIGHT?
#
# The decliners share a profile, but sharing a profile does not make the cut correct. Three checks.
#   1. Does the miss shrink for the arms that were cut, or does the compression simply move them?
#   2. Where do the archetypes land, and does the decline concentrate in the group already known to
#      underperform?
#   3. Does the decline hold up out of sample - if an arm is cut in one season, is the cut vindicated
#      by the NEXT season's actual whiff rate?

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 210); MDIR <- "data/statcast_model"
D <- readRDS(file.path(MDIR, "smean_movers.rds")); setDT(D)

## ---- 1. does the miss shrink where the cut happens ----
cat("=== 1. DOES THE CUT REDUCE THE ERROR? ===\n\n")
D[, bucket := cut(mv, c(-Inf,-20,-10,-5,5,10,20,Inf),
                  labels = c("cut > 20","cut 10-20","cut 5-10","little change",
                             "raise 5-10","raise 10-20","raise > 20"))]
B <- D[, .(seasons = .N, swings = sum(nsw), actual = weighted.mean(act, nsw),
           before = weighted.mean(s0, nsw), after = weighted.mean(s1, nsw),
           rmse0 = sqrt(weighted.mean(e0^2, nsw)), rmse1 = sqrt(weighted.mean(e1^2, nsw)),
           better = mean(abs(e1) < abs(e0))), by = bucket][order(bucket)]
print(B[, .(bucket, seasons, swings, actual = round(actual,1), before = round(before,1),
            after = round(after,1), miss_before = round(actual-before,2),
            miss_after = round(actual-after,2), rmse_before = round(rmse0,2),
            rmse_after = round(rmse1,2), `% improved` = round(100*better))], row.names = FALSE)

## ---- 2. archetype membership among the cut ----
cat("\n=== 2. WHO GETS CUT, BY ARCHETYPE ===\n\n")
TAGS <- c("A2_wide","A4_broad","A5_seam","A6_extreme","A7_trad","A8_mismatch")
LAB <- c(A2_wide = "Matched axis, wide", A4_broad = "Broad (axis + separation)",
         A5_seam = "Seam-shifted matched", A6_extreme = "Extreme seam shift",
         A7_trad = "Traditional, low seam dev", A8_mismatch = "Mismatch (underperforms)")
A <- rbindlist(lapply(TAGS, function(t) {
  X <- D[get(t) == TRUE]
  data.table(group = LAB[[t]], seasons = nrow(X), move = mean(X$mv),
             `% cut` = 100*mean(X$mv < 0), miss_before = weighted.mean(X$e0, X$nsw),
             miss_after = weighted.mean(X$e1, X$nsw))
}))
A <- rbind(A, data.table(group = "everyone else", seasons = D[tag == "", .N],
                         move = D[tag == "", mean(mv)], `% cut` = 100*D[tag == "", mean(mv < 0)],
                         miss_before = D[tag == "", weighted.mean(e0, nsw)],
                         miss_after = D[tag == "", weighted.mean(e1, nsw)]))
print(A[, .(group, seasons, move = round(move,1), `% cut` = round(`% cut`),
            miss_before = round(miss_before,2), miss_after = round(miss_after,2))], row.names = FALSE)

## ---- 3. out of sample: does the next season agree with the cut ----
cat("\n=== 3. IS THE CUT VINDICATED BY THE FOLLOWING SEASON? ===\n\n")
N <- D[, .(id, nm, season, mv, nsw, act, s0, s1)]
N2 <- copy(N)[, season := season - 1][, .(id, season, nxt_act = act, nxt_n = nsw,
                                          nxt_s0 = s0, nxt_s1 = s1)]
J <- merge(N, N2, by = c("id","season"))[nsw >= 75 & nxt_n >= 75]
cat(sprintf("%d consecutive-season pairs, %d arms.\n\n", nrow(J), uniqueN(J$id)))
J[, grp := fifelse(mv < -10, "cut > 10", fifelse(mv > 10, "raised > 10", "little change"))]
K <- J[, .(pairs = .N, this_yr = weighted.mean(act, nsw), next_yr = weighted.mean(nxt_act, nxt_n),
           next_pred_before = weighted.mean(nxt_s0, nxt_n),
           next_pred_after = weighted.mean(nxt_s1, nxt_n)), by = grp]
K[, `:=`(miss_before = next_yr - next_pred_before, miss_after = next_yr - next_pred_after)]
print(K[order(-pairs), .(grp, pairs, this_yr = round(this_yr,1), next_yr = round(next_yr,1),
                         next_miss_before = round(miss_before,2),
                         next_miss_after = round(miss_after,2))], row.names = FALSE)
cat(sprintf("\n  correlation of this year's move with next year's actual whiff: %+.3f\n",
            J[, cor(mv, nxt_act)]))
cat(sprintf("  correlation of this year's move with next year's MISS under the old model: %+.3f\n",
            J[, cor(mv, nxt_act - nxt_s0)]))
