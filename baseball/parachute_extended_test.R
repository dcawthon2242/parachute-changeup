#!/usr/bin/env Rscript

# THE PARACHUTE BIN, TESTED ON THREE SEASONS IT WAS NEVER FITTED ON.
#
# Everything about this bin was found in 2023-2026. Statcast's arm angle and active-spin
# leaderboards both reach back to 2020, so 2020-2022 is untouched data measured the same way,
# and the prediction to beat was written down before it was scraped:
#
#     ground balls   +6.9 pp above a shape-and-location model   (p = .0018 in 2023-2026)
#     whiffs         no reliable edge                            (p = .41)
#     run value      no reliable edge                            (p = .72)
#
# Every threshold is frozen at the value it had before the new seasons existed: arm angle
# >= 44 degrees, spin-axis gap <= 10 degrees, active-spin gap within 10 points, 60 swings
# minimum. Nothing is re-tuned. The 44-degree line in particular is NOT recomputed as the top
# quartile of the new population, because that would let the definition drift.

suppressPackageStartupMessages({ library(data.table); library(lightgbm); library(ggplot2) })
set.seed(1); options(width = 205)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")
ARC <- 44; AXC <- 10; EFC <- 0.10; MINSW <- 60L
RES <- file.path(MDIR, "parachute_extended_resid.rds")

FEAT <- c("release_speed","release_spin_rate","release_extension","release_pos_x",
          "release_pos_z","ax","az","speed_diff","ax_diff","az_diff","plate_x","plate_z",
          "plate_x_in","plate_x_arm","z_rel_bot","z_rel_top","VAA","HAA","HAA_in",
          "stand_R","throws_R","same_hand","balls","strikes")

oof <- function(D, y, obj, tag) {
  K <- 4; fold <- sample(rep(1:K, length.out = nrow(D))); p <- rep(NA_real_, nrow(D))
  for (f in 1:K) {
    tr <- D[fold != f]; n <- nrow(tr); vi <- sample(n, floor(.12*n))
    dtr <- lgb.Dataset(as.matrix(tr[-vi, ..FEAT]), label = y[fold != f][-vi])
    dva <- lgb.Dataset.create.valid(dtr, as.matrix(tr[vi, ..FEAT]), label = y[fold != f][vi])
    m <- lgb.train(params = list(objective = obj,
                   metric = if (obj == "binary") "binary_logloss" else "l2",
                   learning_rate = .06, num_leaves = 31, min_data_in_leaf = 300,
                   feature_fraction = .8, bagging_fraction = .8, bagging_freq = 1),
                   data = dtr, nrounds = 1500, valids = list(v = dva),
                   early_stopping_rounds = 50, verbose = -1)
    p[fold == f] <- predict(m, as.matrix(D[fold == f, ..FEAT]))
  }
  cat(sprintf("  %-9s n=%8s  R2=%.4f\n", tag, format(nrow(D), big.mark=","),
              1 - var(y - p)/var(y)))
  y - p
}

if (!file.exists(RES) || nzchar(Sys.getenv("REFIT"))) {
  CH <- readRDS(file.path(MDIR, "parachute_extended.rds"))
  CH <- CH[is.finite(axis_diff) & stats::complete.cases(CH[, ..FEAT])]
  cat("=== out-of-fold models, 2020-2026 pooled; axis gap and arm angle both excluded ===\n")
  CH[, rv_res := oof(CH, CH$rv, "regression", "run value")]
  SW <- CH[is_swing == TRUE]; SW[, wh_res := oof(SW, SW$whiff, "binary", "whiff")]
  BP <- CH[is_bip == TRUE]; BP[, gb := as.integer(bb_type == "ground_ball")]
  BP[, gb_res := oof(BP, BP$gb, "binary", "grounder")]
  saveRDS(list(CH = CH[, .(pitcher, player_name, season, axis_diff, az_diff, speed_diff,
                           arm_angle, release_spin_rate, rv_res, is_swing, is_bip)],
               SW = SW[, .(pitcher, season, wh_res)],
               BP = BP[, .(pitcher, season, gb_res, launch_speed)]), RES)
}
L <- readRDS(RES)
S <- L$CH[, .(np = .N, nsw = sum(is_swing), axis = mean(axis_diff), kill = -mean(az_diff),
              velo_sep = -mean(speed_diff), arm = mean(arm_angle, na.rm = TRUE),
              spin = mean(release_spin_rate, na.rm = TRUE), rv100 = 100*mean(rv_res)),
          by = .(pitcher, player_name, season)]
S <- merge(S, L$SW[, .(wh = 100*mean(wh_res)), by = .(pitcher, season)], by = c("pitcher","season"))
S <- merge(S, L$BP[, .(gb = 100*mean(gb_res)), by = .(pitcher, season)], by = c("pitcher","season"))

# active spin for the changeup and for the pitcher's primary fastball
AS   <- readRDS(file.path(MDIR, "active_spin_long.rds"))
prim <- fread(file.path(MDIR, "parachute_extended_fbtype.csv"))
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, as_ch = active_spin)],
           by = c("pitcher","season"), all.x = TRUE)
S <- merge(S, prim, by = c("pitcher","season"), all.x = TRUE)
S <- merge(S, AS[, .(pitcher, season, fb_type = pitch_type, as_fb = active_spin)],
           by = c("pitcher","season","fb_type"), all.x = TRUE)
S[, as_gap := as_ch - as_fb]
S <- S[nsw >= MINSW & is.finite(axis) & is.finite(arm)]
S[, era := fifelse(season <= 2022L, "2020-2022 (new, never fitted)", "2023-2026 (discovery)")]
S[, para := axis <= AXC & arm >= ARC & is.finite(as_gap) & abs(as_gap) <= EFC]

cat(sprintf("\n%d pitcher-seasons; %d with measured active spin on both pitches\n",
            nrow(S), sum(is.finite(S$as_gap))))
print(S[, .(seasons = .N, with_as = sum(is.finite(as_gap)), in_bin = sum(para)),
        by = season][order(season)], row.names = FALSE)

A <- S[is.finite(as_gap)]
cat("\n=== THE REPLICATION ===\n")
out <- rbindlist(lapply(c("2020-2022 (new, never fitted)","2023-2026 (discovery)","ALL"), function(e) {
  q <- if (e == "ALL") A else A[era == e]
  rbindlist(lapply(c("gb","wh","rv100"), function(v) {
    t <- t.test(q[para == TRUE][[v]], q[para == FALSE][[v]])
    data.table(era = e, outcome = v, n_bin = sum(q$para), diff = diff(rev(t$estimate)),
               lo = t$conf.int[1], hi = t$conf.int[2], p = t$p.value) })) }))
print(out[, .(era, outcome, n_bin, diff = round(diff,2), ci = sprintf("[%+.2f, %+.2f]", lo, hi),
              p = round(p,4))], row.names = FALSE)

cat("\n=== who is in the bin in the new seasons ===\n")
print(A[para == TRUE & season <= 2022, .(player_name, season, pitches = np,
        arm = round(arm,1), axis = round(axis,1), as_gap = round(as_gap,2),
        gb_over = round(gb,1), whiff_over = round(wh,1), rv100 = round(rv100,2))][order(-gb_over)],
      row.names = FALSE)

cat("\n=== threshold grid on the NEW seasons only ===\n")
G <- CJ(ax = c(8,10,12,15,20), ar = c(40,42,44,46,48))
N <- A[season <= 2022]
G[, c("n","d","p") := { z <- mapply(function(a,r) {
    i <- N$axis <= a & N$arm >= r & abs(N$as_gap) <= EFC
    if (sum(i) < 5) return(c(sum(i), NA, NA))
    t <- t.test(N$gb[i], N$gb[!i]); c(sum(i), diff(rev(t$estimate)), t$p.value) }, ax, ar)
  list(z[1,], z[2,], z[3,]) }]
print(dcast(G, ax ~ ar, value.var = "d")[, lapply(.SD, function(x) round(x,2))], row.names = FALSE)
cat(sprintf("  %d of %d cells positive, %d reach p<.05\n",
            sum(G$d > 0, na.rm=TRUE), sum(!is.na(G$d)), sum(G$p < .05, na.rm=TRUE)))

fwrite(S[order(-para, -gb)], file.path(AST, "ext_parachute_extended.csv"))

## ---- figure --------------------------------------------------------------------------
TARGET <- 6.93   # the pre-registered value, published before 2020-2022 existed
P <- out[outcome == "gb"][era != "ALL"]
P[, era := factor(era, levels = c("2023-2026 (discovery)","2020-2022 (new, never fitted)"))]
gg <- ggplot(P, aes(era, diff, fill = era)) +
  geom_hline(yintercept = 0, linewidth = .4) +
  geom_hline(yintercept = TARGET, linetype = "dashed", colour = "#c0392b", linewidth = .5) +
  annotate("text", x = Inf, y = TARGET + 0.55, hjust = 1.03, size = 3.1, colour = "#c0392b",
           fontface = "bold", label = "pre-registered target: +6.93 pp") +
  geom_col(width = .55, colour = "black", linewidth = .3, show.legend = FALSE) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .12, linewidth = .5) +
  geom_text(aes(label = sprintf("%+.2f pp\np = %.3f\n%d seasons in bin", diff, p, n_bin)),
            vjust = -0.35, size = 3.4, lineheight = .95) +
  scale_fill_manual(values = c("#9aa5b1", "#1d7870")) +
  coord_cartesian(ylim = c(min(P$lo) - 1, max(TARGET, max(P$hi)) + 3)) +
  labs(title = "The parachute bin scored on three seasons it was never fitted on: it does not replicate",
       subtitle = paste0("Ground-ball rate above a shape-and-location model that never sees the spin axis or the arm angle. Every threshold is frozen at the value it held before\n",
                         "2020-2022 was scraped: arm angle 44 degrees, axis gap 10 degrees, active-spin gap 10 points, 60 swings minimum. Thirteen brand-new pitcher-seasons\n",
                         "qualify and they land at -0.03 points. The discovery seasons themselves fall from +6.93 to +3.04, not because the model changed - the original eleven\n",
                         "seasons still read +6.16 under the new model - but because including the first half of 2023 reshuffles who qualifies, which is its own verdict on how\n",
                         "stable this bin was."),
       x = NULL, y = "Ground-ball% above model, bin minus rest") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 8.2),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
ggsave(file.path(AST, "fig26_parachute_replication.png"), gg, width = 9.5, height = 6.4, dpi = 150)
cat("\nwrote fig26_parachute_replication.png, ext_parachute_extended.csv\n")
