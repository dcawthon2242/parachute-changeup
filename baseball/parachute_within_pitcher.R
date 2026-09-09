#!/usr/bin/env Rscript

# PARACHUTE CHANGEUPS, WITHIN PITCHER.
#
# Every between-pitcher test of this idea is confounded, and the confounds all run the same
# way: a pitcher whose changeup spins like his fastball is a different pitcher from one whose
# doesn't - different spin rate, so different Hawk-Eye axis precision; different arsenal;
# different role. Differencing a pitcher against himself removes all of it. If matching the
# fastball's spin axis does anything, a pitcher who moves his changeup axis TOWARD his
# fastball should get better, and the same man a year earlier is the control.
#
# Pre-registered: outcomes are residuals from shape-and-location models that exclude axis gap,
# so a grip change that also changes movement is already priced in and only the axis part is
# left. Direction of interest is NEGATIVE - shrinking the axis gap should raise performance.
# A backward placebo is run alongside: this season's axis change must not explain last
# season's performance change. If the placebo fires, the result is noise.

suppressPackageStartupMessages({ library(data.table); library(lightgbm); library(ggplot2) })
set.seed(1); options(width = 205)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")
RES <- file.path(MDIR, "parachute_within.rds")

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
                   data = dtr, nrounds = 1500, valids = list(val = dva),
                   early_stopping_rounds = 50, verbose = -1)
    p[fold == f] <- predict(m, as.matrix(D[fold == f, ..FEAT]))
  }
  cat(sprintf("  %-9s n=%7d  R2=%.4f\n", tag, nrow(D), 1 - var(y - p)/var(y)))
  y - p
}

if (!file.exists(RES) || nzchar(Sys.getenv("REFIT"))) {
  d <- readRDS(file.path(MDIR, "parachute_rv.rds"))
  sw   <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip","hit_into_play",
            "foul_bunt","missed_bunt","bunt_foul_tip")
  miss <- c("swinging_strike","swinging_strike_blocked","foul_tip","missed_bunt")
  CH <- d[pitch_type == "CH" & is.finite(axis_diff)]
  CH <- CH[stats::complete.cases(CH[, ..FEAT])]
  CH[, `:=`(is_swing = description %in% sw, is_bip = !is.na(launch_speed) & bb_type != "",
            whiff = as.integer(description %in% miss))]
  cat("=== out-of-fold models, axis gap excluded from every one ===\n")
  CH[, rv_res := oof(CH, CH$rv, "regression", "run value")]
  SWD <- CH[is_swing == TRUE]; SWD[, wh_res := oof(SWD, SWD$whiff, "binary", "whiff")]
  BIP <- CH[is_bip == TRUE]; BIP[, gb := as.integer(bb_type == "ground_ball")]
  BIP[, gb_res := oof(BIP, BIP$gb, "binary", "grounder")]
  saveRDS(list(CH = CH[, .(pitcher, player_name, season, axis_diff, az_diff, ax_diff,
                           speed_diff, release_spin_rate, rv_res, is_swing, is_bip)],
               SWD = SWD[, .(pitcher, season, wh_res)],
               BIP = BIP[, .(pitcher, season, gb_res, launch_speed,
                             xw = estimated_woba_using_speedangle)]), RES)
}
L <- readRDS(RES); CH <- L$CH; SWD <- L$SWD; BIP <- L$BIP

S <- CH[, .(np = .N, nsw = sum(is_swing), axis = mean(axis_diff), kill = -mean(az_diff),
            velo_sep = -mean(speed_diff), spin = mean(release_spin_rate, na.rm = TRUE),
            rv100 = 100*mean(rv_res)), by = .(pitcher, player_name, season)]
S <- merge(S, SWD[, .(wh = 100*mean(wh_res)), by = .(pitcher, season)], by = c("pitcher","season"))
S <- merge(S, BIP[, .(nbip = .N, gb = 100*mean(gb_res), ev = mean(launch_speed, na.rm=TRUE),
                      xw = mean(xw, na.rm = TRUE)), by = .(pitcher, season)],
           by = c("pitcher","season"))
MINSW <- 50L; S <- S[nsw >= MINSW]

## ---- consecutive-season pairs, each pitcher his own control ---------------------
# Each row is season t carrying season t-1 alongside it. Written out longhand because the
# obvious one-line merge silently puts the later season in the earlier season's slot.
PREV <- S[, .(pitcher, season = season + 1L, axis_p = axis, kill_p = kill, velo_p = velo_sep,
              rv_p = rv100, wh_p = wh, gb_p = gb, ev_p = ev, xw_p = xw, spin_p = spin)]
P <- merge(S, PREV, by = c("pitcher","season"))
# d_axis is prev minus current, so POSITIVE means the gap CLOSED toward the fastball.
# Every performance delta is current minus prev, so POSITIVE means the pitcher improved.
P[, `:=`(d_axis = axis_p - axis, d_kill = kill - kill_p, d_velo = velo_sep - velo_p,
         d_rv = rv100 - rv_p, d_wh = wh - wh_p, d_gb = gb - gb_p, d_ev = ev - ev_p,
         d_xw = xw - xw_p, d_spin = spin - spin_p)]
cat(sprintf("\n%d consecutive-season pairs, %d pitchers, min %d swings each side\n",
            nrow(P), uniqueN(P$pitcher), MINSW))
cat(sprintf("axis-gap change: SD %.1f deg, 10th/90th pct %+.1f / %+.1f\n",
            sd(P$d_axis), quantile(P$d_axis,.1), quantile(P$d_axis,.9)))

sp <- function(x, y) { z <- suppressWarnings(cor.test(x, y, method = "spearman"))
                       list(r = unname(z$estimate), p = z$p.value, n = sum(is.finite(x*y))) }
# d_axis is (last year - this year), so POSITIVE = the gap shrank = moved toward the fastball.
# A positive correlation with a change in performance is therefore the parachute prediction.
cat("\n=== within-pitcher: did closing the axis gap improve anything? ===\n")
cat("    (d_axis > 0 means the changeup moved TOWARD the fastball's spin axis)\n")
OUT <- rbindlist(lapply(list(c("d_rv","Run value /100 over model"), c("d_wh","Whiff% over model"),
                             c("d_gb","Ground-ball% over model"), c("d_ev","Exit velo (neg=better)"),
                             c("d_xw","xwOBAcon (neg=better)")), function(v) {
  z <- sp(P$d_axis, P[[v[1]]]); data.table(outcome = v[2], r = z$r, p = z$p, n = z$n) }))
OUT[, sig := fifelse(p < .05, "*", "")]
print(OUT[, .(outcome, r = round(r,3), p = round(p,4), n, sig)], row.names = FALSE)

## ---- controls: grip changes move more than the axis -----------------------------
cat("\n=== after removing co-moving changes in drop, velocity separation and spin ===\n")
rz <- function(v) residuals(lm(as.formula(sprintf("%s ~ d_kill + d_velo + d_spin", v)), data = P))
ra <- residuals(lm(d_axis ~ d_kill + d_velo + d_spin, data = P))
for (v in c("d_rv","d_wh","d_gb")) { z <- sp(ra, rz(v))
  cat(sprintf("  %-6s partial r=%+.3f  p=%.3f\n", v, z$r, z$p)) }

## ---- placebo: this year's axis change must not explain last year's results ------
cat("\n=== backward placebo (must be null) ===\n")
B <- merge(P[, .(pitcher, season, d_axis)],
           P[, .(pitcher, season = season + 1L, d_rv_prev = d_rv, d_wh_prev = d_wh,
                 d_gb_prev = d_gb)], by = c("pitcher","season"))
for (v in c("d_rv_prev","d_wh_prev","d_gb_prev")) { z <- sp(B$d_axis, B[[v]])
  cat(sprintf("  %-11s r=%+.3f  p=%.3f  n=%d\n", v, z$r, z$p, z$n)) }

## ---- event study: pitchers who crossed the bin boundary -------------------------
P[, grp := fifelse(axis_p > 15 & axis <= 15, "moved INTO parachute",
           fifelse(axis_p <= 15 & axis > 15, "moved OUT of parachute",
           fifelse(axis_p <= 15 & axis <= 15, "stayed parachute", "never parachute")))]
cat("\n=== event study: change in performance by boundary crossing ===\n")
E <- P[, .(n = .N, d_axis = round(mean(d_axis),1), d_rv = round(mean(d_rv),3),
           d_wh = round(mean(d_wh),2), d_gb = round(mean(d_gb),2)),
       by = grp][order(-d_axis)]
print(E, row.names = FALSE)
ino <- P[grp == "moved INTO parachute"]; oth <- P[grp == "never parachute"]
for (v in c("d_rv","d_wh","d_gb")) { t <- t.test(ino[[v]], oth[[v]])
  cat(sprintf("  INTO vs never, %-5s  %+.3f vs %+.3f   p=%.3f\n", v,
              mean(ino[[v]]), mean(oth[[v]]), t$p.value)) }

cat("\n=== biggest axis-gap closers ===\n")
print(P[order(-d_axis)][1:12, .(player_name, yr = sprintf("%d->%d", season-1L, season),
        axis = sprintf("%.0f->%.0f", axis_p, axis), d_rv = round(d_rv,2),
        d_wh = round(d_wh,2), d_gb = round(d_gb,2))], row.names = FALSE)
fwrite(P[order(-d_axis)], file.path(AST, "ext_parachute_within_pitcher.csv"))

## ---- figure ---------------------------------------------------------------------
LAB <- c(d_rv = "Run value /100 over model", d_wh = "Whiff% over model",
         d_gb = "Ground-ball% over model")
M <- melt(P, id.vars = c("player_name","d_axis"), measure.vars = names(LAB))
M[, variable := factor(LAB[as.character(variable)], levels = LAB)]
ann <- M[, { z <- sp(d_axis, value)
             .(lab = sprintf("r = %+.3f   p = %.2f   n = %d", z$r, z$p, z$n)) }, by = variable]
gg <- ggplot(M, aes(d_axis, value)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = .35) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = .35) +
  geom_point(alpha = .32, size = 1.25, colour = "#1d3557") +
  geom_smooth(method = "lm", formula = y ~ x, colour = "#c0392b", fill = "#c0392b",
              alpha = .13, linewidth = .85) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab), hjust = -0.06, vjust = 1.5,
            size = 3.1, fontface = "bold", inherit.aes = FALSE) +
  facet_wrap(~ variable, scales = "free_y", nrow = 1) +
  labs(title = "The parachute changeup, tested against the only fair control there is: the same pitcher a year earlier",
       subtitle = paste0("One point per consecutive pitcher-season pair, min ", MINSW, " changeup swings on each side. The x-axis is how far the changeup's spin axis moved TOWARD the primary\n",
                         "fastball between seasons, so points on the right are pitchers who became more parachute-like. Every outcome is a residual from a shape-and-location model that\n",
                         "never sees the axis gap, so a grip change that also altered drop or velocity separation is already accounted for and only the spin-match part remains. The\n",
                         "parachute claim predicts upward slopes in all three panels."),
       x = "Axis gap closed vs primary fastball, season over season (degrees)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12.5),
        plot.subtitle = element_text(size = 8.1), panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 9.8),
        panel.spacing.x = unit(15, "pt"))
ggsave(file.path(AST, "fig17_parachute_within_pitcher.png"), gg, width = 13.5, height = 5.6, dpi = 150)
cat("\nwrote fig17_parachute_within_pitcher.png + ext_parachute_within_pitcher.csv\n")
