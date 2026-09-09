#!/usr/bin/env Rscript

# THE CHASE TEST, RUN ON MLB WHERE EVERY BIN CONDITION IS MEASURED RATHER THAN INFERRED.
#
# D1 produced a positive chase effect confined to the small axis-gap subset, growing from +1.8 to
# +4.1 points as the residual gate tightened, on 8 to 19 pitcher-seasons. That is the right shape
# but far too little of it, and three of the four bin conditions there are substitutes: efficiency
# is inferred from the lift law, arm slot is a height-informed estimate, and the spin axis is
# movement-derived rather than measured.
#
# MLB has all four directly. Active spin comes from the Savant leaderboards, arm angle is tracked by
# Hawk-Eye, and the spin axis is measured off the ball's rotation instead of reconstructed from the
# break. So this is the clean version of the same test, and the honest one: if the effect is real it
# should be clearer here, and if it is absent here the D1 cell was noise.
#
# The zone is the batter-specific one via z_rel_bot / z_rel_top rather than a fixed band, and the
# distance outside it is handed to the model explicitly - chase is mostly a function of how far the
# pitch missed, and the residual only means anything once that is properly absorbed.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(11); options(width = 205)
MDIR <- "data/statcast_model"
CACHE <- file.path(MDIR, "mlb_chase_resid.rds")
GATE <- { a <- commandArgs(TRUE); if (length(a)) as.integer(a[1]) else 40L }

F <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(F)
F <- F[pitch_type == "CH" & is.finite(plate_x) & is.finite(plate_z) & is.finite(z_rel_bot) &
       is.finite(z_rel_top) & is.finite(release_speed) & is.finite(ax) & is.finite(az) &
       is.finite(speed_diff) & is.finite(is_swing)]
L <- "p_throws" %in% names(F)
F[, lh := if (L) p_throws == "L" else throws_R == 0]

# Outside the zone horizontally, or below the bottom, or above the top. z_rel_bot and z_rel_top are
# signed distances to those edges, so their orientation is checked rather than assumed.
cat(sprintf("z_rel_bot median %.2f, z_rel_top median %.2f\n",
            median(F$z_rel_bot), median(F$z_rel_top)))
F[, `:=`(dx = pmax(abs(plate_x) - 0.83, 0),
         dz = pmax(pmax(-z_rel_bot, 0), pmax(z_rel_top, 0)))]
F[, ooz := dx > 0 | dz > 0]
F[, dist_out := sqrt(dx^2 + dz^2)]
O <- F[ooz == TRUE]
cat(sprintf("out-of-zone MLB changeups: %s from %d pitcher-seasons, chase rate %.3f\n",
            format(nrow(O), big.mark = ","), uniqueN(paste(O$pitcher, O$season)), mean(O$is_swing)))

O[, `:=`(tj_x0 = fifelse(lh, -release_pos_x, release_pos_x),
         tj_ax = fifelse(lh, -ax, ax),
         tj_ax_diff = fifelse(lh, -ax_diff, ax_diff),
         tj_px = fifelse(lh, -plate_x, plate_x),
         tj_haa = fifelse(lh, -HAA, HAA),
         tj_sax = fifelse(lh, -sax, sax))]
# sax/cax are the sine and cosine of the changeup's own spin axis. They belong in the model, and
# leaving them out would rig the later test: the bin condition is an axis GAP, and a model blind to
# axis entirely would leave axis-shaped signal in the residual for trivial reasons. Encoding the
# angle as a sine/cosine pair rather than degrees avoids the wrap at 360.
TJ  <- c("release_speed","release_spin_rate","release_extension","tj_ax","az","tj_x0",
         "release_pos_z","speed_diff","tj_ax_diff","az_diff","tj_sax","cax")
LOC <- c("tj_px","plate_z","dist_out","dx","dz","VAA","tj_haa","same_hand")
for (v in c("balls","strikes")) if (v %in% names(O)) LOC <- c(LOC, v)

if (!file.exists(CACHE)) {
  fitm <- function(FEAT, tag) {
    K <- 4; fold <- sample(rep(1:K, length.out = nrow(O))); p <- rep(NA_real_, nrow(O))
    for (f in 1:K) {
      tr <- O[fold != f]; n <- nrow(tr); vi <- sample(n, floor(.12*n))
      dtr <- lgb.Dataset(as.matrix(tr[-vi, ..FEAT]), label = tr$is_swing[-vi])
      dva <- lgb.Dataset.create.valid(dtr, as.matrix(tr[vi, ..FEAT]), label = tr$is_swing[vi])
      m <- lgb.train(params = list(objective = "binary", metric = "binary_logloss",
                     learning_rate = .06, num_leaves = 31, min_data_in_leaf = 300,
                     feature_fraction = .8, bagging_fraction = .8, bagging_freq = 1),
                     data = dtr, nrounds = 1500, valids = list(v = dva),
                     early_stopping_rounds = 50, verbose = -1)
      p[fold == f] <- predict(m, as.matrix(O[fold == f, ..FEAT]))
    }
    cat(sprintf("  %-30s %2d features  R2=%.4f\n", tag, length(FEAT),
                1 - var(O$is_swing - p)/var(O$is_swing)))
    O$is_swing - p
  }
  cat("\n=== out-of-fold chase models on out-of-zone MLB changeups ===\n")
  O[, r_tj  := fitm(TJ, "tjStuff+ (shape only)")]
  O[, r_all := fitm(c(TJ, LOC), "tjStuff+ + location")]
  saveRDS(O[, .(pitcher, season, r_tj, r_all)], CACHE)
} else cat("\n(using cached MLB chase fits)\n")
R <- readRDS(CACHE); setDT(R)

## ---- the bin, every condition measured --------------------------------------------------------
AS <- readRDS(file.path(MDIR, "active_spin_long.rds")); setDT(AS)
S <- F[, .(axis = mean(axis_diff, na.rm = TRUE), arm = mean(arm_angle, na.rm = TRUE)),
       by = .(pitcher, season)]
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, as_ch = active_spin)], by = c("pitcher","season"))
S <- merge(S, AS[pitch_type == "FF", .(pitcher, season, as_fb = active_spin)], by = c("pitcher","season"))
CH <- R[, .(c_all = 100*mean(r_all), c_tj = 100*mean(r_tj), nooz = .N), by = .(pitcher, season)]
S <- merge(S, CH, by = c("pitcher","season"))[nooz >= GATE & is.finite(arm)]
cat(sprintf("\ngate: %d+ out-of-zone changeups -> %d pitcher-seasons with all four conditions measured\n",
            GATE, nrow(S)))

row <- function(lab, i) {
  x <- S$c_all[i]; y <- S$c_all[!i]
  if (length(x) < 3 || length(y) < 3)
    return(data.table(cut = lab, n = length(x), bin = NA_real_, rest = NA_real_,
                      diff = NA_real_, p = NA_real_, se = NA_real_))
  t <- t.test(x, y)
  data.table(cut = lab, n = length(x), bin = round(mean(x),2), rest = round(mean(y),2),
             diff = round(mean(x)-mean(y),2), p = round(t$p.value,4),
             se = round(sd(x)/sqrt(length(x)),2))
}
for (EFF in c(.85, .90)) {
  cat(sprintf("\n=== MLB, efficiency floor %.2f, target = chase above model (location controlled) ===\n", EFF))
  print(rbind(
    row("eff only",                S[, as_ch >= EFF & as_fb >= EFF]),
    row("eff + arm 44",            S[, as_ch >= EFF & as_fb >= EFF & arm >= 44]),
    row("eff + arm 44 + axis<=10", S[, as_ch >= EFF & as_fb >= EFF & arm >= 44 & axis <= 10]),
    row("eff + axis<=10 (no arm)", S[, as_ch >= EFF & as_fb >= EFF & axis <= 10]),
    row("axis<=10 only",           S[, axis <= 10])), row.names = FALSE)
}

cat("\n=== who is in the .85 + arm + axis bin ===\n")
NM <- unique(readRDS(file.path(MDIR, "parachute_rv.rds"))[, .(pitcher, player_name)])
B <- merge(S[as_ch >= .85 & as_fb >= .85 & arm >= 44 & axis <= 10], NM, by = "pitcher", all.x = TRUE)
if (nrow(B)) print(B[order(-c_all), .(
  Pitcher = trimws(paste(sub(".*,\\s*","",player_name), sub(",.*","",player_name))),
  Season = season, OOZ_CH = nooz, ArmAngle = round(arm,1), AxisGap = round(axis,1),
  CH_active = round(as_ch,2), FF_active = round(as_fb,2),
  Chase_above = round(c_all,2))], row.names = FALSE)
