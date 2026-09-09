#!/usr/bin/env Rscript

# CHASE AS THE TARGET INSTEAD OF WHIFF.
#
# A parachute changeup is supposed to work by looking like the fastball out of the hand, which is a
# claim about the hitter's SWING DECISION rather than about his ability to hit the ball once he has
# committed. Whiff conflates the two: a pitch can be correctly identified and still missed because
# it is hard to square up, and it can be badly misread yet fouled off. Chase - swinging at a pitch
# that is not a strike - isolates the decision, so it is the closer measurement of being fooled.
#
# It is also a far better-powered test here. There are 470k out-of-zone D1 changeups against 195k
# swings, and the per-pitcher denominators are correspondingly larger, which is what the whiff
# analysis kept running out of.
#
# The critical control is location. Chase rate is mostly a function of how far out of the zone the
# pitch finished, so a model without location would rank command rather than deception, and would
# reward pitchers who miss narrowly. Location therefore goes in as a feature and the residual
# answers the question that matters: given a changeup of this shape finishing in this exact spot,
# did hitters offer at it more than they should have?

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(11); options(width = 205)
MDIR <- "data/statcast_model"
CACHE <- file.path(MDIR, "ncaa_chase_resid.rds")

D <- readRDS(file.path(MDIR, "ncaa_d1_pitches.rds"))
D[, L := PitcherThrows == "Left"]
FF <- D[pt == "Four-Seam" & is.finite(RelSpeed),
        .(nff = .N, ff_speed = mean(RelSpeed), ff_ax = mean(ax, na.rm = TRUE),
          ff_az = mean(az, na.rm = TRUE), ff_axis = mean(SpinAxis, na.rm = TRUE)),
        by = .(PitcherId, season)][nff >= 30]
C <- merge(D[pt == "Changeup"], FF, by = c("PitcherId","season"))
C <- C[is.finite(RelSpeed) & is.finite(ax) & is.finite(az) & is.finite(SpinAxis) &
       is.finite(px) & is.finite(pz) & is.finite(SpinRate) & is.finite(Extension) &
       is.finite(VertApprAngle) & is.finite(HorzApprAngle)]
C[, `:=`(speed_diff = RelSpeed - ff_speed, ax_diff = ax - ff_ax, az_diff = az - ff_az,
         axis_diff = pmin(abs(SpinAxis - ff_axis), 360 - abs(SpinAxis - ff_axis)))]
C[, `:=`(tj_x0 = fifelse(L, -RelSide, RelSide), tj_ax = fifelse(L, -ax, ax),
         tj_ax_diff = fifelse(L, -ax_diff, ax_diff), tj_px = fifelse(L, -px, px),
         tj_haa = fifelse(L, -HorzApprAngle, HorzApprAngle),
         tj_axis = fifelse(L, (360 - SpinAxis) %% 360, SpinAxis),
         same_hand = as.integer((PitcherThrows == "Left") == (BatterSide == "Left")))]

# Zone edges are the conventional ones. Distance outside is given to the model explicitly rather
# than left to be inferred from px and pz, because it is the single strongest driver of chase and
# the residual is only trustworthy if it is absorbed properly.
C[, ooz := abs(px) > 0.83 | pz < 1.5 | pz > 3.5]
C[, `:=`(dx = pmax(abs(px) - 0.83, 0), dz = pmax(pmax(1.5 - pz, pz - 3.5), 0))]
C[, dist_out := sqrt(dx^2 + dz^2)]
O <- C[ooz == TRUE]
cat(sprintf("out-of-zone D1 changeups: %s from %d pitcher-seasons, chase rate %.3f\n",
            format(nrow(O), big.mark = ","), uniqueN(paste(O$PitcherId, O$season)), mean(O$swing)))

TJ  <- c("RelSpeed","SpinRate","Extension","tj_ax","az","tj_x0","RelHeight","tj_axis",
         "speed_diff","tj_ax_diff","az_diff")
LOC <- c("tj_px","pz","dist_out","dx","dz","VertApprAngle","tj_haa","same_hand","Balls","Strikes")

if (!file.exists(CACHE)) {
  fit <- function(FEAT, tag) {
    K <- 4; fold <- sample(rep(1:K, length.out = nrow(O))); p <- rep(NA_real_, nrow(O))
    for (f in 1:K) {
      tr <- O[fold != f]; n <- nrow(tr); vi <- sample(n, floor(.12*n))
      dtr <- lgb.Dataset(as.matrix(tr[-vi, ..FEAT]), label = tr$swing[-vi])
      dva <- lgb.Dataset.create.valid(dtr, as.matrix(tr[vi, ..FEAT]), label = tr$swing[vi])
      m <- lgb.train(params = list(objective = "binary", metric = "binary_logloss",
                     learning_rate = .06, num_leaves = 31, min_data_in_leaf = 300,
                     feature_fraction = .8, bagging_fraction = .8, bagging_freq = 1),
                     data = dtr, nrounds = 1500, valids = list(v = dva),
                     early_stopping_rounds = 50, verbose = -1)
      p[fold == f] <- predict(m, as.matrix(O[fold == f, ..FEAT]))
    }
    cat(sprintf("  %-30s %2d features  R2=%.4f\n", tag, length(FEAT),
                1 - var(O$swing - p)/var(O$swing)))
    O$swing - p
  }
  cat("\n=== out-of-fold chase models on out-of-zone D1 changeups ===\n")
  O[, r_tj  := fit(TJ, "tjStuff+ v3.0 (shape only)")]
  O[, r_all := fit(c(TJ, LOC), "tjStuff+ v3.0 + location")]
  saveRDS(O[, .(PitcherId, Pitcher, season, axis_diff, speed_diff, r_tj, r_all)], CACHE)
} else cat("\n(using cached chase fits)\n")
R <- readRDS(CACHE); setDT(R)

## ---- test the bins ---------------------------------------------------------------------------
ARM <- as.data.table(readRDS(file.path(MDIR, "ncaa_armangle.rds")))
PR  <- as.data.table(readRDS(file.path(MDIR, "ncaa_spineff_pairs.rds")))
PR  <- merge(PR, ARM[, .(PitcherId, season, arm_hat)], by = c("PitcherId","season"))
cat(sprintf("\n%d pitcher-seasons carry an efficiency pair and an arm estimate\n", nrow(PR)))

# r_all is the headline: location is controlled, so what is left is shape-driven deception.
for (EFF in c(.85, .90)) {
  cat(sprintf("\n=== efficiency floor %.2f, target = chase above model (location controlled) ===\n", EFF))
  res <- rbindlist(lapply(c(15, 25, 40, 60), function(g) {
    S <- R[!is.na(PitcherId), .(c_all = 100*mean(r_all, na.rm = TRUE),
                                c_tj = 100*mean(r_tj, na.rm = TRUE), n = .N),
           by = .(PitcherId, season)][n >= g]
    Dd <- merge(PR, S, by = c("PitcherId","season"))
    Dd[, `:=`(bin  = eff_ch >= EFF & eff_ff >= EFF & arm_hat >= 44,
              bina = eff_ch >= EFF & eff_ff >= EFF & arm_hat >= 44 & axis_gap <= 10,
              binx = eff_ch >= EFF & eff_ff >= EFF & axis_gap <= 10)]
    row <- function(lab, i) {
      x <- Dd$c_all[i]; y <- Dd$c_all[!i]
      if (length(x) < 3 || length(y) < 3)
        return(data.table(gate = g, cut = lab, n = length(x), bin = NA_real_, rest = NA_real_,
                          diff = NA_real_, p = NA_real_, se = NA_real_))
      t <- t.test(x, y)
      data.table(gate = g, cut = lab, n = length(x), bin = round(mean(x),2),
                 rest = round(mean(y),2), diff = round(mean(x)-mean(y),2),
                 p = round(t$p.value,4), se = round(sd(x)/sqrt(length(x)),2))
    }
    rbind(row("eff + arm 44", Dd$bin),
          row("eff + arm 44 + axis<=10", Dd$bina),
          row("eff + axis<=10 (no arm)", Dd$binx))
  }))
  print(res, row.names = FALSE)
}

cat("\n=== is chase a more stable measure than whiff at the pitcher level? ===\n")
# Split-half reliability decides how much of any bin difference could be real. A measure that does
# not correlate with itself across a random split of a pitcher's own pitches cannot correlate with
# anything else either, so this bounds every result above.
SP <- R[!is.na(PitcherId)][, half := sample(rep_len(1:2, .N)), by = .(PitcherId, season)]
H <- dcast(SP[, .(v = mean(r_all), n = .N), by = .(PitcherId, season, half)],
           PitcherId + season ~ half, value.var = c("v","n"))
H <- H[n_1 >= 50 & n_2 >= 50]
cat(sprintf("  chase residual, split-half r = %+.3f over %d pitcher-seasons (100+ ooz changeups)\n",
            cor(H$v_1, H$v_2), nrow(H)))
W <- readRDS(file.path(MDIR, "ncaa_whiff_resid.rds")); setDT(W)
SW <- W[!is.na(PitcherId)][, half := sample(rep_len(1:2, .N)), by = .(PitcherId, season)]
HW <- dcast(SW[, .(v = mean(r_tj), n = .N), by = .(PitcherId, season, half)],
            PitcherId + season ~ half, value.var = c("v","n"))
HW <- HW[n_1 >= 20 & n_2 >= 20]
cat(sprintf("  whiff residual, split-half r = %+.3f over %d pitcher-seasons (40+ swings)\n",
            cor(HW$v_1, HW$v_2), nrow(HW)))
