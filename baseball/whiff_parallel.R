#!/usr/bin/env Rscript

# THE SAME TEST, RUN ON WHIFFS.
#
# Chase asks whether the hitter was fooled into offering. Whiff asks whether he was beaten once he
# had committed. They are different failures and there is no reason a look-alike cue has to produce
# both, so this is a genuine second look rather than a re-slicing of the first.
#
# D1's whiff residual already exists and is built the right way - location controlled and with the
# pitch's own spin axis among the features. MLB's cached whiff residuals are shape-only, so a
# matching model is fit here instead of reusing them. Comparing a location-controlled residual in
# one league against a shape-only one in the other would guarantee a difference that had nothing to
# do with baseball.
#
# The comparison is deliberately the SAME set of cells as the chase run, decided in advance, so this
# is not another pass through the grid looking for a survivor.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(11); options(width = 205)
MDIR <- "data/statcast_model"; CACHE <- file.path(MDIR, "mlb_whiff_locaware.rds")
GATE <- 40L; EFF <- .85

F <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(F)
F <- F[pitch_type == "CH" & is.finite(plate_x) & is.finite(plate_z) & is.finite(release_speed) &
       is.finite(ax) & is.finite(az) & is.finite(speed_diff) & is.finite(sax) & is.finite(cax)]
F[, lh := if ("p_throws" %in% names(F)) p_throws == "L" else throws_R == 0]
W <- F[is_swing == 1 & is.finite(whiff)]
cat(sprintf("MLB changeup swings: %s from %d pitcher-seasons, whiff rate %.3f\n",
            format(nrow(W), big.mark = ","), uniqueN(paste(W$pitcher, W$season)), mean(W$whiff)))

W[, `:=`(tj_x0 = fifelse(lh, -release_pos_x, release_pos_x), tj_ax = fifelse(lh, -ax, ax),
         tj_ax_diff = fifelse(lh, -ax_diff, ax_diff), tj_px = fifelse(lh, -plate_x, plate_x),
         tj_haa = fifelse(lh, -HAA, HAA), tj_sax = fifelse(lh, -sax, sax))]
TJ  <- c("release_speed","release_spin_rate","release_extension","tj_ax","az","tj_x0",
         "release_pos_z","speed_diff","tj_ax_diff","az_diff","tj_sax","cax")
LOC <- c("tj_px","plate_z","VAA","tj_haa","same_hand","z_rel_bot","z_rel_top")
for (v in c("balls","strikes")) if (v %in% names(W)) LOC <- c(LOC, v)

if (!file.exists(CACHE)) {
  FEAT <- c(TJ, LOC); K <- 4
  fold <- sample(rep(1:K, length.out = nrow(W))); p <- rep(NA_real_, nrow(W))
  for (f in 1:K) {
    tr <- W[fold != f]; n <- nrow(tr); vi <- sample(n, floor(.12*n))
    dtr <- lgb.Dataset(as.matrix(tr[-vi, ..FEAT]), label = tr$whiff[-vi])
    dva <- lgb.Dataset.create.valid(dtr, as.matrix(tr[vi, ..FEAT]), label = tr$whiff[vi])
    m <- lgb.train(params = list(objective = "binary", metric = "binary_logloss",
                   learning_rate = .06, num_leaves = 31, min_data_in_leaf = 300,
                   feature_fraction = .8, bagging_fraction = .8, bagging_freq = 1),
                   data = dtr, nrounds = 1500, valids = list(v = dva),
                   early_stopping_rounds = 50, verbose = -1)
    p[fold == f] <- predict(m, as.matrix(W[fold == f, ..FEAT]))
  }
  cat(sprintf("  whiff model, %d features, out-of-fold R2 = %.4f\n", length(FEAT),
              1 - var(W$whiff - p)/var(W$whiff)))
  saveRDS(W[, .(pitcher, season, r_all = whiff - p)], CACHE)
} else cat("(using cached MLB whiff fits)\n")
RW <- readRDS(CACHE); setDT(RW)

## ---- assemble both leagues ---------------------------------------------------------------------
AS <- readRDS(file.path(MDIR, "active_spin_long.rds")); setDT(AS)
S <- F[, .(axis = mean(axis_diff, na.rm = TRUE), arm = mean(arm_angle, na.rm = TRUE),
           vs = -mean(speed_diff, na.rm = TRUE)), by = .(pitcher, season)]
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, ec = active_spin)], by = c("pitcher","season"))
S <- merge(S, AS[pitch_type == "FF", .(pitcher, season, ef = active_spin)], by = c("pitcher","season"))
S <- merge(S, RW[, .(y = 100*mean(r_all), n = .N), by = .(pitcher, season)], by = c("pitcher","season"))
S <- S[n >= GATE & is.finite(arm), .(id = pitcher, y, axis, arm, ec, ef, vs)]

ARM <- as.data.table(readRDS(file.path(MDIR, "ncaa_armangle.rds")))
PR  <- as.data.table(readRDS(file.path(MDIR, "ncaa_spineff_pairs.rds")))
PR  <- merge(PR, ARM[, .(PitcherId, season, arm_hat)], by = c("PitcherId","season"))
setnames(PR, "axis_gap", "axis", skip_absent = TRUE)
DW  <- readRDS(file.path(MDIR, "ncaa_whiff_resid.rds")); setDT(DW)
D <- merge(PR, DW[!is.na(PitcherId), .(y = 100*mean(r_all, na.rm = TRUE), n = .N),
                  by = .(PitcherId, season)], by = c("PitcherId","season"))[n >= GATE]
D <- D[, .(id = PitcherId, y, axis, arm = arm_hat, ec = eff_ch, ef = eff_ff, vs = velo_sep)]
cat(sprintf("\ngate %d swings: D1 pool %d pitcher-seasons, MLB pool %d\n", GATE, nrow(D), nrow(S)))

## ---- the pre-specified cells --------------------------------------------------------------------
cat("\n=== whiff above model, same cells as the chase run ===\n")
for (ax in c(10, 15, 20)) for (sl in c(FALSE, TRUE)) {
  o <- sapply(list(D, S), function(X) { thr <- as.numeric(quantile(X$arm, 2/3))
    i <- X[, ec >= EFF & ef >= EFF & axis <= ax & (if (sl) arm >= thr else TRUE)]
    if (sum(i) < 4) return(c(sum(i), NA, NA))
    c(sum(i), mean(X$y[i]) - mean(X$y[!i]), t.test(X$y[i], X$y[!i])$p.value) })
  cat(sprintf("  axis<=%2d slot=%-9s  D1: n=%3d %+.2f (p=%.3f)   MLB: n=%3d %+.2f (p=%.3f)\n",
              ax, ifelse(sl, "top third", "any"), o[1,1], o[2,1], o[3,1], o[1,2], o[2,2], o[3,2])) }

cat("\n=== continuous: axis match against whiff above model, arm and efficiency partialled out ===\n")
for (lg in c("D1","MLB")) { X <- if (lg == "D1") D else S
  r0 <- cor.test(-X$axis, X$y)
  ra <- resid(lm(I(-axis) ~ arm + ec + ef, X)); ry <- resid(lm(y ~ arm + ec + ef, X))
  rp <- cor.test(ra, ry)
  cat(sprintf("  %-3s n=%d | raw r=%+.3f (p=%.3f) | partial r=%+.3f (p=%.3f)\n",
              lg, nrow(X), r0$estimate, r0$p.value, rp$estimate, rp$p.value)) }

cat("\n=== does velo separation interact with bin membership on whiffs? ===\n")
for (ax in c(10, 15, 20)) for (lg in c("D1","MLB")) { X <- copy(if (lg == "D1") D else S)
  X[, inb := ec >= EFF & ef >= EFF & axis <= ax]
  m <- summary(lm(y ~ inb*scale(vs), X))$coefficients
  cat(sprintf("  %-3s axis<=%2d  bin %+.2f (p=%.3f) | bin x velo %+.2f (p=%.3f)\n", lg, ax,
              m["inbTRUE","Estimate"], m["inbTRUE","Pr(>|t|)"],
              m[nrow(m),"Estimate"], m[nrow(m),"Pr(>|t|)"])) }

cat("\n=== disjoint axis bands, arm slot and efficiency held ===\n")
for (lg in c("D1","MLB")) { X <- if (lg == "D1") D else S
  thr <- as.numeric(quantile(X$arm, 2/3)); P <- X[ec >= EFF & ef >= EFF & arm >= thr]
  cat(sprintf("  %s (n = %d):\n", lg, nrow(P)))
  print(P[, .(n = .N, whiff_above = round(mean(y),2), se = round(sd(y)/sqrt(.N),2)),
          by = .(band = cut(axis, c(-1,10,20,30,45,360), c("0-10","10-20","20-30","30-45","45+")))
         ][order(band)], row.names = FALSE) }
