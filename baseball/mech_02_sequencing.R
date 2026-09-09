#!/usr/bin/env Rscript

# M2: DOES THE ADVANTAGE CONCENTRATE ON CHANGEUPS THROWN AFTER A FOUR-SEAMER?
#
# This is the strongest of the three mechanism tests because the pitcher is held fixed. Within a
# single pitcher-season, some changeups follow a four-seamer and some do not. If the pitch works by
# being mistaken for the fastball, the hitter who just saw the fastball is the one who should be
# fooled, and the bin's advantage should be larger on those pitches. A pitcher-quality confound
# cannot produce that: "this arm is good" applies equally to both halves of his own season.
#
# The estimator is a difference-in-differences with pitcher-season fixed effects, so the bin main
# effect - the whole +1.87 - is absorbed and contributes nothing. Only the interaction identifies.
# Standard errors are clustered on the pitcher.
#
# MLB runs 2023-2026 here rather than 2020-2026, because the pitch stream with sequencing keys only
# reaches back to 2023. D1 covers its full three seasons. The unit is the pitch, not the season, so
# the loss is smaller than it sounds.

# bit64 is mandatory here: see the note at the top of spec_lock.R.
suppressPackageStartupMessages({ library(data.table); library(bit64); library(lightgbm) })
options(width = 200); set.seed(13)
MDIR <- "data/statcast_model"
SPEC <- readRDS(file.path(MDIR, "locked_spec.rds"))
L <- SPEC$data; setDT(L); L[, id := as.character(id)]

## ---- MLB: refit the whiff residual carrying the sequencing keys ---------------------------------
MO <- file.path(MDIR, "mlb_whiff_resid_seq.rds")
if (!file.exists(MO)) {
  X <- readRDS(file.path(MDIR, "parachute_rv.rds")); setDT(X)
  # The previous pitch is found by matching on the pitch one lower in the same plate appearance.
  # Only four pitch types live in this file, so a first-pitch changeup and a changeup after a
  # slider both fail to match. They are told apart by pitch_number: anything past the first pitch
  # that fails to match followed some pitch this file does not carry, which is a breaking ball.
  P <- X[, .(game_pk, at_bat_number, pn1 = pitch_number, prev_pt = pitch_type)]
  C <- X[pitch_type == "CH"][, pn1 := pitch_number - 1L]
  C <- merge(C, P, by = c("game_pk","at_bat_number","pn1"), all.x = TRUE)
  C[, has_prev := pitch_number > 1L]
  C[, prev_ff := has_prev & !is.na(prev_pt) & prev_pt == "FF"]

  C[, is_swing := description %in% c("swinging_strike","swinging_strike_blocked","foul",
                                     "foul_tip","hit_into_play","foul_bunt","missed_bunt")]
  C[, whiff := as.integer(description %in% c("swinging_strike","swinging_strike_blocked",
                                             "foul_tip","missed_bunt"))]
  C[, lh := p_throws == "L"]
  C[, `:=`(tj_ax = fifelse(lh, -ax, ax), tj_ax_diff = fifelse(lh, -ax_diff, ax_diff),
           tj_x0 = fifelse(lh, -release_pos_x, release_pos_x),
           tj_px = fifelse(lh, -plate_x, plate_x), tj_haa = fifelse(lh, -HAA, HAA),
           tj_sax = fifelse(lh, -sax, sax))]
  FEAT <- c("release_speed","release_spin_rate","release_extension","tj_ax","az","tj_x0",
            "release_pos_z","tj_sax","cax","speed_diff","tj_ax_diff","az_diff",
            "tj_px","plate_z","VAA","tj_haa","same_hand","balls","strikes")
  SW <- C[is_swing == TRUE & complete.cases(C[, ..FEAT])]
  cat(sprintf("MLB changeup swings 2023-2026: %s  whiff rate %.3f\n",
              format(nrow(SW), big.mark = ","), mean(SW$whiff)))
  K <- 4; fold <- sample(rep(1:K, length.out = nrow(SW))); p <- rep(NA_real_, nrow(SW))
  for (f in 1:K) {
    tr <- SW[fold != f]; n <- nrow(tr); vi <- sample(n, floor(.12*n))
    dtr <- lgb.Dataset(as.matrix(tr[-vi, ..FEAT]), label = tr$whiff[-vi])
    dva <- lgb.Dataset.create.valid(dtr, as.matrix(tr[vi, ..FEAT]), label = tr$whiff[vi])
    m <- lgb.train(params = list(objective = "binary", metric = "binary_logloss",
                   learning_rate = .06, num_leaves = 31, min_data_in_leaf = 300,
                   feature_fraction = .8, bagging_fraction = .8, bagging_freq = 1),
                   data = dtr, nrounds = 1500, valids = list(v = dva),
                   early_stopping_rounds = 50, verbose = -1)
    p[fold == f] <- predict(m, as.matrix(SW[fold == f, ..FEAT]))
  }
  SW[, r_all := whiff - p]
  cat(sprintf("  out-of-fold R2 = %.4f\n", 1 - var(SW$r_all)/var(SW$whiff)))
  saveRDS(SW[, .(pitcher, season, game_pk, at_bat_number, pitch_number, has_prev, prev_pt,
                 prev_ff, whiff, r_all)], MO)
} else cat("(using cached MLB sequenced residuals)\n")
MS <- readRDS(MO); setDT(MS)

# Sanity: does this refit land in the same place as the residual the locked spec was built on?
OLD <- readRDS(file.path(MDIR, "mlb_whiff_locaware.rds")); setDT(OLD)
a <- OLD[, .(old = 100*mean(r_all), n = .N), by = .(pitcher, season)][n >= 40]
b <- MS[, .(new = 100*mean(r_all), n2 = .N), by = .(pitcher, season)][n2 >= 40]
m <- merge(a, b, by = c("pitcher","season"))
cat(sprintf("  agreement with the locked-spec MLB residual: %d shared seasons, r = %.3f\n",
            nrow(m), cor(m$old, m$new)))

## ---- assemble the pooled pitch-level table ------------------------------------------------------
MS[, `:=`(league = "MLB", id = as.character(pitcher))]
DS <- readRDS(file.path(MDIR, "ncaa_whiff_resid_seq.rds")); setDT(DS)
DS[, `:=`(league = "D1", id = as.character(PitcherId), has_prev = !is.na(prev_pt),
          prev_ff = !is.na(prev_pt) & prev_pt == "Four-Seam")]
K <- rbind(MS[, .(league, id, season, has_prev, prev_ff, whiff, r_all)],
           DS[, .(league, id, season, has_prev, prev_ff, whiff, r_all)])
K <- merge(K, L[, .(league, id, season, bin)], by = c("league","id","season"))
K <- K[has_prev == TRUE & !is.na(id)]
K[, g := paste(league, id, season)]

cat(sprintf("\npooled changeup swings with a known predecessor: %s\n",
            format(nrow(K), big.mark = ",")))
print(K[, .(swings = .N, after_ff = sum(prev_ff), arms = uniqueN(id)),
        by = .(league, bin)][order(league, -bin)], row.names = FALSE)

## ---- the raw 2x2 before any modelling ------------------------------------------------------------
cat("\n=== residual whiff points, by bin and by what came before ===\n")
T2 <- K[, .(n = .N, r = round(100*mean(r_all),2)), by = .(league, bin, prev_ff)]
print(dcast(T2, league + bin ~ prev_ff, value.var = c("r","n")), row.names = FALSE)
cat("\n  the DiD is the (after-FF minus after-other) gap in the bin, less the same gap outside it\n")
for (lgv in c("MLB","D1")) {
  z <- K[league == lgv]
  gap <- function(b) { s <- z[bin == b]
    100*(mean(s[prev_ff == TRUE]$r_all) - mean(s[prev_ff == FALSE]$r_all)) }
  cat(sprintf("  %-4s  in-bin gap %+6.2f   out-of-bin gap %+6.2f   DiD %+6.2f\n",
              lgv, gap(TRUE), gap(FALSE), gap(TRUE) - gap(FALSE)))
}

## ---- fixed-effects DiD ---------------------------------------------------------------------------
# Everything is demeaned within pitcher-season, which is algebraically the same as putting a dummy
# on every pitcher-season but does not require inverting a 3,000-column design matrix. Because bin
# does not vary within a pitcher-season, its main effect vanishes under the transform - which is
# the point.
K[, `:=`(y = 100*r_all, x = as.numeric(prev_ff))]
K[, `:=`(yd = y - mean(y), xd = x - mean(x)), by = g]
K[, bx := as.numeric(bin) * xd]
K[, lgd := as.numeric(league == "MLB") * xd]

clus <- function(fit, cl) {
  u <- residuals(fit); Xm <- model.matrix(fit); bread <- solve(crossprod(Xm))
  meat <- crossprod(rowsum(Xm * u, cl))
  nc <- length(unique(cl)); k <- ncol(Xm)
  sqrt(diag(bread %*% meat %*% bread) * (nc/(nc-1)) * ((nrow(Xm)-1)/(nrow(Xm)-k)))
}
fit <- lm(yd ~ 0 + xd + lgd + bx, data = K)
se <- clus(fit, K$id); co <- coef(fit)
cat("\n=== M2: pooled DiD with pitcher-season fixed effects, SEs clustered on pitcher ===\n")
for (r in names(co)) {
  t <- co[r]/se[r]
  cat(sprintf("  %-5s %+7.3f  se %5.3f  t %+5.2f  p %.4f\n", r, co[r], se[r], t,
              2*pnorm(-abs(t)))) }
cat(sprintf("\n  bx is the estimate: %+.2f whiff points of extra advantage for bin changeups\n", co["bx"]))
cat(sprintf("  thrown after a four-seamer, over and above the same contrast in every other arm.\n"))
# A null is only informative if the test could have seen something. At 80 percent power a two-sided
# test needs roughly 2.8 standard errors, so this is the smallest interaction the design could have
# distinguished from zero - worth stating alongside the estimate.
cat(sprintf("  smallest DiD this design could detect at 80%% power: %.1f whiff points\n",
            2.8*se["bx"]))

cat("\n=== the same DiD run separately in each league ===\n")
for (lgv in c("MLB","D1")) {
  z <- K[league == lgv]
  f2 <- lm(yd ~ 0 + xd + bx, data = z); s2 <- clus(f2, z$id)
  cat(sprintf("  %-4s  n=%s (bin %s)  DiD %+6.2f  se %.2f  p %.3f\n", lgv,
              format(nrow(z), big.mark = ","), format(sum(z$bin), big.mark = ","),
              coef(f2)["bx"], s2["bx"], 2*pnorm(-abs(coef(f2)["bx"]/s2["bx"]))))
}
saveRDS(K, file.path(MDIR, "mech_sequencing.rds"))
cat("\nwrote mech_sequencing.rds\n")
