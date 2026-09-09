#!/usr/bin/env Rscript

# IS THE MISS-DISTANCE GAP SIGNIFICANT?
#
# The season-level t-test said +0.38 inches, p = 0.21, on 8 seasons. That test throws away the
# ~900 individual swings and treats a 187-swing season the same as a 47-swing one, so it is worth
# asking the question with the pitches themselves before calling it null.
#
# Three versions, in increasing order of how much I would trust them:
#
#   raw pitch level      every swing, but the "n" is fake - swings from one arm are not
#                        independent draws, so the naive p here is far too small
#   residual             miss distance above a model with the same features the whiff model uses,
#                        including location. This is the apples-to-apples analog of the locked
#                        whiff residual rather than a raw rate
#   arm level            one number per pitcher, and a permutation that shuffles the recipe label
#                        across arms. With 7 arms this is the only inference that is honest
#
# miss_grade_data carries no plate location, so it is joined to oof_whiff_resid on the pitch keys
# to recover plate_x/plate_z/VAA/HAA. Without location in the model the residual would just be
# rediscovering that these pitchers bury the changeup.

suppressPackageStartupMessages({ library(data.table); library(bit64); library(lightgbm) })
options(width = 200); set.seed(41)
MDIR <- "data/statcast_model"
OUT <- file.path(MDIR, "miss_resid_ch.rds")

L <- readRDS(file.path(MDIR, "locked_spec.rds"))$data; setDT(L)
L <- L[league == "MLB", .(id, season, bin, nsw, axis, arm, ec, ef)]

if (!file.exists(OUT)) {
  MG <- readRDS(file.path(MDIR, "miss_grade_data.rds")); setDT(MG)
  MG <- MG[pitch_type == "CH" & is.finite(miss_distance)]
  LOC <- readRDS(file.path(MDIR, "oof_whiff_resid.rds")); setDT(LOC)
  MG <- merge(MG, LOC[, .(game_pk, at_bat_number, pitch_number, plate_x, plate_z, VAA, HAA)],
              by = c("game_pk","at_bat_number","pitch_number"))
  TH <- unique(readRDS(file.path(MDIR, "parachute_ff.rds"))[, .(pitcher, season, p_throws, stand)])
  setDT(TH)
  MG <- merge(MG, unique(TH[, .(pitcher, season, p_throws)]), by = c("pitcher","season"))
  MG[, lh := p_throws == "L"]
  MG[, `:=`(tj_ax = fifelse(lh, -ax, ax), tj_ax_diff = fifelse(lh, -ax_diff, ax_diff),
            tj_x0 = fifelse(lh, -release_pos_x, release_pos_x),
            tj_px = fifelse(lh, -plate_x, plate_x), tj_haa = fifelse(lh, -HAA, HAA),
            tj_sax = fifelse(lh, -sax, sax))]
  FEAT <- c("release_speed","release_spin_rate","release_extension","tj_ax","az","tj_x0",
            "release_pos_z","tj_sax","cax","speed_diff","tj_ax_diff","az_diff",
            "tj_px","plate_z","VAA","tj_haa","balls","strikes")
  MG <- MG[complete.cases(MG[, ..FEAT])]
  cat(sprintf("changeup swings with miss distance and location, 2023-2026: %s\n",
              format(nrow(MG), big.mark = ",")))
  K <- 4; fold <- sample(rep(1:K, length.out = nrow(MG))); p <- rep(NA_real_, nrow(MG))
  for (f in 1:K) {
    tr <- MG[fold != f]; n <- nrow(tr); vi <- sample(n, floor(.12*n))
    dtr <- lgb.Dataset(as.matrix(tr[-vi, ..FEAT]), label = tr$miss_distance[-vi])
    dva <- lgb.Dataset.create.valid(dtr, as.matrix(tr[vi, ..FEAT]), label = tr$miss_distance[vi])
    m <- lgb.train(params = list(objective = "regression", metric = "l2",
                   learning_rate = .06, num_leaves = 31, min_data_in_leaf = 300,
                   feature_fraction = .8, bagging_fraction = .8, bagging_freq = 1),
                   data = dtr, nrounds = 1500, valids = list(v = dva),
                   early_stopping_rounds = 50, verbose = -1)
    p[fold == f] <- predict(m, as.matrix(MG[fold == f, ..FEAT]))
  }
  MG[, mres := miss_distance - p]
  cat(sprintf("  out-of-fold R2 on miss distance = %.4f\n",
              1 - var(MG$mres)/var(MG$miss_distance)))
  saveRDS(MG[, .(id = as.character(pitcher), season, miss_distance, mres, is_whiff)], OUT)
} else cat("(using cached miss residuals)\n")

M <- readRDS(OUT); setDT(M)
M <- merge(M, L, by = c("id","season"))
cat(sprintf("\njoined to the locked MLB pool: %s swings, %d seasons, recipe %d seasons / %d arms\n",
            format(nrow(M), big.mark = ","), uniqueN(M[, .(id, season)]),
            uniqueN(M[bin == TRUE, .(id, season)]), uniqueN(M[bin == TRUE]$id)))

## ---- 1. raw pitch level (the p-value that lies) -------------------------------------------------
A <- M[bin == TRUE]; B <- M[bin == FALSE]
t1 <- t.test(A$miss_distance, B$miss_distance)
t2 <- t.test(A$mres, B$mres)
cat("\n=== pitch level, treating every swing as independent (overstates significance) ===\n")
cat(sprintf("  raw miss   bin %.2f in  out %.2f in  gap %+.2f  p %.4f  (n = %s vs %s)\n",
            t1$estimate[1], t1$estimate[2], diff(rev(t1$estimate)), t1$p.value,
            format(nrow(A), big.mark = ","), format(nrow(B), big.mark = ",")))
cat(sprintf("  residual   bin %+.2f  out %+.2f  gap %+.2f  p %.4f\n",
            t2$estimate[1], t2$estimate[2], diff(rev(t2$estimate)), t2$p.value))

## ---- 2. cluster-robust on the pitcher ------------------------------------------------------------
clus_se <- function(y, x, cl) {
  fit <- lm(y ~ x); u <- residuals(fit); Xm <- model.matrix(fit)
  bread <- solve(crossprod(Xm)); meat <- crossprod(rowsum(Xm * u, cl))
  nc <- length(unique(cl))
  se <- sqrt(diag(bread %*% meat %*% bread) * (nc/(nc-1)))
  c(est = coef(fit)[2], se = se[2], nc = nc)
}
cr <- clus_se(M$mres, as.numeric(M$bin), M$id)
cat("\n=== pitch level, standard errors clustered on the pitcher ===\n")
cat(sprintf("  residual gap %+.2f  se %.2f  t %+.2f  p %.3f  (%d clusters)\n",
            cr[1], cr[2], cr[1]/cr[2], 2*pnorm(-abs(cr[1]/cr[2])), cr[3]))
cat("  with 7 arms inside the bin, even this understates the uncertainty.\n")

## ---- 3. arm level, and a permutation that respects it --------------------------------------------
ARM <- M[, .(miss = mean(miss_distance), mres = mean(mres), n = .N), by = .(id, bin)]
a <- ARM[bin == TRUE]; b <- ARM[bin == FALSE]
cat("\n=== arm level, one number per pitcher ===\n")
cat(sprintf("  raw miss   bin %.2f (%d arms)  out %.2f (%d arms)  gap %+.2f  p %.3f\n",
            mean(a$miss), nrow(a), mean(b$miss), nrow(b), mean(a$miss) - mean(b$miss),
            t.test(a$miss, b$miss)$p.value))
cat(sprintf("  residual   bin %+.2f  out %+.2f  gap %+.2f  p %.3f\n",
            mean(a$mres), mean(b$mres), mean(a$mres) - mean(b$mres),
            t.test(a$mres, b$mres)$p.value))

obs <- mean(a$mres) - mean(b$mres)
k <- nrow(a); pool <- ARM$mres
perm <- replicate(20000, { i <- sample(length(pool), k)
  mean(pool[i]) - mean(pool[-i]) })
cat(sprintf("\n  permutation, 20000 reshuffles of which %d arms carry the label:\n", k))
cat(sprintf("  observed %+.2f | null mean %+.2f sd %.2f | two-sided p = %.4f\n",
            obs, mean(perm), sd(perm), mean(abs(perm) >= abs(obs))))

## ---- how much would it take? ---------------------------------------------------------------------
sd_arm <- sd(ARM$mres)
need <- ceiling((2.8 * sd_arm / obs)^2)
cat(sprintf("\n  arm-level residual sd = %.2f inches. To call a %+.2f gap significant at 80%% power\n",
            sd_arm, obs))
cat(sprintf("  this design needs roughly %d bin arms; it has %d.\n", need, k))
