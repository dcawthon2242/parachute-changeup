#!/usr/bin/env Rscript

# THE AXIS x ARM-SLOT INTERACTION, RE-TESTED UNDER THE FROZEN SPEC.
#
# parachute_slot_audit.R found the one result in this project that survived its own audit: the
# spin-axis gap only tracks whiffs for pitchers who throw from over the top. Interaction p = .002.
# That number was produced before three things changed:
#
#   1. the four-seam-primary filter, which removed sinker-primary arms for whom the four-seam
#      anchor was the wrong reference in the first place
#   2. the season-level velocity adjustment, which strips the velo-separation channel out of the
#      residual so a slot effect cannot ride on velocity
#   3. arm-clustered inference. p = .002 counted Kikuchi's four seasons as four observations
#
# It was also MLB-only. D1 has an arm-slot estimate and an axis gap, so the interaction is a real
# out-of-sample replication that this project never actually ran - only the weaker question of
# whether slot adds anything INSIDE the locked bin, which is a different test and came back flat.
#
# Both leagues are z-scored within league before interacting. D1's arm angle is a regression
# estimate with 68% of MLB's spread, so a raw-degree coefficient is not comparable across them;
# in z units the coefficient means "change in the axis slope per 1 SD of slot" in both.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
set.seed(19); options(width = 200)
MDIR <- "data/statcast_model"

L <- readRDS(file.path(MDIR, "locked_spec.rds")); P <- copy(L$data); setDT(P)
P[, `:=`(za = as.numeric(scale(axis)), zr = as.numeric(scale(arm)),
         zv = as.numeric(scale(vs)), ze = as.numeric(scale(ec))), by = league]

# Cluster-robust on the pitcher. One arm is one cluster no matter how many seasons it contributes.
crob <- function(f, D, key) {
  D <- D[complete.cases(D[, c(all.vars(f), key), with = FALSE])]
  m <- lm(f, data = D); u <- residuals(m); X <- model.matrix(m)
  cl <- D[[key]]
  nc <- uniqueN(cl); b <- solve(crossprod(X))
  V <- b %*% crossprod(rowsum(X * u, cl)) %*% b * (nc/(nc-1)) * ((nrow(X)-1)/(nrow(X)-ncol(X)))
  k <- grep(":", names(coef(m)), value = TRUE)[1]
  est <- coef(m)[k]; se <- sqrt(diag(V))[k]
  list(k = k, est = est, se = se, t = est/se,
       p_naive = summary(m)$coefficients[k, 4],
       p_clus = 2*pt(-abs(est/se), nc - 1), nc = nc, n = nrow(X))
}

cat("=== the pool (frozen spec, four-seam primary, velocity-adjusted residual) ===\n")
print(P[, .(seasons = .N, arms = uniqueN(id), axis_sd = round(sd(axis),1),
            arm_sd = round(sd(arm),1)), by = league], row.names = FALSE)

## ---- 1. the interaction, per league ---------------------------------------------------------
cat("\n=== 1. y ~ axis * arm, both z-scored within league ===\n")
cat("   a NEGATIVE interaction is the claimed effect: the axis slope gets steeper as slot rises\n\n")
SPECS <- list("bare"                  = y ~ za * zr,
              "+ velo sep"            = y ~ za * zr + zv,
              "+ velo sep, changeup active spin" = y ~ za * zr + zv + ze)
for (lg in c("MLB","D1")) {
  cat(sprintf("  ---- %s ----\n", lg))
  for (nm in names(SPECS)) {
    r <- crob(SPECS[[nm]], P[league == lg], "id")
    cat(sprintf("    %-34s beta %+.4f  se %.4f  t %+.2f | p naive %.4f | p arm-clustered %.4f\n",
                nm, r$est, r$se, r$t, r$p_naive, r$p_clus))
  }
  cat(sprintf("    (%d seasons, %d arms)\n", P[league == lg, .N], P[league == lg, uniqueN(id)]))
}

## ---- 2. stratified: where does the axis slope live? -----------------------------------------
cat("\n=== 2. Spearman(axis, residual) by slot tercile, one point per pitcher-season ===\n")
P[, s3 := cut(arm, quantile(arm, c(0,1/3,2/3,1)), labels = c("low","mid","HIGH"),
              include.lowest = TRUE), by = league]
ST <- P[, { z <- suppressWarnings(cor.test(axis, y, method = "spearman", exact = FALSE))
            .(n = .N, arms = uniqueN(id), r = round(unname(z$estimate),3),
              p = round(z$p.value,4)) }, by = .(league, s3)]
print(ST[order(league, s3)], row.names = FALSE)

# The same thing collapsed to one row per arm, which is the honest n.
cat("\n   collapsed to one row per arm:\n")
A <- P[, .(axis = mean(axis), y = mean(y), arm = mean(arm)), by = .(league, id)]
A[, s3 := cut(arm, quantile(arm, c(0,1/3,2/3,1)), labels = c("low","mid","HIGH"),
              include.lowest = TRUE), by = league]
print(A[, { z <- suppressWarnings(cor.test(axis, y, method = "spearman", exact = FALSE))
            .(arms = .N, r = round(unname(z$estimate),3), p = round(z$p.value,4)) },
        by = .(league, s3)][order(league, s3)], row.names = FALSE)

## ---- 3. within pitcher, the test that decides ------------------------------------------------
cat("\n=== 3. within pitcher: does an arm that CLOSES its axis gap gain more when it is high slot? ===\n")
W <- merge(P[, .(league, id, season, axis, y, arm, vs)],
           P[, .(league, id, season = season + 1L, axis_p = axis, y_p = y, arm_p = arm, vs_p = vs)],
           by = c("league","id","season"))
W[, `:=`(d_axis = axis_p - axis, d_y = y - y_p, d_vs = vs - vs_p, arm_avg = (arm + arm_p)/2)]
W[, `:=`(zd = as.numeric(scale(d_axis)), zam = as.numeric(scale(arm_avg))), by = league]
for (lg in c("MLB","D1")) {
  Wl <- W[league == lg]
  if (nrow(Wl) < 30) { cat(sprintf("  %-4s only %d season pairs, skipped\n", lg, nrow(Wl))); next }
  r <- crob(d_y ~ zd * zam + d_vs, Wl, "id")
  cat(sprintf("  %-4s %d pairs / %d arms   d_axis x slot beta %+.4f  se %.4f | p clustered %.3f\n",
              lg, nrow(Wl), uniqueN(Wl$id), r$est, r$se, r$p_clus))
  hi <- quantile(Wl$arm_avg, 2/3)
  for (g in c(TRUE, FALSE)) {
    z <- suppressWarnings(cor.test(Wl[(arm_avg >= hi) == g]$d_axis,
                                   Wl[(arm_avg >= hi) == g]$d_y, method = "spearman", exact = FALSE))
    cat(sprintf("         %-9s n=%3d  r = %+.3f  p = %.3f\n",
                if (g) "top third" else "rest", sum((Wl$arm_avg >= hi) == g), z$estimate, z$p.value))
  }
}

## ---- 4. permutation on arms -------------------------------------------------------------------
# Shuffling the residual across ARMS (not seasons) inside a league keeps the repeat structure and
# destroys only the association, which is what the clustered p is approximating analytically.
cat("\n=== 4. permutation: shuffle the residual across arms within league, 4000 draws ===\n")
for (lg in c("MLB","D1")) {
  D <- P[league == lg]
  obs <- coef(lm(y ~ za * zr, D))["za:zr"]
  ids <- unique(D$id)
  nul <- replicate(4000, {
    map <- data.table(id = ids, yy = D[, .(m = mean(y)), by = id][sample(.N)]$m)
    Dp <- merge(D[, .(id, za, zr)], map, by = "id")
    coef(lm(yy ~ za * zr, Dp))["za:zr"] })
  cat(sprintf("  %-4s observed %+.4f | null mean %+.4f sd %.4f | one-sided p = %.4f\n",
              lg, obs, mean(nul), sd(nul), mean(nul <= obs)))
}

## ---- 5. what the original number was ----------------------------------------------------------
cat("\n=== 5. reproducing the original p = .002 on its own terms, for comparison ===\n")
if (file.exists(file.path(MDIR, "parachute_within.rds"))) {
  LW <- readRDS(file.path(MDIR, "parachute_within.rds"))
  S <- LW$CH[, .(np = .N, nsw = sum(is_swing), axis = mean(axis_diff), kill = -mean(az_diff),
                 velo_sep = -mean(speed_diff), spin = mean(release_spin_rate, na.rm = TRUE)),
             by = .(pitcher, season)]
  S <- merge(S, LW$SWD[, .(wh = 100*mean(wh_res)), by = .(pitcher, season)], by = c("pitcher","season"))
  aa <- fread(file.path(MDIR, "arm_angle_tunnel.csv"))[pitch_type == "CH", .(pitcher, season, arm = mean_arm)]
  S <- merge(S, unique(aa, by = c("pitcher","season")), by = c("pitcher","season"))
  S <- S[nsw >= 60 & is.finite(axis) & is.finite(arm)]
  S[, id := as.character(pitcher)]
  m0 <- summary(lm(wh ~ axis * arm, S))$coefficients["axis:arm", ]
  r0 <- crob(wh ~ axis * arm, S, "id")
  cat(sprintf("  original sample (n=%d seasons, %d arms, 60+ swings, no ff filter, raw residual)\n",
              nrow(S), uniqueN(S$id)))
  cat(sprintf("    as published, season-level p            : beta %+.5f  p = %.4f\n", m0[1], m0[4]))
  cat(sprintf("    same model, arm-clustered               : beta %+.5f  p = %.4f (%d arms)\n",
              r0$est, r0$p_clus, r0$nc))
}
