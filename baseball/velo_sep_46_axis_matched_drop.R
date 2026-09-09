#!/usr/bin/env Rscript

# INSIDE THE AXIS-MATCHED POPULATION, WHAT IS NECESSARY FOR WHIFF RESIDUAL?
#
# Two questions, both restricted to changeups that already spin like the fastball.
#
#   1. Necessity. Does an axis-matched changeup need seam shift, or velocity separation, or both,
#      before it beats its whiff model? A trait is necessary in the useful sense if the group that
#      lacks it sits at zero.
#
#   2. Mechanism. Is the operative thing DOWNWARD movement the hitter cannot pick up - the pitch
#      falling further than its own spin says it should? That is a specific, computable quantity,
#      not the same as raw drop and not the same as the unsigned seam deviation used so far.
#
# The distinction in (2) matters because the whiff model already carries the pitch's movement. Raw
# drop is priced. What is not priced is the DISCREPANCY between where the spin says the ball should
# go and where it actually goes, because the model is never given the spin direction. So if a hitter
# is reading spin and getting betrayed by the seams, the residual is where that shows up.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 225); MDIR <- "data/statcast_model"

P <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(P)
P <- P[is.finite(ax) & is.finite(az) & is.finite(sax) & is.finite(cax) & is.finite(ax_diff) &
       is.finite(az_diff) & is.finite(release_speed) & is.finite(release_extension)]
P[, `:=`(id = as.character(pitcher), lh = p_throws == "L")]

# mirror lefties, strip gravity so what is left is the magnus (spin-driven) component
P[, `:=`(mx = fifelse(lh, -ax, ax), mz = az + 32.174)]
P[, `:=`(fx = fifelse(lh, -(ax - ax_diff), ax - ax_diff), fz = (az - az_diff) + 32.174)]
P[, spin_axis := (atan2(sax, cax)*180/pi) %% 360]
P[, spin_axis := fifelse(lh, (360 - spin_axis) %% 360, spin_axis)]
wrap <- function(d) ((d + 180) %% 360) - 180

# flight time, so accelerations can be quoted as inches of break
P[, tf := (60.5 - release_extension)/(release_speed*1.467)]
P[, inch := 0.5 * tf^2 * 12]

P[, phi_obs := atan2(mx, mz)*180/pi]
P[, phi_spin_raw := 180 - spin_axis]
OFF <- P[, median(wrap(phi_obs - phi_spin_raw), na.rm = TRUE)]   # calibration offset, as before
P[, phi_spin := phi_spin_raw + OFF]
P[, dev := wrap(phi_obs - phi_spin)]
P[, mag := sqrt(mx^2 + mz^2)]

# where the spin says the ball should go, at the same total magnus magnitude
P[, `:=`(pred_z = mag*cos(phi_spin*pi/180), pred_x = mag*sin(phi_spin*pi/180))]
P[, `:=`(vdev = (mz - pred_z)*inch,      # negative -> falls MORE than its spin predicts
         hdev = (mx - pred_x)*inch)]
P[, `:=`(ivb = mz*inch, hb = mx*inch, f_ivb = fz*inch, f_hb = fx*inch)]
P[, `:=`(ivb_gap = f_ivb - ivb, sep = -speed_diff)]

S <- P[, .(np = .N, nm = player_name[1], nsw2 = sum(is_swing == 1),
           axis = mean(axis_diff), sep = mean(sep), absdev = mean(abs(dev)), sgndev = mean(dev),
           vdev = mean(vdev), hdev = mean(hdev), ivb = mean(ivb), hb = mean(hb),
           ivb_gap = mean(ivb_gap), vaa = mean(VAA, na.rm = TRUE), slot = mean(arm_angle),
           spin = mean(release_spin_rate, na.rm = TRUE)), by = .(id, season)]

M <- readRDS(file.path(MDIR, "archetype_roster.rds")); setDT(M)
D <- merge(M[league == "MLB", .(id, season, nsw, y, ych, yrv, ssw, sswsign)], S,
           by = c("id","season"))
D <- D[nsw >= 75 & is.finite(y) & is.finite(vdev)]

cat(sprintf("Full population: %d seasons, %d arms.\n", nrow(D), uniqueN(D$id)))
cat(sprintf("Calibration offset applied to the spin-implied direction: %+.1f deg\n\n", OFF))

A <- D[axis <= 16.3]
cat(sprintf("AXIS-MATCHED (gap <= 16.3 deg): %d seasons, %d arms.\n\n", nrow(A), uniqueN(A$id)))

clus <- function(dat, f, k, wt = "nsw") {
  environment(f) <- environment()
  dat <- dat[complete.cases(dat[, c(all.vars(f), "id", wt), with = FALSE])]
  w <- as.numeric(dat[[wt]]); m <- lm(f, dat, weights = w)
  u <- residuals(m)*sqrt(w); X <- model.matrix(m)*sqrt(w); nc <- uniqueN(dat$id)
  keep <- !is.na(coef(m)); X <- X[, keep, drop = FALSE]
  b <- solve(crossprod(X)); V <- b %*% crossprod(rowsum(X*u, dat$id)) %*% b * (nc/(nc-1))
  j <- which(names(coef(m))[keep] == k)
  e <- unname(coef(m)[keep][j]); s <- unname(sqrt(diag(V))[j]); t <- qt(.975, nc-1)
  c(est = e, se = s, lo = e-t*s, hi = e+t*s, p = 2*pt(-abs(e/s), nc-1))
}
st <- function(r, d = 2) paste0(formatC(r["est"], width = d+5, digits = d, format = "f", flag = "+"),
  ifelse(r["p"] < .01, " **", ifelse(r["p"] < .05, " * ", "   ")))
ci <- function(r, d = 2) sprintf("[%+.*f,%+.*f]", d, r["lo"], d, r["hi"])

## =============================================================================================
cat("=== 1. WHAT DOES THE SIGNED DEVIATION PHYSICALLY MEAN? ===\n\n")
cat("  before using it, establish which way is which. correlations across all seasons:\n\n")
for (v in c("vdev","hdev","ivb","hb")) {
  cat(sprintf("  signed deviation vs %-32s r = %+.3f\n",
              c(vdev = "LIFT beyond spin (in, + = rides)", hdev = "run beyond spin (in)",
                ivb = "vertical break (in)", hb = "horizontal break (in)")[[v]],
              D[, cor(sgndev, get(v), use = "complete.obs")]))
}
cat(sprintf("\n  lift-beyond-spin ranges %+.1f to %+.1f inches, mean %+.2f\n",
            min(D$vdev), max(D$vdev), mean(D$vdev)))
cat("  POSITIVE means the pitch holds its plane better than its spin direction predicts.\n")
cat("  NEGATIVE means it finishes lower than its spin predicts - the hypothesis under test.\n")
cat(sprintf("  share of all changeup seasons that fall below their spin prediction: %.0f%%\n",
            100*mean(D$vdev < 0)))

## =============================================================================================
cat("\n=== 2. NECESSITY: THE FOUR CORNERS INSIDE THE AXIS-MATCHED GROUP ===\n\n")
A[, `:=`(hs = sep >= 8.7, hd = absdev >= 7.8)]
CN <- A[, .(seasons = .N, arms = uniqueN(id), sep = round(mean(sep),1),
            `seam deg` = round(mean(absdev),1), `drop beyond spin` = round(mean(vdev),2),
            `whiff above model` = round(weighted.mean(y, nsw),2),
            `chase above` = round(weighted.mean(ych, nsw),2),
            `rv above` = round(weighted.mean(yrv, nsw),2)),
        by = .(separation = fifelse(hs, "8.7+", "under"), seam = fifelse(hd, "7.8+", "under"))]
setorder(CN, -separation, -seam)
print(CN, row.names = FALSE)
cat("\n  each corner tested against zero:\n\n")
for (i in 1:nrow(CN)) {
  s <- CN$separation[i]; d <- CN$seam[i]
  X <- A[fifelse(hs, "8.7+", "under") == s & fifelse(hd, "7.8+", "under") == d]
  if (nrow(X) < 8) next
  r <- clus(X, y ~ 1, "(Intercept)")
  cat(sprintf("  separation %-6s seam %-6s  %s  %s  p = %.4f  (n = %d)\n",
              s, d, st(r), ci(r), r["p"], nrow(X)))
}

## =============================================================================================
cat("\n=== 3. WHICH MOVEMENT QUANTITY CARRIES IT? ===\n\n")
cat("  all standardised within the axis-matched group, one at a time, arm-clustered.\n")
cat("  'drop beyond spin' is signed so that MORE DROP THAN SPIN PREDICTS is positive.\n\n")
A[, `:=`(z_sep = scale(sep)[,1], z_abs = scale(absdev)[,1], z_sgn = scale(sgndev)[,1],
         z_drop = scale(-vdev)[,1], z_run = scale(hdev)[,1], z_ivb = scale(-ivb)[,1],
         z_gap = scale(ivb_gap)[,1], z_vaa = scale(-vaa)[,1], z_slot = scale(slot)[,1])]
LAB <- c(z_sep = "Velocity separation", z_abs = "Seam deviation (unsigned)",
         z_sgn = "Seam deviation (signed)", z_drop = "Drop beyond spin prediction",
         z_run = "Arm-side run beyond spin", z_ivb = "Raw drop (less vertical break)",
         z_gap = "Vertical break gap vs fastball", z_vaa = "Steeper approach angle",
         z_slot = "Arm slot")
cat(sprintf("  %-34s %10s %18s %9s\n", "trait, per 1 SD", "whiff pts", "95% interval", "p"))
for (v in names(LAB)) {
  r <- clus(A, as.formula(paste("y ~", v)), v)
  cat(sprintf("  %-34s %10s %18s %9.4f\n", LAB[[v]], st(r), ci(r), r["p"]))
}

## =============================================================================================
cat("\n=== 4. HORSE RACE AMONG THE MOVEMENT TERMS ===\n\n")
RACE <- list(
  "drop beyond spin vs unsigned seam"  = y ~ z_drop + z_abs,
  "drop beyond spin vs raw drop"       = y ~ z_drop + z_ivb,
  "drop beyond spin vs separation"     = y ~ z_drop + z_sep,
  "drop, run, separation"              = y ~ z_drop + z_run + z_sep,
  "everything movement-side"           = y ~ z_drop + z_run + z_sep + z_ivb + z_slot)
KEYS <- c("z_drop","z_abs","z_ivb","z_sep","z_run","z_slot")
cat(sprintf("  %-36s %11s %11s %11s %11s %11s %11s\n", "specification", "drop>spin",
            "seam abs", "raw drop", "separation", "run>spin", "slot"))
for (nm in names(RACE)) {
  f <- RACE[[nm]]; vars <- all.vars(f)[-1]
  cells <- sapply(KEYS, function(v) if (v %in% vars) st(clus(A, f, v)) else "      -    ")
  cat(sprintf("  %-36s %11s %11s %11s %11s %11s %11s\n", nm, cells[1], cells[2], cells[3],
              cells[4], cells[5], cells[6]))
}

## =============================================================================================
cat("\n=== 5. DOSE RESPONSE ON DROP BEYOND SPIN, INSIDE AXIS-MATCHED ===\n\n")
A[, dq := cut(-vdev, quantile(-vdev, 0:4/4), include.lowest = TRUE,
              labels = c("rises vs spin","neutral","drops","drops hard"))]
print(A[, .(seasons = .N, arms = uniqueN(id),
            `drop beyond spin (in)` = round(mean(-vdev),2), `seam deg` = round(mean(absdev),1),
            `separation` = round(mean(sep),1), `raw vert break (in)` = round(mean(ivb),1),
            `whiff above model` = round(weighted.mean(y, nsw),2),
            `chase above` = round(weighted.mean(ych, nsw),2),
            `rv above` = round(weighted.mean(yrv, nsw),2)), by = dq][order(dq)], row.names = FALSE)

cat("\n  and the same split OUTSIDE the axis-matched group, as a control:\n\n")
B <- D[axis > 16.3]
B[, dq := cut(-vdev, quantile(-vdev, 0:4/4), include.lowest = TRUE,
              labels = c("rises vs spin","neutral","drops","drops hard"))]
print(B[, .(seasons = .N, `drop beyond spin (in)` = round(mean(-vdev),2),
            `whiff above model` = round(weighted.mean(y, nsw),2)), by = dq][order(dq)],
      row.names = FALSE)
Bz <- copy(B)[, z_drop := scale(-vdev)[,1]]
r <- clus(Bz, y ~ z_drop, "z_drop")
cat(sprintf("\n  drop-beyond-spin effect outside the axis match: %s  %s  p = %.4f\n",
            st(r), ci(r), r["p"]))
r <- clus(copy(D)[, `:=`(z_drop = scale(-vdev)[,1], am = axis <= 16.3)], y ~ z_drop*am, "z_drop:amTRUE")
cat(sprintf("  interaction, drop-beyond-spin x axis matched:   %s  %s  p = %.4f\n",
            st(r), ci(r), r["p"]))

## =============================================================================================
cat("\n=== 6. IS IT NECESSARY? LOOKING FROM THE OUTCOME BACK ===\n\n")
# necessity runs the other direction from prediction: among the axis-matched arms that DID beat the
# model, how many had each trait? a trait present in nearly all of them is a candidate necessity; a
# trait present in half is not.
W <- A[y >= 2]
cat(sprintf("  axis-matched seasons beating the whiff model by 2+ points: %d of %d (%.0f%%)\n\n",
            nrow(W), nrow(A), 100*nrow(W)/nrow(A)))
for (nm in c("separation 8.7+","seam deviation 7.8+","drops more than its spin predicts",
             "any two of the three","all three")) {
  f <- switch(nm,
    "separation 8.7+" = W$sep >= 8.7,
    "seam deviation 7.8+" = W$absdev >= 7.8,
    "drops more than its spin predicts" = W$vdev < 0,
    "any two of the three" = (W$sep >= 8.7) + (W$absdev >= 7.8) + (W$vdev < 0) >= 2,
    (W$sep >= 8.7) + (W$absdev >= 7.8) + (W$vdev < 0) == 3)
  g <- switch(nm,
    "separation 8.7+" = A$sep >= 8.7,
    "seam deviation 7.8+" = A$absdev >= 7.8,
    "drops more than its spin predicts" = A$vdev < 0,
    "any two of the three" = (A$sep >= 8.7) + (A$absdev >= 7.8) + (A$vdev < 0) >= 2,
    (A$sep >= 8.7) + (A$absdev >= 7.8) + (A$vdev < 0) == 3)
  cat(sprintf("  %-36s present in %3.0f%% of them, %3.0f%% of all axis-matched\n",
              nm, 100*mean(f), 100*mean(g)))
}
## =============================================================================================
cat("\n=== 7. SEPARATION PAYS A RESIDUAL ONLY WHEN THE AXIS MATCHES ===\n\n")
# across the whole population separation buys raw whiff and no residual, because the model already
# prices the changeup's own speed. inside the axis match it buys a residual. that is an interaction
# claim and it should be tested as one rather than inferred from two subgroup runs.
F <- copy(D)[, `:=`(z_sep = scale(sep)[,1], z_abs = scale(absdev)[,1], am = axis <= 16.3)]
for (k in c("z_sep","z_abs")) {
  r <- clus(F, as.formula(paste0("y ~ ", k, "*am")), paste0(k, ":amTRUE"))
  m <- clus(F, as.formula(paste0("y ~ ", k, "*am")), k)
  cat(sprintf("  %-26s main effect (unmatched) %s   interaction with axis match %s  p = %.4f\n",
              c(z_sep = "Velocity separation", z_abs = "Seam deviation")[[k]],
              st(m), st(r), r["p"]))
}
cat("\n  separation effect on whiff residual, measured inside each population:\n\n")
for (a in c(TRUE, FALSE)) {
  X <- copy(D[(axis <= 16.3) == a])[, z_sep := scale(sep)[,1]]
  r <- clus(X, y ~ z_sep, "z_sep")
  cat(sprintf("  %-22s %s  %s  p = %.4f  (n = %d)\n",
              if (a) "axis matched" else "axis not matched", st(r), ci(r), r["p"], nrow(X)))
}
cat("\n  the conjunction against what the linear terms would predict for it:\n\n")
A2 <- copy(A)[, both := sep >= 8.7 & absdev >= 7.8]
r1 <- clus(A2, y ~ both, "bothTRUE")
r2 <- clus(A2, y ~ both + z_sep + z_abs, "bothTRUE")
cat(sprintf("  both gates, alone                    %s  %s\n", st(r1), ci(r1)))
cat(sprintf("  both gates, net of the linear terms  %s  %s\n", st(r2), ci(r2)))
saveRDS(list(all = D, matched = A), file.path(MDIR, "axis_matched_drop.rds"))
