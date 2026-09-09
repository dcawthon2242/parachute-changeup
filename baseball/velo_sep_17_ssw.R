#!/usr/bin/env Rscript

# SEPARATING SEAM-SHIFT EFFECTS FROM AXIS MATCHING.
#
# Active spin and seam-shifted wake are different things and the parachute definition only ever
# controlled the first. Active spin asks how much of the total spin is transverse - how much of it
# COULD produce Magnus movement. Seam-shifted wake is a DIRECTIONAL failure: the ball moves somewhere
# the spin axis says it should not, because the seam orientation steers the wake. A pitch can have
# ordinary active spin and still be heavily seam-shifted.
#
# The signature is therefore an angle, not a magnitude: the gap between where the ball actually
# moved and where its measured spin axis predicts it should have moved. Statcast reports spin_axis
# from Hawk-Eye imaging of the ball, independent of the trajectory, so the two are separately
# measured and their disagreement is the seam-shift proxy. This is the standard "axis deviation"
# construction.
#
# Two conventions have to be handled. Statcast reports spin_axis so that 180 degrees is pure
# backspin, and movement must be taken from the MAGNUS acceleration, meaning gravity is added back
# out of az first - otherwise every pitch points nearly straight down and the angle is meaningless.
# Rather than trust a hand-derived mapping, the sign and offset are calibrated empirically below by
# picking the form that makes the population distribution tightest.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 200); MDIR <- "data/statcast_model"
P <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(P)
P <- P[is.finite(ax) & is.finite(az) & is.finite(sax) & is.finite(cax)]
P[, lh := p_throws == "L"]
# Mirror lefties so both hands share a frame, then take the movement angle from Magnus acceleration.
P[, `:=`(mx = fifelse(lh, -ax, ax), mz = az + 32.174)]
P[, spin_axis := (atan2(sax, cax) * 180/pi) %% 360]
P[, spin_axis := fifelse(lh, (360 - spin_axis) %% 360, spin_axis)]
P[, theta_mov := atan2(mx, mz) * 180/pi]

wrap <- function(d) ((d + 180) %% 360) - 180
cand <- list(
  "mov + (axis-180)" = quote(wrap(theta_mov + (spin_axis - 180))),
  "mov - (axis-180)" = quote(wrap(theta_mov - (spin_axis - 180))),
  "mov + (180-axis)" = quote(wrap(theta_mov + (180 - spin_axis))))
cat("=== calibrating the spin-axis to movement-direction mapping ===\n")
for (k in names(cand)) {
  d <- P[, eval(cand[[k]])]
  cat(sprintf("  %-18s median %+7.1f deg | IQR width %5.1f | sd %5.1f\n",
              k, median(d, na.rm = TRUE), IQR(d, na.rm = TRUE), sd(d, na.rm = TRUE)))
}
P[, dev_raw := wrap(theta_mov + (spin_axis - 180))]
off <- median(P$dev_raw, na.rm = TRUE)
P[, dev := wrap(dev_raw - off)]
cat(sprintf("\n  using 'mov + (axis-180)', centred on the population median of %+.1f degrees.\n", off))
cat("  dev is now DEGREES OF SEAM-SHIFT DEVIATION: how far the ball moved from where its\n")
cat("  imaged spin axis said it would, relative to the typical changeup.\n")

S <- P[, .(nsw = sum(is_swing), dev = mean(dev, na.rm = TRUE),
           absdev = mean(abs(dev), na.rm = TRUE),
           theta = mean(theta_mov), spinax = mean(spin_axis),
           ivb = mean(mz)/32.174*12, hb = mean(mx)/32.174*12,
           spin = mean(release_spin_rate, na.rm = TRUE),
           velo = mean(release_speed)), by = .(pitcher, season)][nsw >= 40]
S[, id := as.character(pitcher)]
cat(sprintf("\n%d pitcher-seasons. deviation: median %+.1f, 10th pct %+.1f, 90th pct %+.1f degrees\n",
            nrow(S), median(S$dev), quantile(S$dev,.10), quantile(S$dev,.90)))

## ---- where do the 13 members sit? ------------------------------------------------------------------
L <- readRDS(file.path(MDIR, "locked_spec.rds"))$data; setDT(L); M <- L[league == "MLB"]
NM <- unique(readRDS(file.path(MDIR, "parachute_ff.rds"))[, .(pitcher, player_name)]); setDT(NM)
NM <- unique(NM, by = "pitcher")[, .(id = as.character(pitcher), nm = sub(",.*", "", player_name))]
M <- merge(M, NM, by = "id", all.x = TRUE)
M[, pa := axis < 10 & arm >= arm_thr]
H <- merge(M[vs >= quantile(M$vs, 2/3)], S[, .(id, season, dev, absdev, ivb, hb, spin, velo)],
           by = c("id","season"))
S[, pct := frank(dev)/.N]
H <- merge(H, S[, .(id, season, dev_pct = pct)], by = c("id","season"))
cat("\n=== the 13 members, ranked by seam-shift deviation ===\n")
cat("    dev_pct is the percentile against all changeups; near 0 or 1 means unusual\n\n")
print(H[pa == TRUE][order(-abs(dev)), .(nm, season, swings = nsw, whiff_over = round(y,1),
        ch_spin = round(ec,3), seam_dev = round(dev,1), pctile = round(dev_pct,2),
        ivb = round(ivb,1), hb = round(hb,1), spin = round(spin))], row.names = FALSE)
cat(sprintf("\n  member mean |deviation| %.1f deg vs %.1f for other high-sep seasons\n",
            H[pa == TRUE, mean(abs(dev))], H[pa == FALSE, mean(abs(dev))]))

## ---- does deviation explain the effect, or is it separate? -------------------------------------------
crob <- function(D, f, k) {
  D <- D[complete.cases(D[, c(all.vars(f), "id"), with = FALSE])]
  m <- lm(f, D); u <- residuals(m); X <- model.matrix(m); nc <- uniqueN(D$id)
  if (!k %in% colnames(X)) return(c(NA,NA,NA))
  b <- solve(crossprod(X)); V <- b %*% crossprod(rowsum(X*u, D$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[k]); s <- unname(sqrt(diag(V))[k]); c(e, s, 2*pt(-abs(e/s), nc-1))
}
cat("\n=== is seam-shift deviation its own driver of overperformance? ===\n")
H[, zdev := as.numeric(scale(abs(dev)))]
for (spec in list(list("|deviation| alone", y ~ zdev + vs, "zdev"),
                  list("axis+arm alone", y ~ pa + vs, "paTRUE"),
                  list("both together (pa)", y ~ pa + zdev + vs, "paTRUE"),
                  list("both together (dev)", y ~ pa + zdev + vs, "zdev"),
                  list("also controlling active spin", y ~ pa + zdev + ec + vs, "paTRUE"))) {
  r <- crob(H, spec[[2]], spec[[3]])
  cat(sprintf("  %-30s %-8s %+.3f (se %.3f) p = %.4f\n", spec[[1]], spec[[3]], r[1], r[2], r[3]))
}
cat("\n=== dropping the most seam-shifted members and refitting ===\n")
for (q in c(1.00, 0.90, 0.80, 0.70)) {
  cut <- quantile(abs(H$dev), q); D <- H[abs(dev) <= cut]
  r <- crob(D, y ~ pa + vs, "paTRUE")
  cat(sprintf("  keep |dev| below the %.0fth pct (%.1f deg): %+.3f (se %.3f) p = %.4f | %d members\n",
              100*q, cut, r[1], r[2], r[3], sum(D$pa)))
}
saveRDS(S, file.path(MDIR, "ch_seam_deviation.rds"))
cat("\nwrote ch_seam_deviation.rds (pitcher-season changeup seam-shift deviation)\n")
