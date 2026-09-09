#!/usr/bin/env Rscript

# D1 REPLICATION: DOES AXIS MATCH PAY OFF MORE WHEN VELOCITY SEPARATION IS LARGE?
#
# MLB gave axis match +1.13 whiff points per SD inside the top third of separation and nothing
# outside it, but the interaction that would establish the conditionality came back p = .145. D1
# has 733 arms against MLB's 442, so it is the higher-powered test, and it did not generate the
# hypothesis.
#
# The D1 whiff residual is the right comparison: its model already carries the pitch's own spin
# axis, plate location, approach angles and count, so an axis effect there is not an omitted
# feature the way it was against the fastball-blind model.
#
# ONE CONSTRUCT PROBLEM GOVERNS EVERYTHING BELOW. TrackMan's changeup spin axis sits within about
# five degrees of a deterministic function of the break, so the D1 axis gap is close to a
# movement-direction gap rather than imaged rotation. Recomputed that way, only 0.9 percent of MLB
# seasons land inside ten degrees against 10.2 percent measured. So a null in D1 is ambiguous
# between "the effect is not real" and "D1 cannot measure the thing." To separate those, section 3
# rebuilds the MLB axis gap from ax/az the way TrackMan effectively reports it and reruns the MLB
# test on that degraded construct. If MLB survives degradation and D1 still nulls, the effect is
# suspect. If MLB dies under degradation too, D1 was never able to see it and the null is empty.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 200); set.seed(67); MDIR <- "data/statcast_model"

crob <- function(D, f, k) {
  D <- D[complete.cases(D[, c(all.vars(f), "id"), with = FALSE])]
  m <- lm(f, D); u <- residuals(m); X <- model.matrix(m); nc <- uniqueN(D$id)
  if (!k %in% colnames(X)) return(list(est = NA, se = NA, p = NA, n = nrow(D), arms = nc))
  b <- solve(crossprod(X)); V <- b %*% crossprod(rowsum(X*u, D$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[k]); s <- unname(sqrt(diag(V))[k])
  list(est = e, se = s, p = 2*pt(-abs(e/s), nc-1), n = nrow(D), arms = nc)
}
# One routine for both leagues so nothing differs but the data.
report <- function(D, tag) {
  D <- copy(D)[is.finite(y) & is.finite(axis) & is.finite(vs)]
  D[, `:=`(axis_match = as.numeric(scale(-axis)), zsep = as.numeric(scale(vs)))]
  i <- crob(D, y ~ axis_match * zsep, "axis_match:zsep")
  cat(sprintf("\n  %s  (%d seasons, %d arms)\n", tag, nrow(D), uniqueN(D$id)))
  cat(sprintf("    interaction axis_match x separation : %+.3f  se %.3f  p = %.4f\n",
              i$est, i$se, i$p))
  thr <- quantile(D$vs, 2/3); D[, hi := vs >= thr]
  for (h in c(TRUE, FALSE)) {
    S <- copy(D[hi == h]); S[, axis_match := as.numeric(scale(axis_match))]
    r <- crob(S, y ~ axis_match + vs, "axis_match")
    cat(sprintf("    %-22s (>= %.1f mph) : %+.3f  se %.3f  p = %.4f   n=%d arms=%d\n",
                if (h) "top third of sep" else "rest", thr, r$est, r$se, r$p, r$n, r$arms))
  }
  # Smallest interaction this sample could have caught, two-sided at 80 percent power.
  cat(sprintf("    minimum detectable interaction at 80%% power: %.3f\n", 2.8 * i$se))
  invisible(i)
}

## ---- 1. D1 whiff -----------------------------------------------------------------------------------
L <- readRDS(file.path(MDIR, "locked_spec.rds"))$data; setDT(L)
D1 <- L[league == "D1" & !is.na(id) & id != "0" & id != ""]
cat("=== 1. D1 whiff residual (model already carries spin axis, location, approach angles, count) ===")
report(D1, "D1 whiff")

## ---- 2. D1 chase -----------------------------------------------------------------------------------
CH <- readRDS(file.path(MDIR, "ncaa_chase_resid.rds")); setDT(CH)
cat("\n\n=== 2. D1 chase residual ===\n")
cat("   cols:", paste(names(CH), collapse = ", "), "\n")
rc <- grep("^r_all$|^r_tj$", names(CH), value = TRUE)[1]
idc <- grep("PitcherId", names(CH), value = TRUE)[1]
if (!is.na(rc) && !is.na(idc)) {
  CS <- CH[, .(nz = .N, ych = 100*mean(get(rc))), by = .(id = as.character(get(idc)), season)][nz >= 40]
  D1c <- merge(D1[, .(id, season, axis, vs)], CS, by = c("id","season"))
  setnames(D1c, "ych", "y"); report(D1c, "D1 chase")
} else cat("   (no usable chase residual column)\n")

## ---- 3. is the MLB signal visible through a TrackMan-style axis? --------------------------------------
cat("\n\n=== 3. construct check: rebuild the MLB axis gap as a movement-direction gap ===\n")
P <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(P)
# The file is changeup rows only; the four-seam anchor lives in ax_diff/az_diff as CH minus FF, so
# the fastball movement is recovered by subtraction. Lefties mirrored so both hands share a frame.
# az is raw vertical acceleration and is dominated by gravity, which points both pitches nearly
# straight down and collapses the gap to near zero. Movement direction is the MAGNUS component, so
# gravity is added back out before taking the angle.
P <- P[is.finite(ax) & is.finite(az) & is.finite(ax_diff) & is.finite(az_diff)]
P[, `:=`(chx = fifelse(p_throws == "L", -ax, ax), chz = az + 32.174,
         ffx = fifelse(p_throws == "L", -(ax - ax_diff), ax - ax_diff),
         ffz = (az - az_diff) + 32.174)]
A <- P[, .(n = .N, ca = atan2(mean(chx), mean(chz))*180/pi,
           fa = atan2(mean(ffx), mean(ffz))*180/pi), by = .(pitcher, season)][n >= 30]
A[, mv_axis := abs(((ca - fa + 180) %% 360) - 180)]
A[, id := as.character(pitcher)]
MLB <- merge(L[league == "MLB", .(id, season, y, axis, vs)], A[, .(id, season, mv_axis)],
             by = c("id","season"))
cat(sprintf("   measured axis gap vs movement-direction gap: r = %+.3f (n = %d)\n",
            cor(MLB$axis, MLB$mv_axis, use = "complete.obs"), nrow(MLB)))
cat(sprintf("   inside 10 degrees: measured %.1f%% of seasons, movement-direction %.1f%%\n",
            100*mean(MLB$axis < 10, na.rm = TRUE), 100*mean(MLB$mv_axis < 10, na.rm = TRUE)))
cat("\n   MLB on its own MEASURED axis (the construct the hypothesis was built on):")
report(MLB, "MLB measured axis")
cat("\n   MLB on the DEGRADED movement-direction axis (what D1 can actually see):")
report(copy(MLB)[, axis := mv_axis], "MLB movement axis")
