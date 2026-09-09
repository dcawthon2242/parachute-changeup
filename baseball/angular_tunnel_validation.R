#!/usr/bin/env Rscript

# Does the angular tunnel metric predict anything, tested at scale on three outcomes:
# the swing-decision run value, chase rate, and miss distance.
#
# The earlier whiff test was underpowered by construction: it aggregated to
# pitcher x pitch-type, leaving 80-134 units, and then asked for a Spearman
# correlation. Here the unit is the pitch, which is ~96k pairs, with pitcher fixed
# effects absorbing arm quality and cluster-robust standard errors by pitcher.
#
# WHY THESE THREE OUTCOMES ARE BETTER TARGETS THAN WHIFF RATE
#   sd_rv_exp     already residualised on location, movement, velocity and count by the
#                 event cascade, so "where it finished" cannot drive it. Sequencing is
#                 one of the few things left that could.
#   chase         a pure decision variable. If tunneling does anything, it should show up
#                 in whether the hitter commits to a ball, not in how he contacts it.
#   miss_distance only exists on whiffs (22.6% of swings), so it is conditional on already
#                 having missed and measures how badly. Smaller sample, cleaner mechanism.
#
# EXPECTED SIGNS if tunneling helps the pitcher: a longer tunnel (higher break fraction)
# should push sd_rv_exp DOWN (negative coefficient, hitter loses), chase UP, and miss
# distance UP.
#
# PRIMARY SPECIFICATION, fixed before looking: tau = 0.05 deg, combined position+depth
# cue (brk_any_005). Everything else in the sweep is reported as robustness, and the
# multiplicity that implies is stated with the results rather than hidden.

suppressPackageStartupMessages(library(data.table))
options(width = 220)
PRIMARY <- "brk_any_005"
TYPES <- c("ST","SL","CU","KC","CH","FS")

d <- readRDS("data/statcast_model/angular_tunnel_2026.rds")
sd_rv <- readRDS("data/statcast_model/swing_decision_rv_2026.rds")[
  , .(game_pk, at_bat_number, pitch_number, sd_rv_exp, sd_rv_argmax)]
d <- merge(d, sd_rv, by = c("game_pk","at_bat_number","pitch_number"), all.x = TRUE)

# strike-zone geometry and velocity, needed as controls for chase and miss distance
cf <- list.files("data/statcast_2026/chunks", pattern = "csv$", full.names = TRUE)
ex <- rbindlist(lapply(cf, function(f) fread(f, select = c(
  "game_pk","at_bat_number","pitch_number","sz_top","sz_bot","release_speed","bat_speed"),
  showProgress = FALSE)))
ex <- unique(ex, by = c("game_pk","at_bat_number","pitch_number"))
d <- merge(d, ex, by = c("game_pk","at_bat_number","pitch_number"), all.x = TRUE)
cat(sprintf("pairs: %s | sd_rv attached %.1f%%\n", format(nrow(d), big.mark = ","),
            100*mean(!is.na(d$sd_rv_exp))))

# distance outside the rulebook zone, 0 when the pitch is in the zone
d[, dx := pmax(abs(plate_x) - (0.83 + 0.12), 0)]
d[, dz := pmax(plate_z - sz_top, sz_bot - plate_z, 0)]
d[, out_zone := (dx > 0 | dz > 0)]
d[, zdist := sqrt(dx^2 + dz^2)]
d[, chase := swing]

## ---- fixed-effects OLS with cluster-robust SE ---------------------------
# Pitcher fixed effects are applied by within-demeaning, which is exact for a single
# factor and avoids building a 1000-column design matrix.
fe_fit <- function(dt, yvar, xvar, controls = character()) {
  dt <- dt[is.finite(get(yvar)) & is.finite(get(xvar))]
  for (cv in controls) dt <- dt[is.finite(get(cv))]
  if (uniqueN(dt$pitcher) < 20 || nrow(dt) < 400) return(NULL)
  vars <- c(yvar, xvar, controls)
  # Everything must be double BEFORE demeaning. data.table assigns in place, so a
  # logical outcome like chase would have its demeaned values coerced back to
  # TRUE/FALSE, and integer controls like balls/strikes would be truncated.
  dm <- dt[, lapply(.SD, as.numeric), .SDcols = vars]
  dm[, pitcher := dt$pitcher]
  for (v in vars) dm[, (v) := get(v) - mean(get(v)), by = pitcher]
  X <- cbind(1, as.matrix(dm[, c(xvar, controls), with = FALSE]))
  y <- dm[[yvar]]
  XtXi <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
  if (is.null(XtXi)) return(NULL)
  b <- XtXi %*% crossprod(X, y)
  e <- as.vector(y - X %*% b)
  g <- split(seq_len(nrow(X)), dm$pitcher)
  meat <- matrix(0, ncol(X), ncol(X))
  for (ix in g) { u <- crossprod(X[ix, , drop = FALSE], e[ix]); meat <- meat + u %*% t(u) }
  V <- XtXi %*% meat %*% XtXi
  nc <- length(g)
  V <- V * (nc/(nc - 1))
  data.table(n = nrow(dt), arms = nc, beta = b[2], se = sqrt(V[2,2]),
             t = b[2]/sqrt(V[2,2]), p = 2*pt(-abs(b[2]/sqrt(V[2,2])), nc - 1),
             sd_x = sd(dt[[xvar]]))
}

LOC <- c("plate_x","plate_z","zdist","balls","strikes","release_speed")

run_block <- function(xvar) {
  rbindlist(lapply(c("ALL", TYPES), function(ty) {
    s <- if (ty == "ALL") d else d[pitch_type == ty]
    rbindlist(list(
      cbind(outcome = "sd_rv_exp",     ptype = ty, fe_fit(s, "sd_rv_exp", xvar, "strikes")),
      cbind(outcome = "sd_rv_argmax",  ptype = ty, fe_fit(s, "sd_rv_argmax", xvar, "strikes")),
      cbind(outcome = "chase",         ptype = ty, fe_fit(s[out_zone == TRUE], "chase", xvar, LOC)),
      cbind(outcome = "miss_distance", ptype = ty, fe_fit(s[whiff == TRUE], "miss_distance", xvar, LOC))
    ), fill = TRUE)
  }), fill = TRUE)
}

cat("\n\n================ PRIMARY SPECIFICATION:", PRIMARY, "================\n")
cat("beta is the outcome change per unit of break fraction; per_sd rescales it to a\n")
cat("one-standard-deviation move in the metric, which is the interpretable size.\n")
cat("Wanted signs: sd_rv negative, chase positive, miss_distance positive.\n\n")
pr <- run_block(PRIMARY)
pr[, per_sd := beta*sd_x]
print(pr[, .(outcome, ptype, n, arms, beta = round(beta,5), se = round(se,5),
             per_sd = round(per_sd,5), t = round(t,2), p = round(p,4))])

cat("\n\n================ THRESHOLD SWEEP ================\n")
cat("t-statistics only. 7 thresholds x 2 cue variants x 4 outcomes x 7 groupings = 392\n")
cat("tests, so at alpha=0.05 about 20 false positives are expected by chance alone.\n")
allx <- grep("^brk_(any|pos)_", names(d), value = TRUE)
sw <- rbindlist(lapply(allx, function(v) { r <- run_block(v); if (nrow(r)) r[, metric := v]; r }),
                fill = TRUE)
for (oc in c("sd_rv_exp","chase","miss_distance")) {
  cat("\n---", oc, "--- (t-stats; sign convention as above)\n")
  print(dcast(sw[outcome == oc], metric ~ ptype, value.var = "t")[
    , lapply(.SD, function(z) if (is.numeric(z)) round(z,2) else z)])
}

cat("\n\n================ HOW MANY SURVIVE MULTIPLICITY ================\n")
sw[, q := p.adjust(p, method = "BH")]
cat(sprintf("tests run: %d | raw p<0.05: %d | BH q<0.05: %d | q<0.10: %d\n",
            nrow(sw), sum(sw$p < 0.05, na.rm = TRUE),
            sum(sw$q < 0.05, na.rm = TRUE), sum(sw$q < 0.10, na.rm = TRUE)))
srv <- sw[q < 0.10][order(q)]
if (nrow(srv)) {
  cat("\nsurvivors at q<0.10, with the sign we wanted marked:\n")
  srv[, wanted := fifelse(outcome %in% c("sd_rv_exp","sd_rv_argmax"), beta < 0, beta > 0)]
  print(srv[, .(metric, outcome, ptype, n, arms, beta = round(beta,5),
                per_sd = round(beta*sd_x,5), t = round(t,2),
                p = signif(p,2), q = signif(q,2), right_sign = wanted)])
} else cat("\nnothing survives BH correction at q<0.10.\n")

cat("\n\n================ BENCHMARK: THE OLD WEIGHTED-DISTANCE METRIC ================\n")
cat("Same rows, same specification, so the comparison is apples to apples.\n")
d[, tunnel_neg := -tunnel_old]      # sign-flipped so 'higher = tighter tunnel', as brk
ob <- run_block("tunnel_neg")
ob[, per_sd := beta*sd_x]
print(ob[, .(outcome, ptype, n, arms, beta = round(beta,6), per_sd = round(per_sd,5),
             t = round(t,2), p = round(p,4))])

saveRDS(list(primary = pr, sweep = sw, old = ob),
        "data/statcast_model/angular_tunnel_validation.rds")
cat("\nwrote data/statcast_model/angular_tunnel_validation.rds\n")
