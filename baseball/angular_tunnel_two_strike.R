#!/usr/bin/env Rscript

# The chase effect, restricted to two-strike counts.
#
# Two strikes is the interesting place to look because the decision problem changes
# shape. The hitter can no longer take a close pitch, so he expands the zone and swings
# at things he would let go in a neutral count. That cuts both ways for a tunneling
# metric: protection could swamp any perceptual effect, because he is swinging at
# marginal pitches regardless of what they looked like early; or it could amplify it,
# because he has to commit earlier and on less information.
#
# The full count grid is reported alongside 0-2/1-2/2-2 so it is possible to see whether
# two strikes is genuinely special or just where the sample happens to be.
#
# Specification is the hardest one from angular_tunnel_chase_robust.R -- location,
# plate separation between the pair, the previous pitch's location, and fixed effects --
# so nothing here is a weaker test than the headline result.

suppressPackageStartupMessages(library(data.table))
options(width = 200)
PRIMARY <- "brk_any_005"
TS <- c("0-2","1-2","2-2")

d <- readRDS("data/statcast_model/angular_tunnel_2026.rds")
cf <- list.files("data/statcast_2026/chunks", pattern = "csv$", full.names = TRUE)
ex <- rbindlist(lapply(cf, function(f) fread(f, select = c(
  "game_pk","at_bat_number","pitch_number","sz_top","sz_bot","release_speed",
  "plate_x","plate_z"), showProgress = FALSE)))
ex <- unique(ex, by = c("game_pk","at_bat_number","pitch_number"))
prev <- ex[, .(game_pk, at_bat_number, pitch_number = pitch_number + 1L,
               pp_x = plate_x, pp_z = plate_z)]
ex[, c("plate_x","plate_z") := NULL]
d <- merge(d, ex, by = c("game_pk","at_bat_number","pitch_number"), all.x = TRUE)
d <- merge(d, prev, by = c("game_pk","at_bat_number","pitch_number"), all.x = TRUE)
sdr <- readRDS("data/statcast_model/swing_decision_rv_2026.rds")[
  , .(game_pk, at_bat_number, pitch_number, sd_rv_exp)]
d <- merge(d, sdr, by = c("game_pk","at_bat_number","pitch_number"), all.x = TRUE)

d[, plate_sep := sqrt((plate_x - pp_x)^2 + (plate_z - pp_z)^2)]
d[, zdist := sqrt(pmax(abs(plate_x) - 0.95, 0)^2 +
                  pmax(plate_z - sz_top, sz_bot - plate_z, 0)^2)]
d[, `:=`(out_zone = zdist > 0, chase = as.numeric(swing),
         count = paste(balls, strikes, sep = "-"))]
d[, pit_count := paste(pitcher, count)]
d[, two_strike := count %in% TS]

fe_fit <- function(dt, yvar, xvar, controls, fe) {
  vars <- c(yvar, xvar, controls)
  dt <- dt[stats::complete.cases(dt[, vars, with = FALSE])]
  if (nrow(dt) < 300 || uniqueN(dt$pitcher) < 15) return(NULL)
  dm <- dt[, lapply(.SD, as.numeric), .SDcols = vars]
  dm[, `:=`(grp = dt[[fe]], cl = dt$pitcher)]
  for (v in vars) dm[, (v) := get(v) - mean(get(v)), by = grp]
  X <- cbind(1, as.matrix(dm[, c(xvar, controls), with = FALSE])); y <- dm[[yvar]]
  XtXi <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
  if (is.null(XtXi)) return(NULL)
  b <- XtXi %*% crossprod(X, y); e <- as.vector(y - X %*% b)
  meat <- matrix(0, ncol(X), ncol(X))
  for (ix in split(seq_len(nrow(X)), dm$cl)) {
    u <- crossprod(X[ix, , drop = FALSE], e[ix]); meat <- meat + u %*% t(u) }
  V <- XtXi %*% meat %*% XtXi; nc <- uniqueN(dm$cl); V <- V*(nc/(nc-1))
  data.table(n = nrow(dt), arms = nc, beta = b[2], se = sqrt(V[2,2]),
             t = b[2]/sqrt(V[2,2]), p = 2*pt(-abs(b[2]/sqrt(V[2,2])), nc-1),
             per_sd = b[2]*sd(dt[[xvar]]))
}

CTL <- c("plate_x","plate_z","zdist","release_speed","plate_sep","pp_x","pp_z")
OZ <- d[out_zone == TRUE & is.finite(zdist)]

## ---- 1. the three counts asked for, each on its own ---------------------
cat("=== CHASE IN TWO-STRIKE COUNTS ===\n")
cat("metric", PRIMARY, "| pitcher FE + location + plate separation + previous location\n")
cat("per_sd_pp is the change in chase probability, in percentage points, for a 1 SD\n")
cat("longer tunnel. base is the chase rate in that cell.\n\n")
one <- rbindlist(lapply(TS, function(cc) {
  s <- OZ[count == cc]
  r <- fe_fit(s, "chase", PRIMARY, CTL, "pitcher")
  if (is.null(r)) return(NULL)
  cbind(count = cc, base = 100*mean(s$chase), r) }), fill = TRUE)
print(one[, .(count, n, arms, base = round(base,1), beta = round(beta,4),
              se = round(se,4), per_sd_pp = round(100*per_sd,2),
              t = round(t,2), p = signif(p,3))])

cat("\npooled across 0-2/1-2/2-2, with pitcher x count fixed effects:\n")
pl <- fe_fit(OZ[two_strike == TRUE], "chase", PRIMARY, CTL, "pit_count")
cat(sprintf("  n=%s arms=%d base=%.1f%%  per_sd=%+.2f pp  t=%.2f  p=%.2g\n",
            format(pl$n, big.mark = ","), pl$arms,
            100*mean(OZ[two_strike == TRUE]$chase), 100*pl$per_sd, pl$t, pl$p))
cat("\nfor contrast, everything that is not a two-strike count:\n")
nt <- fe_fit(OZ[two_strike == FALSE], "chase", PRIMARY, CTL, "pit_count")
cat(sprintf("  n=%s arms=%d base=%.1f%%  per_sd=%+.2f pp  t=%.2f  p=%.2g\n",
            format(nt$n, big.mark = ","), nt$arms,
            100*mean(OZ[two_strike == FALSE]$chase), 100*nt$per_sd, nt$t, nt$p))

## ---- 2. the whole count grid, so two strikes can be put in context ------
cat("\n\n=== EVERY COUNT, SAME SPECIFICATION ===\n")
grid <- rbindlist(lapply(sort(unique(OZ$count)), function(cc) {
  s <- OZ[count == cc]; r <- fe_fit(s, "chase", PRIMARY, CTL, "pitcher")
  if (is.null(r)) return(data.table(count = cc, n = nrow(s), base = 100*mean(s$chase)))
  cbind(count = cc, base = 100*mean(s$chase), r) }), fill = TRUE)
grid[, strikes := as.integer(substr(count, 3, 3))]
print(grid[order(strikes, count), .(count, n, arms, base = round(base,1),
      per_sd_pp = round(100*per_sd,2), t = round(t,2), p = signif(p,3))])

## ---- 3. by pitch type inside two strikes -------------------------------
cat("\n\n=== BY PITCH TYPE, TWO-STRIKE COUNTS ONLY ===\n")
d[, tunnel_neg := -tunnel_old]; OZ <- d[out_zone == TRUE & is.finite(zdist)]
bt <- rbindlist(lapply(c("ST","SL","CU","KC","CH","FS"), function(ty) {
  s <- OZ[two_strike == TRUE & pitch_type == ty]
  a <- fe_fit(s, "chase", PRIMARY, CTL, "pit_count")
  o <- fe_fit(s, "chase", "tunnel_neg", CTL, "pit_count")
  if (is.null(a)) return(NULL)
  data.table(ptype = ty, n = a$n, base = 100*mean(s$chase),
             ang_per_sd_pp = 100*a$per_sd, ang_t = a$t,
             old_per_sd_pp = if (is.null(o)) NA_real_ else 100*o$per_sd,
             old_t = if (is.null(o)) NA_real_ else o$t) }), fill = TRUE)
print(bt[, lapply(.SD, function(z) if (is.numeric(z)) round(z,2) else z)])

## ---- 4. does the acuity optimum move with two strikes? -----------------
cat("\n\n=== THRESHOLD SWEEP, TWO-STRIKE vs NOT ===\n")
cat("If the hitter is committing earlier with two strikes, the threshold that matters\n")
cat("could shift. t-stats and effect sizes for both subsets.\n\n")
swp <- rbindlist(lapply(grep("^brk_any_", names(d), value = TRUE), function(v) {
  a <- fe_fit(OZ[two_strike == TRUE],  "chase", v, CTL, "pit_count")
  b <- fe_fit(OZ[two_strike == FALSE], "chase", v, CTL, "pit_count")
  if (is.null(a) || is.null(b)) return(NULL)
  data.table(tau_deg = as.numeric(sub("brk_any_", "", v))/100,
             ts_per_sd_pp = 100*a$per_sd, ts_t = a$t,
             non_per_sd_pp = 100*b$per_sd, non_t = b$t) }), fill = TRUE)
print(swp[order(tau_deg), lapply(.SD, function(z) round(z,2))])

## ---- 5. do the other two outcomes come alive with two strikes? ---------
cat("\n\n=== SWING-DECISION RV AND MISS DISTANCE, TWO-STRIKE COUNTS ===\n")
cat("Two strikes is where strikeouts live, so if sd_rv was going to work anywhere\n")
cat("it should be here. Wanted signs: sd_rv negative, miss distance positive.\n\n")
oth <- rbindlist(list(
  cbind(outcome = "sd_rv_exp", subset = "two-strike",
        fe_fit(d[two_strike == TRUE], "sd_rv_exp", PRIMARY, c("plate_x","plate_z","release_speed","plate_sep","pp_x","pp_z"), "pit_count")),
  cbind(outcome = "sd_rv_exp", subset = "other counts",
        fe_fit(d[two_strike == FALSE], "sd_rv_exp", PRIMARY, c("plate_x","plate_z","release_speed","plate_sep","pp_x","pp_z"), "pit_count")),
  cbind(outcome = "miss_distance", subset = "two-strike",
        fe_fit(d[two_strike == TRUE & whiff == TRUE], "miss_distance", PRIMARY, c("plate_x","plate_z","release_speed","plate_sep","pp_x","pp_z"), "pit_count")),
  cbind(outcome = "miss_distance", subset = "other counts",
        fe_fit(d[two_strike == FALSE & whiff == TRUE], "miss_distance", PRIMARY, c("plate_x","plate_z","release_speed","plate_sep","pp_x","pp_z"), "pit_count"))
), fill = TRUE)
print(oth[, .(outcome, subset, n, arms, beta = round(beta,4),
              per_sd = round(per_sd,4), t = round(t,2), p = signif(p,3))])
