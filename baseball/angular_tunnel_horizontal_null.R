#!/usr/bin/env Rscript

# Why the metric does nothing on sliders and sweepers, and whether that is a real null or
# just an underpowered one.
#
# A confidence interval settles the difference. If the sweeper interval excludes the
# changeup effect size, the two pitch classes genuinely behave differently. If it contains
# it, the honest statement is "cannot tell" rather than "no effect".
#
# Also reported: the minimum effect each cell could have detected at 80% power, and the
# distribution of WHEN each pitch type becomes distinguishable, which is the mechanism
# candidate -- a sweeper diverges sideways early and unambiguously, so there may be no
# window of genuine ambiguity for the metric to grade.

suppressPackageStartupMessages(library(data.table))
options(width = 200)
PRIMARY <- "brk_any_005"
GRP <- list(Slider = "SL", Sweeper = "ST", Curveball = c("CU","KC"),
            Changeup = "CH", Splitter = "FS")

d <- readRDS("data/statcast_model/angular_tunnel_2026.rds")
cf <- list.files("data/statcast_2026/chunks", pattern = "csv$", full.names = TRUE)
sh <- unique(rbindlist(lapply(cf, function(f) fread(f, select = c(
  "game_pk","at_bat_number","pitch_number","sz_top","sz_bot","release_speed","game_type"),
  showProgress = FALSE)))[game_type == "R"], by = c("game_pk","at_bat_number","pitch_number"))
sh[, game_type := NULL]
d <- merge(d, sh, by = c("game_pk","at_bat_number","pitch_number"), all.x = TRUE)
d[, zdist := sqrt(pmax(abs(plate_x) - 0.95, 0)^2 + pmax(plate_z - sz_top, sz_bot - plate_z, 0)^2)]
d[, `:=`(out_zone = zdist > 0 & is.finite(zdist), chase = as.numeric(swing))]
d[, grp := NA_character_]
for (g in names(GRP)) d[pitch_type %in% GRP[[g]], grp := g]
d <- d[!is.na(grp)]
d[, brk := get(PRIMARY)]

fe_fit <- function(dt, yvar, xvar, controls, fe) {
  vars <- c(yvar, xvar, controls)
  dt <- dt[stats::complete.cases(dt[, vars, with = FALSE])]
  dm <- dt[, lapply(.SD, as.numeric), .SDcols = vars]
  dm[, `:=`(g2 = dt[[fe]], cl = dt$pitcher)]
  for (v in vars) dm[, (v) := get(v) - mean(get(v)), by = g2]
  X <- cbind(1, as.matrix(dm[, c(xvar, controls), with = FALSE])); y <- dm[[yvar]]
  XtXi <- solve(crossprod(X)); b <- XtXi %*% crossprod(X, y)
  e <- as.vector(y - X %*% b)
  meat <- matrix(0, ncol(X), ncol(X))
  for (ix in split(seq_len(nrow(X)), dm$cl)) {
    u <- crossprod(X[ix, , drop = FALSE], e[ix]); meat <- meat + u %*% t(u) }
  V <- XtXi %*% meat %*% XtXi; nc <- uniqueN(dm$cl); V <- V*(nc/(nc-1))
  sdx <- sd(dt[[xvar]])
  # everything expressed per SD of the metric, in percentage points
  data.table(n = nrow(dt), arms = nc, sd_x = sdx,
             est = 100*b[2]*sdx, se = 100*sqrt(V[2,2])*sdx)
}

LOC <- c("plate_x","plate_z","zdist","release_speed","balls","strikes")
OZ <- d[out_zone == TRUE]

cat("=== CHASE EFFECT WITH INTERVALS, PER SD OF THE TUNNEL METRIC ===\n")
cat("est and bounds in percentage points of chase. mde80 is the smallest effect the cell\n")
cat("could have detected at 80% power; anything below it is invisible to this sample.\n\n")
r <- rbindlist(lapply(names(GRP), function(g) {
  a <- fe_fit(OZ[grp == g], "chase", "brk", LOC, "pitcher")
  cbind(group = g, a) }), fill = TRUE)
r[, `:=`(lo = est - 1.96*se, hi = est + 1.96*se, mde80 = 2.80*se)]
print(r[, .(group, n, arms, sd_metric = round(sd_x,3), est = round(est,2),
            se = round(se,2), lo = round(lo,2), hi = round(hi,2),
            mde80 = round(mde80,2))], row.names = FALSE)

ch <- r[group == "Changeup"]$est; fs <- r[group == "Splitter"]$est
cat(sprintf("\nchangeup effect = %+.2f pp, splitter = %+.2f pp\n", ch, fs))
for (g in c("Slider","Sweeper")) {
  x <- r[group == g]
  cat(sprintf("  %s CI [%+.2f, %+.2f] -- excludes the changeup effect? %s | the splitter effect? %s\n",
              g, x$lo, x$hi, ifelse(ch < x$lo | ch > x$hi, "YES", "no"),
              ifelse(fs < x$lo | fs > x$hi, "YES", "no")))
}

cat("\n\n=== IS IT A VARIANCE PROBLEM? HOW MUCH THE METRIC MOVES BY PITCH TYPE ===\n")
cat("If the metric barely varies for a pitch type there is nothing for it to grade.\n\n")
v <- d[, .(n = .N, mean_brk = mean(brk), sd_brk = sd(brk),
           p10 = quantile(brk,.10), p90 = quantile(brk,.90),
           share_zero = 100*mean(brk <= 0.02),
           mean_overlap = mean(overlap_frac)), by = grp]
print(v[order(-sd_brk), lapply(.SD, function(z) if (is.numeric(z)) round(z,3) else z)],
      row.names = FALSE)

cat("\n\n=== WHEN EACH PITCH TYPE BECOMES DISTINGUISHABLE ===\n")
cat("share of pairs already separable at each point in the decision window (tau=0.05).\n")
cat("A pitch that is separable immediately has no ambiguity for the hitter to lose.\n\n")
qs <- d[, .(n = .N,
            sep_immediately = 100*mean(brk <= 0.05),
            by_quarter = 100*mean(brk <= 0.25),
            by_half = 100*mean(brk <= 0.50),
            survives_to_commit = 100*mean(brk >= 0.95)), by = grp]
print(qs[order(-sep_immediately), lapply(.SD, function(z) if (is.numeric(z)) round(z,1) else z)],
      row.names = FALSE)

cat("\n\n=== HOW THE OLD METRIC DOES ON THE SAME CELLS, FOR CONTRAST ===\n")
d[, tunnel_neg := -tunnel_old]; OZ <- d[out_zone == TRUE]
ro <- rbindlist(lapply(names(GRP), function(g) {
  a <- fe_fit(OZ[grp == g], "chase", "tunnel_neg", LOC, "pitcher")
  cbind(group = g, a) }), fill = TRUE)
ro[, `:=`(lo = est - 1.96*se, hi = est + 1.96*se)]
print(merge(r[, .(group, angular = round(est,2), ang_t = round(est/se,2))],
            ro[, .(group, old = round(est,2), old_t = round(est/se,2))], by = "group"),
      row.names = FALSE)
