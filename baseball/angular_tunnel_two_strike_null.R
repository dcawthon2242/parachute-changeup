#!/usr/bin/env Rscript

# The slider/sweeper null, restricted to two strikes with 3-2 excluded.
#
# 3-2 is dropped because it is a different decision problem: a ball ends the plate
# appearance, so the hitter protects harder than in any other count, and it showed the
# largest chase effect of any cell earlier. Keeping it in would let one unusual count carry
# the two-strike result.
#
# The point of the exercise is the same as the full-sample version: an interval, not just a
# p-value. Restricting to three counts costs roughly 60% of the sample, so the minimum
# detectable effect grows by about 1.6x and a null that was informative on the full sample
# may stop being informative here. That has to be stated rather than glossed.

suppressPackageStartupMessages(library(data.table))
options(width = 200)
PRIMARY <- "brk_any_005"
TS <- c("0-2","1-2","2-2")
GRP <- list(Slider = "SL", Sweeper = "ST", Curveball = c("CU","KC"),
            Changeup = "CH", Splitter = "FS")

d <- readRDS("data/statcast_model/angular_tunnel_2026.rds")
cf <- list.files("data/statcast_2026/chunks", pattern = "csv$", full.names = TRUE)
sh <- unique(rbindlist(lapply(cf, function(f) fread(f, select = c(
  "game_pk","at_bat_number","pitch_number","sz_top","sz_bot","release_speed",
  "plate_x","plate_z","game_type"), showProgress = FALSE)))[game_type == "R"],
  by = c("game_pk","at_bat_number","pitch_number"))
prev <- sh[, .(game_pk, at_bat_number, pitch_number = pitch_number + 1L,
               pp_x = plate_x, pp_z = plate_z)]
sh[, c("game_type","plate_x","plate_z") := NULL]
K <- c("game_pk","at_bat_number","pitch_number")
d <- merge(d, sh, by = K, all.x = TRUE)
d <- merge(d, prev, by = K, all.x = TRUE)
d[, plate_sep := sqrt((plate_x - pp_x)^2 + (plate_z - pp_z)^2)]
d[, zdist := sqrt(pmax(abs(plate_x) - 0.95, 0)^2 + pmax(plate_z - sz_top, sz_bot - plate_z, 0)^2)]
d[, `:=`(out_zone = zdist > 0 & is.finite(zdist), chase = as.numeric(swing),
         count = paste(balls, strikes, sep = "-"))]
d[, pit_count := paste(pitcher, count)]
d[, grp := NA_character_]
for (g in names(GRP)) d[pitch_type %in% GRP[[g]], grp := g]
d <- d[!is.na(grp)]
d[, `:=`(brk = get(PRIMARY), tunnel_neg = -tunnel_old)]

fe_fit <- function(dt, yvar, xvar, controls, fe) {
  vars <- c(yvar, xvar, controls)
  dt <- dt[stats::complete.cases(dt[, vars, with = FALSE])]
  if (nrow(dt) < 300 || uniqueN(dt$pitcher) < 15) return(NULL)
  dm <- dt[, lapply(.SD, as.numeric), .SDcols = vars]
  dm[, `:=`(g2 = dt[[fe]], cl = dt$pitcher)]
  for (v in vars) dm[, (v) := get(v) - mean(get(v)), by = g2]
  X <- cbind(1, as.matrix(dm[, c(xvar, controls), with = FALSE])); y <- dm[[yvar]]
  XtXi <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
  if (is.null(XtXi)) return(NULL)
  b <- XtXi %*% crossprod(X, y); e <- as.vector(y - X %*% b)
  meat <- matrix(0, ncol(X), ncol(X))
  for (ix in split(seq_len(nrow(X)), dm$cl)) {
    u <- crossprod(X[ix, , drop = FALSE], e[ix]); meat <- meat + u %*% t(u) }
  V <- XtXi %*% meat %*% XtXi; nc <- uniqueN(dm$cl); V <- V*(nc/(nc-1))
  sdx <- sd(dt[[xvar]])
  data.table(n = nrow(dt), arms = nc, sd_x = sdx,
             est = 100*b[2]*sdx, se = 100*sqrt(V[2,2])*sdx)
}

LOC <- c("plate_x","plate_z","zdist","release_speed","plate_sep","pp_x","pp_z")
OZ <- d[out_zone == TRUE]
TSD <- OZ[count %in% TS]

report <- function(dt, lab, fe) {
  r <- rbindlist(lapply(names(GRP), function(g) {
    a <- fe_fit(dt[grp == g], "chase", "brk", LOC, fe)
    o <- fe_fit(dt[grp == g], "chase", "tunnel_neg", LOC, fe)
    if (is.null(a)) return(NULL)
    cbind(group = g, base = 100*mean(dt[grp == g]$chase), a,
          old_est = if (is.null(o)) NA_real_ else o$est,
          old_t = if (is.null(o)) NA_real_ else o$est/o$se) }), fill = TRUE)
  r[, `:=`(t = est/se, lo = est - 1.96*se, hi = est + 1.96*se, mde80 = 2.80*se)]
  cat(sprintf("\n=== %s ===\n", lab))
  print(r[, .(group, n, arms, base = round(base,1), est = round(est,2), se = round(se,2),
              t = round(t,2), lo = round(lo,2), hi = round(hi,2), mde80 = round(mde80,2),
              old = round(old_est,2), old_t = round(old_t,2))], row.names = FALSE)
  r
}

cat("Effects are percentage points of chase per 1 SD of the tunnel metric.\n")
cat("mde80 = smallest effect detectable at 80% power in that cell.\n")
r_ts  <- report(TSD, "TWO STRIKES, 3-2 EXCLUDED (0-2 / 1-2 / 2-2), pitcher x count FE", "pit_count")
r_all <- report(OZ,  "ALL COUNTS, for reference", "pitcher")

## ---- is the two-strike null informative? -------------------------------
cat("\n\n=== IS THE TWO-STRIKE SLIDER/SWEEPER NULL INFORMATIVE? ===\n")
ch <- r_ts[group == "Changeup"]$est; fs <- r_ts[group == "Splitter"]$est
cat(sprintf("two-strike changeup effect %+.2f pp, splitter %+.2f pp\n\n", ch, fs))
for (g in c("Slider","Sweeper")) {
  x <- r_ts[group == g]
  cat(sprintf("%s: est %+.2f, CI [%+.2f, %+.2f], mde80 %.2f\n", g, x$est, x$lo, x$hi, x$mde80))
  cat(sprintf("   excludes the changeup effect? %s | the splitter effect? %s\n",
              ifelse(ch < x$lo | ch > x$hi, "YES", "NO - cannot rule it out"),
              ifelse(fs < x$lo | fs > x$hi, "YES", "NO - cannot rule it out")))
  cat(sprintf("   powered to see a changeup-sized effect? %s\n",
              ifelse(abs(ch) > x$mde80, "yes", "NO - underpowered for that size")))
}

## ---- pooled horizontal vs pooled offspeed, and a formal interaction ----
cat("\n\n=== POOLED HORIZONTAL (SL+ST) vs OFFSPEED (CH+FS), TWO STRIKES EX 3-2 ===\n")
TSD[, fam := fifelse(grp %in% c("Slider","Sweeper"), "horizontal",
              fifelse(grp %in% c("Changeup","Splitter"), "offspeed", "curveball"))]
for (f in c("horizontal","offspeed","curveball")) {
  a <- fe_fit(TSD[fam == f], "chase", "brk", LOC, "pit_count")
  cat(sprintf("  %-11s n=%-6s %+.2f pp  se %.2f  t=%.2f  CI [%+.2f, %+.2f]\n", f,
              format(a$n, big.mark=","), a$est, a$se, a$est/a$se,
              a$est - 1.96*a$se, a$est + 1.96*a$se))
}
cat("\nformal test of the difference: interaction of the metric with family\n")
sub <- TSD[fam %in% c("horizontal","offspeed")]
sub[, off := as.numeric(fam == "offspeed")]
sub[, brk_off := brk * off]
a <- fe_fit(sub, "chase", "brk_off", c("brk", LOC), "pit_count")
cat(sprintf("  interaction term: %+.2f pp per SD  se %.2f  t=%.2f  p=%.4f\n",
            a$est, a$se, a$est/a$se, 2*pt(-abs(a$est/a$se), a$arms - 1)))
cat("  (positive and significant means offspeed genuinely responds more than horizontal)\n")
