#!/usr/bin/env Rscript

# The chase result is the one place the angular metric clearly beats the old one, so it
# gets stress-tested before being believed.
#
# THE CONFOUND. Unlike the old path-to-location ratio, the break fraction is NOT divided
# by plate separation. So "long tunnel" and "the two pitches finished in nearly the same
# place" are entangled, and a breaking ball landing where a fastball just was could draw
# swings because the hitter expects that location, with nothing perceptual about it. If
# the effect is really about tunneling it has to survive controlling for where the pair
# finished.
#
# Ladder of specifications, each strictly harder than the last:
#   1  pitcher FE, location controls                       (the headline spec)
#   2  + plate separation between the two pitches
#   3  + the previous pitch's own location
#   4  pitcher x count FE, so counts are compared like for like
#   5  + previous pitch location, on 3-0/3-1/2-0 excluded (removes auto-take counts)
#   6  same as 4 but the outcome is restricted to pitches well outside the zone

suppressPackageStartupMessages(library(data.table))
options(width = 200)
PRIMARY <- "brk_any_005"

d <- readRDS("data/statcast_model/angular_tunnel_2026.rds")
cf <- list.files("data/statcast_2026/chunks", pattern = "csv$", full.names = TRUE)
ex <- rbindlist(lapply(cf, function(f) fread(f, select = c(
  "game_pk","at_bat_number","pitch_number","sz_top","sz_bot","release_speed",
  "plate_x","plate_z"), showProgress = FALSE)))
ex <- unique(ex, by = c("game_pk","at_bat_number","pitch_number"))
# The previous pitch is the setup fastball, which is NOT among the pair rows in d, so its
# location has to come from the full pitch table.
prev <- ex[, .(game_pk, at_bat_number, pitch_number = pitch_number + 1L,
               pp_x = plate_x, pp_z = plate_z)]
ex[, c("plate_x","plate_z") := NULL]
d <- merge(d, ex, by = c("game_pk","at_bat_number","pitch_number"), all.x = TRUE)
d <- merge(d, prev, by = c("game_pk","at_bat_number","pitch_number"), all.x = TRUE)
d[, plate_sep := sqrt((plate_x - pp_x)^2 + (plate_z - pp_z)^2)]

d[, dx := pmax(abs(plate_x) - 0.95, 0)]
d[, dz := pmax(plate_z - sz_top, sz_bot - plate_z, 0)]
d[, zdist := sqrt(dx^2 + dz^2)]
d[, out_zone := zdist > 0]
d[, chase := as.numeric(swing)]
d[, count := paste(balls, strikes, sep = "-")]

fe_fit <- function(dt, yvar, xvar, controls, fe) {
  vars <- c(yvar, xvar, controls)
  dt <- dt[stats::complete.cases(dt[, vars, with = FALSE])]
  if (nrow(dt) < 400) return(NULL)
  dm <- dt[, lapply(.SD, as.numeric), .SDcols = vars]
  dm[, grp := dt[[fe]]]
  dm[, cl := dt$pitcher]
  for (v in vars) dm[, (v) := get(v) - mean(get(v)), by = grp]
  X <- cbind(1, as.matrix(dm[, c(xvar, controls), with = FALSE])); y <- dm[[yvar]]
  XtXi <- solve(crossprod(X)); b <- XtXi %*% crossprod(X, y)
  e <- as.vector(y - X %*% b)
  meat <- matrix(0, ncol(X), ncol(X))
  for (ix in split(seq_len(nrow(X)), dm$cl)) {
    u <- crossprod(X[ix, , drop = FALSE], e[ix]); meat <- meat + u %*% t(u) }
  V <- XtXi %*% meat %*% XtXi; nc <- uniqueN(dm$cl); V <- V*(nc/(nc-1))
  data.table(n = nrow(dt), arms = nc, beta = b[2], se = sqrt(V[2,2]),
             t = b[2]/sqrt(V[2,2]), p = 2*pt(-abs(b[2]/sqrt(V[2,2])), nc-1),
             per_sd = b[2]*sd(dt[[xvar]]))
}

d[, pit_count := paste(pitcher, count)]
LOC <- c("plate_x","plate_z","zdist","release_speed")
OZ <- d[out_zone == TRUE]

specs <- list(
  list(nm = "1 pitcher FE + location",        dt = OZ, ctl = c(LOC,"balls","strikes"), fe = "pitcher"),
  list(nm = "2 + plate separation",           dt = OZ, ctl = c(LOC,"balls","strikes","plate_sep"), fe = "pitcher"),
  list(nm = "3 + previous pitch location",    dt = OZ, ctl = c(LOC,"balls","strikes","plate_sep","pp_x","pp_z"), fe = "pitcher"),
  list(nm = "4 pitcher x count FE",           dt = OZ, ctl = c(LOC,"plate_sep","pp_x","pp_z"), fe = "pit_count"),
  list(nm = "5 drop 3-0/3-1/2-0 counts",      dt = OZ[!count %in% c("3-0","3-1","2-0")],
       ctl = c(LOC,"plate_sep","pp_x","pp_z"), fe = "pit_count"),
  list(nm = "6 well outside zone (>3in)",     dt = OZ[zdist > 0.25],
       ctl = c(LOC,"plate_sep","pp_x","pp_z"), fe = "pit_count")
)

cat("=== CHASE, ALL SECONDARY PITCHES AFTER A FASTBALL, metric =", PRIMARY, "===\n")
cat("per_sd is the change in chase probability for a 1 SD longer tunnel.\n\n")
res <- rbindlist(lapply(specs, function(s) {
  r <- fe_fit(s$dt, "chase", PRIMARY, s$ctl, s$fe); if (!is.null(r)) r[, spec := s$nm]; r }), fill = TRUE)
print(res[, .(spec, n, arms, beta = round(beta,4), se = round(se,4),
              per_sd_pp = round(100*per_sd,2), t = round(t,2), p = signif(p,3))])

cat("\n=== SAME LADDER, THE OLD WEIGHTED-DISTANCE METRIC ===\n")
d[, tunnel_neg := -tunnel_old]; OZ <- d[out_zone == TRUE]
specs[[1]]$dt <- specs[[2]]$dt <- specs[[3]]$dt <- specs[[4]]$dt <- OZ
specs[[5]]$dt <- OZ[!count %in% c("3-0","3-1","2-0")]; specs[[6]]$dt <- OZ[zdist > 0.25]
ro <- rbindlist(lapply(specs, function(s) {
  r <- fe_fit(s$dt, "chase", "tunnel_neg", s$ctl, s$fe); if (!is.null(r)) r[, spec := s$nm]; r }), fill = TRUE)
print(ro[, .(spec, n, arms, beta = round(beta,5),
             per_sd_pp = round(100*per_sd,2), t = round(t,2), p = signif(p,3))])

cat("\n=== BY PITCH TYPE UNDER THE HARDEST SPEC (4) ===\n")
hard <- c(LOC,"plate_sep","pp_x","pp_z")
bt <- rbindlist(lapply(c("ST","SL","CU","KC","CH","FS"), function(ty) {
  a <- fe_fit(OZ[pitch_type == ty], "chase", PRIMARY, hard, "pit_count")
  b <- fe_fit(OZ[pitch_type == ty], "chase", "tunnel_neg", hard, "pit_count")
  if (is.null(a)) return(NULL)
  data.table(ptype = ty, n = a$n, ang_per_sd_pp = 100*a$per_sd, ang_t = a$t,
             old_per_sd_pp = if (is.null(b)) NA else 100*b$per_sd,
             old_t = if (is.null(b)) NA else b$t) }), fill = TRUE)
print(bt[, lapply(.SD, function(z) if (is.numeric(z)) round(z,2) else z)])

cat("\n=== THRESHOLD SWEEP UNDER THE HARDEST SPEC (4), pooled ===\n")
cat("This is where tau earns its keep: if the acuity idea is real there should be an\n")
cat("interior optimum rather than a monotone drift.\n\n")
tw <- rbindlist(lapply(grep("^brk_any_", names(d), value = TRUE), function(v) {
  r <- fe_fit(OZ, "chase", v, hard, "pit_count"); if (!is.null(r)) r[, metric := v]; r }), fill = TRUE)
tw[, tau_deg := as.numeric(sub("brk_any_", "", metric))/100]
print(tw[, .(tau_deg, n, beta = round(beta,4), per_sd_pp = round(100*per_sd,2),
             t = round(t,2), p = signif(p,3))][order(tau_deg)])

cat("\n=== DOES THE DEPTH CUE ADD ANYTHING OVER POSITION ALONE? (spec 4) ===\n")
cmp <- rbindlist(lapply(c("brk_pos_005","brk_any_005","brk_pos_010","brk_any_010"), function(v) {
  r <- fe_fit(OZ, "chase", v, hard, "pit_count"); if (!is.null(r)) r[, metric := v]; r }), fill = TRUE)
print(cmp[, .(metric, per_sd_pp = round(100*per_sd,2), t = round(t,2), p = signif(p,3))])
