#!/usr/bin/env Rscript

# The angular tunnel effect on offspeed: changeups and splitters.
#
# These were the strongest cells in the by-type breakdowns, so this is the place the
# hypothesis has the best chance. The same battery the breaking balls got is applied here:
# the six-step robustness ladder, the count grid, the acuity sweep, the grading model, and
# leaderboards.
#
# There is a mechanical reason to expect offspeed to do better. The metric's second cue is
# angular SIZE difference, which is a depth cue -- and depth is exactly what separates a
# changeup from a fastball. A sweeper diverges sideways, which the positional cue catches
# early and which no amount of acuity tuning will hide. A changeup stays on the fastball's
# line and arrives late, so the only early evidence is a subtle difference in how fast the
# ball is growing on the retina. If the angular formulation is doing real work rather than
# repackaging movement, offspeed is where it should show.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(11); options(width = 215)
PRIMARY <- "brk_any_005"
TS <- c("0-2","1-2","2-2")
GROUPS <- list(Changeup = "CH", Splitter = "FS")

d <- readRDS("data/statcast_model/angular_tunnel_2026.rds")
cf <- list.files("data/statcast_2026/chunks", pattern = "csv$", full.names = TRUE)
NEED <- c("game_pk","at_bat_number","pitch_number","pitcher","pitch_type","p_throws",
          "release_speed","release_spin_rate","release_extension","ax","az","spin_axis",
          "release_pos_x","release_pos_z","sz_top","sz_bot","delta_run_exp",
          "plate_x","plate_z","game_type","vx0","vy0","vz0","ay")
full <- unique(rbindlist(lapply(cf, function(f) fread(f, select = NEED, showProgress = FALSE)))[
  game_type == "R"], by = c("game_pk","at_bat_number","pitch_number"))
prev <- full[, .(game_pk, at_bat_number, pitch_number = pitch_number + 1L,
                 pp_x = plate_x, pp_z = plate_z)]
fb <- full[pitch_type %in% c("FF","SI","FC") & is.finite(release_speed),
           .(n = .N, fb_speed = mean(release_speed), fb_ax = mean(ax), fb_az = mean(az)),
           by = .(pitcher, pitch_type)]
setorder(fb, pitcher, -n); fb <- unique(fb, by = "pitcher")[, .(pitcher, fb_speed, fb_ax, fb_az)]

KEY <- c("game_pk","at_bat_number","pitch_number","pitcher")
newc <- setdiff(names(full), c(names(d), "game_type","pitch_type","p_throws"))
d <- merge(d, full[, c(KEY, newc), with = FALSE], by = KEY, all.x = TRUE)
d <- merge(d, prev, by = c("game_pk","at_bat_number","pitch_number"), all.x = TRUE)
d <- merge(d, fb, by = "pitcher", all.x = TRUE)
sdr <- readRDS("data/statcast_model/swing_decision_rv_2026.rds")[
  , .(game_pk, at_bat_number, pitch_number, sd_rv_exp)]
d <- merge(d, sdr, by = c("game_pk","at_bat_number","pitch_number"), all.x = TRUE)

d[, `:=`(speed_diff = release_speed - fb_speed, ax_diff = ax - fb_ax, az_diff = az - fb_az)]
d[, L := p_throws == "L"]
d[, `:=`(tj_x0 = fifelse(L, -release_pos_x, release_pos_x), tj_ax = fifelse(L, -ax, ax),
         tj_ax_diff = fifelse(L, -ax_diff, ax_diff),
         tj_axis = fifelse(L, (360 - spin_axis) %% 360, spin_axis),
         plate_x_arm = fifelse(L, -plate_x, plate_x))]
yf <- 17/12; y0 <- 50
d[, vy_f := -sqrt(pmax(vy0^2 - 2*ay*(y0 - yf), 0))][, tt := (vy_f - vy0)/ay]
d[, VAA := -atan((vz0 + az*tt)/vy_f)*180/pi]
d[, HAA_in := fifelse(L, -1, 1)*(-atan((vx0 + ax*tt)/vy_f)*180/pi)]
d[, `:=`(z_rel_bot = plate_z - sz_bot, z_rel_top = plate_z - sz_top,
         stand_R = as.integer(stand == "R"))]
d[, plate_sep := sqrt((plate_x - pp_x)^2 + (plate_z - pp_z)^2)]
d[, zdist := sqrt(pmax(abs(plate_x) - 0.95, 0)^2 + pmax(plate_z - sz_top, sz_bot - plate_z, 0)^2)]
d[, `:=`(out_zone = zdist > 0 & is.finite(zdist), chase = as.numeric(swing),
         whiff_n = as.numeric(whiff), count = paste(balls, strikes, sep = "-"))]
d[, `:=`(pit_count = paste(pitcher, count), two_strike = count %in% TS,
         tunnel_neg = -tunnel_old)]
d[, grp := fifelse(pitch_type == "CH", "Changeup",
            fifelse(pitch_type == "FS", "Splitter", NA_character_))]

fe_fit <- function(dt, yvar, xvar, controls, fe) {
  vars <- c(yvar, xvar, controls)
  dt <- dt[stats::complete.cases(dt[, vars, with = FALSE])]
  if (nrow(dt) < 300 || uniqueN(dt$pitcher) < 15) return(NULL)
  dm <- dt[, lapply(.SD, as.numeric), .SDcols = vars]
  dm[, `:=`(grp2 = dt[[fe]], cl = dt$pitcher)]
  for (v in vars) dm[, (v) := get(v) - mean(get(v)), by = grp2]
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

LOC <- c("plate_x","plate_z","zdist","release_speed")
OZ <- d[out_zone == TRUE & grp %in% c("Changeup","Splitter")]

## ---- 1. the robustness ladder -------------------------------------------
cat("=== 1. CHASE ROBUSTNESS LADDER, OFFSPEED ===\n")
cat("per_sd_pp = change in chase probability, percentage points, per 1 SD longer tunnel.\n")
lad <- list(
  list(nm = "1 pitcher FE + location",     ctl = c(LOC,"balls","strikes"), fe = "pitcher", f = function(x) x),
  list(nm = "2 + plate separation",        ctl = c(LOC,"balls","strikes","plate_sep"), fe = "pitcher", f = function(x) x),
  list(nm = "3 + previous pitch location", ctl = c(LOC,"balls","strikes","plate_sep","pp_x","pp_z"), fe = "pitcher", f = function(x) x),
  list(nm = "4 pitcher x count FE",        ctl = c(LOC,"plate_sep","pp_x","pp_z"), fe = "pit_count", f = function(x) x),
  list(nm = "5 drop 3-0/3-1/2-0",          ctl = c(LOC,"plate_sep","pp_x","pp_z"), fe = "pit_count",
       f = function(x) x[!count %in% c("3-0","3-1","2-0")]),
  list(nm = "6 well outside zone (>3in)",  ctl = c(LOC,"plate_sep","pp_x","pp_z"), fe = "pit_count",
       f = function(x) x[zdist > 0.25]))
for (g in c("Changeup","Splitter","POOLED")) {
  s0 <- if (g == "POOLED") OZ else OZ[grp == g]
  cat(sprintf("\n--- %s (base chase %.1f%%) ---\n", g, 100*mean(s0$chase)))
  r <- rbindlist(lapply(lad, function(L) {
    a <- fe_fit(L$f(s0), "chase", PRIMARY, L$ctl, L$fe)
    o <- fe_fit(L$f(s0), "chase", "tunnel_neg", L$ctl, L$fe)
    if (is.null(a)) return(NULL)
    data.table(spec = L$nm, n = a$n, ang_pp = 100*a$per_sd, ang_t = a$t, ang_p = a$p,
               old_pp = if (is.null(o)) NA_real_ else 100*o$per_sd,
               old_t = if (is.null(o)) NA_real_ else o$t) }), fill = TRUE)
  print(r[, .(spec, n, ang_pp = round(ang_pp,2), ang_t = round(ang_t,2),
              ang_p = signif(ang_p,3), old_pp = round(old_pp,2), old_t = round(old_t,2))],
        row.names = FALSE)
}

## ---- 2. counts ----------------------------------------------------------
cat("\n\n=== 2. BY COUNT, OFFSPEED POOLED (spec 1, pitcher FE) ===\n")
gr <- rbindlist(lapply(sort(unique(OZ$count)), function(cc) {
  s <- OZ[count == cc]; r <- fe_fit(s, "chase", PRIMARY, c(LOC,"plate_sep","pp_x","pp_z"), "pitcher")
  if (is.null(r)) return(data.table(count = cc, n = nrow(s), base = 100*mean(s$chase)))
  cbind(count = cc, base = 100*mean(s$chase), r) }), fill = TRUE)
gr[, strikes2 := as.integer(substr(count,3,3))]
print(gr[order(strikes2, count), .(count, n, base = round(base,1),
      per_sd_pp = round(100*per_sd,2), t = round(t,2), p = signif(p,3))], row.names = FALSE)
cat("\ntwo-strike vs rest, pitcher x count FE:\n")
for (lab in c(TRUE, FALSE)) {
  r <- fe_fit(OZ[two_strike == lab], "chase", PRIMARY, c(LOC,"plate_sep","pp_x","pp_z"), "pit_count")
  cat(sprintf("  %-14s n=%s base=%.1f%%  %+.2f pp  t=%.2f  p=%.2g\n",
              if (lab) "two-strike" else "other counts", format(r$n, big.mark=","),
              100*mean(OZ[two_strike == lab]$chase), 100*r$per_sd, r$t, r$p))
}

## ---- 3. acuity sweep, and whether the depth cue is what is working ------
cat("\n\n=== 3. ACUITY SWEEP, OFFSPEED (spec 4) ===\n")
cat("pos = positional cue only; any = position OR angular-size depth cue.\n")
cat("The gap between them is the depth cue's contribution.\n\n")
swp <- rbindlist(lapply(c(1,2,5,10,20,40,80), function(k) {
  va <- sprintf("brk_any_%03d", k); vp <- sprintf("brk_pos_%03d", k)
  a <- fe_fit(OZ, "chase", va, c(LOC,"plate_sep","pp_x","pp_z"), "pit_count")
  p <- fe_fit(OZ, "chase", vp, c(LOC,"plate_sep","pp_x","pp_z"), "pit_count")
  if (is.null(a)) return(NULL)
  data.table(tau_deg = k/100, any_pp = 100*a$per_sd, any_t = a$t,
             pos_pp = 100*p$per_sd, pos_t = p$t, depth_gain_t = a$t - p$t) }), fill = TRUE)
print(swp[, lapply(.SD, function(z) round(z,2))], row.names = FALSE)

## ---- 4. the other two outcomes -----------------------------------------
cat("\n\n=== 4. SWING-DECISION RV AND MISS DISTANCE, OFFSPEED ===\n")
cat("Wanted signs: sd_rv negative, miss distance positive.\n\n")
CTL2 <- c("plate_x","plate_z","release_speed","plate_sep","pp_x","pp_z")
oth <- rbindlist(lapply(c("Changeup","Splitter"), function(g) {
  s <- d[grp == g]
  rbindlist(list(
    cbind(grp = g, outcome = "sd_rv_exp", fe_fit(s, "sd_rv_exp", PRIMARY, CTL2, "pit_count")),
    cbind(grp = g, outcome = "miss_distance", fe_fit(s[whiff == TRUE], "miss_distance", PRIMARY, CTL2, "pit_count")),
    cbind(grp = g, outcome = "whiff | swing", fe_fit(s[swing == TRUE], "whiff_n", PRIMARY, CTL2, "pit_count"))
  ), fill = TRUE) }), fill = TRUE)
print(oth[, .(grp, outcome, n, arms, beta = round(beta,4), per_sd = round(per_sd,4),
              t = round(t,2), p = signif(p,3))], row.names = FALSE)

## ---- 5. grading model --------------------------------------------------
TJ  <- c("release_speed","release_spin_rate","release_extension","tj_ax","az","tj_x0",
         "release_pos_z","tj_axis","speed_diff","tj_ax_diff","az_diff")
LOCF <- c("plate_x_arm","plate_z","z_rel_bot","z_rel_top","VAA","HAA_in","stand_R","balls","strikes")
SETS <- list(A_tjstuff = TJ, B_plus_loc = c(TJ, LOCF),
             C_plus_angular = c(TJ, LOCF, "brk_any_005","overlap_frac"),
             D_plus_old = c(TJ, LOCF, "tunnel_old"))
auc <- function(y,p){ r <- rank(p); n1 <- sum(y==1); n0 <- sum(y==0)
  if (n1==0||n0==0) return(NA_real_); (sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0) }
oofit <- function(D, target, FEAT, binary = TRUE) {
  D <- D[stats::complete.cases(D[, c(FEAT, target), with = FALSE])]
  if (nrow(D) < 1200) return(NULL)
  K <- 4; fold <- sample(rep(1:K, length.out = nrow(D))); p <- rep(NA_real_, nrow(D))
  for (f in 1:K) {
    tr <- D[fold != f]; vi <- sample(nrow(tr), floor(.12*nrow(tr)))
    dtr <- lgb.Dataset(as.matrix(tr[-vi, ..FEAT]), label = tr[[target]][-vi])
    dva <- lgb.Dataset.create.valid(dtr, as.matrix(tr[vi, ..FEAT]), label = tr[[target]][vi])
    m <- lgb.train(params = list(objective = if (binary) "binary" else "regression",
                   metric = if (binary) "binary_logloss" else "l2", learning_rate = .06,
                   num_leaves = 31, min_data_in_leaf = 200, feature_fraction = .8,
                   bagging_fraction = .8, bagging_freq = 1,
                   num_threads = parallel::detectCores()),
                   data = dtr, nrounds = 1500, valids = list(v = dva),
                   early_stopping_rounds = 50, verbose = -1)
    p[fold == f] <- predict(m, as.matrix(D[fold == f, ..FEAT]))
  }
  y <- D[[target]]
  data.table(n = nrow(D), r2 = 1 - var(y-p)/var(y), auc = if (binary) auc(y,p) else NA_real_)
}
cat("\n\n=== 5. GRADING MODEL, OFFSPEED (out-of-fold, 4-fold) ===\n")
TG <- list(list(nm="whiff | swing", tg="whiff_n", sub=quote(swing==TRUE), bin=TRUE),
           list(nm="chase | out of zone", tg="chase", sub=quote(out_zone==TRUE), bin=TRUE),
           list(nm="run value", tg="delta_run_exp", sub=quote(rep(TRUE,.N)), bin=FALSE))
for (g in c("Changeup","Splitter")) for (T_ in TG) {
  S <- d[grp == g][eval(T_$sub)]
  b <- rbindlist(lapply(names(SETS), function(sn) {
    r <- oofit(S, T_$tg, SETS[[sn]], T_$bin); if (is.null(r)) return(NULL); cbind(set = sn, r)
  }), fill = TRUE)
  if (!nrow(b)) next
  base <- b[set == "B_plus_loc"]
  cat(sprintf("\n--- %s | %s (n=%s) ---\n", g, T_$nm, format(max(b$n), big.mark=",")))
  print(b[, .(set, r2 = round(r2,5), d_r2 = round(r2 - base$r2,5),
              auc = round(auc,5), d_auc = round(auc - base$auc,5))], row.names = FALSE)
}

## ---- 6. leaderboards ---------------------------------------------------
cat("\n\n=== 6. TOP 10 FASTBALL-TO-OFFSPEED TUNNELS, 2026 (min 25 pairs) ===\n")
for (g in c("Changeup","Splitter")) {
  s <- d[grp == g]
  a <- s[, .(pairs = .N, fb = names(sort(table(p_pitch_type), decreasing = TRUE))[1],
             brk = mean(brk_any_005), overlap = mean(overlap_frac),
             chase = 100*mean(swing[out_zone]), n_oz = sum(out_zone),
             whiff = 100*mean(whiff[swing]), n_sw = sum(swing), old = mean(tunnel_old)),
         by = .(pitcher, player_name, p_throws)][pairs >= 25]
  setorder(a, old); a[, old_rk := .I]; setorder(a, -brk); a[, rk := .I]
  a[, pct := round(100*(1 - (rk-1)/.N))]
  cat(sprintf("\n--- %s (%d qualified | league brk %.3f, chase %.1f%%, whiff %.1f%%) ---\n",
              g, nrow(a), mean(a$brk), 100*mean(s$swing[s$out_zone]), 100*mean(s$whiff[s$swing])))
  print(a[1:10, .(rk, pitcher = player_name, T = p_throws, fb, pairs, brk = round(brk,3),
                  pct, overlap = round(overlap,3), chase = round(chase,1), n_oz,
                  whiff = round(whiff,1), n_sw, old_rk)], row.names = FALSE)
  cat(sprintf("    rank agreement with old metric: Spearman %.3f\n",
              cor(a$rk, a$old_rk, method = "spearman")))
}
