#!/usr/bin/env Rscript

# PARACHUTE CHANGEUPS: TWO PRE-SPECIFIED HYPOTHESES THAT THE MISS-DISTANCE TEST COULD NOT SEE.
#
# The direct test came back null: parachute changeups sit +0.029 inches of miss distance above
# a shape-and-location model, p = 0.49. Rather than re-cut the bin until something appears,
# this keeps the bin FROZEN exactly as it was - mean axis gap <= 15 deg from the primary
# fastball, IVB kill >= 13.9, min 60 swings - and changes only the OUTCOME, because there are
# two specific ways the earlier test could have been blind rather than the pitch being ordinary.
#
# H1  WRONG OUTCOME. Miss distance conditions on a swing and scores contact as zero. A pitch
#     whose job is weak contact instead of whiffs is invisible to it and would look average
#     while suppressing runs. Test run value per 100 pitches, xwOBA on contact, ground-ball
#     rate and exit velocity, each against a shape-and-location model of the same quantity.
#
# H2  WRONG UNIT. A changeup that mirrors the fastball's spin may pay off on the FASTBALL, by
#     making it harder to sit on. No per-pitch model of the changeup can see that. Test
#     whether a pitcher's four-seam and sinker overperform when he throws a parachute changeup.
#
# Both were written down before running. Whatever they return is the answer.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(1); options(width = 205)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")
FASTBALLS <- c("FF","SI","FC"); yf <- 17/12; g <- 32.174
CACHE <- file.path(MDIR, "parachute_rv.rds")

if (!file.exists(CACHE) || nzchar(Sys.getenv("REFIT"))) {
COLS <- c("game_pk","at_bat_number","pitch_number","game_date","game_type","pitcher",
          "player_name","pitch_type","p_throws","stand","balls","strikes","release_speed",
          "release_spin_rate","release_extension","release_pos_x","release_pos_y",
          "release_pos_z","spin_axis","plate_x","plate_z","sz_bot","sz_top","vx0","vy0","vz0",
          "ax","ay","az","description","events","bb_type","launch_speed","launch_angle",
          "estimated_woba_using_speedangle","delta_run_exp","delta_pitcher_run_exp")
d <- rbindlist(lapply(2023:2026, function(y) {
  k <- fread(sprintf("data/statcast_%d/statcast_%d_all.csv", y, y), select = COLS,
             showProgress = FALSE); k[, season := y]; k }))
d <- d[game_type == "R" & !is.na(vx0) & pitch_type != "" & !is.na(release_pos_y)]
d[season == 2023L, game_date := as.Date(game_date)]
d <- d[season != 2023L | game_date >= as.Date("2023-07-14")]
d <- unique(d, by = c("game_pk","at_bat_number","pitch_number"))

# Primary fastball anchor per pitcher-season, FF > SI > FC, as used throughout.
cmean <- function(a) { r <- a*pi/180; (atan2(mean(sin(r)), mean(cos(r)))*180/pi) %% 360 }
circd <- function(a,b) { z <- abs(a-b) %% 360; pmin(z, 360-z) }
fb <- d[pitch_type %in% FASTBALLS, .(nfb = .N, fb_velo = mean(release_speed, na.rm=TRUE),
        fb_ax = mean(ax, na.rm=TRUE), fb_az = mean(az, na.rm=TRUE),
        fb_axis = cmean(spin_axis[!is.na(spin_axis)])), by = .(pitcher, season, pitch_type)][nfb >= 50]
fb[, rk := match(pitch_type, FASTBALLS)]
fb <- fb[order(pitcher, season, rk)][, .SD[1], by = .(pitcher, season)][
  , .(pitcher, season, fb_velo, fb_ax, fb_az, fb_axis)]
d <- merge(d, fb, by = c("pitcher","season"), all.x = TRUE)
d[, `:=`(speed_diff = release_speed - fb_velo, ax_diff = ax - fb_ax, az_diff = az - fb_az,
         axis_diff = circd(spin_axis, fb_axis))]

d[, vyf := -sqrt(pmax(vy0^2 - 2*ay*(50-yf), 0))][, tf := (vyf - vy0)/ay]
d[, `:=`(VAA = atan2(vz0 + az*tf, vyf)*180/pi, HAA = atan2(vx0 + ax*tf, vyf)*180/pi,
         below_zone = as.integer(plate_z < sz_bot), z_rel_bot = plate_z - sz_bot,
         z_rel_top = plate_z - sz_top, sax = sin(spin_axis*pi/180), cax = cos(spin_axis*pi/180))]
d[, `:=`(bs = fifelse(stand == "R", -1, 1), ps = fifelse(p_throws == "R", -1, 1))]
d[, `:=`(plate_x_in = bs*plate_x, HAA_in = bs*HAA, plate_x_arm = ps*plate_x,
         stand_R = as.integer(stand == "R"), throws_R = as.integer(p_throws == "R"),
         same_hand = as.integer(stand == p_throws))]
# Run value from the pitcher's side: positive = good for the pitcher.
d[, rv := -delta_run_exp]
d <- d[is.finite(rv)]
saveRDS(d[pitch_type %in% c("CH", FASTBALLS)], CACHE)
} else d <- readRDS(CACHE)

BASE <- c("release_speed","release_spin_rate","release_extension","release_pos_x",
          "release_pos_z","sax","cax","ax","az","speed_diff","ax_diff","az_diff")
LOC  <- c("plate_x","plate_z","below_zone","VAA","HAA","z_rel_bot","plate_x_in",
          "plate_x_arm","HAA_in","z_rel_top","stand_R","throws_R","same_hand")
CNT  <- c("balls","strikes")

# axis_diff is excluded from the features: it defines the bin under test.
oof <- function(D, tag) {
  FEAT <- c(BASE, LOC, CNT); K <- 4
  D[, fold := sample(rep(1:K, length.out = .N))]
  p <- rep(NA_real_, nrow(D))
  for (f in 1:K) {
    tr <- D[fold != f]; n <- nrow(tr); vi <- sample(n, floor(.12*n))
    dtr <- lgb.Dataset(as.matrix(tr[-vi, ..FEAT]), label = tr$rv[-vi])
    dva <- lgb.Dataset.create.valid(dtr, as.matrix(tr[vi, ..FEAT]), label = tr$rv[vi])
    m <- lgb.train(params = list(objective="regression", metric="l2", learning_rate=.06,
                   num_leaves=31, min_data_in_leaf=400, feature_fraction=.8,
                   bagging_fraction=.8, bagging_freq=1), data = dtr, nrounds = 1200,
                   valids = list(val = dva), early_stopping_rounds = 50, verbose = -1)
    p[D$fold == f] <- predict(m, as.matrix(D[fold == f, ..FEAT]))
  }
  cat(sprintf("  %s run-value model: n=%d  R2=%.4f\n", tag, nrow(D), 1 - var(D$rv - p)/var(D$rv)))
  D$rv - p
}

CH <- d[pitch_type == "CH" & stats::complete.cases(d[pitch_type == "CH", c(BASE, LOC), with = FALSE])]
cat("=== out-of-fold run-value models (shape + location + count, axis_diff excluded) ===\n")
CH[, rv_res := oof(CH, "changeup")]

## ---- the FROZEN bin ------------------------------------------------------------
sw <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip","hit_into_play",
        "foul_bunt","missed_bunt","bunt_foul_tip")
CH[, `:=`(is_swing = description %in% sw,
          is_bip = !is.na(launch_speed) & bb_type != "",
          gb = bb_type == "ground_ball")]
S <- CH[, .(np = .N, nsw = sum(is_swing), axis = mean(axis_diff, na.rm = TRUE),
            ivb_kill = -mean(az_diff, na.rm = TRUE), spin = mean(release_spin_rate, na.rm = TRUE),
            rv100 = 100*mean(rv), rv_res100 = 100*mean(rv_res),
            xwobacon = mean(estimated_woba_using_speedangle[is_bip], na.rm = TRUE),
            gbpct = 100*mean(gb[is_bip], na.rm = TRUE),
            ev = mean(launch_speed[is_bip], na.rm = TRUE)),
        by = .(pitcher, player_name, season)]
S <- S[nsw >= 60 & is.finite(axis) & is.finite(ivb_kill)]
S[, para := axis <= 15 & ivb_kill >= 13.9]      # FROZEN, identical to the miss-distance test
cat(sprintf("\nfrozen bin: %d parachute pitcher-seasons of %d\n", sum(S$para), nrow(S)))

tst <- function(v, lab, digits = 3) {
  a <- S[para == TRUE][[v]]; b <- S[para == FALSE][[v]]
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  t <- t.test(a, b)
  cat(sprintf("  %-26s parachute %8.*f   rest %8.*f   diff %+.*f   p=%.3f\n",
              lab, digits, mean(a), digits, mean(b), digits, mean(a)-mean(b), t$p.value))
}
cat("\n=== H1: is the pitch effective in a way miss distance cannot see? ===\n")
tst("rv100",     "Run value /100 (raw)")
tst("rv_res100", "Run value /100 vs model")
tst("xwobacon",  "xwOBA on contact")
tst("gbpct",     "Ground-ball %", 1)
tst("ev",        "Exit velocity (mph)", 1)

## ---- H2: does it make the fastball better? --------------------------------------
FBd <- d[pitch_type %in% c("FF","SI") &
         stats::complete.cases(d[pitch_type %in% c("FF","SI"), c(BASE, LOC), with = FALSE])]
FBd[, rv_res := oof(FBd, "fastball")]
FB <- FBd[, .(nfb = .N, fb_rv_res100 = 100*mean(rv_res)), by = .(pitcher, season)]
H2 <- merge(S[, .(pitcher, season, para, np)], FB, by = c("pitcher","season"))[nfb >= 300]
cat(sprintf("\n=== H2: does the fastball play up? (%d pitcher-seasons, %d with a parachute CH) ===\n",
            nrow(H2), sum(H2$para)))
t2 <- t.test(H2[para == TRUE]$fb_rv_res100, H2[para == FALSE]$fb_rv_res100)
cat(sprintf("  FF/SI run value /100 vs model:  parachute %+.3f   rest %+.3f   diff %+.3f   p=%.3f\n",
            mean(H2[para == TRUE]$fb_rv_res100), mean(H2[para == FALSE]$fb_rv_res100),
            diff(rev(t2$estimate)), t2$p.value))

fwrite(S[order(-rv_res100)], file.path(AST, "ext_parachute_run_value.csv"))
cat("\n=== the frozen bin, by run value over model ===\n")
print(S[para == TRUE][order(-rv_res100), .(player_name, season, np, axis = round(axis,1),
        rv100 = round(rv100,2), rv_res100 = round(rv_res100,2),
        xwobacon = round(xwobacon,3), gbpct = round(gbpct,1))], row.names = FALSE)
cat("\nwrote ext_parachute_run_value.csv\n")
