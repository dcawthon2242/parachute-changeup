#!/usr/bin/env Rscript

# Do kick changes beat their STUFF-model expectation?
# Build a LightGBM stuff model (the same family used by Stuff+/PitchingBot) on pitch
# SHAPE ONLY -- velo, movement, spin, spin axis, release point, extension, and the
# diffs from the pitcher's own fastball -- with NO location, count, or pitch-type
# label. Take out-of-fold predictions, then compare kick changes vs the rest on:
#   (1) whiff over expected   (higher = beats stuff)
#   (2) run value over expected (lower = beats stuff; RV is runs allowed)

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
dir <- file.path("data","statcast_2026"); set.seed(42)

dt <- fread(file.path(dir,"statcast_2026_all.csv"), showProgress=FALSE, select=c(
  "player_name","pitcher","pitch_type","game_type","description","delta_run_exp",
  "release_speed","release_spin_rate","spin_axis","pfx_x","pfx_z",
  "release_pos_x","release_pos_z","release_extension",
  "vx0","vy0","vz0","ax","ay","az"))
dt <- dt[game_type=="R" & !is.na(release_speed) & !is.na(pfx_x) & !is.na(release_spin_rate)]

# VAA (vertical approach angle) at the plate -- a core stuff feature
yf <- 17/12
dt[, t := (-vy0 - sqrt(vy0^2 - 2*ay*(release_extension*0+50-yf)))/ay]
dt[, vy_f := vy0 + ay*t]; dt[, vz_f := vz0 + az*t]
dt[, vaa := -atan2(vz_f, -vy_f)*180/pi]

# fastball anchor per pitcher (FF>SI>FC), then shape diffs vs the fastball
fb <- dt[pitch_type %in% c("FF","SI","FC"),
  .(nfb=.N, fb_velo=mean(release_speed), fb_x=mean(pfx_x), fb_z=mean(pfx_z),
    fb_spin=mean(release_spin_rate)),
  by=.(pitcher, pitch_type)]
fb[, rank := fifelse(pitch_type=="FF",1,fifelse(pitch_type=="SI",2,3))]
fb <- fb[nfb>=50][order(pitcher, rank)][, .SD[1], by=pitcher][, .(pitcher, fb_velo, fb_x, fb_z, fb_spin)]
dt <- merge(dt, fb, by="pitcher", all.x=TRUE)
dt[, `:=`(velo_diff=fb_velo-release_speed, hb_diff=pfx_x-fb_x, ivb_diff=pfx_z-fb_z,
          spin_diff=fb_spin-release_spin_rate,
          sax=sin(spin_axis*pi/180), cax=cos(spin_axis*pi/180))]

feat <- c("release_speed","pfx_x","pfx_z","release_spin_rate","sax","cax",
          "release_pos_x","release_pos_z","release_extension","vaa",
          "velo_diff","hb_diff","ivb_diff","spin_diff")

whiff <- c("swinging_strike","swinging_strike_blocked","foul_tip","missed_bunt")
swing <- c(whiff,"foul","hit_into_play","foul_bunt","bunt_foul_tip")
dt[, is_whiff := as.integer(description %in% whiff)]
dt[, is_swing := description %in% swing]

# out-of-fold GBM predictions
oof <- function(d, target, feats, obj, nr=400){
  d <- d[complete.cases(d[, ..feats]) & !is.na(get(target))]
  X <- as.matrix(d[, ..feats]); y <- d[[target]]
  fold <- sample(rep(1:5, length.out=nrow(d)))
  pred <- numeric(nrow(d))
  for(k in 1:5){
    tr <- fold!=k; te <- fold==k
    m <- lgb.train(list(objective=obj, learning_rate=0.05, num_leaves=31,
      min_data_in_leaf=200, feature_fraction=0.8, bagging_fraction=0.8, bagging_freq=1, verbose=-1),
      lgb.Dataset(X[tr,], label=y[tr]), nrounds=nr)
    pred[te] <- predict(m, X[te,])
  }
  d[, exp := pred]; d[]
}

# tag kick changes (same fingerprint as before)
kn <- dt[pitch_type=="CH", .(n=.N, spin=mean(release_spin_rate), ivb=mean(pfx_z)*12), by=player_name][
  n>=15 & spin<=1500 & ivb<=3]$player_name
report <- function(d, val, better){
  d[, grp := fifelse(pitch_type!="CH","non-CH",
             fifelse(player_name %in% kn,"kick change","other CH"))]
  s <- d[, .(n=.N, actual=mean(get(val)), expected=mean(exp), resid=mean(get(val)-exp)), by=grp]
  s[, beats_stuff := if(better=="high") resid>0 else resid<0]
  s[order(grp)]
}

cat("=== STUFF MODEL 1: WHIFF over expected (among swings; higher resid = beats stuff) ===\n")
sw <- oof(dt[is_swing==TRUE], "is_whiff", feat, "binary")
r1 <- report(sw, "is_whiff", "high")
print(r1[, .(grp, swings=n, whiff=round(actual,3), xwhiff=round(expected,3),
  whiff_over_exp=round(resid,3), beats_stuff)])

cat("\n=== STUFF MODEL 2: RUN VALUE over expected (all pitches; RV=runs allowed, lower=better) ===\n")
rv <- oof(dt, "delta_run_exp", feat, "regression")
r2 <- report(rv, "delta_run_exp", "low")
print(r2[, .(grp, pitches=n, rv100=round(actual*100,2), xrv100=round(expected*100,2),
  rv_over_exp100=round(resid*100,2), beats_stuff)])

cat("\n=== Best individual kick changes vs their stuff expectation (whiff, min 25 swings) ===\n")
ind <- sw[player_name %in% kn & pitch_type=="CH", .(sw=.N, whiff=mean(is_whiff), xwhiff=mean(exp),
  over=mean(is_whiff-exp)), by=player_name][sw>=25][order(-over)]
print(head(ind[, .(player_name, swings=sw, whiff=round(whiff,3), xwhiff=round(xwhiff,3),
  whiff_over_exp=round(over,3))], 15))
cat("\nWorst:\n")
print(tail(ind[, .(player_name, swings=sw, whiff=round(whiff,3), xwhiff=round(xwhiff,3),
  whiff_over_exp=round(over,3))], 6))
