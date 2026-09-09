#!/usr/bin/env Rscript

# tjStuff+ v3.0 — Thomas Nestico's methodology and hyperparameters
# Source: Modelling tjStuff+ v3.0 (Medium / tjstats.ca) and
#         github.com/tnestico/tjstuff_plus/tj_stuff_plus_v3.ipynb
#
# Train 2020-2022 (his original window). Score 2025-2026.
# Target is count-conditional average run value from his lookup table,
# not the pitch's realized delta_run_exp.

suppressPackageStartupMessages({
  library(data.table)
  library(lightgbm)
})
set.seed(42)
options(width = 210)

FEATS <- c("start_speed","spin_rate","extension","az","ax","x0","z0",
           "speed_diff","az_diff","ax_diff")

# =============================================================================
# 1. Load Statcast and assign Nestico's target
# =============================================================================
cols <- c("game_year","game_type","pitcher","player_name","p_throws","pitch_type",
          "description","events","balls","strikes",
          "release_speed","release_spin_rate","release_extension",
          "release_pos_x","release_pos_z","ax","az")

load_year <- function(yr) {
  f <- file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr))
  if (!file.exists(f)) return(NULL)
  d <- fread(f, showProgress=FALSE, select=cols)
  d[game_type=="R" & pitch_type!="" & balls<=3 & strikes<=2]
}

cat("Loading Statcast...\n")
dt <- rbindlist(lapply(2020:2026, load_year), use.names=TRUE, fill=TRUE)
cat(sprintf("  %s regular-season pitches\n", format(nrow(dt), big.mark=",")))

# Map Statcast description / events onto Nestico's lookup keys
# Event (BIP / walk / K) wins when present; otherwise the pitch description.
des_map <- c(
  ball="ball", blocked_ball="ball", pitchout="ball",
  called_strike="called_strike",
  swinging_strike="swinging_strike", swinging_strike_blocked="swinging_strike",
  foul_tip="swinging_strike", missed_bunt="swinging_strike", bunt_foul_tip="swinging_strike",
  foul="foul", foul_bunt="foul",
  hit_by_pitch="hit_by_pitch"
)
evt_map <- c(
  single="single", double="double", triple="triple", home_run="home_run",
  field_out="field_out", force_out="field_out", grounded_into_double_play="field_out",
  fielders_choice_out="field_out", fielders_choice="field_out", double_play="field_out",
  sac_fly="field_out", sac_bunt="field_out", sac_fly_double_play="field_out",
  triple_play="field_out", other_out="field_out",
  walk="walk", hit_by_pitch="hit_by_pitch", strikeout="strikeout"
)
dt[, des_key := des_map[description]]
dt[, evt_key := evt_map[events]]

RV <- fread("baseball/tjstuff_run_values.csv")
setnames(RV, "delta_run_exp", "rv")
dt <- merge(dt, RV, by.x=c("evt_key","balls","strikes"),
            by.y=c("event","balls","strikes"), all.x=TRUE)
setnames(dt, "rv", "rv_evt")
dt <- merge(dt, RV, by.x=c("des_key","balls","strikes"),
            by.y=c("event","balls","strikes"), all.x=TRUE)
setnames(dt, "rv", "rv_des")
dt[, target := fifelse(is.finite(rv_evt), rv_evt, rv_des)]
dt <- dt[is.finite(target)]
cat(sprintf("  with target: %s\n", format(nrow(dt), big.mark=",")))

# =============================================================================
# 2. Feature engineering (notebook cell 6, verbatim)
# =============================================================================
dt[, `:=`(
  start_speed = release_speed,
  spin_rate   = release_spin_rate,
  extension   = release_extension,
  x0_raw      = release_pos_x,
  z0          = release_pos_z
)]
# LHP: flip ax. RHP: flip x0. Both end up on a righty-equivalent scale.
dt[, ax := fifelse(p_throws=="L", -ax, ax)]
dt[, x0 := fifelse(p_throws=="L", x0_raw, -x0_raw)]

FB <- dt[pitch_type %in% c("SI","FF","FC") & is.finite(start_speed) & is.finite(az) & is.finite(ax),
         .(n=.N, avg_fastball_speed=mean(start_speed),
           avg_fastball_az=mean(az), avg_fastball_ax=mean(ax)),
         by=.(pitcher, game_year, pitch_type)]
setorder(FB, pitcher, game_year, -n, -avg_fastball_speed)
PRIM <- FB[, .SD[1], by=.(pitcher, game_year)]
dt <- merge(dt, PRIM[, .(pitcher, game_year, avg_fastball_speed, avg_fastball_az, avg_fastball_ax)],
            by=c("pitcher","game_year"), all.x=TRUE)
# no-fastball fallback: fastest pitch that year
fbk <- dt[, .(fb_spd=max(start_speed, na.rm=TRUE),
              fb_az=max(az, na.rm=TRUE),
              fb_ax=max(ax, na.rm=TRUE)), by=.(pitcher, game_year)]
dt <- merge(dt, fbk, by=c("pitcher","game_year"))
dt[is.na(avg_fastball_speed), avg_fastball_speed := fb_spd]
dt[is.na(avg_fastball_az),    avg_fastball_az    := fb_az]
dt[is.na(avg_fastball_ax),    avg_fastball_ax    := fb_ax]
dt[, `:=`(
  speed_diff = start_speed - avg_fastball_speed,
  az_diff    = az - avg_fastball_az,
  ax_diff    = abs(ax - avg_fastball_ax)
)]
dt <- dt[is.finite(start_speed) & is.finite(spin_rate) & is.finite(extension) &
         is.finite(az) & is.finite(ax) & is.finite(x0) & is.finite(z0) &
         is.finite(speed_diff) & is.finite(az_diff) & is.finite(ax_diff)]
cat(sprintf("  with features: %s\n", format(nrow(dt), big.mark=",")))

# =============================================================================
# 3. RobustScaler + LightGBM — exact v3.0 hyperparameters
# =============================================================================
train <- dt[game_year %in% 2020:2022]
cat(sprintf("\nTraining pitches (2020-22): %s\n", format(nrow(train), big.mark=",")))

robust_fit <- function(X) {
  med <- apply(X, 2, median, na.rm=TRUE)
  iqr <- apply(X, 2, IQR, na.rm=TRUE)
  iqr[iqr==0] <- 1
  list(med=med, iqr=iqr)
}
robust_apply <- function(X, sc) sweep(sweep(X, 2, sc$med, "-"), 2, sc$iqr, "/")

Xtr <- as.matrix(train[, ..FEATS])
ytr <- train$target
sc  <- robust_fit(Xtr)
Xtr_s <- robust_apply(Xtr, sc)

dtrain <- lgb.Dataset(Xtr_s, label=ytr)
params <- list(
  objective         = "regression",
  metric            = "l2",
  learning_rate     = 0.01,
  num_leaves        = 31,
  max_depth         = -1,
  min_data_in_leaf  = 20,
  bagging_fraction  = 0.8,
  bagging_freq      = 1,
  feature_fraction  = 0.8,
  lambda_l1         = 0.1,
  lambda_l2         = 0.2,
  force_row_wise    = TRUE,
  verbosity         = -1,
  seed              = 42
)
cat("Fitting LGBMRegressor (n_estimators=1000, lr=0.01)...\n")
model <- lgb.train(params, dtrain, nrounds=1000)
imp <- lgb.importance(model)
cat("\nFeature importance (gain):\n")
print(imp)

# =============================================================================
# 4. Score 2025-2026 and convert to Stuff+
#    Stuff+ = 100 - 10 * z(predicted RV)   [higher = better for pitcher]
# =============================================================================
score <- dt[game_year %in% 2025:2026]
Xs <- robust_apply(as.matrix(score[, ..FEATS]), sc)
score[, pred := predict(model, Xs)]
# pitch-level mean/sd on the scored years, as Nestico does on the test year
mu <- mean(score$pred); sdv <- sd(score$pred)
score[, stuff_plus := 100 - 10*(pred - mu)/sdv]
cat(sprintf("\nScored %s pitches (2025-26). pred RV mean %.5f  sd %.5f\n",
            format(nrow(score), big.mark=","), mu, sdv))
cat(sprintf("  → 1 Stuff+ point = %.4f RV/100\n", 100*sdv/10))
cat(sprintf("  pitch-level Stuff+ mean %.1f  sd %.1f\n",
            mean(score$stuff_plus), sd(score$stuff_plus)))

# 20-80 pitch grade by type: 80 at p99.9, 20 at p0.1 of that type
score[, grade := {
  lo <- quantile(stuff_plus, 0.001, na.rm=TRUE)
  hi <- quantile(stuff_plus, 0.999, na.rm=TRUE)
  50 + 30*(stuff_plus - mean(stuff_plus))/(hi-lo)*2
}, by=pitch_type]

# =============================================================================
# 5. Leaderboards
# =============================================================================
NAME <- unique(score[, .(pitcher, player_name)], by="pitcher")

PT <- score[, .(
  pitches=.N,
  sl_velo=mean(start_speed),
  gap=-mean(speed_diff),          # FB minus this pitch
  stuff=mean(stuff_plus),
  grade=mean(grade),
  pred_rv100=100*mean(pred)
), by=.(game_year, pitcher, pitch_type)]
PT <- merge(PT, NAME, by="pitcher")

PIT <- score[, .(
  pitches=.N,
  stuff=mean(stuff_plus),
  pred_rv100=100*mean(pred)
), by=.(game_year, pitcher)]
PIT <- merge(PIT, NAME, by="pitcher")

cat("\n############ 2026 pitcher tjStuff+ (min 400 pitches) ############\n")
print(PIT[game_year==2026 & pitches>=400][order(-stuff)][1:15,
      .(player_name, pitches, stuff=round(stuff,1), pred_rv100=round(pred_rv100,2))])

cat("\n############ 2026 pitch-type leaders (min 150) ############\n")
for (pt in c("FF","SI","FC","SL","ST","CH","FS","CU")) {
  top <- PT[game_year==2026 & pitch_type==pt & pitches>=150][order(-stuff)][1:5]
  if (nrow(top)) {
    cat(sprintf("\n  %s\n", pt))
    print(top[, .(player_name, pitches, velo=round(sl_velo,1), gap=round(gap,1),
                  stuff=round(stuff,1), grade=round(grade,0))], row.names=FALSE)
  }
}

# =============================================================================
# 6. The slider question, now on a real Stuff+ model
# =============================================================================
cat("\n############ Slider Stuff+ vs velo vs gap (2025-26, min 80) ############\n\n")
SL <- PT[pitch_type=="SL" & pitches>=80]
cat(sprintf("  SL pitcher-seasons: %d\n", nrow(SL)))
cat(sprintf("  cor(Stuff+, SL velo) = %+.3f\n", cor(SL$stuff, SL$sl_velo)))
cat(sprintf("  cor(Stuff+, gap)     = %+.3f\n", cor(SL$stuff, SL$gap)))
m <- lm(stuff ~ sl_velo + gap, data=SL)
cat("\n  Stuff+ ~ SL velo + gap:\n")
print(round(summary(m)$coefficients, 4))
cat(sprintf("  R2 = %.3f\n", summary(m)$r.squared))

SL[, vQ := cut(sl_velo, quantile(sl_velo, 0:4/4), include.lowest=TRUE,
               labels=c("Q1 softest","Q2","Q3","Q4 hardest"))]
cat("\n  Stuff+ by slider-velocity quartile:\n")
print(SL[, .(n=.N, sl=round(mean(sl_velo),1), gap=round(mean(gap),1),
             stuff=round(mean(stuff),1)), by=vQ][order(vQ)], row.names=FALSE)

fwrite(PT[game_year==2026],  "data/statcast_2026/tjstuff_plus_2026_pitch.csv")
fwrite(PIT[game_year==2026], "data/statcast_2026/tjstuff_plus_2026_pitcher.csv")
saveRDS(list(model=model, scaler=sc, feats=FEATS, mu=mu, sdv=sdv),
        "data/statcast_2026/tjstuff_plus_v3_model.rds")
cat("\nWrote 2026 pitch-type and pitcher tables, plus the fitted model.\n")
