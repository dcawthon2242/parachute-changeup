#!/usr/bin/env Rscript

# Does the angular tunnel metric earn a place in a pitch grading model?
#
# The grading model is built in the house style of baseball/whiff_tjstuff.R: the
# tjStuff+ v3.0 eleven-feature set, mirrored for left-handers, scored out-of-fold with
# LightGBM. Feature sets are nested so the marginal contribution of the tunnel metric is
# attributable rather than assumed:
#
#   A  tjStuff+ v3.0                        shape only, no location -- Nestico's set
#   B  + location                           plate_x/z, zone-relative height, VAA/HAA, count
#   C  + angular tunnel                     break fraction at tau=0.05 and overlap fraction
#   D  + old weighted-distance tunnel       the same slot, filled with the previous metric
#
# C minus B is the answer to the question. D exists so the comparison is against the
# incumbent rather than against nothing.
#
# Three targets, because "grading" means different things depending on what is being
# grade: whiff on swing (the house standard), chase on out-of-zone pitches (where the
# metric was shown to work), and run value (what a grade ultimately has to move).
#
# Breaking balls only, grouped as the request specifies: Slider, Sweeper, and Curveball
# with knuckle-curve folded in.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(11); options(width = 215)
FASTBALLS <- c("FF","SI","FC")
GROUPS <- list(Slider = "SL", Sweeper = "ST", Curveball = c("CU","KC"))

d <- readRDS("data/statcast_model/angular_tunnel_2026.rds")
cf <- list.files("data/statcast_2026/chunks", pattern = "csv$", full.names = TRUE)
NEED <- c("game_pk","at_bat_number","pitch_number","pitcher","pitch_type","p_throws",
          "release_speed","release_spin_rate","release_extension","ax","az","spin_axis",
          "release_pos_x","release_pos_z","arm_angle","sz_top","sz_bot","delta_run_exp",
          "plate_x","plate_z","game_type","vx0","vy0","vz0","ay")
full <- rbindlist(lapply(cf, function(f) fread(f, select = NEED, showProgress = FALSE)))
full <- unique(full[game_type == "R"], by = c("game_pk","at_bat_number","pitch_number"))

# tjStuff+ differentials are taken against the pitcher's own primary fastball.
fb <- full[pitch_type %in% FASTBALLS & is.finite(release_speed),
           .(n = .N, fb_speed = mean(release_speed), fb_ax = mean(ax), fb_az = mean(az)),
           by = .(pitcher, pitch_type)]
setorder(fb, pitcher, -n)
fb <- unique(fb, by = "pitcher")[, .(pitcher, fb_speed, fb_ax, fb_az)]

KEY <- c("game_pk","at_bat_number","pitch_number","pitcher")
# The angular table already carries release point and the kinematic coefficients, so only
# the genuinely new columns are merged; otherwise the join silently renames ax to ax.x.
newc <- setdiff(names(full), c(names(d), "game_type", "pitch_type", "p_throws"))
d <- merge(d, full[, c(KEY, newc), with = FALSE], by = KEY, all.x = TRUE)
d <- merge(d, fb, by = "pitcher", all.x = TRUE)

d[, `:=`(speed_diff = release_speed - fb_speed, ax_diff = ax - fb_ax, az_diff = az - fb_az)]
d[, L := p_throws == "L"]
d[, `:=`(tj_x0 = fifelse(L, -release_pos_x, release_pos_x),
         tj_ax = fifelse(L, -ax, ax),
         tj_ax_diff = fifelse(L, -ax_diff, ax_diff),
         tj_axis = fifelse(L, (360 - spin_axis) %% 360, spin_axis),
         plate_x_arm = fifelse(L, -plate_x, plate_x))]
yf <- 17/12; y0 <- 50
d[, vy_f := -sqrt(pmax(vy0^2 - 2*ay*(y0 - yf), 0))]
d[, tt := (vy_f - vy0)/ay]
d[, VAA := -atan((vz0 + az*tt)/vy_f)*180/pi]
d[, HAA_in := fifelse(L, -1, 1) * (-atan((vx0 + ax*tt)/vy_f)*180/pi)]
d[, `:=`(z_rel_bot = plate_z - sz_bot, z_rel_top = plate_z - sz_top,
         stand_R = as.integer(stand == "R"))]
d[, zdist := sqrt(pmax(abs(plate_x) - 0.95, 0)^2 + pmax(plate_z - sz_top, sz_bot - plate_z, 0)^2)]
d[, `:=`(out_zone = zdist > 0, chase = as.numeric(swing), whiff_n = as.numeric(whiff))]
d[, grp := fifelse(pitch_type == "SL", "Slider",
            fifelse(pitch_type == "ST", "Sweeper",
            fifelse(pitch_type %in% c("CU","KC"), "Curveball", NA_character_)))]

TJ  <- c("release_speed","release_spin_rate","release_extension","tj_ax","az","tj_x0",
         "release_pos_z","tj_axis","speed_diff","tj_ax_diff","az_diff")
LOC <- c("plate_x_arm","plate_z","z_rel_bot","z_rel_top","VAA","HAA_in","stand_R","balls","strikes")
ANG <- c("brk_any_005","overlap_frac")
OLDT <- c("tunnel_old")
SETS <- list(A_tjstuff = TJ, B_plus_loc = c(TJ, LOC),
             C_plus_angular = c(TJ, LOC, ANG), D_plus_old = c(TJ, LOC, OLDT))

auc <- function(y, p) {
  r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[y == 1]) - n1*(n1+1)/2)/(n1*n0)
}
oof <- function(D, target, FEAT, binary = TRUE) {
  D <- D[stats::complete.cases(D[, c(FEAT, target), with = FALSE])]
  if (nrow(D) < 1500) return(NULL)
  K <- 4; fold <- sample(rep(1:K, length.out = nrow(D))); p <- rep(NA_real_, nrow(D))
  for (f in 1:K) {
    tr <- D[fold != f]; vi <- sample(nrow(tr), floor(.12*nrow(tr)))
    dtr <- lgb.Dataset(as.matrix(tr[-vi, ..FEAT]), label = tr[[target]][-vi])
    dva <- lgb.Dataset.create.valid(dtr, as.matrix(tr[vi, ..FEAT]), label = tr[[target]][vi])
    m <- lgb.train(params = list(objective = if (binary) "binary" else "regression",
                   metric = if (binary) "binary_logloss" else "l2",
                   learning_rate = .06, num_leaves = 31, min_data_in_leaf = 200,
                   feature_fraction = .8, bagging_fraction = .8, bagging_freq = 1,
                   num_threads = parallel::detectCores()),
                   data = dtr, nrounds = 1500, valids = list(v = dva),
                   early_stopping_rounds = 50, verbose = -1)
    p[fold == f] <- predict(m, as.matrix(D[fold == f, ..FEAT]))
  }
  y <- D[[target]]
  ll <- if (binary) -mean(y*log(pmax(p,1e-9)) + (1-y)*log(pmax(1-p,1e-9))) else NA_real_
  data.table(n = nrow(D), r2 = 1 - var(y - p)/var(y),
             auc = if (binary) auc(y, p) else NA_real_, logloss = ll)
}

TARGETS <- list(
  list(nm = "whiff | swing",       tg = "whiff_n",       sub = quote(swing == TRUE),    bin = TRUE),
  list(nm = "chase | out of zone", tg = "chase",         sub = quote(out_zone == TRUE), bin = TRUE),
  list(nm = "run value",           tg = "delta_run_exp", sub = quote(rep(TRUE, .N)),    bin = FALSE))

cat("=== PITCH GRADING MODEL: DOES THE TUNNEL METRIC ADD ANYTHING? ===\n")
cat("Out-of-fold, 4-fold. Nested feature sets; the row that matters is C minus B.\n")
res <- rbindlist(lapply(names(GROUPS), function(g) {
  rbindlist(lapply(TARGETS, function(T_) {
    S <- d[grp == g][eval(T_$sub)]
    rbindlist(lapply(names(SETS), function(sn) {
      r <- oof(S, T_$tg, SETS[[sn]], T_$bin)
      if (is.null(r)) return(NULL)
      cbind(group = g, target = T_$nm, set = sn, r)
    }), fill = TRUE)
  }), fill = TRUE)
}), fill = TRUE)

for (g in names(GROUPS)) for (tn in sapply(TARGETS, `[[`, "nm")) {
  b <- res[group == g & target == tn]
  if (!nrow(b)) next
  base <- b[set == "B_plus_loc"]
  cat(sprintf("\n--- %s | %s  (n=%s) ---\n", g, tn, format(max(b$n), big.mark = ",")))
  b[, d_auc := if (all(is.na(auc))) NA_real_ else auc - base$auc]
  b[, d_ll  := if (all(is.na(logloss))) NA_real_ else logloss - base$logloss]
  b[, d_r2  := r2 - base$r2]
  print(b[, .(set, n, r2 = round(r2,5), d_r2 = round(d_r2,5),
              auc = round(auc,5), d_auc = round(d_auc,5),
              logloss = round(logloss,5), d_ll = round(d_ll,5))], row.names = FALSE)
}
saveRDS(res, "data/statcast_model/angular_tunnel_grading.rds")

## ---- feature importance: where does the tunnel metric rank? -------------
cat("\n\n=== WHERE THE TUNNEL METRIC RANKS AMONG THE FEATURES (chase model) ===\n")
for (g in names(GROUPS)) {
  S <- d[grp == g & out_zone == TRUE]
  FEAT <- SETS$C_plus_angular
  S <- S[stats::complete.cases(S[, c(FEAT,"chase"), with = FALSE])]
  m <- lgb.train(params = list(objective = "binary", learning_rate = .06, num_leaves = 31,
                 min_data_in_leaf = 200, feature_fraction = .8, bagging_fraction = .8,
                 bagging_freq = 1), verbose = -1, nrounds = 400,
                 data = lgb.Dataset(as.matrix(S[, ..FEAT]), label = S$chase))
  imp <- as.data.table(lgb.importance(m))[order(-Gain)]
  imp[, rank := .I]
  cat(sprintf("\n%s (n=%s): brk_any_005 ranks %d/%d with %.1f%% of gain; overlap_frac ranks %d\n",
              g, format(nrow(S), big.mark = ","),
              imp[Feature == "brk_any_005"]$rank, nrow(imp),
              100*imp[Feature == "brk_any_005"]$Gain,
              imp[Feature == "overlap_frac"]$rank))
  print(head(imp[, .(rank, feature = Feature, gain_pct = round(100*Gain,1))], 6), row.names = FALSE)
}
