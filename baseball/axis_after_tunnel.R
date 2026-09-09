#!/usr/bin/env Rscript

# DOES AXIS SIMILARITY SURVIVE PUTTING THE TUNNEL METRIC IN THE MODEL ITSELF?
#
# Partialling path_ratio out of the pitcher-level correlation afterwards is a weak control:
# it only removes a linear, pitcher-averaged version of the tunnel. The real test is to
# give the whiff model path_ratio as a feature so it can price the tunnel per pitch and
# non-linearly, regenerate the out-of-fold residual, and ask whether axis similarity still
# explains what is left.
#
# Nested variants, all trained on breaking balls only so the comparison is internal:
#   M0  base shape features                     (reproduces the current residual)
#   M1  + path_ratio                            (tunnel priced in)
#   M2  + location & approach angle             (location priced in, replacing the post-hoc lm)
#   M3  + path_ratio + location                 (the strongest control)
#   M4  M3 + axis_diff itself                   (is the cue predictive, not just correlated?)
#
# M4 turns the question from "does the correlation survive" into "does knowing the axis gap
# make the model better at calling whiffs", which is the claim that actually matters.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(1)
options(width = 200)
MDIR <- file.path("data","statcast_model")
BASE <- readRDS(file.path(MDIR, "arch_compare.rds"))$FEAT
# Location is meaningless without knowing who is standing where. Raw plate_x cannot
# distinguish inside to a righty from outside to a lefty, so the block carries the pitch
# in three frames - umpire (raw), batter (inside positive) and pitcher (arm side positive)
# - plus both handednesses and the matchup, and the zone edges top and bottom.
LOC  <- c("plate_x","plate_z","below_zone","VAA","HAA","z_rel_bot",
          "plate_x_in","plate_x_arm","HAA_in","z_rel_top","stand_R","throws_R","same_hand")

GRP <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(GRP)) GRP <- "breaking"
GRPS <- if (GRP == "breaking") {
  list(sliders = "SL", sweepers = "ST", curves = c("CU","KC"))
} else {
  list(changeups = "CH", splitters = "FS")
}
cat("### group:", GRP, "\n\n")

d <- readRDS(file.path(MDIR, "miss_grade_data_activespin.rds"))
# Active-spin gap to the primary fastball, so the offspeed run can test the efficiency
# cue alongside the axis cue rather than only the axis cue.
as_long <- readRDS(file.path(MDIR, "active_spin_long.rds")); setDT(as_long)
fbA <- as_long[pitch_type %in% c("FF","SI","FC")][
  , pr := match(pitch_type, c("FF","SI","FC"))][order(pitcher, season, pr)][
  , .SD[1], by = .(pitcher, season)][, .(pitcher, season, fb_active = active_spin)]
d <- merge(d, fbA, by = c("pitcher","season"), all.x = TRUE)
d[, as_gap := active_spin - fb_active]
d[, spin_sim := exp(-(abs(as_gap)/0.10)^2) * exp(-(axis_diff/45)^2)]

d <- d[grp == GRP & !is.na(is_whiff)]

kin <- rbindlist(lapply(2023:2026, function(y)
  fread(sprintf("data/statcast_%d/statcast_%d_all.csv", y, y),
        select = c("game_pk","at_bat_number","pitch_number","plate_x","plate_z","sz_bot",
                   "sz_top","stand","p_throws","vx0","vy0","vz0","ax","ay","az"),
        showProgress = FALSE)))
kin <- unique(kin, by = c("game_pk","at_bat_number","pitch_number"))
yf <- 17/12
kin[, vyf := -sqrt(pmax(vy0^2 - 2*ay*(50-yf), 0))][, tf := (vyf - vy0)/ay]
kin[, `:=`(VAA = atan2(vz0 + az*tf, vyf)*180/pi, HAA = atan2(vx0 + ax*tf, vyf)*180/pi,
           below_zone = as.integer(plate_z < sz_bot), z_rel_bot = plate_z - sz_bot,
           z_rel_top = plate_z - sz_top)]
# Verified against hit-by-pitch: mean plate_x is -1.92 for RHB and +1.98 for LHB, so a
# righty's inside is negative plate_x. Flipping by batter side puts inside positive for
# both. Flipping by pitcher hand puts the pitcher's arm side positive for both.
kin[, `:=`(bs = fifelse(stand == "R", -1, 1), ps = fifelse(p_throws == "R", -1, 1))]
kin[, `:=`(plate_x_in  = bs * plate_x, HAA_in = bs * HAA,
           plate_x_arm = ps * plate_x,
           stand_R = as.integer(stand == "R"), throws_R = as.integer(p_throws == "R"),
           same_hand = as.integer(stand == p_throws))]
d <- merge(d, kin[, c("game_pk","at_bat_number","pitch_number", LOC), with = FALSE],
           by = c("game_pk","at_bat_number","pitch_number"), all.x = TRUE)

D <- d[stats::complete.cases(d[, c(BASE, LOC), with = FALSE]) & is.finite(axis_diff)]
cat(GRP, "pitches:", nrow(D), " path_ratio present on",
    round(100*mean(!is.na(D$path_ratio)),1), "%,  active-spin gap on",
    round(100*mean(is.finite(D$as_gap)),1), "%\n\n")

K <- 4; D[, fold := sample(rep(1:K, length.out = .N))]
oof <- function(FEAT, tag) {
  p <- rep(NA_real_, nrow(D))
  for (f in 1:K) {
    tr <- D[fold != f]; n <- nrow(tr); vi <- sample(n, floor(0.12*n))
    dtr <- lgb.Dataset(as.matrix(tr[-vi, ..FEAT]), label = tr$is_whiff[-vi])
    dva <- lgb.Dataset.create.valid(dtr, as.matrix(tr[vi, ..FEAT]), label = tr$is_whiff[vi])
    mf <- lgb.train(params = list(objective="binary", metric="binary_logloss",
                    learning_rate=0.06, num_leaves=31, min_data_in_leaf=200,
                    feature_fraction=0.8, bagging_fraction=0.8, bagging_freq=1),
                    data = dtr, nrounds = 1200, valids = list(val = dva),
                    early_stopping_rounds = 50, verbose = -1)
    p[D$fold == f] <- predict(mf, as.matrix(D[fold == f, ..FEAT]))
  }
  y <- D$is_whiff
  ll <- -mean(y*log(pmax(p,1e-9)) + (1-y)*log(pmax(1-p,1e-9)))
  r <- rank(p); n1 <- as.numeric(sum(y == 1)); n0 <- as.numeric(sum(y == 0))
  auc <- (sum(r[y == 1]) - n1*(n1+1)/2) / (n1*n0)
  cat(sprintf("  %-3s  %2d feats   logloss %.5f   AUC %.5f\n", tag, length(FEAT), ll, auc))
  list(p = p, ll = ll, auc = auc)
}

cat("=== OUT-OF-FOLD WHIFF MODELS (", GRP, ") ===\n", sep = "")
M <- list(
  M0 = oof(BASE,                              "M0"),
  M1 = oof(c(BASE, "path_ratio"),             "M1"),
  M2 = oof(c(BASE, LOC),                      "M2"),
  M3 = oof(c(BASE, "path_ratio", LOC),        "M3"),
  M4 = oof(c(BASE, "path_ratio", LOC, "axis_diff"), "M4"))

## ---- does axis similarity still explain the leftover? -------------------------
detrend <- function(res) {   # the post-hoc location adjustment used by Fig 11b
  residuals(lm(res ~ poly(plate_x,3)*poly(plate_z,3) + below_zone + VAA + HAA, data = D))
}
sp <- function(x, y) { ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 12) return("       (too few)")
  ct <- suppressWarnings(cor.test(x, y, method="spearman", exact=FALSE))
  sprintf("%+.3f (p=%-8.2g n=%d)", ct$estimate, ct$p.value, length(x)) }

# Cues under test. Breaking balls only ever had the axis cue; offspeed also had the
# active-spin match and the old spin_sim kernel, so all three get carried through.
CUES <- list(axis_sim = quote(-mean(axis_diff)))
if (GRP == "offspeed") CUES <- c(CUES, list(act_sim  = quote(-mean(abs(as_gap), na.rm = TRUE)),
                                            spin_sim = quote(mean(spin_sim, na.rm = TRUE))))
test <- function(res, types, cue) {
  D[, .r := res]
  A <- D[pitch_type %in% types, .(n = .N, cue = eval(CUES[[cue]]), r = 100*mean(.r)),
         by = pitcher][n >= 200]
  sp(A$cue, A$r)
}

lbl <- c(M0="base shape", M1="+ tunnel", M2="+ location", M3="+ tunnel + location",
         M4="+ tunnel + location + axis_diff")
rows <- rbindlist(lapply(names(CUES), function(cue) rbindlist(lapply(names(M), function(nm) {
  res <- D$is_whiff - M[[nm]]$p
  dres <- detrend(res)
  cbind(data.table(cue = cue, model = nm, what = lbl[nm]),
        as.data.table(setNames(lapply(names(GRPS), function(g)
          test(if (nm %in% c("M0","M1")) dres else res, GRPS[[g]], cue)), names(GRPS))))
}))))

cat("\n=== CUE vs THE LEFTOVER RESIDUAL ===\n")
cat("    pitcher level, min 200 pitches; positive = more similar to the fastball helps.\n")
cat("    M0/M1 use the old post-hoc location detrend; M2-M4 have location in the model.\n\n")
for (cu in names(CUES)) { cat("--", cu, "--\n")
  print(rows[rows$cue == cu][, !"cue"], row.names = FALSE); cat("\n") }

cat("\n=== IS THE CUE PREDICTIVE, NOT JUST CORRELATED? (M4 vs M3) ===\n")
cat(sprintf("  logloss  %.5f -> %.5f  (%+.5f)\n", M$M3$ll, M$M4$ll, M$M4$ll - M$M3$ll))
cat(sprintf("  AUC      %.5f -> %.5f  (%+.5f)\n", M$M3$auc, M$M4$auc, M$M4$auc - M$M3$auc))
cat("  Adding the axis gap to a model that already knows shape, tunnel and location either\n",
    "  helps out of fold or it does not. This is the cleanest version of the question.\n", sep="")

out <- sprintf("axis_after_tunnel_%s.rds", GRP)
saveRDS(list(models = rows, ll = sapply(M, `[[`, "ll"), auc = sapply(M, `[[`, "auc")),
        file.path(MDIR, out))

# Per-pitch residuals under the old and the corrected definition, so Figure 11 can be
# rebuilt for any cue without retraining.
res_out <- D[, .(game_pk, at_bat_number, pitch_number, season, pitcher, pitch_type, grp,
                 is_whiff, path_ratio, axis_diff, as_gap, spin_sim)]
res_out[, `:=`(res_shape = D$is_whiff - M$M0$p,            # shape only, no location at all
               res_old   = detrend(D$is_whiff - M$M0$p),   # shape only + post-hoc lm (published)
               res_new   = D$is_whiff - M$M2$p,            # shape + location, NO path_ratio
               res_m3    = D$is_whiff - M$M3$p)]
# res_new deliberately excludes path_ratio: it is one of the cues under test, and a cue
# that is also a model feature is orthogonalised by construction rather than tested.
saveRDS(res_out, file.path(MDIR, sprintf("resid_old_vs_new_%s.rds", GRP)))
cat("\nsaved", out, "and resid_old_vs_new_", GRP, ".rds\n", sep = "")
