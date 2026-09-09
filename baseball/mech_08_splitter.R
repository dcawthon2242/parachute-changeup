#!/usr/bin/env Rscript

# SPLITTER REPLICATION OF THE LOCKED CHANGEUP CELL.
#
# Same gates, same four-seam-primary filter, disjoint pitches. If the combined bin is under 10
# pitcher-seasons the script stops after the count: that is underpowered, not a null.
#
# MLB splitters live in miss_grade / oof_whiff_resid (2023-2026 only). parachute_rv has no FS, so
# the replication cannot reach 2020-2022. Arm angle is not on those files for splitters; it is
# reconstructed from release point using the same regression the D1 arm file was built on, fit
# where measured arm exists (changeup seasons in parachute_ff.rds).
# D1 splitters are scored with the inferred-efficiency model from ncaa_03_build.R.

suppressPackageStartupMessages({ library(data.table); library(bit64); library(lightgbm) })
options(width = 200); set.seed(31)
MDIR <- "data/statcast_model"
LOCK <- readRDS(file.path(MDIR, "locked_spec.rds"))
SPEC <- LOCK$spec
circ <- function(a, b) pmin(abs(a - b), 360 - abs(a - b))

## ---- usage (same sources as spec_lock; missing rows dropped) -----------------------------------
AR <- rbindlist(lapply(2020:2026, function(y) {
  x <- fread(sprintf("data/savant/arsenal_%d.csv", y), showProgress = FALSE)
  setnames(x, names(x)[2], "pitcher"); x[, season := y]
  x[, .(id = as.character(pitcher), season,
        ff_use = suppressWarnings(as.numeric(n_ff)),
        si_use = suppressWarnings(as.numeric(n_si)))] }), use.names = TRUE)
AR[is.na(ff_use), ff_use := 0][is.na(si_use), si_use := 0]
UN <- readRDS(file.path(MDIR, "ncaa_usage.rds")); setDT(UN)
UN[, id := as.character(PitcherId)]
U <- rbind(AR[, .(league = "MLB", id, season, ff_use, si_use)],
           UN[, .(league = "D1", id, season, ff_use, si_use)])

## ---- MLB: season table from existing FS residuals ----------------------------------------------
OOF <- readRDS(file.path(MDIR, "oof_whiff_resid.rds")); setDT(OOF)
MG  <- readRDS(file.path(MDIR, "miss_grade_data.rds")); setDT(MG)
FS  <- merge(OOF[pitch_type == "FS", .(game_pk, at_bat_number, pitch_number, pitcher, season,
                                       is_whiff, wres)],
             MG[pitch_type == "FS", .(game_pk, at_bat_number, pitch_number, axis_diff,
                                      speed_diff, release_pos_x, release_pos_z,
                                      release_extension)],
             by = c("game_pk","at_bat_number","pitch_number"))
AS <- readRDS(file.path(MDIR, "active_spin_long.rds")); setDT(AS)
S <- FS[, .(y_raw = 100*mean(wres), nsw = .N, axis = mean(axis_diff, na.rm = TRUE),
            vs = -mean(speed_diff, na.rm = TRUE),
            rx = mean(abs(release_pos_x), na.rm = TRUE),
            rz = mean(release_pos_z, na.rm = TRUE),
            ext = mean(release_extension, na.rm = TRUE)),
        by = .(id = as.character(pitcher), season)][nsw >= SPEC$swing_gate]
S <- merge(S, AS[pitch_type == "FS", .(id = as.character(pitcher), season, ec = active_spin)],
           by = c("id","season"))
S <- merge(S, AS[pitch_type == "FF", .(id = as.character(pitcher), season, ef = active_spin)],
           by = c("id","season"))
# Release-point arm proxy, fit where measured arm exists so the MLB FS slot is on the same
# scale the D1 file already uses.
CH <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(CH)
MS <- CH[is.finite(arm_angle), .(arm = mean(arm_angle), rx = mean(abs(release_pos_x)),
        rz = mean(release_pos_z), ext = mean(release_extension, na.rm = TRUE)),
        by = .(pitcher, season)]
MS[, proxy := atan2(rz - 4.7, rx) * 180/pi]
fit_arm <- lm(arm ~ proxy + rz + rx + ext, MS)
S[, proxy := atan2(rz - 4.7, rx) * 180/pi]
S[, arm := predict(fit_arm, S)]
S[, league := "MLB"]
cat(sprintf("MLB splitter seasons with 40+ swings, both active spins, arm proxy: %d\n",
            S[is.finite(ec) & is.finite(ef) & is.finite(arm), .N]))

## ---- D1: infer efficiency, fit a splitter residual, pair to four-seam --------------------------
D1C <- file.path(MDIR, "ncaa_fs_whiff_resid.rds")
D <- readRDS(file.path(MDIR, "ncaa_d1_seq_pitches.rds")); setDT(D)
D[, PitcherId := as.character(PitcherId)]

# Efficiency model: same specification as ncaa_03_build.R, trained on MLB measured active spin.
M <- readRDS(file.path(MDIR, "parachute_rv.rds"))[pitch_type %in% c("CH","FF")]
M[, magnus := sqrt(ax^2 + (az + 32.174)^2)]
M[, ratio := magnus / (release_spin_rate * release_speed / 1000)]
PM <- M[is.finite(ratio) & release_spin_rate > 500,
        .(n = .N, ratio = mean(ratio), magnus = mean(magnus), spin = mean(release_spin_rate),
          velo = mean(release_speed), ext = mean(release_extension, na.rm = TRUE)),
        by = .(pitcher, season, pitch_type)]
PM <- merge(PM, AS[, .(pitcher, season, pitch_type, measured = active_spin)],
            by = c("pitcher","season","pitch_type"))[n >= 50 & is.finite(measured) & is.finite(ext)]
FE <- c("ratio","magnus","spin","velo","ext")
mod <- lgb.train(params = list(objective = "regression", learning_rate = .05, num_leaves = 15,
                 min_data_in_leaf = 40, feature_fraction = .9), verbose = -1, nrounds = 400,
                 data = lgb.Dataset(as.matrix(PM[, ..FE]), label = PM$measured))
SH <- D[pt %in% c("Splitter","Four-Seam") & is.finite(ratio) & SpinRate > 500 & is.finite(Extension),
        .(n = .N, ratio = mean(ratio), magnus = mean(magnus), spin = mean(SpinRate),
          velo = mean(RelSpeed), ext = mean(Extension), axis = mean(SpinAxis, na.rm = TRUE)),
        by = .(PitcherId, season, pt)][n >= 30]
SH[, eff := pmin(pmax(predict(mod, as.matrix(SH[, ..FE])), 0), 1)]
PR <- merge(SH[pt == "Splitter", .(PitcherId, season, n_fs = n, ec = eff, axis_fs = axis, velo_fs = velo)],
            SH[pt == "Four-Seam", .(PitcherId, season, ef = eff, axis_ff = axis, velo_ff = velo)],
            by = c("PitcherId","season"))
PR[, `:=`(axis = circ(axis_fs, axis_ff), vs = velo_ff - velo_fs)]

ARM <- as.data.table(readRDS(file.path(MDIR, "ncaa_armangle.rds")))
ARM[, PitcherId := as.character(PitcherId)]
PR <- merge(PR, ARM[, .(PitcherId, season, arm = arm_hat)], by = c("PitcherId","season"))

if (!file.exists(D1C)) {
  D[, L := PitcherThrows == "Left"]
  FF <- D[pt == "Four-Seam" & is.finite(RelSpeed),
          .(nff = .N, ff_speed = mean(RelSpeed), ff_ax = mean(ax, na.rm = TRUE),
            ff_az = mean(az, na.rm = TRUE), ff_axis = mean(SpinAxis, na.rm = TRUE)),
          by = .(PitcherId, season)][nff >= 30]
  C <- merge(D[pt == "Splitter"], FF, by = c("PitcherId","season"))
  C <- C[is.finite(RelSpeed) & is.finite(ax) & is.finite(az) & is.finite(SpinAxis) &
         is.finite(px) & is.finite(pz) & is.finite(SpinRate) & is.finite(Extension) &
         is.finite(VertApprAngle) & is.finite(HorzApprAngle) & swing == TRUE]
  C[, `:=`(speed_diff = RelSpeed - ff_speed, ax_diff = ax - ff_ax, az_diff = az - ff_az)]
  C[, `:=`(tj_x0 = fifelse(L, -RelSide, RelSide), tj_ax = fifelse(L, -ax, ax),
           tj_ax_diff = fifelse(L, -ax_diff, ax_diff), tj_px = fifelse(L, -px, px),
           tj_haa = fifelse(L, -HorzApprAngle, HorzApprAngle),
           tj_axis = fifelse(L, (360 - SpinAxis) %% 360, SpinAxis),
           same_hand = as.integer((PitcherThrows == "Left") == (BatterSide == "Left")))]
  FEAT <- c("RelSpeed","SpinRate","Extension","tj_ax","az","tj_x0","RelHeight","tj_axis",
            "speed_diff","tj_ax_diff","az_diff",
            "tj_px","pz","VertApprAngle","tj_haa","same_hand","Balls","Strikes")
  cat(sprintf("D1 splitter swings for the residual: %s  whiff rate %.3f\n",
              format(nrow(C), big.mark = ","), mean(C$whiff)))
  K <- 4; fold <- sample(rep(1:K, length.out = nrow(C))); p <- rep(NA_real_, nrow(C))
  for (f in 1:K) {
    tr <- C[fold != f]; n <- nrow(tr); vi <- sample(n, floor(.12*n))
    dtr <- lgb.Dataset(as.matrix(tr[-vi, ..FEAT]), label = tr$whiff[-vi])
    dva <- lgb.Dataset.create.valid(dtr, as.matrix(tr[vi, ..FEAT]), label = tr$whiff[vi])
    m <- lgb.train(params = list(objective = "binary", metric = "binary_logloss",
                   learning_rate = .06, num_leaves = 31, min_data_in_leaf = 80,
                   feature_fraction = .8, bagging_fraction = .8, bagging_freq = 1),
                   data = dtr, nrounds = 800, valids = list(v = dva),
                   early_stopping_rounds = 40, verbose = -1)
    p[fold == f] <- predict(m, as.matrix(C[fold == f, ..FEAT]))
  }
  C[, r_all := whiff - p]
  saveRDS(C[, .(PitcherId, season, r_all, speed_diff)], D1C)
} else cat("(using cached D1 splitter residuals)\n")
DR <- readRDS(D1C); setDT(DR); DR[, PitcherId := as.character(PitcherId)]
DS <- DR[, .(y_raw = 100*mean(r_all), nsw = .N), by = .(PitcherId, season)][nsw >= SPEC$swing_gate]
N <- merge(PR, DS, by = c("PitcherId","season"))
N[, `:=`(league = "D1", id = PitcherId)]
cat(sprintf("D1 splitter seasons with 40+ swings, efficiency, arm: %d\n", nrow(N)))

## ---- assemble, same population filter, same gates ----------------------------------------------
COLS <- c("league","id","season","y_raw","nsw","vs","axis","arm","ec","ef")
P <- rbind(S[, ..COLS], N[, ..COLS])
P <- merge(P, U, by = c("league","id","season"), all.x = TRUE)
n_miss <- P[is.na(ff_use) | is.na(si_use), .N]
if (n_miss > 0) cat(sprintf("dropping %d splitter seasons with missing usage\n", n_miss))
P <- P[is.finite(ff_use) & is.finite(si_use) & is.finite(arm) & is.finite(ec) & is.finite(ef)]
if (isTRUE(SPEC$require_ff_primary)) P <- P[ff_use >= si_use]
P[, y := resid(lm(y_raw ~ vs)), by = league]
P[, arm_thr := as.numeric(quantile(arm, SPEC$arm_pctile)), by = league]
P[, bin := ec >= SPEC$eff_min & ef >= SPEC$eff_min & axis <= SPEC$axis_max & arm >= arm_thr]

cat("\n=== splitter analog of the locked cell ===\n")
print(P[, .(pool = .N, arms = uniqueN(id), arm_threshold = round(arm_thr[1],1),
            bin_seasons = sum(bin), bin_arms = uniqueN(id[bin])), by = league], row.names = FALSE)
cat(sprintf("combined bin: %d seasons from %d arms\n", sum(P$bin), uniqueN(P[bin == TRUE]$id)))

if (sum(P$bin) < 10) {
  cat("\nGATE: combined bin is under 10 pitcher-seasons. Underpowered, not a null. Stopping.\n")
  if (sum(P$bin) > 0) print(P[bin == TRUE, .(league, id, season, axis = round(axis,1),
      arm = round(arm,1), ec, ef, nsw, y = round(y,2))], row.names = FALSE)
  saveRDS(list(spec = SPEC, data = P, gated = TRUE), file.path(MDIR, "splitter_spec.rds"))
  quit(save = "no", status = 0)
}

est <- P[, { t <- t.test(y[bin], y[!bin]); d <- as.numeric(diff(rev(t$estimate)))
             .(n = sum(bin), arms = uniqueN(id[bin]), diff = d, se = as.numeric(d/t$statistic),
               p = t$p.value) }, by = league]
cat("\n=== per league ===\n")
print(est[, .(league, n, arms, diff = round(diff,2), se = round(se,2), p = round(p,4))],
      row.names = FALSE)
ok <- est[is.finite(se) & se > 0]
if (nrow(ok) >= 1) {
  w <- 1/ok$se^2; m <- sum(ok$diff*w)/sum(w); se <- sqrt(1/sum(w))
  cat(sprintf("\npooled: %+.2f +/- %.2f | z = %.2f | p = %.4f\n",
              m, se, m/se, 2*pnorm(-abs(m/se))))
  if (nrow(ok) == 2) {
    perm <- replicate(2000, {
      Pp <- copy(P)[, y := sample(y), by = league]
      e <- Pp[, .(d = mean(y[bin]) - mean(y[!bin])), by = league]
      sum(e$d * w)/sum(w) })
    cat(sprintf("permutation p (pooled >= observed): %.4f\n", mean(perm >= m)))
  }
}
saveRDS(list(spec = SPEC, data = P, gated = FALSE), file.path(MDIR, "splitter_spec.rds"))
cat("wrote splitter_spec.rds\n")
