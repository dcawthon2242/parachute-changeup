#!/usr/bin/env Rscript

# DO THE DOWNGRADED CHANGEUPS EARN IT BACK ON CONTACT?
#
# The compression cuts low-separation, flat-slot changeups on whiff. If those pitches are really
# doing sinker work - inducing grounders and weak contact rather than swings and misses - then the
# whiff-based grade is measuring the wrong thing for them and the cut is unfair.
#
# Same population and same stuff-only feature set as the whiff and run value models, so the three
# are directly comparable. Four contact outcomes, all conditional on a ball in play:
#     ground ball, hard hit (95+ mph), weak (85 or less), and xwOBA on contact.
# Plus two per-pitch rates, since a pitch that never gets put in play manages contact by avoiding it.
#
# The sinker benchmark comes from the full four-pitch file, so "sinker-like" is measured against
# actual sinkers rather than asserted.

suppressPackageStartupMessages({ library(data.table); library(lightgbm); library(bit64) })
set.seed(41); options(width = 220); MDIR <- "data/statcast_model"
CACHE <- file.path(MDIR, "contact_oof.rds")

P <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(P)
P <- P[is.finite(release_speed) & is.finite(release_spin_rate) & is.finite(release_extension) &
       is.finite(ax) & is.finite(az) & is.finite(release_pos_x) & is.finite(release_pos_z) &
       is.finite(sax) & is.finite(cax) & is.finite(speed_diff) & is.finite(rv) &
       is.finite(axis_diff) & is.finite(arm_angle) & is.finite(arm_diff)]
P[, lh := p_throws == "L"]
P[, `:=`(tj_ax = fifelse(lh, -ax, ax), tj_x0 = fifelse(lh, -release_pos_x, release_pos_x),
         tj_ax_diff = fifelse(lh, -ax_diff, ax_diff), sep = -speed_diff,
         mx = fifelse(lh, -ax, ax), mz = az + 32.174)]
P[, spin_axis := (atan2(sax, cax) * 180/pi) %% 360]
P[, spin_axis := fifelse(lh, (360 - spin_axis) %% 360, spin_axis)]
wrap <- function(d) ((d + 180) %% 360) - 180
P[, dev := wrap(wrap(atan2(mx, mz)*180/pi + (spin_axis - 180)) -
                median(wrap(atan2(mx, mz)*180/pi + (spin_axis - 180)), na.rm = TRUE))]
P[, absdev := abs(dev)][, id := as.character(pitcher)]

SHAPE <- c("release_speed","release_spin_rate","release_extension","tj_ax","az","tj_x0",
           "release_pos_z","speed_diff","tj_ax_diff","az_diff","axis_diff","dev","absdev",
           "arm_angle","arm_diff")

B <- P[is_bip == 1 & !is.na(bb_type) & bb_type != ""]
B[, `:=`(gb = as.integer(bb_type == "ground_ball"),
         hard = as.integer(launch_speed >= 95),
         weak = as.integer(launch_speed <= 85),
         xw = estimated_woba_using_speedangle)]
cat(sprintf("%d changeups, %d put in play (%.1f%%). exit velo present on %.1f%% of them.\n",
            nrow(P), nrow(B), 100*nrow(B)/nrow(P), 100*mean(is.finite(B$launch_speed))))

grp_oof <- function(D, feats, lab, obj, K = 5) {
  D <- D[is.finite(get(lab))]
  arms <- unique(D$id); fa <- data.table(id = arms, fold = sample(rep(1:K, length.out = length(arms))))
  D <- merge(D, fa, by = "id", sort = FALSE); y <- D[[lab]]; p <- rep(NA_real_, nrow(D))
  for (f in 1:K) {
    tri <- which(D$fold != f); va <- sample(tri, floor(.10*length(tri))); tr <- setdiff(tri, va)
    dtr <- lgb.Dataset(as.matrix(D[tr, ..feats]), label = y[tr])
    dva <- lgb.Dataset.create.valid(dtr, as.matrix(D[va, ..feats]), label = y[va])
    m <- lgb.train(params = list(objective = obj,
                     metric = if (obj == "binary") "binary_logloss" else "l2",
                     learning_rate = .05, num_leaves = 31, min_data_in_leaf = 200,
                     feature_fraction = .8, bagging_fraction = .8, bagging_freq = 1),
                   data = dtr, nrounds = 2000, valids = list(v = dva),
                   early_stopping_rounds = 60, verbose = -1)
    p[D$fold == f] <- predict(m, as.matrix(D[D$fold == f, ..feats]))
  }
  data.table(id = D$id, season = D$season, y = y, p = p)
}

if (!file.exists(CACHE)) {
  cat("fitting contact models with pitcher-grouped folds ...\n")
  saveRDS(list(gb   = grp_oof(B, SHAPE, "gb",   "binary"),
               hard = grp_oof(B[is.finite(launch_speed)], SHAPE, "hard", "binary"),
               weak = grp_oof(B[is.finite(launch_speed)], SHAPE, "weak", "binary"),
               xw   = grp_oof(B[is.finite(xw)], SHAPE, "xw", "regression"),
               bip  = grp_oof(P[is_swing == 1], SHAPE, "is_bip", "binary")), CACHE)
}
CO <- readRDS(CACHE)

TR <- P[, .(velo = mean(release_speed), sep = mean(sep), axis = mean(axis_diff),
            seamdev = mean(absdev), signed = mean(dev), spin = mean(release_spin_rate),
            slot = mean(arm_angle), armgap = mean(arm_diff), ivb = mean(mz)/32.174*12,
            hb = mean(mx)/32.174*12, ext = mean(release_extension), relz = mean(release_pos_z),
            npit = .N), by = .(id, season)]
NM <- unique(P[, .(id, season, nm = player_name)]); TR <- merge(TR, NM, by = c("id","season"))
MV <- readRDS(file.path(MDIR, "smean_movers.rds")); setDT(MV)
TR <- merge(TR, MV[, .(id, season, mv, nsw, whiff_act = act, whiff_pred = s0, e0)],
            by = c("id","season"))

agg <- function(X) X[, .(n = .N, act = mean(y), pred = mean(p)), by = .(id, season)]
for (k in names(CO)) {
  A <- agg(CO[[k]]); setnames(A, c("n","act","pred"), paste0(k, c("_n","_act","_pred")))
  TR <- merge(TR, A, by = c("id","season"), all.x = TRUE)
}
D <- TR[nsw >= 75 & gb_n >= 40]
cat(sprintf("%d pitcher-seasons with 75+ swings and 40+ balls in play, %d arms. median %d BIP.\n\n",
            nrow(D), uniqueN(D$id), median(D$gb_n)))

VARS <- c("sep","slot","axis","seamdev","signed","velo","spin","armgap","ivb","hb","ext","relz")
NICE <- c(sep = "Separation (mph)", slot = "Arm slot (deg)", axis = "Spin-axis gap (deg)",
          seamdev = "Seam deviation (deg)", signed = "Signed seam dev (deg)", velo = "Velocity (mph)",
          spin = "Spin rate (rpm)", armgap = "Arm-angle gap (deg)", ivb = "Ind. vert. break (in)",
          hb = "Horizontal break (in)", ext = "Extension (ft)", relz = "Release height (ft)")
OUT <- c(gb = "Ground ball %", hard = "Hard hit % (95+)", weak = "Weak % (85 or less)",
         xw = "xwOBA on contact", bip = "Balls in play per swing")
SCALE <- c(gb = 100, hard = 100, weak = 100, xw = 1000, bip = 100)
GOOD <- c(gb = 1, hard = -1, weak = 1, xw = -1, bip = 0)   # sign that favours the pitcher

for (k in names(OUT)) {
  D[, (paste0(k, "_r")) := SCALE[[k]] * (get(paste0(k,"_act")) - get(paste0(k,"_pred")))]
  D[, (paste0(k, "_l")) := SCALE[[k]] * get(paste0(k,"_act"))]
}
D[, whiff_r := whiff_act - whiff_pred]      # both already in percentage points

clus <- function(dat, yv, xv, wv) {
  dat <- dat[is.finite(get(yv)) & is.finite(get(xv))]
  if (nrow(dat) < 40) return(c(NA, NA, NA))
  x <- scale(dat[[xv]])[,1]; y <- dat[[yv]]; w <- as.numeric(dat[[wv]])
  m <- lm(y ~ x, weights = w); u <- residuals(m)*sqrt(w); X <- model.matrix(m)*sqrt(w)
  nc <- uniqueN(dat$id); b <- solve(crossprod(X))
  V <- b %*% crossprod(rowsum(X*u, dat$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[2]); s <- unname(sqrt(diag(V))[2]); c(e, s, 2*pt(-abs(e/s), nc-1))
}
st <- function(e, p) if (is.na(e)) "     -   " else
  sprintf("%+7.2f%s", e, ifelse(p < .01, "**", ifelse(p < .05, "* ", "  ")))

## =============================================================================================
cat("=== 1. IS THE MODEL EVEN MISSING ON CONTACT? CALIBRATION BY PREDICTED QUINTILE ===\n\n")
for (k in names(OUT)) {
  wv <- paste0(k, "_n")
  D[, qq := cut(get(paste0(k,"_pred")), quantile(get(paste0(k,"_pred")), 0:5/5),
                include.lowest = TRUE, labels = paste0("Q",1:5))]
  Q <- D[, .(seasons = .N, pred = weighted.mean(SCALE[[k]]*get(paste0(k,"_pred")), get(wv)),
             act = weighted.mean(get(paste0(k,"_l")), get(wv))), by = qq][order(qq)]
  cat(sprintf("  %-26s %s\n", OUT[[k]],
              paste(sprintf("%s %5.1f/%5.1f", Q$qq, Q$pred, Q$act), collapse = "  ")))
}
cat("  (predicted / actual. a well-calibrated model tracks the diagonal.)\n")

## =============================================================================================
cat("\n=== 2. WHICH TRAITS PREDICT CONTACT-QUALITY RESIDUALS? ===\n\n")
cat("  standardised slope per 1 SD of trait, arm-clustered, weighted by balls in play.\n")
cat("  positive is more of the outcome, not necessarily better. ** p<.01, * p<.05\n\n")
cat(sprintf("  %-24s %9s %9s %9s %9s %9s %9s\n", "trait", "whiff", "GB%", "hard%", "weak%",
            "xwOBAcon", "BIP/swing"))
for (v in VARS) {
  cells <- c(st(clus(D, "whiff_r", v, "nsw")[1], clus(D, "whiff_r", v, "nsw")[3]),
             sapply(names(OUT), function(k) { r <- clus(D, paste0(k,"_r"), v, paste0(k,"_n")); st(r[1], r[3]) }))
  cat(sprintf("  %-24s %s\n", NICE[[v]], paste(cells, collapse = " ")))
}

## =============================================================================================
cat("\n=== 3. DO THE WHIFF DECLINERS COMPENSATE ON CONTACT? ===\n\n")
D[, grp := fifelse(mv < -10, "cut > 10", fifelse(mv > 10, "raised > 10", "little change"))]
G <- D[, .(seasons = .N, sep = mean(sep), slot = mean(slot),
           whiff = weighted.mean(whiff_act, nsw), whiff_r = weighted.mean(whiff_r, nsw),
           gb = weighted.mean(gb_l, gb_n), gb_r = weighted.mean(gb_r, gb_n),
           hard = weighted.mean(hard_l, hard_n), hard_r = weighted.mean(hard_r, hard_n),
           weak = weighted.mean(weak_l, weak_n), weak_r = weighted.mean(weak_r, weak_n),
           xw = weighted.mean(xw_l, xw_n), xw_r = weighted.mean(xw_r, xw_n),
           bip = weighted.mean(bip_l, bip_n)), by = grp]
print(G[order(-seasons), .(grp, seasons, sep = round(sep,1), slot = round(slot,0),
                           whiff = round(whiff,1), whiff_miss = round(whiff_r,2),
                           gb = round(gb,1), gb_miss = round(gb_r,2),
                           hard = round(hard,1), hard_miss = round(hard_r,2),
                           weak = round(weak,1), weak_miss = round(weak_r,2),
                           xwcon = round(xw,0), xw_miss = round(xw_r,1),
                           bip_per_sw = round(bip,1))], row.names = FALSE)

cat("\n  same split, but on separation directly (tertiles):\n\n")
D[, st3 := cut(sep, quantile(sep, 0:3/3), include.lowest = TRUE,
               labels = c("T1 least sep","T2","T3 most sep"))]
S <- D[, .(seasons = .N, sep = mean(sep), slot = mean(slot),
           whiff = weighted.mean(whiff_act, nsw),
           gb = weighted.mean(gb_l, gb_n), gb_r = weighted.mean(gb_r, gb_n),
           hard = weighted.mean(hard_l, hard_n), hard_r = weighted.mean(hard_r, hard_n),
           weak = weighted.mean(weak_l, weak_n), xw = weighted.mean(xw_l, xw_n),
           xw_r = weighted.mean(xw_r, xw_n), bip = weighted.mean(bip_l, bip_n)), by = st3][order(st3)]
print(S[, .(st3, seasons, sep = round(sep,1), slot = round(slot,0), whiff = round(whiff,1),
            gb = round(gb,1), gb_miss = round(gb_r,2), hard = round(hard,1),
            hard_miss = round(hard_r,2), weak = round(weak,1), xwcon = round(xw,0),
            xw_miss = round(xw_r,1), bip_per_sw = round(bip,1))], row.names = FALSE)

## =============================================================================================
cat("\n=== 4. THE SINKER BENCHMARK ===\n\n")
A <- readRDS(file.path(MDIR, "parachute_rv.rds")); setDT(A)
A <- A[pitch_type %in% c("SI","CH","FF") & !is.na(bb_type) & bb_type != ""]
A[, `:=`(gb = as.integer(bb_type == "ground_ball"), hard = as.integer(launch_speed >= 95),
         weak = as.integer(launch_speed <= 85))]
RA <- readRDS(file.path(MDIR, "parachute_rv.rds")); setDT(RA)
RA <- RA[pitch_type %in% c("SI","CH","FF") & is.finite(rv)]
# raw run value per pitch type is confounded by count: changeups skew to two strikes, where the run
# value swings are wider. reweight every pitch type to the same count distribution before comparing.
RA[, cnt := paste0(balls, "-", strikes)]
CW <- RA[, .(w = .N/nrow(RA)), by = cnt]
PC <- merge(RA[, .(n = .N, m = mean(rv)), by = .(pitch_type, cnt)], CW, by = "cnt")
PC <- PC[, .(pitches = sum(n), rv100 = 100*mean(rv_raw <- NA), rv100_raw = NA_real_), by = pitch_type]
PC <- merge(RA[, .(pitches = .N, rv100_raw = 100*mean(rv)), by = pitch_type],
            merge(RA[, .(m = mean(rv)), by = .(pitch_type, cnt)], CW, by = "cnt")[
              , .(rv100 = 100*sum(m*w)/sum(w)), by = pitch_type], by = "pitch_type")
BQ <- A[, .(bip = .N, gb = 100*mean(gb), hard = 100*mean(hard, na.rm = TRUE),
            weak = 100*mean(weak, na.rm = TRUE), ev = mean(launch_speed, na.rm = TRUE),
            la = mean(launch_angle, na.rm = TRUE),
            xw = 1000*mean(estimated_woba_using_speedangle, na.rm = TRUE)), by = pitch_type]
BM <- merge(PC, BQ, by = "pitch_type")
LB <- c(SI = "Sinker", CH = "Changeup (all)", FF = "Four-seam")
BM[, pitch := LB[pitch_type]]
# the changeup tertiles on the same footing
CT <- D[, .(pitch = as.character(st3), pitches = sum(npit), rv100 = NA_real_,
            bip = sum(gb_n), gb = weighted.mean(gb_l, gb_n), hard = weighted.mean(hard_l, hard_n),
            weak = weighted.mean(weak_l, weak_n), ev = NA_real_, la = NA_real_,
            xw = weighted.mean(xw_l, xw_n)), by = st3][, st3 := NULL]
CT[, pitch := paste("  changeup,", pitch)]
print(rbind(BM[order(pitch_type), .(pitch, pitches, bip, gb = round(gb,1), hard = round(hard,1),
                   weak = round(weak,1), ev = round(ev,1), la = round(la,1), xwcon = round(xw,0),
                   rv100_raw = round(rv100_raw,2), rv100_evencount = round(rv100,2))],
            CT[, .(pitch, pitches, bip, gb = round(gb,1), hard = round(hard,1), weak = round(weak,1),
                   ev = round(ev,1), la = round(la,1), xwcon = round(xw,0),
                   rv100_raw = NA_real_, rv100_evencount = NA_real_)], fill = TRUE), row.names = FALSE)

cat("\n  same-pitcher check: arms in the study who also threw 200+ sinkers, their own changeup\n")
cat("  against their own sinker on balls in play.\n\n")
SI <- A[pitch_type == "SI"][, id := as.character(pitcher)]
SIS <- SI[, .(si_bip = .N, si_gb = 100*mean(gb), si_hard = 100*mean(hard, na.rm = TRUE),
              si_xw = 1000*mean(estimated_woba_using_speedangle, na.rm = TRUE)), by = .(id, season)]
SP <- merge(D, SIS[si_bip >= 50], by = c("id","season"))
SP[, s3 := cut(sep, quantile(sep, 0:3/3), include.lowest = TRUE,
               labels = c("T1 least sep","T2","T3 most sep"))]
PR <- SP[, .(seasons = .N, ch_gb = weighted.mean(gb_l, gb_n), si_gb = weighted.mean(si_gb, si_bip),
             ch_hard = weighted.mean(hard_l, hard_n), si_hard = weighted.mean(si_hard, si_bip),
             ch_xw = weighted.mean(xw_l, xw_n), si_xw = weighted.mean(si_xw, si_bip)), by = s3][order(s3)]
print(PR[, .(s3, seasons, ch_gb = round(ch_gb,1), si_gb = round(si_gb,1),
             ch_hard = round(ch_hard,1), si_hard = round(si_hard,1),
             ch_xwcon = round(ch_xw,0), si_xwcon = round(si_xw,0))], row.names = FALSE)
cat(sprintf("\n  %d pitcher-seasons have both. paired difference in xwOBAcon (changeup minus sinker): %+.0f, p %.4f\n",
            nrow(SP), SP[, mean(xw_l - si_xw)], t.test(SP$xw_l, SP$si_xw, paired = TRUE)$p.value))

## =============================================================================================
cat("\n=== 5. WHERE THE RUN VALUE ACTUALLY COMES FROM ===\n\n")
# split each pitcher-season's run value into the part from balls in play and the part from everything
# else, then see which piece separation is buying.
RV <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(RV)
RV[, id := as.character(pitcher)]
RS <- RV[is.finite(rv), .(np = .N, rv100 = 100*mean(rv),
                          rv_bip = 100*sum(rv*(is_bip == 1), na.rm = TRUE)/.N,
                          rv_oth = 100*sum(rv*(is_bip != 1 | is.na(is_bip)), na.rm = TRUE)/.N,
                          bip_share = mean(is_bip == 1, na.rm = TRUE)), by = .(id, season)]
D2 <- merge(D, RS, by = c("id","season"))
T2 <- D2[, .(seasons = .N, rv100 = weighted.mean(rv100, np),
             from_bip = weighted.mean(rv_bip, np), from_rest = weighted.mean(rv_oth, np),
             bip_rate = 100*weighted.mean(bip_share, np)), by = st3][order(st3)]
print(T2[, .(st3, seasons, rv_per100 = round(rv100,2), from_balls_in_play = round(from_bip,2),
             from_everything_else = round(from_rest,2), bip_rate = round(bip_rate,1))],
      row.names = FALSE)
cat("\n  run value is pitcher-positive. slopes per 1 SD of separation, arm-clustered:\n")
for (yv in c("rv100","rv_bip","rv_oth")) {
  r <- clus(D2, yv, "sep", "np")
  cat(sprintf("    %-24s %s\n", yv, st(r[1], r[3])))
}
## =============================================================================================
cat("\n=== 6. IS ANY CONTACT OVERPERFORMANCE A REPEATABLE SKILL? ===\n\n")
# a miss the model could learn has to persist. whiff overperformance does; test whether the contact
# residuals do, using the same year-over-year design.
cat(sprintf("  %-24s %8s %8s %10s %10s\n", "residual", "seasons", "sd", "yr to yr r", "p"))
for (k in c("whiff", names(OUT))) {
  yv <- paste0(k, "_r"); wv <- if (k == "whiff") "nsw" else paste0(k, "_n")
  X <- D[is.finite(get(yv)), .(id, season, r = get(yv), w = get(wv))]
  Y <- copy(X)[, season := season - 1][, .(id, season, r2 = r, w2 = w)]
  J <- merge(X, Y, by = c("id","season"))
  if (nrow(J) < 30) next
  ct <- cor.test(J$r, J$r2)
  cat(sprintf("  %-24s %8d %8.2f %10s %10.4f\n", if (k == "whiff") "Whiff" else OUT[[k]],
              nrow(J), sd(X$r), sprintf("%+.3f", ct$estimate), ct$p.value))
}

cat("\n  are the group contact residuals different from zero at all?\n")
cat("  (arm-clustered t-test of the residual against zero, weighted by balls in play)\n\n")
zt <- function(dat, yv, wv) {
  dat <- dat[is.finite(get(yv))]
  w <- as.numeric(dat[[wv]]); y <- dat[[yv]]
  m <- lm(y ~ 1, weights = w); u <- residuals(m)*sqrt(w); X <- model.matrix(m)*sqrt(w)
  nc <- uniqueN(dat$id); b <- solve(crossprod(X))
  V <- b %*% crossprod(rowsum(X*u, dat$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[1]); s <- unname(sqrt(diag(V))[1]); c(e, 2*pt(-abs(e/s), nc-1))
}
cat(sprintf("  %-16s %s\n", "group",
            paste(sprintf("%16s", c("whiff", OUT[names(OUT)])), collapse = "")))
for (g in c("cut > 10","little change","raised > 10")) {
  X <- D[grp == g]
  cells <- sapply(c("whiff", names(OUT)), function(k) {
    r <- zt(X, paste0(k, "_r"), if (k == "whiff") "nsw" else paste0(k, "_n"))
    sprintf("%12.2f%s", r[1], ifelse(r[2] < .01, "**", ifelse(r[2] < .05, "* ", "  ")))
  })
  cat(sprintf("  %-16s %s\n", g, paste(cells, collapse = "")))
}
saveRDS(D, file.path(MDIR, "contact_season.rds"))
