#!/usr/bin/env Rscript

# THE REMAINING GAP, AND WHETHER FEATURE ENGINEERING CLOSES IT.
#
# Adding shape features to the stuff model absorbed roughly a third of the matched-axis miss and
# almost none of the extreme seam-shift miss. That pattern is what you get when the signal lives in
# a narrow interaction: boosted trees find main effects easily and conjunctions of three continuous
# variables only when the cell is large enough to survive the leaf minimum. The cell in question is
# ten pitcher-seasons.
#
# So the question becomes whether the model needs the conjunction handed to it. Four specifications,
# each a strict superset of the last:
#
#   SHAPE    the marginal shape features from the previous run
#   INTER    + explicit products, so the tree can split on the conjunction in one cut
#   FLAG     + the archetype indicator itself, thresholds FROZEN on 2020-2023
#   ORACLE   + the indicator with thresholds fit on all years - not honest, included as a ceiling
#
# FLAG is the fair test and ORACLE bounds it. If FLAG lands near ORACLE the conjunction is genuinely
# learnable from frozen thresholds; if it lands near SHAPE then the search that found those
# thresholds was fitting noise and the whole line of work is weaker than it looks.
#
# Everything is evaluated on 2024-2026 after training on 2020-2023, plus pitcher-grouped folds for
# the in-sample view. Calibration is reported by decile as well as in aggregate, because a feature
# that fixes one tail while bending the middle is not an improvement.

suppressPackageStartupMessages({ library(data.table); library(lightgbm); library(bit64) })
set.seed(707); options(width = 215); MDIR <- "data/statcast_model"

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
P[, dev_raw := wrap(atan2(mx, mz) * 180/pi + (spin_axis - 180))]
P[, dev := wrap(dev_raw - median(dev_raw, na.rm = TRUE))][, absdev := abs(dev)]
P[, id := as.character(pitcher)]

# season-mean versions, because the archetypes are season-level properties and a single pitch's
# axis reading is far noisier than the season average the definitions were built on
SM <- P[, .(m_axis = mean(axis_diff), m_dev = mean(dev), m_abs = mean(absdev),
            m_sep = mean(sep), m_arm = mean(arm_angle), m_ad = mean(arm_diff),
            nsw = sum(is_swing)), by = .(id, season)]
P <- merge(P, SM, by = c("id","season"), sort = FALSE)

# explicit conjunctions. centred so the products are interpretable and not collinear with mains.
ctr <- function(x) x - mean(x)
P[, `:=`(x_axis_dev = ctr(m_axis) * ctr(m_abs),
         x_axis_sep = ctr(m_axis) * ctr(m_sep),
         x_dev_sep  = ctr(m_abs)  * ctr(m_sep),
         x_axis_arm = ctr(m_axis) * ctr(m_arm),
         x_triple   = ctr(m_axis) * ctr(m_abs) * ctr(m_sep))]

## thresholds frozen on the training era only ----------------------------------------------------
TRS <- SM[nsw >= 75 & season <= 2023]
q <- list(ax = quantile(TRS$m_axis, .25), vs = quantile(TRS$m_sep, .60),
          d25 = quantile(TRS$m_abs, .25), d75 = quantile(TRS$m_abs, .75),
          s75 = quantile(TRS$m_dev, .75))
ALLS <- SM[nsw >= 75]
qa <- list(ax = quantile(ALLS$m_axis, .25), vs = quantile(ALLS$m_sep, .60),
           d25 = quantile(ALLS$m_abs, .25), d75 = quantile(ALLS$m_abs, .75),
           s75 = quantile(ALLS$m_dev, .75))
cat("thresholds  frozen on 2020-23:  axis<=%.1f sep>=%.1f dev>=%.1f dev75>=%.1f sgn75>=%.1f\n")
cat(sprintf("            frozen 2020-23 : axis<=%5.1f  sep>=%4.1f  dev>=%4.1f  dev75>=%5.1f  sgn75>=%5.1f\n",
            q$ax, q$vs, q$d25, q$d75, q$s75))
cat(sprintf("            all years      : axis<=%5.1f  sep>=%4.1f  dev>=%4.1f  dev75>=%5.1f  sgn75>=%5.1f\n\n",
            qa$ax, qa$vs, qa$d25, qa$d75, qa$s75))
P[, `:=`(f_seam = as.integer(m_axis <= q$ax  & m_sep >= q$vs & m_abs >= q$d25),
         f_extr = as.integer(m_axis <= q$ax  & m_abs >= q$d75 & m_dev >= q$s75),
         o_seam = as.integer(m_axis <= qa$ax & m_sep >= qa$vs & m_abs >= qa$d25),
         o_extr = as.integer(m_axis <= qa$ax & m_abs >= qa$d75 & m_dev >= qa$s75))]

SHAPE <- c("release_speed","release_spin_rate","release_extension","tj_ax","az","tj_x0",
           "release_pos_z","speed_diff","tj_ax_diff","az_diff","axis_diff","dev","absdev",
           "arm_angle","arm_diff")
INTER <- c(SHAPE, "x_axis_dev","x_axis_sep","x_dev_sep","x_axis_arm","x_triple")
FLAG  <- c(INTER, "f_seam","f_extr")
ORACLE<- c(INTER, "o_seam","o_extr")
SETS <- list(SHAPE = SHAPE, INTER = INTER, FLAG = FLAG, ORACLE = ORACLE)

M <- readRDS(file.path(MDIR, "archetype_roster.rds")); setDT(M)
TAGS <- c("A2_wide","A4_broad","A5_seam","A6_extreme","A7_trad","A8_mismatch")
LAB <- c(A2_wide = "matched axis (optimized)", A4_broad = "axis + separation",
         A5_seam = "seam-shifted", A6_extreme = "extreme seam shift",
         A7_trad = "traditional (low seam dev)", A8_mismatch = "mismatch")
KEY <- M[, c("id","season","nm","nsw", TAGS), with = FALSE]

W <- P[is_swing == 1 & is.finite(whiff)]
tr <- W[season <= 2023]; te <- W[season >= 2024]
cat(sprintf("train %s swings (2020-23), test %s swings (2024-26), %d test arms\n\n",
            format(nrow(tr), big.mark=","), format(nrow(te), big.mark=","), uniqueN(te$id)))

fit1 <- function(feats, dtrain, dtest) {
  vi <- sample(nrow(dtrain), floor(.1*nrow(dtrain)))
  d1 <- lgb.Dataset(as.matrix(dtrain[-vi, ..feats]), label = dtrain$whiff[-vi])
  d2 <- lgb.Dataset.create.valid(d1, as.matrix(dtrain[vi, ..feats]), label = dtrain$whiff[vi])
  m <- lgb.train(params = list(objective = "binary", metric = "binary_logloss",
                   learning_rate = .06, num_leaves = 31, min_data_in_leaf = 300,
                   feature_fraction = .8, bagging_fraction = .8, bagging_freq = 1),
                 data = d1, nrounds = 2000, valids = list(v = d2),
                 early_stopping_rounds = 60, verbose = -1)
  list(m = m, p = predict(m, as.matrix(dtest[, ..feats])))
}
ll <- function(y,p){p<-pmin(pmax(p,1e-9),1-1e-9); -mean(y*log(p)+(1-y)*log(1-p))}
crob <- function(D, f, k) {
  environment(f) <- environment()
  D <- D[complete.cases(D[, c(all.vars(f), "id", "n"), with = FALSE])]
  wt <- as.numeric(D$n); m <- lm(f, D, weights = wt)
  u <- residuals(m)*sqrt(wt); X <- model.matrix(m)*sqrt(wt)
  nc <- uniqueN(D$id); b <- solve(crossprod(X))
  V <- b %*% crossprod(rowsum(X*u, D$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[k]); se <- unname(sqrt(diag(V))[k]); c(e, se, 2*pt(-abs(e/se), nc-1))
}

FITS <- lapply(names(SETS), function(s) { cat(sprintf("  fitting %s ...\n", s)); fit1(SETS[[s]], tr, te) })
names(FITS) <- names(SETS)

## =============================================================================================
cat("\n=== 1. FORWARD TEST: 2024-2026, TRAINED ON 2020-2023 ===\n\n")
base <- readRDS(file.path(MDIR, "absorb_temporal.rds"))
agg <- function(X) X[, .(n = .N, act = mean(y), pred = mean(p)), by = .(id, season)]
BASE <- merge(agg(base$STUFF), KEY, by = c("id","season"))[n >= 40]
cat(sprintf("  %-8s %9s %8s %12s %11s\n", "features", "logloss", "vs stuff", "season RMSE", "vs stuff"))
b_ll <- ll(base$STUFF$y, base$STUFF$p); b_rm <- sqrt(mean((100*(BASE$act-BASE$pred))^2))
cat(sprintf("  %-8s %9.5f %8s %12.3f %11s\n", "STUFF", b_ll, "-", b_rm, "-"))
SEAS <- list()
for (s in names(SETS)) {
  T2 <- data.table(id = te$id, season = te$season, y = te$whiff, p = FITS[[s]]$p)
  A <- merge(agg(T2), KEY, by = c("id","season"))[n >= 40]; SEAS[[s]] <- A
  cat(sprintf("  %-8s %9.5f %+8.2f%% %12.3f %+10.1f%%\n", s, ll(T2$y, T2$p),
              100*(ll(T2$y,T2$p)-b_ll)/b_ll, sqrt(mean((100*(A$act-A$pred))^2)),
              100*(sqrt(mean((100*(A$act-A$pred))^2))-b_rm)/b_rm))
}

## =============================================================================================
cat("\n=== 2. DOES THE ARCHETYPE GAP CLOSE ON HELD-OUT SEASONS? ===\n\n")
GG <- rbindlist(lapply(TAGS, function(t) {
  o <- list(archetype = LAB[[t]])
  B <- copy(BASE)[, `:=`(r = 100*(act-pred), gg = get(t))]
  if (B[gg == TRUE, .N] < 4) return(NULL)
  o$n <- B[gg == TRUE, .N]; e <- crob(B, r ~ gg, "ggTRUE"); o$STUFF <- e[1]; o$p_stuff <- e[3]
  for (s in names(SETS)) {
    A <- copy(SEAS[[s]])[, `:=`(r = 100*(act-pred), gg = get(t))]
    e <- crob(A, r ~ gg, "ggTRUE"); o[[s]] <- e[1]; o[[paste0("p_",s)]] <- e[3]
  }
  as.data.table(o)
}))
print(GG[, .(archetype, seasons = n, stuff = round(STUFF,2), shape = round(SHAPE,2),
             inter = round(INTER,2), flag = round(FLAG,2), oracle = round(ORACLE,2),
             p_flag = round(p_FLAG,4),
             `closed%` = round(100*(1 - FLAG/STUFF)))], row.names = FALSE)
cat("\n  Positive means the group still beats the model after that feature set. 'closed%' is how\n")
cat("  much of the original stuff-model gap the FLAG specification removes on held-out seasons.\n")

## =============================================================================================
cat("\n=== 3. WHERE DOES THE GAIN COME FROM, AND DOES ANYTHING GET WORSE? ===\n\n")
A0 <- copy(BASE); A1 <- copy(SEAS$FLAG)
CMP <- merge(A0[, .(id, season, n, act, p0 = pred)], A1[, .(id, season, p1 = pred)], by = c("id","season"))
CMP[, `:=`(e0 = 100*(act-p0), e1 = 100*(act-p1), dec = cut(100*p0, quantile(100*p0, 0:5/5),
                                                            include.lowest = TRUE, dig.lab = 3))]
D <- CMP[, .(seasons = .N, pred = round(100*mean(p0),1), actual = round(100*mean(act),1),
             bias_stuff = round(mean(e0),2), bias_flag = round(mean(e1),2),
             rmse_stuff = round(sqrt(mean(e0^2)),2), rmse_flag = round(sqrt(mean(e1^2)),2)),
         by = dec][order(dec)]
D[, better := fifelse(rmse_flag < rmse_stuff, "yes", "no")]
cat("  by quintile of the stuff model's predicted whiff rate:\n")
print(D, row.names = FALSE)
cat(sprintf("\n  seasons where the grade moved by more than 2 whiff points: %d of %d (%.0f%%)\n",
            CMP[abs(100*(p1-p0)) > 2, .N], nrow(CMP), 100*CMP[abs(100*(p1-p0)) > 2, .N]/nrow(CMP)))
cat(sprintf("  largest upgrades and downgrades the FLAG model applies:\n"))
CMP2 <- merge(CMP, unique(M[, .(id, season, nm)]), by = c("id","season"))
CMP2[, mv := 100*(p1-p0)]
print(rbind(CMP2[order(-mv)][1:6], CMP2[order(mv)][1:6])[
  , .(pitcher = nm, season, swings = n, stuff = round(100*p0,1), flag = round(100*p1,1),
      moved = round(mv,2), actual = round(100*act,1))], row.names = FALSE)

## =============================================================================================
cat("\n=== 4. WHICH FEATURES THE MODEL ACTUALLY USES ===\n\n")
imp <- lgb.importance(FITS$FLAG$m)
setDT(imp); imp[, share := 100*Gain/sum(Gain)]
print(imp[, .(feature = Feature, `gain%` = round(share,1), splits = Frequency)][1:22], row.names = FALSE)
