#!/usr/bin/env Rscript

# TWO LOOSE ENDS.
#
# A. The interaction features that helped were built from SEASON MEANS, not per-pitch readings. That
#    conflates two different improvements: averaging a noisy per-pitch measurement over a season, and
#    the conjunction itself. A specification with the season means but no products separates them. If
#    SMEAN captures most of the gain then the honest story is "give the model a stable arsenal
#    summary", which is a much cheaper recommendation than "engineer these five products".
#
# B. What the underperforming group physically is. It has been named by its rule, not described, and
#    the rule does not say whether these are hard changeups, seam-shifted ones, or both.

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
SM <- P[, .(m_axis = mean(axis_diff), m_dev = mean(dev), m_abs = mean(absdev), m_sep = mean(sep),
            m_arm = mean(arm_angle), m_ad = mean(arm_diff), m_velo = mean(release_speed),
            m_spin = mean(release_spin_rate), nsw = sum(is_swing)), by = .(id, season)]
P <- merge(P, SM, by = c("id","season"), sort = FALSE)
ctr <- function(x) x - mean(x)
P[, `:=`(x_axis_dev = ctr(m_axis)*ctr(m_abs), x_axis_sep = ctr(m_axis)*ctr(m_sep),
         x_dev_sep = ctr(m_abs)*ctr(m_sep), x_axis_arm = ctr(m_axis)*ctr(m_arm),
         x_triple = ctr(m_axis)*ctr(m_abs)*ctr(m_sep))]

SHAPE <- c("release_speed","release_spin_rate","release_extension","tj_ax","az","tj_x0",
           "release_pos_z","speed_diff","tj_ax_diff","az_diff","axis_diff","dev","absdev",
           "arm_angle","arm_diff")
SMEAN <- c(SHAPE, "m_axis","m_dev","m_abs","m_sep","m_arm","m_ad")
INTER <- c(SMEAN, "x_axis_dev","x_axis_sep","x_dev_sep","x_axis_arm","x_triple")
PRODONLY <- c(SHAPE, "x_axis_dev","x_axis_sep","x_dev_sep","x_axis_arm","x_triple")
SETS <- list(SHAPE = SHAPE, SMEAN = SMEAN, PRODONLY = PRODONLY, INTER = INTER)

M <- readRDS(file.path(MDIR, "archetype_roster.rds")); setDT(M)
TAGS <- c("A2_wide","A4_broad","A5_seam","A6_extreme","A7_trad","A8_mismatch")
LAB <- c(A2_wide = "matched axis (optimized)", A4_broad = "axis + separation",
         A5_seam = "seam-shifted", A6_extreme = "extreme seam shift",
         A7_trad = "traditional (low seam dev)", A8_mismatch = "mismatch")
KEY <- M[, c("id","season","nm","nsw", TAGS), with = FALSE]

W <- P[is_swing == 1 & is.finite(whiff)]
tr <- W[season <= 2023]; te <- W[season >= 2024]
fit1 <- function(feats) {
  vi <- sample(nrow(tr), floor(.1*nrow(tr)))
  d1 <- lgb.Dataset(as.matrix(tr[-vi, ..feats]), label = tr$whiff[-vi])
  d2 <- lgb.Dataset.create.valid(d1, as.matrix(tr[vi, ..feats]), label = tr$whiff[vi])
  m <- lgb.train(params = list(objective = "binary", metric = "binary_logloss",
                   learning_rate = .06, num_leaves = 31, min_data_in_leaf = 300,
                   feature_fraction = .8, bagging_fraction = .8, bagging_freq = 1),
                 data = d1, nrounds = 2000, valids = list(v = d2),
                 early_stopping_rounds = 60, verbose = -1)
  predict(m, as.matrix(te[, ..feats]))
}
ll <- function(y,p){p<-pmin(pmax(p,1e-9),1-1e-9); -mean(y*log(p)+(1-y)*log(1-p))}
agg <- function(X) X[, .(n = .N, act = mean(y), pred = mean(p)), by = .(id, season)]
crob <- function(D, f, k) {
  environment(f) <- environment()
  D <- D[complete.cases(D[, c(all.vars(f), "id", "n"), with = FALSE])]
  wt <- as.numeric(D$n); m <- lm(f, D, weights = wt)
  u <- residuals(m)*sqrt(wt); X <- model.matrix(m)*sqrt(wt)
  nc <- uniqueN(D$id); b <- solve(crossprod(X))
  V <- b %*% crossprod(rowsum(X*u, D$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[k]); se <- unname(sqrt(diag(V))[k]); c(e, se, 2*pt(-abs(e/se), nc-1))
}

cat("=== A. AVERAGING VERSUS INTERACTION ===\n\n")
base <- readRDS(file.path(MDIR, "absorb_temporal.rds"))
B <- merge(agg(base$STUFF), KEY, by = c("id","season"))[n >= 40]
b_ll <- ll(base$STUFF$y, base$STUFF$p); b_rm <- sqrt(mean((100*(B$act-B$pred))^2))
DESC <- c(SHAPE = "per-pitch shape only", SMEAN = "+ season-mean arsenal summary",
          PRODONLY = "per-pitch shape + products", INTER = "season means + products")
cat(sprintf("  %-9s %-32s %9s %10s %13s %10s %14s\n", "features", "", "logloss", "vs stuff",
            "season RMSE", "vs stuff", "seam-shift gap"))
cat(sprintf("  %-9s %-32s %9.5f %10s %13.3f %10s %8.2f\n", "STUFF", "the baseline", b_ll, "-", b_rm, "-",
            crob(copy(B)[, `:=`(r = 100*(act-pred), gg = A5_seam)], r ~ gg, "ggTRUE")[1]))
RES <- list()
for (s in names(SETS)) {
  T2 <- data.table(id = te$id, season = te$season, y = te$whiff, p = fit1(SETS[[s]]))
  A <- merge(agg(T2), KEY, by = c("id","season"))[n >= 40]; RES[[s]] <- A
  g <- crob(copy(A)[, `:=`(r = 100*(act-pred), gg = A5_seam)], r ~ gg, "ggTRUE")
  cat(sprintf("  %-9s %-32s %9.5f %+9.2f%% %13.3f %+9.1f%% %8.2f (p=%.4f)\n", s, DESC[[s]],
              ll(T2$y,T2$p), 100*(ll(T2$y,T2$p)-b_ll)/b_ll,
              sqrt(mean((100*(A$act-A$pred))^2)),
              100*(sqrt(mean((100*(A$act-A$pred))^2))-b_rm)/b_rm, g[1], g[3]))
}

cat("\n  every archetype under the best specification:\n")
GG <- rbindlist(lapply(TAGS, function(t) {
  o <- list(archetype = LAB[[t]])
  b <- copy(B)[, `:=`(r = 100*(act-pred), gg = get(t))]
  if (b[gg == TRUE, .N] < 4) return(NULL)
  o$n <- b[gg == TRUE, .N]; o$stuff <- crob(b, r ~ gg, "ggTRUE")[1]
  for (s in names(SETS)) {
    a <- copy(RES[[s]])[, `:=`(r = 100*(act-pred), gg = get(t))]
    e <- crob(a, r ~ gg, "ggTRUE"); o[[s]] <- e[1]; o[[paste0("p_",s)]] <- e[3]
  }
  as.data.table(o)
}))
print(GG[, .(archetype, seasons = n, stuff = round(stuff,2), shape = round(SHAPE,2),
             smean = round(SMEAN,2), prodonly = round(PRODONLY,2), inter = round(INTER,2),
             p_inter = round(p_INTER,4), `closed%` = round(100*(1-INTER/stuff)))], row.names = FALSE)

## =============================================================================================
cat("\n=== B. WHAT THE UNDERPERFORMING GROUP PHYSICALLY IS ===\n\n")
S <- merge(SM[nsw >= 75], KEY, by = c("id","season"))
pop <- S[, .(m_velo = mean(m_velo), m_sep = mean(m_sep), m_axis = mean(m_axis),
             m_abs = mean(m_abs), m_dev = mean(m_dev), m_spin = mean(m_spin),
             m_arm = mean(m_arm), m_ad = mean(m_ad))]
CH <- rbindlist(c(list(cbind(group = "all changeups", seasons = nrow(S), pop)),
  lapply(TAGS, function(t) {
    D <- S[get(t) == TRUE]
    cbind(group = LAB[[t]], seasons = nrow(D),
          D[, .(m_velo = mean(m_velo), m_sep = mean(m_sep), m_axis = mean(m_axis),
                m_abs = mean(m_abs), m_dev = mean(m_dev), m_spin = mean(m_spin),
                m_arm = mean(m_arm), m_ad = mean(m_ad))])
  })))
print(CH[, .(group, seasons, velo = round(m_velo,1), separation = round(m_sep,1),
             axis_gap = round(m_axis,1), seam_dev = round(m_abs,1), signed = round(m_dev,1),
             spin = round(m_spin), slot = round(m_arm,1), arm_gap = round(m_ad,1))],
      row.names = FALSE)
cat("\n  seam deviation of the mismatch group versus everyone else:\n")
cat(sprintf("    mismatch  %.1f deg   |  all changeups  %.1f deg   |  seam-shifted set  %.1f deg\n",
            S[A8_mismatch == TRUE, mean(m_abs)], S[, mean(m_abs)], S[A5_seam == TRUE, mean(m_abs)]))
cat(sprintf("    share of the mismatch group above the population median seam deviation: %.0f%%\n",
            100*S[A8_mismatch == TRUE, mean(m_abs > median(S$m_abs))]))
