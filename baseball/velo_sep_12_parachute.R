#!/usr/bin/env Rscript

# DO THE PARACHUTE TRAITS MATTER *INSIDE* THE BIG-SEPARATION POPULATION?
#
# The frozen parachute cell - high active spin on both pitches, a spin-axis gap under 10 degrees,
# a top-third arm slot - was tested against the whole league and came back small and fragile. It
# was never tested conditional on velocity separation, which is a different question: given that a
# pitcher already takes a lot off the changeup, does matching the fastball's spin add anything?
#
# There is a plausible story for why it might. Separation buys you the swing; spin match is about
# whether the pitch looks like a fastball until it is too late. Those are complementary in
# principle and could be substitutes in practice, where a big enough velocity gap makes the
# disguise irrelevant.
#
# THIS IS A SUBGROUP ANALYSIS INSIDE AN ALREADY-SEARCHED SPACE and is treated as exploratory. The
# grid is fixed before looking: four components (changeup active spin, four-seam active spin, axis
# gap, arm slot) plus the assembled bin, against three residuals, in two separation strata. That is
# 30 tests and Bonferroni is reported for all of them. Nothing here is confirmatory.
#
# All three residuals are fit on the same anchor table with folds grouped by pitcher, so whiff,
# chase and run value are on the same footing rather than inherited from three earlier models.

suppressPackageStartupMessages({ library(data.table); library(lightgbm); library(bit64) })
set.seed(67); options(width = 210)
MDIR <- "data/statcast_model"; CACHE <- file.path(MDIR, "anchor_resid3.rds")
C <- readRDS(file.path(MDIR, "anchor_model_table2.rds")); setDT(C)

BASE <- c("release_speed","release_spin_rate","release_extension","tj_ax","az","tj_x0",
          "release_pos_z","tj_ax_diff","az_diff","sep_season")
run <- function(D, label, obj) {
  arms <- unique(D$pitcher); fa <- sample(rep(1:4, length.out = length(arms)))
  fold <- fa[match(D$pitcher, arms)]; y <- D[[label]]; p <- rep(NA_real_, nrow(D))
  for (f in 1:4) {
    tri <- which(fold != f); n <- length(tri); vi <- sample(n, floor(.12*n))
    dtr <- lgb.Dataset(as.matrix(D[tri[-vi], ..BASE]), label = y[tri[-vi]])
    dva <- lgb.Dataset.create.valid(dtr, as.matrix(D[tri[vi], ..BASE]), label = y[tri[vi]])
    m <- lgb.train(params = list(objective = obj,
                   metric = if (obj == "binary") "binary_logloss" else "l2",
                   learning_rate = .06, num_leaves = 31, min_data_in_leaf = 300,
                   feature_fraction = .8, bagging_fraction = .8, bagging_freq = 1),
                   data = dtr, nrounds = 1500, valids = list(v = dva),
                   early_stopping_rounds = 50, verbose = -1)
    p[fold == f] <- predict(m, as.matrix(D[fold == f, ..BASE]))
  }
  y - p
}
if (!file.exists(CACHE)) {
  W <- C[is_swing == 1]; Z <- C[ooz == 1]
  set.seed(67); W[, r := run(W, "whiff", "binary")]
  set.seed(67); Z[, r := run(Z, "is_swing", "binary")]
  set.seed(67); C[, r := run(C, "rv", "regression")]
  saveRDS(list(
    w = W[, .(nsw = .N, wh = 100*mean(r)), by = .(pitcher, season)],
    z = Z[, .(nz  = .N, ch = 100*mean(r)), by = .(pitcher, season)],
    v = C[, .(np  = .N, rv = 100*mean(r), sep = mean(sep_game)), by = .(pitcher, season)]), CACHE)
} else cat("(using cached residuals)\n")

L <- readRDS(CACHE)
S <- merge(merge(L$w, L$z, by = c("pitcher","season")), L$v, by = c("pitcher","season"))
S <- S[nsw >= 40 & nz >= 40]
S[, id := as.character(pitcher)]

# Frozen parachute components, from the locked spec so the definitions are not re-litigated here.
P <- readRDS(file.path(MDIR, "locked_spec.rds"))$data; setDT(P)
# y is the project's own whiff residual: velocity-adjusted, and from a model that DOES carry plate
# location and the pitch's own spin-axis components. Carried through as a control outcome, because
# the three residuals above come from a model with neither. A spin-axis trait is guaranteed to
# correlate with a residual from a model that cannot see spin axis - that is not a finding, it is
# the definition of an omitted feature. Only a trait that also moves `locked` is telling us
# something the properly specified model missed.
P <- P[league == "MLB", .(id, season, axis, arm, ec, ef, bin, locked = y)]
S <- merge(S, P, by = c("id","season"))
cat(sprintf("\n%d pitcher-seasons with all three residuals and the parachute components, %d arms\n",
            nrow(S), uniqueN(S$id)))

## ---- strata ----------------------------------------------------------------------------------------
thr <- quantile(S$sep, 2/3)
S[, hi := sep >= thr]
cat(sprintf("high-separation stratum = top third, %.1f mph or more: %d seasons (%d arms)\n",
            thr, sum(S$hi), uniqueN(S[hi == TRUE]$id)))
cat(sprintf("  mean separation: high %.1f mph, rest %.1f mph\n",
            S[hi == TRUE, mean(sep)], S[hi == FALSE, mean(sep)]))
cat("\n  do the parachute traits even look different in the high-separation group?\n")
print(S[, .(seasons = .N, ch_active_spin = round(mean(ec),3), ff_active_spin = round(mean(ef),3),
            axis_gap = round(mean(axis),1), arm_slot = round(mean(arm),1),
            in_bin = sum(bin)), by = .(group = fifelse(hi, "top third", "rest"))], row.names = FALSE)

## ---- the grid ------------------------------------------------------------------------------------
crob <- function(D, f, k) {
  D <- D[complete.cases(D[, c(all.vars(f), "id"), with = FALSE])]
  m <- lm(f, D); u <- residuals(m); X <- model.matrix(m); nc <- uniqueN(D$id)
  b <- solve(crossprod(X)); V <- b %*% crossprod(rowsum(X*u, D$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[k]); s <- unname(sqrt(diag(V))[k])
  list(est = e, se = s, p = 2*pt(-abs(e/s), nc-1), n = nrow(D), arms = nc)
}
# Components z-scored within stratum so one coefficient scale covers all four; axis is negated so
# that for every trait, higher = more parachute-like.
COMP <- c(ch_spin = "ec", ff_spin = "ef", axis_match = "axis", arm_slot = "arm")
G <- rbindlist(lapply(c(TRUE, FALSE), function(h) {
  D <- copy(S[hi == h])
  for (nm in names(COMP)) D[, (nm) := as.numeric(scale(get(COMP[[nm]]) * if (nm == "axis_match") -1 else 1))]
  rbindlist(lapply(c("wh","ch","rv","locked"), function(o) {
    a <- rbindlist(lapply(names(COMP), function(nm) {
      r <- crob(D, as.formula(paste(o, "~", nm, "+ sep")), nm)
      data.table(stratum = if (h) "top third" else "rest", outcome = o, trait = nm,
                 est = r$est, se = r$se, p = r$p, n = r$n) }))
    r <- crob(D, as.formula(paste(o, "~ bin + sep")), "binTRUE")
    rbind(a, data.table(stratum = if (h) "top third" else "rest", outcome = o, trait = "LOCKED BIN",
                        est = r$est, se = r$se, p = r$p, n = r$n)) }))
}))
G[, `:=`(bonf = p < .05/.N, raw = p < .05)]
OUTN <- c(wh = "whiff", ch = "chase", rv = "run value", locked = "whiff (loc+spin model)")
cat(sprintf("\n=== %d tests. Coefficients are points of residual per 1 SD of the trait ===\n", nrow(G)))
cat("   every trait oriented so higher = more parachute-like (axis gap negated)\n")
for (st in c("top third","rest")) {
  cat(sprintf("\n  ############ separation: %s\n", st))
  X <- G[stratum == st]
  print(X[order(outcome, p), .(outcome = OUTN[outcome], trait, est = round(est,3),
          se = round(se,3), p = round(p,4), raw, bonf)], row.names = FALSE)
}
cat(sprintf("\n  %d of %d clear p<.05 (chance alone gives %.1f); %d clear Bonferroni at p<%.4f\n",
            sum(G$raw), nrow(G), .05*nrow(G), sum(G$bonf), .05/nrow(G)))

## ---- does anything INTERACT with separation? --------------------------------------------------------
cat("\n=== interactions: does a parachute trait pay off differently when separation is large? ===\n")
D <- copy(S); D[, zsep := as.numeric(scale(sep))]
for (nm in names(COMP)) D[, (nm) := as.numeric(scale(get(COMP[[nm]]) * if (nm == "axis_match") -1 else 1))]
IT <- rbindlist(lapply(c("wh","ch","rv","locked"), function(o)
  rbindlist(lapply(c(names(COMP), "bin"), function(nm) {
    k <- if (nm == "bin") "binTRUE:zsep" else paste0(nm, ":zsep")
    f <- as.formula(paste(o, "~", nm, "* zsep"))
    r <- crob(D, f, k)
    data.table(outcome = OUTN[o], trait = nm, est = round(r$est,3), se = round(r$se,3),
               p = round(r$p,4)) }))))
print(IT[order(p)], row.names = FALSE)
cat(sprintf("  %d of %d interactions clear p<.05 (chance gives %.1f)\n",
            sum(IT$p < .05), nrow(IT), .05*nrow(IT)))
fwrite(G, file.path(MDIR, "article_assets", "ext_parachute_in_sep.csv"))
cat("\nwrote ext_parachute_in_sep.csv\n")
