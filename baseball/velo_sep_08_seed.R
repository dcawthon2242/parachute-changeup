#!/usr/bin/env Rscript

# IS THE ANCHOR IMPROVEMENT REAL, OR ONE LUCKY FOLD DRAW?
#
# The outing anchor beat the season constant by 0.000244 of RMSE and a third of a percent of
# season RMSE on a single fold assignment. Differences that small are exactly the size of fold
# noise, so the comparison is repeated across seeds with the three anchors sharing an identical
# fold assignment within each seed - paired, so the seed-to-seed variation cancels out of the
# difference even though it dominates the levels.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
options(width = 205)
MDIR <- "data/statcast_model"
C <- readRDS(file.path(MDIR, "anchor_model_table.rds")); setDT(C)
W <- C[is_swing == 1]

BASE <- c("release_speed","release_spin_rate","release_extension","tj_ax","az","tj_x0",
          "release_pos_z","tj_ax_diff","az_diff")
SETS <- list(season = c(BASE,"sep_season"), game = c(BASE,"sep_game"), recent = c(BASE,"sep_recent"))
SEEDS <- c(53, 101, 202, 303, 404)

one <- function(D, feats, fold, label, obj) {
  y <- D[[label]]; p <- rep(NA_real_, nrow(D))
  for (f in unique(fold)) {
    tri <- which(fold != f); n <- length(tri); vi <- sample(n, floor(.12*n))
    dtr <- lgb.Dataset(as.matrix(D[tri[-vi], ..feats]), label = y[tri[-vi]])
    dva <- lgb.Dataset.create.valid(dtr, as.matrix(D[tri[vi], ..feats]), label = y[tri[vi]])
    m <- lgb.train(params = list(objective = obj, metric = "binary_logloss", learning_rate = .06,
                   num_leaves = 31, min_data_in_leaf = 300, feature_fraction = .8,
                   bagging_fraction = .8, bagging_freq = 1), data = dtr, nrounds = 1500,
                   valids = list(v = dva), early_stopping_rounds = 50, verbose = -1)
    p[fold == f] <- predict(m, as.matrix(D[fold == f, ..feats]))
  }
  S <- data.table(k = paste(D$pitcher, D$season), y = y, p = p)[, .(n = .N, ay = mean(y), ap = mean(p)), by = k]
  c(rmse = sqrt(mean((y - p)^2)), r2 = 1 - var(y - p)/var(y),
    srmse = sqrt(weighted.mean((S$ay - S$ap)^2, S$n)))
}

R <- rbindlist(lapply(SEEDS, function(sd0) {
  set.seed(sd0)
  arms <- unique(W$pitcher); fa <- sample(rep(1:4, length.out = length(arms)))
  fold <- fa[match(W$pitcher, arms)]
  rbindlist(lapply(names(SETS), function(nm) {
    set.seed(sd0); v <- one(W, SETS[[nm]], fold, "whiff", "binary")
    data.table(seed = sd0, anchor = nm, rmse = v[["rmse"]], r2 = v[["r2"]], srmse = v[["srmse"]]) }))
}))

cat("=== whiff, five paired fold assignments ===\n")
print(dcast(R, seed ~ anchor, value.var = "srmse")[, lapply(.SD, function(x)
       if (is.numeric(x) && max(x) < 1) round(x, 5) else x)], row.names = FALSE)

cat("\n=== paired differences against the season anchor (negative = the anchor helps) ===\n")
P <- merge(R[anchor != "season"], R[anchor == "season", .(seed, s_rmse = rmse, s_r2 = r2, s_srmse = srmse)],
           by = "seed")
P[, `:=`(d_rmse = rmse - s_rmse, d_srmse = 100*(srmse - s_srmse)/s_srmse, d_r2 = r2 - s_r2)]
print(P[, .(mean_d_rmse = round(mean(d_rmse), 6), sd = round(sd(d_rmse), 6),
            mean_d_season_rmse_pct = round(mean(d_srmse), 3), sd_pct = round(sd(d_srmse), 3),
            mean_d_r2 = round(mean(d_r2), 5), wins = sum(d_rmse < 0), n = .N), by = anchor],
      row.names = FALSE)
for (a in c("game","recent")) {
  x <- P[anchor == a]
  t1 <- t.test(x$d_rmse); t2 <- t.test(x$d_srmse)
  cat(sprintf("  %-7s paired t on RMSE: p = %.4f | on season RMSE: p = %.4f\n", a, t1$p.value, t2$p.value))
}
cat(sprintf("\n  for scale, seed-to-seed spread of the season anchor's own season RMSE: %.5f\n",
            sd(R[anchor == "season"]$srmse)))
fwrite(R, file.path(MDIR, "article_assets", "ext_anchor_seeds.csv"))
