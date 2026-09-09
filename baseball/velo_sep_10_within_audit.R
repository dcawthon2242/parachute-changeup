#!/usr/bin/env Rscript

# THE BETWEEN-OUTING TERM WON. AUDIT IT BEFORE BELIEVING IT.
#
# The prediction going in was that within-outing repeatability is the real skill and game-to-game
# drift is a workload story that should not matter. The data said the reverse: with both terms in
# the model, between-outing drift survives at p = .005 on whiff and p < .0001 on chase while
# within-outing repeatability fades. That is a surprise, and surprises in this project have a poor
# track record, so:
#
#   A  SCALE. The two terms have different spreads, so raw per-mph coefficients are not comparable.
#      Restate per standard deviation of each.
#   B  ESTIMATION NOISE. sd_between is computed from a handful of game means, so it is mechanically
#      noisier for arms with few outings, and noisier estimates attenuate differently. Control for
#      the number of outings and re-check.
#   C  ROLE. A reliever's outings are one inning apart in leverage and a starter's are not. Check
#      that this is not starters versus relievers wearing a consistency costume.
#   D  SEASON SHAPE. If drift were really about fatigue, it should look like a trend across the
#      year rather than scatter. A pitcher's changeup-gap trend is separable from his scatter.

suppressPackageStartupMessages({ library(data.table) })
options(width = 205)
MDIR <- "data/statcast_model"
C <- readRDS(file.path(MDIR, "anchor_model_table2.rds")); setDT(C)
S <- readRDS(file.path(MDIR, "velo_sep_within_seasons.rds")); setDT(S)

G <- C[, .(games = uniqueN(game_pk), np = .N), by = .(pitcher, season)]
GM <- C[, .(mg = mean(sep_game), ng = .N), by = .(pitcher, season, game_pk)][ng >= 5]
GM[, gi := seq_len(.N), by = .(pitcher, season)]
TR <- GM[, .(trend = if (.N >= 5) coef(lm(mg ~ gi))[2] else NA_real_,
             scatter = if (.N >= 5) sd(residuals(lm(mg ~ gi))) else NA_real_), by = .(pitcher, season)]
S <- merge(merge(S, G, by = c("pitcher","season")), TR, by = c("pitcher","season"))
S[, ch_per_game := np/games]

crob <- function(f, D, k) {
  D <- D[complete.cases(D[, c(all.vars(f), "id"), with = FALSE])]
  m <- lm(f, D); u <- residuals(m); X <- model.matrix(m); nc <- uniqueN(D$id)
  b <- solve(crossprod(X)); V <- b %*% crossprod(rowsum(X*u, D$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[k]); s <- unname(sqrt(diag(V))[k])
  c(est = e, se = s, p = 2*pt(-abs(e/s), nc-1))
}
show <- function(lab, r, unit = "") cat(sprintf("    %-34s %+.3f (se %.3f) | p %.4f%s\n",
                                                lab, r["est"], r["se"], r["p"], unit))

cat(sprintf("%d pitcher-seasons, %d arms\n", nrow(S), uniqueN(S$id)))
cat(sprintf("outings per season: median %d (IQR %d-%d) | changeups per outing: median %.1f\n",
            median(S$games), quantile(S$games,.25), quantile(S$games,.75), median(S$ch_per_game)))

## ---- A. same scale -----------------------------------------------------------------------------
cat("\n=== A. per standard deviation of each term, not per mph ===\n")
S[, `:=`(zw = as.numeric(scale(sd_within)), zb = as.numeric(scale(sd_between)))]
cat(sprintf("  spread: sd_within sd = %.3f mph | sd_between sd = %.3f mph\n",
            sd(S$sd_within), sd(S$sd_between)))
for (v in c("wh_aware","ch_aware")) {
  cat(sprintf("  -- %s, both terms in --\n", v))
  f <- as.formula(paste(v, "~ zw + zb + sep"))
  show("within-outing, per 1 SD", crob(f, S, "zw"), " points")
  show("between-outing, per 1 SD", crob(f, S, "zb"), " points")
}

## ---- B and C. estimation noise and role ----------------------------------------------------------
cat("\n=== B/C. controls: number of outings, workload, and starter-vs-reliever shape ===\n")
for (v in c("wh_aware","ch_aware")) {
  cat(sprintf("  -- %s, between-outing term per 1 SD --\n", v))
  show("bare",                       crob(as.formula(paste(v,"~ zw + zb + sep")), S, "zb"))
  show("+ log outings",              crob(as.formula(paste(v,"~ zw + zb + sep + log(games)")), S, "zb"))
  show("+ log outings, ch per outing", crob(as.formula(paste(v,"~ zw + zb + sep + log(games) + ch_per_game")), S, "zb"))
  show("starters only (>=4 ch/outing)",
       crob(as.formula(paste(v,"~ zw + zb + sep")), S[ch_per_game >= 4], "zb"))
  show("relievers only (<4 ch/outing)",
       crob(as.formula(paste(v,"~ zw + zb + sep")), S[ch_per_game < 4], "zb"))
}

## ---- D. is drift a trend or scatter? -----------------------------------------------------------
cat("\n=== D. drift decomposed into a season trend and scatter around it ===\n")
cat(sprintf("  median |trend| = %.4f mph per outing | median scatter = %.3f mph\n",
            median(abs(S$trend), na.rm = TRUE), median(S$scatter, na.rm = TRUE)))
cat(sprintf("  cor(scatter, sd_between) = %+.3f - scatter is nearly all of the between term\n",
            cor(S$scatter, S$sd_between, use = "complete.obs")))
S[, `:=`(zt = as.numeric(scale(abs(trend))), zs = as.numeric(scale(scatter)))]
for (v in c("wh_aware","ch_aware")) {
  cat(sprintf("  -- %s --\n", v))
  f <- as.formula(paste(v, "~ zt + zs + zw + sep"))
  show("season trend magnitude, per 1 SD", crob(f, S, "zt"))
  show("outing-to-outing scatter, per 1 SD", crob(f, S, "zs"))
}

## ---- practical size --------------------------------------------------------------------------------
cat("\n=== practical size of the between-outing term ===\n")
q <- quantile(S$sd_between, c(.1,.9))
for (v in c("wh_aware","ch_aware")) {
  b <- crob(as.formula(paste(v,"~ zw + zb + sep")), S, "zb")
  cat(sprintf("  %-9s p10 to p90 of drift (%.2f to %.2f mph) is worth %.2f points; residual SD is %.2f\n",
              v, q[1], q[2], abs(b["est"])*(q[2]-q[1])/sd(S$sd_between), sd(S[[v]])))
}
print(S[, .(seasons = .N, drift = round(mean(sd_between),2),
            whiff_over = round(mean(wh_aware),2), chase_over = round(mean(ch_aware),2),
            se_ch = round(sd(ch_aware)/sqrt(.N),2)),
        by = .(bin = cut(sd_between, quantile(sd_between, 0:5/5),
               labels = c("steadiest","2","3","4","most drift"), include.lowest = TRUE))][order(bin)],
      row.names = FALSE)
