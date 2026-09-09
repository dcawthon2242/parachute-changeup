#!/usr/bin/env Rscript

# High velo gap + poor tunnel: do any of these pairs upset timing and still lose?

suppressPackageStartupMessages({ library(data.table) })
options(width = 220)

P <- fread("data/statcast_2026/velo_gap_tunnel_pairs.csv")
cat(sprintf("Pairs: %d\n\n", nrow(P)))

# Residualize if the file is the older write without the _a columns
need <- c("xrv_a","tscore_a","whiff_a","tun_a","com_a","vg_a")
have <- names(P)
if (!"xrv_a" %in% have) {
  P[, `:=`(
    xrv_a    = residuals(lm(xrv ~ factor(pitch_type) + use)),
    tscore_a = residuals(lm(tscore ~ factor(pitch_type) + use)),
    whiff_a  = residuals(lm(whiff ~ factor(pitch_type) + use)),
    tun_a    = residuals(lm(tunnel_dist ~ factor(pitch_type) + use)),
    com_a    = residuals(lm(commit_dist ~ factor(pitch_type) + use)),
    vg_a     = residuals(lm(velo_gap ~ factor(pitch_type) + use))
  )]
} else {
  if (!"tun_a" %in% have) P[, tun_a := get(if ("tunnel_dist_a" %in% have) "tunnel_dist_a" else "tun_a")]
  if (!"com_a" %in% have) P[, com_a := get(if ("commit_dist_a" %in% have) "commit_dist_a" else "com_a")]
  if (!"vg_a"  %in% have && "velo_gap_a" %in% have) P[, vg_a := velo_gap_a]
}

# 2x2 on type-adjusted gap and type-adjusted spatial tunnel
P[, hi_gap := velo_gap >= quantile(velo_gap, 2/3)]
P[, lo_tun := tun_a    >= quantile(tun_a,    2/3)]   # high residual = worse tunnel
P[, cell := paste0(fifelse(hi_gap, "big gap", "small gap"), " + ",
                   fifelse(lo_tun, "wide tunnel", "tight tunnel"))]

cat("############ 1. The 2x2, type-adjusted outcomes ############\n\n")
CELL <- P[, .(
  pairs     = .N,
  velo_gap  = mean(velo_gap),
  tunnel    = mean(tunnel_dist),
  commit    = mean(commit_dist),
  tdev      = mean(tdev),
  tscore    = mean(tscore),
  tscore_a  = mean(tscore_a),
  whiff     = mean(whiff),
  whiff_a   = mean(whiff_a),
  xrv       = mean(xrv),
  xrv_a     = mean(xrv_a),
  xrv_wt    = weighted.mean(xrv, pitches)
), by=cell]
setorder(CELL, -xrv_a)
print(CELL[, lapply(.SD, function(x) if (is.numeric(x)) round(x,3) else x)], row.names=FALSE)

# Same 2x2 using time-based commit (the tunnel that actually absorbs velo gap)
P[, lo_com := com_a >= quantile(com_a, 2/3)]
P[, cell2 := paste0(fifelse(hi_gap, "big gap", "small gap"), " + ",
                    fifelse(lo_com, "late reveal", "early overlay"))]
cat("\n  Same 2x2, but the tunnel axis is time-based commit residual:\n")
C2 <- P[, .(pairs=.N, tscore_a=mean(tscore_a), whiff_a=mean(whiff_a),
            xrv_a=mean(xrv_a), tdev=mean(tdev)), by=cell2]
print(C2[order(-xrv_a)][, lapply(.SD, function(x) if (is.numeric(x)) round(x,3) else x)],
      row.names=FALSE)

# =============================================================================
# 2. The specific examples: big gap, wide tunnel, good timing, bad xRV
# =============================================================================
cat("\n############ 2. Pairs that upset timing and still lose ############\n\n")

# Target cell
Q <- P[hi_gap == TRUE & lo_tun == TRUE]
cat(sprintf("  Big-gap + wide-tunnel cell: %d pairs\n", nrow(Q)))
cat(sprintf("  Of those, timing residual > 0 (better than type-expected): %d\n",
            nrow(Q[tscore_a > 0])))
cat(sprintf("  ...and xRV residual < 0 (worse than type-expected): %d\n",
            nrow(Q[tscore_a > 0 & xrv_a < 0])))
cat(sprintf("  ...and xRV residual below -0.4 (clearly bad): %d\n\n",
            nrow(Q[tscore_a > 0 & xrv_a < -0.4])))

# How unusual is the timing-vs-performance split in this cell vs others?
P[, win_time := tscore_a > 0]
P[, lose_xrv := xrv_a < 0]
cat("  Share of pairs that win timing AND lose xRV, by cell:\n")
print(P[, .(pairs=.N, win_time=round(100*mean(win_time),1),
            lose_xrv=round(100*mean(lose_xrv),1),
            both=round(100*mean(win_time & lose_xrv),1)), by=cell][order(-both)],
      row.names=FALSE)

SHOW <- c("player_name","fb_type","pitch_type","use","pitches","swings",
          "velo_gap","tunnel_dist","commit_dist","plate_dist","tdev",
          "tscore","tscore_a","whiff","xrv","xrv_a")

# Require enough pitches that xRV is not a coin flip
EX <- Q[tscore_a > 0 & xrv_a < -0.4 & pitches >= 200]
cat(sprintf("\n  Examples (tscore_a > 0, xrv_a < -0.4, 200+ pitches): %d pairs\n", nrow(EX)))
cat("  Sorted by worst xRV residual:\n")
print(EX[order(xrv_a), ..SHOW][, lapply(.SD, function(x) if (is.numeric(x)) round(x,2) else x)],
      row.names=FALSE)

cat("\n  Contrast: same cell, win timing AND win xRV (the ones the tradeoff did not punish):\n")
OK <- Q[tscore_a > 0 & xrv_a > 0.4 & pitches >= 200]
print(OK[order(-xrv_a)][1:12, ..SHOW][
  , lapply(.SD, function(x) if (is.numeric(x)) round(x,2) else x)], row.names=FALSE)

# =============================================================================
# 3. Is the "win timing, lose xRV" group actually different on tunnel?
# =============================================================================
cat("\n############ 3. Within big-gap pairs: does a worse tunnel predict the split? ############\n\n")
BG <- P[hi_gap == TRUE]
BG[, split := fifelse(tscore_a > 0 & xrv_a < 0, "win timing, lose xRV",
               fifelse(tscore_a > 0 & xrv_a > 0, "win both",
               fifelse(tscore_a < 0 & xrv_a < 0, "lose both", "lose timing, win xRV")))]
print(BG[, .(pairs=.N,
             velo_gap=round(mean(velo_gap),1),
             tunnel=round(mean(tunnel_dist),2),
             tun_a=round(mean(tun_a),2),
             commit=round(mean(commit_dist),2),
             com_a=round(mean(com_a),2),
             tscore_a=round(mean(tscore_a),3),
             xrv_a=round(mean(xrv_a),3),
             whiff_a=round(mean(whiff_a),2)), by=split][order(-tscore_a)],
      row.names=FALSE)

cat("\n  Logistic: among big-gap pairs that win timing, does a wider tunnel predict losing xRV?\n")
W <- BG[tscore_a > 0]
m <- glm(I(xrv_a < 0) ~ tun_a + com_a + velo_gap + factor(pitch_type) + use, data=W, family=binomial)
co <- as.data.table(summary(m)$coefficients, keep.rownames="term")
setnames(co, c("term","b","se","z","p"))
print(co[!grepl("factor|Intercept", term)][
  , .(term, b=round(b,3), z=round(z,2), p=round(p,4))], row.names=FALSE)
cat(sprintf("  n = %d timing-winning big-gap pairs; %d of them lose on xRV\n",
            nrow(W), sum(W$xrv_a < 0)))

# =============================================================================
# 4. Pitch-type mix of the examples, so we know what "poor tunnel" looks like
# =============================================================================
cat("\n############ 4. What pitch types show up in the losing examples? ############\n\n")
print(EX[, .(pairs=.N, mean_gap=round(mean(velo_gap),1),
             mean_tun=round(mean(tunnel_dist),1),
             mean_ts_a=round(mean(tscore_a),2),
             mean_xrv_a=round(mean(xrv_a),2)), by=.(fb_type, pitch_type)][order(-pairs)],
      row.names=FALSE)

fwrite(EX, "data/statcast_2026/poor_tunnel_high_gap_examples.csv")
cat("\nWrote data/statcast_2026/poor_tunnel_high_gap_examples.csv\n")
