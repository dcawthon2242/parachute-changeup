#!/usr/bin/env Rscript

# Do big-gap + tight-tunnel pairs actually outperform — and do they
# "overperform xRV" (realized beating expected), or just post a higher xRV?

suppressPackageStartupMessages({ library(data.table) })
options(width = 210)

P <- fread("data/statcast_2026/velo_gap_tunnel_pairs.csv")
P[, tun_a := tunnel_dist_a]
P[, com_a := commit_dist_a]
P[, hi_gap := velo_gap >= quantile(velo_gap, 2/3)]
P[, tight  := tun_a    <= quantile(tun_a,    1/3)]
P[, wide   := tun_a    >= quantile(tun_a,    2/3)]
P[, cell := fifelse(hi_gap & tight, "big gap + tight tunnel",
            fifelse(hi_gap & wide,  "big gap + wide tunnel",
            fifelse(hi_gap,         "big gap + mid tunnel",
            fifelse(tight,          "small gap + tight tunnel",
            fifelse(wide,           "small gap + wide tunnel",
                                    "small gap + mid tunnel")))))]
P[, over := rv - xrv]   # realized minus expected; + = overperformed xRV

cat(sprintf("Pairs: %d\n\n", nrow(P)))

cat("############ 1. Cell means ############\n\n")
CELL <- P[, .(
  pairs   = .N,
  pitchers= uniqueN(pitcher),
  pitches = sum(pitches),
  velo_gap= mean(velo_gap),
  tunnel  = mean(tunnel_dist),
  tscore_a= mean(tscore_a),
  whiff_a = mean(whiff_a),
  xrv     = mean(xrv),
  xrv_a   = mean(xrv_a),
  xrv_wt  = weighted.mean(xrv, pitches),
  xrv_a_wt= weighted.mean(xrv_a, pitches),
  rv      = mean(rv),
  over    = mean(over),
  over_wt = weighted.mean(over, pitches),
  pct_pos = 100*mean(xrv_a > 0)
), by=cell]
setorder(CELL, -xrv_a)
print(CELL[, lapply(.SD, function(x) if (is.numeric(x)) round(x,3) else x)], row.names=FALSE)

# Head-to-head: the two cells in the question
A <- P[hi_gap == TRUE & tight == TRUE]
B <- P[hi_gap == TRUE & wide  == TRUE]
C <- P[!(hi_gap == TRUE & tight == TRUE)]

cat("\n############ 2. Head-to-head tests ############\n\n")
tt <- function(x, y, lab) {
  t <- t.test(x, y)
  cat(sprintf("  %-40s  %.3f vs %.3f   diff = %+.3f   p = %.3f   n = %d / %d\n",
              lab, mean(x), mean(y), mean(x)-mean(y), t$p.value, length(x), length(y)))
}
cat("  big-gap + tight  vs  big-gap + wide:\n")
tt(A$xrv_a,    B$xrv_a,    "xRV residual")
tt(A$xrv,      B$xrv,      "raw xRV")
tt(A$tscore_a, B$tscore_a, "timing residual")
tt(A$whiff_a,  B$whiff_a,  "whiff residual")
tt(A$over,     B$over,     "realized minus xRV")
tt(A$rv,       B$rv,       "realized RV")

cat("\n  big-gap + tight  vs  everyone else:\n")
tt(A$xrv_a,    C$xrv_a,    "xRV residual")
tt(A$over,     C$over,     "realized minus xRV")

cat("\n  Share with xRV residual > 0:\n")
cat(sprintf("     big gap + tight  %.1f%%   (n = %d)\n", 100*mean(A$xrv_a>0), nrow(A)))
cat(sprintf("     big gap + wide   %.1f%%   (n = %d)\n", 100*mean(B$xrv_a>0), nrow(B)))
cat(sprintf("     everyone else    %.1f%%   (n = %d)\n", 100*mean(C$xrv_a>0), nrow(C)))

cat("\n  Share that overperform xRV (realized > expected):\n")
cat(sprintf("     big gap + tight  %.1f%%\n", 100*mean(A$over>0)))
cat(sprintf("     big gap + wide   %.1f%%\n", 100*mean(B$over>0)))
cat(sprintf("     everyone else    %.1f%%\n", 100*mean(C$over>0)))

# =============================================================================
# 3. Continuous version: among big-gap pairs, does tighter tunnel raise xRV?
# =============================================================================
cat("\n############ 3. Continuous: tunnel residual among big-gap pairs ############\n\n")
BG <- P[hi_gap == TRUE]
m1 <- lm(xrv_a ~ tun_a + velo_gap + use + factor(pitch_type), data=BG)
m2 <- lm(over  ~ tun_a + velo_gap + use + factor(pitch_type), data=BG)
co <- function(m, lab) {
  s <- summary(m)$coefficients
  cat(sprintf("  -- %s | n = %d, R2 = %.3f --\n", lab, nobs(m), summary(m)$r.squared))
  for (v in c("tun_a","velo_gap","use")) {
    cat(sprintf("     %-10s  b = %+.4f   t = %+.2f   p = %.3f\n",
                v, s[v,1], s[v,3], s[v,4]))
  }
}
co(m1, "xRV residual")
co(m2, "realized minus xRV")

# terciles of tunnel within big-gap only
BG[, tq := cut(tun_a, quantile(tun_a, 0:3/3), include.lowest=TRUE,
               labels=c("tightest third","middle","widest third"))]
cat("\n  Big-gap pairs only, by their own tunnel tercile:\n")
print(BG[, .(pairs=.N, tunnel=round(mean(tunnel_dist),2),
             tscore_a=round(mean(tscore_a),3),
             xrv_a=round(mean(xrv_a),3),
             xrv_wt=round(weighted.mean(xrv, pitches),3),
             over=round(mean(over),3),
             pct_pos=round(100*mean(xrv_a>0),1)), by=tq][order(tq)], row.names=FALSE)

# =============================================================================
# 4. Pitcher level: a pitcher "has the combo" if any pair qualifies
# =============================================================================
cat("\n############ 4. Pitcher-level (any qualifying pair) ############\n\n")
# pitcher total xRV, not pair-level (avoids double counting)
PIT <- P[, .(
  n_pairs = .N,
  has_combo = any(hi_gap & tight),
  has_wide  = any(hi_gap & wide),
  pitches   = sum(pitches),
  xrv_a     = weighted.mean(xrv_a, pitches),
  xrv       = weighted.mean(xrv, pitches),
  over      = weighted.mean(over, pitches)
), by=pitcher]
cat(sprintf("  Pitchers with at least one big-gap + tight pair: %d / %d\n",
            sum(PIT$has_combo), nrow(PIT)))
tt(PIT[has_combo==TRUE]$xrv_a, PIT[has_combo==FALSE]$xrv_a, "pitcher xRV residual")
tt(PIT[has_combo==TRUE]$over,  PIT[has_combo==FALSE]$over,  "pitcher realized-xRV")
cat(sprintf("  Share with xRV residual > 0:  combo %.1f%%   others %.1f%%\n",
            100*mean(PIT[has_combo==TRUE]$xrv_a>0),
            100*mean(PIT[has_combo==FALSE]$xrv_a>0)))

# exclusive groups so combo and wide-gap don't overlap
PIT[, grp := fifelse(has_combo, "has big+tight",
             fifelse(has_wide,  "has big+wide only", "neither"))]
cat("\n  Exclusive pitcher groups:\n")
print(PIT[, .(pitchers=.N,
              xrv_a=round(mean(xrv_a),3),
              xrv_wt=round(weighted.mean(xrv, pitches),3),
              over=round(mean(over),3),
              pct_pos=round(100*mean(xrv_a>0),1)), by=grp][order(-xrv_a)],
      row.names=FALSE)

# =============================================================================
# 5. Magnitude in runs
# =============================================================================
cat("\n############ 5. Practical magnitude ############\n\n")
cat(sprintf("  Pair-level xRV residual: tight %.3f vs wide %.3f  = %+.3f runs / 100 pitches\n",
            mean(A$xrv_a), mean(B$xrv_a), mean(A$xrv_a)-mean(B$xrv_a)))
cat(sprintf("  On a 200-pitch secondary: %+.2f runs over a season\n",
            2*(mean(A$xrv_a)-mean(B$xrv_a))))
cat(sprintf("  Weighted xRV (raw):      tight %.3f vs wide %.3f  = %+.3f\n",
            weighted.mean(A$xrv, A$pitches), weighted.mean(B$xrv, B$pitches),
            weighted.mean(A$xrv, A$pitches)-weighted.mean(B$xrv, B$pitches)))
cat(sprintf("  Realized minus xRV:      tight %.3f vs wide %.3f  = %+.3f  (luck channel)\n",
            mean(A$over), mean(B$over), mean(A$over)-mean(B$over)))
