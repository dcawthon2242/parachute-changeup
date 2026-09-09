#!/usr/bin/env Rscript

# Decompose the whiff MISS into a timing axis (y = depth toward the pitcher) and a
# plane axis (x = side), using Statcast's intercept fields, then tie each axis back
# to the earlier findings (velo kill vs spin-look similarity vs kick changes).
#
# NOTE on the data feed: Statcast publishes scalar miss_distance plus the swing
# INTERCEPT point relative to the batter in x (side) and y (depth). There is no
# published vertical (z) intercept, so the "2D plane" here is represented by the
# horizontal (x) axis; the y axis is the timing/depth dimension.

suppressPackageStartupMessages({ library(data.table) })
dir <- file.path("data","statcast_2026")

dt <- fread(file.path(dir,"statcast_2026_all.csv"), showProgress=FALSE, select=c(
  "pitcher","player_name","pitch_type","game_type","description","miss_distance",
  "intercept_ball_minus_batter_pos_x_inches","intercept_ball_minus_batter_pos_y_inches"))
setnames(dt, 7:8, c("side","depth"))
whiff <- c("swinging_strike","swinging_strike_blocked","foul_tip","missed_bunt")
w <- dt[game_type=="R" & description %in% whiff & is.finite(miss_distance) & is.finite(side) & is.finite(depth)]

cat("=== Pitch-level: how much of MISS distance is timing (depth/y) vs plane (side/x)? ===\n")
cat(sprintf("  cor(miss, depth [timing/y]) = %+.3f\n", cor(w$miss_distance, w$depth)))
cat(sprintf("  cor(miss, side  [plane /x]) = %+.3f\n", cor(w$miss_distance, w$side)))
f <- lm(scale(miss_distance) ~ scale(depth) + scale(side), data=w)
cat("  standardized betas (share of the miss each axis carries):\n"); print(round(coef(f),3))
cat(sprintf("  model R2 = %.3f\n", summary(f)$r.squared))

# aggregate per pitcher x pitch type (signed axis means on whiffs)
ag <- w[, .(nwh=.N, miss=mean(miss_distance), depth=mean(depth), side=mean(side)),
        by=.(pitcher, pitch_type)]
sec <- fread(file.path(dir,"spin_lookalike_2026.csv"), showProgress=FALSE)
m <- merge(sec, ag, by=c("pitcher","pitch_type"))[nwh>=15]

pc <- function(x,y){ok<-is.finite(x)&is.finite(y); if(sum(ok)<15) return(c(NA,NA,sum(ok)))
  ct<-cor.test(x[ok],y[ok]); c(ct$estimate, ct$p.value, sum(ok))}

cat("\n=== Does the VELO kill push the miss onto the TIMING (depth) axis, and spin onto the PLANE (side) axis? ===\n")
for (grp in list(CH="CH", FS="FS", breaking=c("SL","ST","CU","KC"))) {
  d <- m[pitch_type %in% grp]; lab <- if(length(grp)==1) grp else "breaking"
  cat(sprintf("  -- %s (n=%d) --\n", lab, nrow(d)))
  for (ax in c("depth","side")) {
    axlab <- if(ax=="depth") "timing/y" else "plane/x"
    for (v in c("velo_diff","spin_diff","axis_diff","eff_diff")) {
      r <- pc(d[[v]], d[[ax]])
      cat(sprintf("     %-6s[%-8s] ~ %-9s r=%+.2f p=%.2g\n", ax, axlab, v, r[1], r[2]))
    }
  }
}

cat("\n=== Which axis carries the overall miss for each pitch family? ===\n")
print(m[, .(n=.N, miss=round(mean(miss),1), depth=round(mean(depth),1),
  side=round(mean(side),1),
  cor_miss_depth=round(cor(miss,depth,use="complete.obs"),2),
  cor_miss_side =round(cor(miss,side, use="complete.obs"),2),
  cor_miss_velo =round(cor(miss,velo_diff,use="complete.obs"),2)),
  by=.(fam=fifelse(pitch_type=="CH","CH",fifelse(pitch_type=="FS","FS",
       fifelse(pitch_type %in% c("SL","ST","CU","KC"),"breaking","other"))))][order(fam)])

cat("\n=== Kick changes: do they miss on TIMING (depth) or on PLANE (side)? ===\n")
ch <- m[pitch_type=="CH"]
ch[, kick := spin<=1600 & ivb<=2]
print(ch[, .(n=.N, spin=round(mean(spin)), ivb=round(mean(ivb),1), Dvelo=round(mean(velo_diff),1),
  miss=round(mean(miss),1), timing_depth=round(mean(depth),1), plane_side=round(mean(side),1)),
  by=kick])

fwrite(m, file.path(dir,"miss_decomp_2026.csv"))
cat(sprintf("\nWrote %s\n", file.path(dir,"miss_decomp_2026.csv")))
