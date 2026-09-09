#!/usr/bin/env Rscript

# Two questions:
#   (1) Do miss distance (on whiffs) and xwOBAcon (on balls in play) relate to each other?
#       They live on mutually exclusive pitch subsets, so a per-pitch correlation does not
#       exist by construction. Section 2 links them through the shared timing axis instead;
#       Section 3 correlates them as aggregates, per pitcher and per pitcher x pitch type.
#   (2) Is the underlying swing-timing signal a repeatable skill, and whose skill is it?
#       Section 4 measures split-half (within 2026) and year-over-year (2025 -> 2026)
#       repeatability for every timing metric, separately for pitchers and for batters.

suppressPackageStartupMessages({ library(data.table); library(splines) })
set.seed(7)

cols <- c("game_year","game_type","player_name","pitcher","batter","stand","pitch_type",
          "description","bb_type","balls","strikes","plate_x","plate_z","sz_top","sz_bot",
          "release_speed","launch_speed","estimated_woba_using_speedangle","delta_run_exp",
          "miss_distance","bat_speed","swing_length",
          "intercept_ball_minus_batter_pos_y_inches")

dt <- rbindlist(lapply(c(2025, 2026), function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress = FALSE, select = cols)))
setnames(dt, "intercept_ball_minus_batter_pos_y_inches", "depth")
dt <- dt[game_type == "R" & pitch_type != "" & balls <= 3 & strikes <= 2]

# =============================================================================
# 1. Coverage: why there is no per-pitch correlation to compute
# =============================================================================
cat("############ 1. Where each field exists ############\n")
cov <- dt[description %in% c("swinging_strike","swinging_strike_blocked","foul_tip",
                             "foul","hit_into_play"),
  .(pitches=.N,
    miss_dist_pct = round(100*mean(is.finite(miss_distance)),1),
    xwobacon_pct  = round(100*mean(is.finite(estimated_woba_using_speedangle)),1),
    depth_pct     = round(100*mean(is.finite(depth)),1)), by=description][order(-pitches)]
print(cov)
cat("\n  miss_distance and xwOBAcon never co-occur on the same pitch, so a per-pitch\n")
cat("  correlation is undefined. Everything below is either mediated by the timing axis\n")
cat("  (Section 2) or computed on aggregates (Section 3).\n")

cat("\n  miss_distance distribution on whiffs (inches):\n")
md <- dt[is.finite(miss_distance), miss_distance]
print(round(quantile(md, c(.01,.05,.25,.5,.75,.95,.99)), 2))
cat(sprintf("  n=%s  mean=%.2f  sd=%.2f\n", format(length(md), big.mark=","), mean(md), sd(md)))

# ---- timing axis ----------------------------------------------------------
SWING <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip","hit_into_play")
sw <- dt[description %in% SWING & is.finite(depth) & is.finite(plate_x) &
         is.finite(plate_z) & is.finite(release_speed)]
FB<-c("FF","SI","FC"); BR<-c("SL","ST","CU","KC","SV","CS"); OS<-c("CH","FS","FO")
sw[, pgrp := fifelse(pitch_type %in% FB,"FB", fifelse(pitch_type %in% BR,"BR",
             fifelse(pitch_type %in% OS,"OS",NA_character_)))]
sw <- sw[!is.na(pgrp)]
sw[, px_bat := fifelse(stand=="R", -plate_x, plate_x)]
sw[, pz_rel := (plate_z - sz_bot)/pmax(sz_top - sz_bot, 0.1)]
fit <- lm(depth ~ ns(px_bat,5)*pgrp + ns(pz_rel,5) + ns(release_speed,4) + stand, data=sw)
sw[, r1 := residuals(fit)]
sw[, nb := .N, by=batter]; sw <- sw[nb >= 200]
sw[, tdev := r1 - mean(r1), by=batter]
sw[, whiff := description %in% c("swinging_strike","swinging_strike_blocked","foul_tip")]
sw[, bipf  := description == "hit_into_play" & bb_type != "" &
              is.finite(estimated_woba_using_speedangle)]
sw[, xw := estimated_woba_using_speedangle]

# =============================================================================
# 2. The shared timing axis: does it drive both quantities the same way?
# =============================================================================
cat("\n############ 2. Both quantities against the same timing axis ############\n")
BR6 <- c(-Inf,-15,-12,-9,-6,-3,0,3,6,9,12,15,Inf)
LB <- c("<=-15","-15..-12","-12..-9","-9..-6","-6..-3","-3..0",
        "0..+3","+3..+6","+6..+9","+9..+12","+12..+15",">=+15")
sw[, bin := cut(tdev, BR6, labels=LB, right=FALSE)]
t2 <- sw[, .(
  whiffs=sum(whiff & is.finite(miss_distance)),
  miss_dist=round(mean(miss_distance[whiff], na.rm=TRUE),2),
  bip=sum(bipf), xwobacon=round(mean(xw[bipf], na.rm=TRUE),3)
), by=bin][order(bin)]
print(t2)

wh <- sw[whiff & is.finite(miss_distance)]
bp <- sw[bipf == TRUE]
cat(sprintf("\n  Within whiffs (n=%s):  cor(tdev, miss_distance)      = %+.3f\n",
  format(nrow(wh), big.mark=","), cor(wh$tdev, wh$miss_distance)))
cat(sprintf("                          cor(|tdev|, miss_distance)    = %+.3f\n",
  cor(abs(wh$tdev), wh$miss_distance)))
cat(sprintf("  Within balls in play (n=%s): cor(tdev, xwOBAcon)   = %+.3f\n",
  format(nrow(bp), big.mark=","), cor(bp$tdev, bp$xw)))
cat(sprintf("                          cor(|tdev|, xwOBAcon)        = %+.3f\n",
  cor(abs(bp$tdev), bp$xw)))

# =============================================================================
# 3. Aggregate correlation between miss distance and xwOBAcon allowed
# =============================================================================
cat("\n############ 3. Aggregate: mean miss distance vs xwOBAcon allowed ############\n")
mk <- function(d, by, minw, minb) {
  a <- d[, .(
    whiffs=sum(whiff & is.finite(miss_distance)),
    miss=mean(miss_distance[whiff], na.rm=TRUE),
    se_miss=sd(miss_distance[whiff], na.rm=TRUE)/sqrt(pmax(sum(whiff & is.finite(miss_distance)),1)),
    bip=sum(bipf),
    xwobacon=mean(xw[bipf], na.rm=TRUE),
    se_xw=sd(xw[bipf], na.rm=TRUE)/sqrt(pmax(sum(bipf),1)),
    whiffpct=100*mean(whiff), tdev=mean(tdev),
    late6=100*mean(tdev <= -6), early12=100*mean(tdev >= 12)
  ), by=by]
  a[whiffs >= minw & bip >= minb]
}
s26 <- sw[game_year == 2026]
rep_pair <- function(d, lab) {
  ct <- cor.test(d$miss, d$xwobacon)
  cat(sprintf("  %-34s n=%4d | r = %+.3f (95%% CI %+.3f to %+.3f, p = %.3g)\n",
    lab, nrow(d), ct$estimate, ct$conf.int[1], ct$conf.int[2], ct$p.value))
  invisible(ct$estimate)
}
P  <- mk(s26, "pitcher", 60, 60)
PT <- mk(s26, c("pitcher","pitch_type"), 40, 40)
B  <- mk(s26, "batter", 60, 60)
r_p  <- rep_pair(P,  "per pitcher (2026)")
r_pt <- rep_pair(PT, "per pitcher x pitch type (2026)")
r_b  <- rep_pair(B,  "per batter (2026)")

cat("\n  Is the pitcher x pitch-type result just pitch-type composition? Centering both\n")
cat("  variables within pitch type isolates the pitcher-specific component:\n")
PTc <- copy(PT)[, `:=`(miss_c = miss - mean(miss), xw_c = xwobacon - mean(xwobacon)), by=pitch_type]
ctp <- cor.test(PTc$miss_c, PTc$xw_c)
cat(sprintf("    within pitch type              n=%4d | r = %+.3f (95%% CI %+.3f to %+.3f, p = %.3g)\n",
  nrow(PTc), ctp$estimate, ctp$conf.int[1], ctp$conf.int[2], ctp$p.value))
cat("\n  And the pitch-type averages themselves (the between-pitch-type pattern):\n")
print(PT[, .(units=.N, miss=round(mean(miss),2), xwobacon=round(mean(xwobacon),3),
             whiffpct=round(mean(whiffpct),1)), by=pitch_type][order(-miss)])

cat("\n  Same, controlling for how often the pitch is put in play at all:\n")
for (nm in list(list(P,"pitcher"), list(PT,"pitcher x pitch type"), list(B,"batter"))) {
  d <- nm[[1]]
  f <- lm(scale(xwobacon) ~ scale(miss) + scale(whiffpct), data=d)
  co <- summary(f)$coefficients
  cat(sprintf("    %-22s miss beta %+.3f (p=%.2g) | whiff%% beta %+.3f (p=%.2g) | R2=%.3f\n",
    nm[[2]], co[2,1], co[2,4], co[3,1], co[3,4], summary(f)$r.squared))
}

# =============================================================================
# 4. Repeatability: split-half within 2026, and 2025 -> 2026
# =============================================================================
cat("\n############ 4. Repeatability of each timing metric ############\n")
cat("  Split-half correlations are on half-samples; the Spearman-Brown column projects\n")
cat("  them to full-season reliability. Year-over-year is the stricter skill test.\n")

cat("\n  CAUTION: tdev is demeaned per batter across both seasons by construction, which forces\n")
cat("  each batter's 2025 and 2026 means to offset. Batter panels therefore use the raw\n")
cat("  location/velo-adjusted residual (praw) instead, which retains between-batter level.\n")

sw[, half := sample(rep_len(1:2, .N)), by=.(pitcher, pitch_type)]
metr <- function(d, by) d[, .(
  n=.N,
  whiffs=sum(whiff & is.finite(miss_distance)),
  bip=sum(bipf),
  miss=mean(miss_distance[whiff], na.rm=TRUE),
  tdev=mean(tdev),
  late6=100*mean(tdev <= -6),
  early12=100*mean(tdev >= 12),
  praw=mean(r1),
  late6r=100*mean(r1 <= -6),
  early12r=100*mean(r1 >= 12),
  whiffpct=100*mean(whiff),
  xwobacon=mean(xw[bipf], na.rm=TRUE)
), by=by]

MS_P <- c("miss","tdev","late6","early12","whiffpct","xwobacon")
MS_B <- c("miss","praw","late6r","early12r","whiffpct","xwobacon")
sb <- function(r) if (is.na(r)) NA else 2*r/(1+r)

# When the unit key contains pitch_type, raw cross-unit correlation is inflated by stable
# pitch-type differences ("curveballs miss big in both years"). ctr centers each metric
# within pitch type to isolate the pitcher-specific component.
ctr <- function(m, v) {
  a <- m[[paste0(v,"_a")]]; b <- m[[paste0(v,"_b")]]
  ok <- is.finite(a) & is.finite(b)
  d <- data.table(a=a, b=b, pt=m$pitch_type)[ok]
  d[, .(a = a - mean(a), b = b - mean(b)), by=pt]
}
panel <- function(A1, B1, by, lab, minn, minw, minb, ms=MS_P) {
  m <- merge(A1[n>=minn & whiffs>=minw & bip>=minb],
             B1[n>=minn & whiffs>=minw & bip>=minb], by=by, suffixes=c("_a","_b"))
  if (nrow(m) < 30) { cat(sprintf("\n  -- %s: too few units (%d) --\n", lab, nrow(m))); return(invisible()) }
  bypt <- "pitch_type" %in% by
  cat(sprintf("\n  -- %s (n=%d units) --\n", lab, nrow(m)))
  for (v in ms) {
    x <- m[[paste0(v,"_a")]]; y <- m[[paste0(v,"_b")]]
    ok <- is.finite(x) & is.finite(y)
    r <- if (sum(ok) >= 25) cor(x[ok], y[ok]) else NA_real_
    extra <- ""
    if (bypt && !is.na(r)) {
      cd <- ctr(m, v)
      extra <- sprintf(" | within-pitch-type %+.3f", cor(cd$a, cd$b))
    }
    cat(sprintf("     %-9s r = %+.3f   (S-B full-season %s)%s\n",
      v, r, if (is.na(r)) "  n/a" else sprintf("%+.3f", sb(r)), extra))
  }
}
cat("\n=== SPLIT-HALF WITHIN 2026 ===")
s <- sw[game_year==2026]
panel(metr(s[half==1], "pitcher"), metr(s[half==2], "pitcher"),
      "pitcher", "PITCHER, split-half 2026", 150, 25, 25)
panel(metr(s[half==1], c("pitcher","pitch_type")), metr(s[half==2], c("pitcher","pitch_type")),
      c("pitcher","pitch_type"), "PITCHER x PITCH TYPE, split-half 2026", 60, 15, 15)
s[, bhalf := sample(rep_len(1:2, .N)), by=batter]
panel(metr(s[bhalf==1], "batter"), metr(s[bhalf==2], "batter"),
      "batter", "BATTER, split-half 2026", 150, 25, 25, MS_B)

cat("\n=== YEAR OVER YEAR, 2025 -> 2026 (no Spearman-Brown needed) ===")
panel(metr(sw[game_year==2025], "pitcher"), metr(sw[game_year==2026], "pitcher"),
      "pitcher", "PITCHER, 2025 vs 2026", 300, 50, 50)
panel(metr(sw[game_year==2025], c("pitcher","pitch_type")),
      metr(sw[game_year==2026], c("pitcher","pitch_type")),
      c("pitcher","pitch_type"), "PITCHER x PITCH TYPE, 2025 vs 2026", 120, 30, 30)
panel(metr(sw[game_year==2025], "batter"), metr(sw[game_year==2026], "batter"),
      "batter", "BATTER, 2025 vs 2026", 300, 50, 50, MS_B)

# =============================================================================
# 5. Disattenuated version of the Section 3 correlation
# =============================================================================
cat("\n############ 5. Section 3 correlation corrected for measurement error ############\n")
relof <- function(A1, B1, by, v, minn, minw, minb) {
  m <- merge(A1[n>=minn & whiffs>=minw & bip>=minb], B1[n>=minn & whiffs>=minw & bip>=minb],
             by=by, suffixes=c("_a","_b"))
  x <- m[[paste0(v,"_a")]]; y <- m[[paste0(v,"_b")]]
  ok <- is.finite(x) & is.finite(y); sb(cor(x[ok], y[ok]))
}
for (cfg in list(
  list(by="pitcher", lab="pitcher", minn=150, minw=25, minb=25, r=r_p),
  list(by=c("pitcher","pitch_type"), lab="pitcher x pitch type", minn=60, minw=15, minb=15, r=r_pt),
  list(by="batter", lab="batter", minn=150, minw=25, minb=25, r=r_b))) {
  hcol <- if (identical(cfg$by, "batter")) "bhalf" else "half"
  A1 <- metr(s[get(hcol)==1], cfg$by); B1 <- metr(s[get(hcol)==2], cfg$by)
  rm_ <- relof(A1, B1, cfg$by, "miss", cfg$minn, cfg$minw, cfg$minb)
  rx_ <- relof(A1, B1, cfg$by, "xwobacon", cfg$minn, cfg$minw, cfg$minb)
  cat(sprintf("  %-22s observed r %+.3f | reliability: miss %.2f, xwOBAcon %.2f | disattenuated r %+.3f\n",
    cfg$lab, cfg$r, rm_, rx_, cfg$r/sqrt(rm_*rx_)))
}
cat("\n  Disattenuated r is the correlation the two would show if both were measured\n")
cat("  without error. It is an upper bound, and it inflates fast when reliability is low.\n")
