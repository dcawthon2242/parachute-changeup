#!/usr/bin/env Rscript

# HOW WELL DOES THE LEAGUE-WIDE NULL HOLD?
#
# Correlation coefficients are not comparable across levels of aggregation here. At the pitch
# level the outcome is a 0/1 whiff, so almost all of its variance is binomial noise unrelated
# to anything, which crushes r toward zero no matter how real the underlying effect is. At the
# season level that noise averages out and the same effect reads as a larger r. The quantity
# that IS comparable is the slope in percentage points of residual per mph of separation, so
# everything below is reported that way.
#
# Five ways of asking whether the null holds:
#   levels     per-pitch vs per-season slope. They should agree if the null is real.
#   seasons    seven independent years. A stable effect does not flip sign annually.
#   strata     by changeup count. Rising slope with sample size would signal attenuation.
#   within     a pitcher's own year-over-year change. Kills all between-pitcher confounding.
#   shape      quadratic and high-separation tail, in case the effect is not linear.
#
# Then the reverse question: given this null, how surprising is the +1.25 slope in the bin?

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
options(width = 200)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")

L  <- readRDS(file.path(MDIR, "parachute_extended_resid.rds"))
CH <- L$CH; SW <- L$SW
stopifnot(identical(CH[is_swing == TRUE]$pitcher, SW$pitcher))
P <- cbind(CH[is_swing == TRUE, .(pitcher, season, speed_diff)], SW[, .(wh_res)])
P[, velo_sep := -speed_diff]

S <- fread(file.path(AST, "ext_parachute_extended.csv"))
S <- S[is.finite(as_gap) & is.finite(arm) & is.finite(velo_sep) & is.finite(wh)]
S[, core := axis <= 10 & abs(as_gap) <= .10 & arm >= 44]

# Slope of residual (in percentage points) on velocity separation, with a 95% interval.
slope <- function(y, x, w = NULL, lab, scale = 1) {
  ok <- is.finite(x) & is.finite(y); y <- y[ok]*scale; x <- x[ok]; w <- if (is.null(w)) NULL else w[ok]
  m <- if (is.null(w)) lm(y ~ x) else lm(y ~ x, weights = w)
  ci <- confint(m)["x", ]; cf <- summary(m)$coefficients["x", ]
  cat(sprintf("  %-44s n=%7d  %+.3f pp/mph  [%+.3f,%+.3f]  p=%.3f\n",
              lab, length(x), cf[1], ci[1], ci[2], cf[4]))
  invisible(cf[1])
}

cat("=== 1. the same effect measured at two levels of aggregation ===\n")
slope(P$wh_res, P$velo_sep, NULL, "per pitch, every changeup swing", 100)
slope(S$wh, S$velo_sep, S$nsw, "per pitcher-season, weighted by swings")
slope(S$wh, S$velo_sep, NULL, "per pitcher-season, unweighted")

cat("\n=== 2. one estimate per season ===\n")
YS <- rbindlist(lapply(sort(unique(S$season)), function(y) {
  D <- S[season == y]; m <- lm(wh ~ velo_sep, D, weights = nsw); ci <- confint(m)["velo_sep",]
  data.table(season = y, n = nrow(D), b = coef(m)[2], lo = ci[1], hi = ci[2],
             p = summary(m)$coefficients["velo_sep",4]) }))
print(YS[, .(season, n, slope = round(b,3), lo = round(lo,3), hi = round(hi,3), p = round(p,3))],
      row.names = FALSE)
cat(sprintf("  %d of 7 seasons positive; %d reach p<.05; sign flips %d times\n",
            sum(YS$b > 0), sum(YS$p < .05), sum(diff(sign(YS$b)) != 0)))

cat("\n=== 3. by changeup volume (does the slope grow as noise falls?) ===\n")
for (m in c(60, 150, 250, 400, 600))
  slope(S[np >= m]$wh, S[np >= m]$velo_sep, S[np >= m]$nsw, sprintf("seasons with %d+ changeups", m))

cat("\n=== 4. within pitcher: his own year-over-year change ===\n")
# Consecutive seasons only, so the comparison is the same arm a year apart rather than two
# different pitchers. This removes every stable pitcher trait at once.
PR <- S[, .(pitcher, season = season + 1L, v0 = velo_sep, w0 = wh, n0 = nsw)]
W  <- merge(S[, .(pitcher, season, velo_sep, wh, nsw)], PR, by = c("pitcher","season"))
W[, `:=`(dv = velo_sep - v0, dw = wh - w0, wt = pmin(nsw, n0))]
slope(W$dw, W$dv, W$wt, "change in residual on change in separation")
cat(sprintf("  %d pitcher-season pairs, %d distinct pitchers, sd of the velo change is %.2f mph\n",
            nrow(W), uniqueN(W$pitcher), sd(W$dv)))
Wc <- W[pitcher %in% S[core == TRUE]$pitcher]
if (nrow(Wc) >= 15) slope(Wc$dw, Wc$dv, Wc$wt, "same, restricted to bin pitchers")

cat("\n=== 5. is it linear? ===\n")
m2 <- lm(wh ~ poly(velo_sep, 2), S, weights = nsw)
cat(sprintf("  quadratic term p=%.3f   full-model R2=%.5f\n",
            summary(m2)$coefficients[3,4], summary(m2)$r.squared))
for (q in c(.80, .90, .95)) {
  th <- quantile(S$velo_sep, q); D <- S[velo_sep >= th]
  cat(sprintf("  top %2.0f%% by separation (>= %.1f mph): mean residual %+.2f pp  (rest %+.2f)  p=%.3f\n",
              100*(1-q), th, mean(D$wh), mean(S[velo_sep < th]$wh),
              t.test(D$wh, S[velo_sep < th]$wh)$p.value))
}

## ---- how unusual is the bin, given this null? ------------------------------------------
cat("\n=== 6. null distribution: 26-season subsets drawn to match the bin ===\n")
C <- S[core == TRUE]; robs <- cor(C$velo_sep, C$wh); bobs <- coef(lm(wh ~ velo_sep, C))[2]
# Matching on changeup count matters: the bin skews to well-sampled seasons, whose residuals
# are less noisy and therefore correlate more readily with anything.
pool <- S[core == FALSE]
set.seed(11)
draw <- replicate(20000, {
  i <- sapply(C$np, function(k) { cand <- which(abs(pool$np - k) <= 0.25*k)
                                  if (!length(cand)) which.min(abs(pool$np - k)) else sample(cand, 1) })
  d <- pool[i]; c(cor(d$velo_sep, d$wh), coef(lm(wh ~ velo_sep, d))[2]) })
cat(sprintf("  observed  r=%+.3f  slope=%+.2f pp/mph\n", robs, bobs))
cat(sprintf("  matched null r: mean %+.3f, 95%% of draws below %+.3f, share >= observed %.2f%%\n",
            mean(draw[1,]), quantile(draw[1,], .95), 100*mean(draw[1,] >= robs)))
cat(sprintf("  matched null slope: mean %+.3f, 95%% below %+.3f, share >= observed %.2f%%\n",
            mean(draw[2,]), quantile(draw[2,], .95), 100*mean(draw[2,] >= bobs)))

## ---- figure -------------------------------------------------------------------------------
YS[, lab := as.character(season)]
A <- rbind(
  data.table(grp = "By season", lab = YS$lab, b = YS$b, lo = YS$lo, hi = YS$hi),
  data.table(grp = "By volume", lab = c("60+","150+","250+","400+","600+"),
             b = sapply(c(60,150,250,400,600), function(m) coef(lm(wh~velo_sep, S[np>=m], weights=nsw))[2]),
             lo = sapply(c(60,150,250,400,600), function(m) confint(lm(wh~velo_sep, S[np>=m], weights=nsw))["velo_sep",1]),
             hi = sapply(c(60,150,250,400,600), function(m) confint(lm(wh~velo_sep, S[np>=m], weights=nsw))["velo_sep",2])),
  data.table(grp = "Reference", lab = c("All seasons","Within pitcher","Core bin"),
             b  = c(coef(lm(wh~velo_sep,S,weights=nsw))[2], coef(lm(dw~dv,W,weights=wt))[2], bobs),
             lo = c(confint(lm(wh~velo_sep,S,weights=nsw))["velo_sep",1], confint(lm(dw~dv,W,weights=wt))["dv",1],
                    confint(lm(wh~velo_sep,C))["velo_sep",1]),
             hi = c(confint(lm(wh~velo_sep,S,weights=nsw))["velo_sep",2], confint(lm(dw~dv,W,weights=wt))["dv",2],
                    confint(lm(wh~velo_sep,C))["velo_sep",2])))
A[, grp := factor(grp, levels = c("By season","By volume","Reference"))]
A[, lab := factor(lab, levels = rev(lab))]
A[, hit := fifelse(lo > 0 | hi < 0, "clears zero", "consistent with zero")]

gg <- ggplot(A, aes(b, lab, colour = hit)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey45") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = .28, linewidth = .7) +
  geom_point(size = 2.5) +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = c("clears zero" = "#1d7870", "consistent with zero" = "grey55"), name = NULL) +
  labs(title = "Across all changeups, velocity separation buys almost nothing the model has not already priced in",
       subtitle = paste0("Slope of whiff residual on velocity separation, in percentage points per mph, with 95 percent intervals. Zero is the expected answer: velocity separation is one\n",
                         "of the 24 features the model trains on, so anything left over is the model failing. All seven seasons sit on zero and the sign flips three times. The estimate\n",
                         "does drift upward as the volume floor rises, clearing zero at the 250 and 400 changeup cuts, but the intervals overlap heavily and the drift most likely\n",
                         "reflects which pitchers survive each cut rather than a growing effect. The decisive row is 'within pitcher': when the same arm changes its own separation\n",
                         "from one year to the next it gains nothing, and that interval rules out anything close to the Core bin's slope. The bin remains the lone exception."),
       x = "Whiff residual per mph of additional separation (percentage points)", y = NULL,
       caption = "Source: Statcast 2020-2026 - 1,226 changeup pitcher-seasons - 230,716 swings - out-of-fold LightGBM residuals") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12.2), plot.subtitle = element_text(size = 8.2),
        strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold", size = 9),
        panel.grid.minor = element_blank(), legend.position = "top",
        plot.caption = element_text(size = 7.5, colour = "grey40"))
ggsave(file.path(AST, "fig30_velo_null_robustness.png"), gg, width = 10.5, height = 6.8, dpi = 150)
cat("\nwrote fig30_velo_null_robustness.png\n")
