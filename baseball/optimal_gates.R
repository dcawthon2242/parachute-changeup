#!/usr/bin/env Rscript

# OPTIMAL GATES ON THE FOUR VARIABLES - AND WHAT THAT OPTIMUM IS ACTUALLY WORTH.
#
# Searching a four-dimensional threshold grid for the subgroup with the largest whiff
# overperformance is guaranteed to return something impressive. With thousands of candidate
# bins and ~1,100 pitcher-seasons, the winner is the best draw from thousands of correlated
# noisy estimates, so its effect size is biased upward and its p-value is meaningless as
# reported. Two calibrations are therefore run alongside the search:
#
#   PERMUTATION NULL   shuffle the whiff residual across pitcher-seasons, destroying any real
#                      relationship, and run the identical search. Whatever it finds is what
#                      this procedure extracts from pure noise. The real optimum has to beat
#                      that distribution, not zero.
#   SPLIT-HALF         choose thresholds on one random half and score them on the untouched
#                      half. The gap between the two is the overfitting, measured directly, and
#                      the held-out number is the honest estimate of what the gate is worth.
#
# Residuals come from the arm-aware model (shape + location + arm angle), because the search
# includes an arm-slot threshold and scoring an arm gate against an arm-blind model would
# rediscover the missing-feature artifact documented earlier in this project.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
set.seed(7); options(width = 205)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")

R  <- readRDS(file.path(MDIR, "whiff_tjstuff.rds"))
AS <- readRDS(file.path(MDIR, "active_spin_long.rds"))
S <- R[, .(nsw = .N, axis = mean(axis_diff), arm = mean(arm_angle), velo_sep = -mean(speed_diff),
           w = 100*mean(r4)), by = .(pitcher, player_name, season)]
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, as_ch = active_spin)], by = c("pitcher","season"))
S <- merge(S, AS[pitch_type == "FF", .(pitcher, season, as_fb = active_spin)], by = c("pitcher","season"))
S <- S[nsw >= 60]
S[, name := trimws(paste(sub(".*,\\s*","",player_name), sub(",.*","",player_name)))]
cat(sprintf("population: %d pitcher-seasons\n", nrow(S)))

ARM <- seq(25, 55, 2.5); AXIS <- seq(6, 30, 2)
AFB <- seq(.85, .96, .02); ACH <- seq(.85, .96, .02)
GRID <- CJ(arm = ARM, axis = AXIS, as_fb = AFB, as_ch = ACH)
cat(sprintf("candidate gate combinations: %s   minimum bin size: 20\n\n", format(nrow(GRID), big.mark=",")))

# Returns the best combination by mean-difference, subject to a floor on bin size. Written to
# take the outcome as an argument so the permutation null reuses the exact same search.
search <- function(D, y, minn = 20L) {
  best <- list(v = -Inf); tot <- sum(y); n <- length(y)
  for (a in ARM) { i1 <- D$arm >= a
    for (x in AXIS) { i2 <- i1 & D$axis <= x
      if (sum(i2) < minn) next
      for (fb in AFB) { i3 <- i2 & D$as_fb >= fb
        if (sum(i3) < minn) next
        for (ch in ACH) { i <- i3 & D$as_ch >= ch; k <- sum(i)
          if (k < minn) next
          m <- sum(y[i])/k; v <- m - (tot - sum(y[i]))/(n - k)
          if (v > best$v) best <- list(v = v, arm = a, axis = x, as_fb = fb, as_ch = ch, n = k)
        } } } }
  best
}

## ---- 1. the optimum on the real data --------------------------------------------------------
B <- search(S, S$w)
S[, opt := arm >= B$arm & axis <= B$axis & as_fb >= B$as_fb & as_ch >= B$as_ch]
t0 <- t.test(S$w[S$opt], S$w[!S$opt])
cat("=== best gates found by exhaustive search ===\n")
cat(sprintf("  arm slot >= %.1f deg | axis gap <= %.0f deg | FF active spin >= %.2f | CH active spin >= %.2f\n",
            B$arm, B$axis, B$as_fb, B$as_ch))
cat(sprintf("  %d seasons, whiff above model %+.2f pp [%+.2f, %+.2f], nominal p = %.4f\n\n",
            B$n, B$v, t0$conf.int[1], t0$conf.int[2], t0$p.value))

## ---- 2. what the same search finds in pure noise ---------------------------------------------
NPERM <- 300
cat(sprintf("=== permutation null: same search on %d shuffles of the outcome ===\n", NPERM))
null <- vapply(seq_len(NPERM), function(i) search(S, sample(S$w))$v, numeric(1))
pval <- (1 + sum(null >= B$v)) / (NPERM + 1)
cat(sprintf("  best effect found in shuffled data: median %+.2f, 90th pct %+.2f, max %+.2f\n",
            median(null), quantile(null,.90), max(null)))
cat(sprintf("  real optimum %+.2f  ->  search-corrected p = %.4f\n\n", B$v, pval))

## ---- 3. how much survives out of sample -------------------------------------------------------
NSPLIT <- 300
cat(sprintf("=== split-half: gates chosen on one half, scored on the other (%d splits) ===\n", NSPLIT))
sp <- rbindlist(lapply(seq_len(NSPLIT), function(i) {
  h <- sample(nrow(S), floor(nrow(S)/2)); A <- S[h]; Bh <- S[-h]
  b <- search(A, A$w, minn = 10L)
  if (is.infinite(b$v)) return(NULL)
  j <- Bh$arm >= b$arm & Bh$axis <= b$axis & Bh$as_fb >= b$as_fb & Bh$as_ch >= b$as_ch
  if (sum(j) < 3 || sum(!j) < 3) return(NULL)
  data.table(train = b$v, test = mean(Bh$w[j]) - mean(Bh$w[!j]), n_test = sum(j),
             arm = b$arm, axis = b$axis, as_fb = b$as_fb, as_ch = b$as_ch) }))
cat(sprintf("  in-sample (chosen on this half):  %+.2f pp\n", mean(sp$train)))
cat(sprintf("  held out (scored on other half):  %+.2f pp  [%+.2f, %+.2f]\n",
            mean(sp$test), quantile(sp$test,.05), quantile(sp$test,.95)))
cat(sprintf("  shrinkage: %.0f%% of the in-sample effect does not replicate\n",
            100*(1 - mean(sp$test)/mean(sp$train))))
cat(sprintf("  held-out effect positive in %.0f%% of splits\n", 100*mean(sp$test > 0)))
cat("\n  thresholds the search picks across splits (median, and how much they wander):\n")
for (v in c("arm","axis","as_fb","as_ch"))
  cat(sprintf("    %-6s median %6.2f   5th-95th %.2f to %.2f\n", v, median(sp[[v]]),
              quantile(sp[[v]],.05), quantile(sp[[v]],.95)))

## ---- 4. a few fixed reference gates, for comparison --------------------------------------------
cat("\n=== reference gates, no searching involved ===\n")
refs <- list("all four above league average" = quote(arm > 37.9 & axis < 22.4 & as_fb > .904 & as_ch > .900),
             "spin floors + axis <= 10, no arm gate" = quote(as_fb > .904 & as_ch > .900 & axis <= 10),
             "spin floors + axis <= 10 + arm >= 44" = quote(as_fb > .904 & as_ch > .900 & axis <= 10 & arm >= 44),
             "searched optimum" = quote(opt))
print(rbindlist(lapply(names(refs), function(k) { i <- S[, eval(refs[[k]])]
  t <- t.test(S$w[i], S$w[!i])
  data.table(gate = k, seasons = sum(i), whiff = round(diff(rev(t$estimate)),2),
             ci = sprintf("[%+.2f, %+.2f]", t$conf.int[1], t$conf.int[2]),
             nominal_p = round(t$p.value,4)) })), row.names = FALSE)

cat("\n=== the searched-optimum roster ===\n")
print(S[opt == TRUE][order(-w), .(Pitcher = name, Season = season, Swings = nsw,
      Arm = round(arm,1), AxisGap = round(axis,1), ActFF = round(as_fb,2), ActCH = round(as_ch,2),
      VeloSep = round(velo_sep,1), WhiffAbove = round(w,1))], row.names = FALSE)

## ---- figure ---------------------------------------------------------------------------------
D <- rbind(data.table(what = "Best bin found in shuffled data\n(300 permutations)", v = null),
           data.table(what = "Best bin found on one half,\nscored on the other (300 splits)", v = sp$test))
gg <- ggplot(D, aes(v, fill = what)) +
  geom_histogram(bins = 40, alpha = .75, colour = NA, position = "identity") +
  geom_vline(xintercept = B$v, colour = "#b4632a", linewidth = .9) +
  annotate("text", x = B$v, y = Inf, vjust = 1.6, hjust = -0.05, size = 3, colour = "#b4632a",
           label = sprintf("searched optimum  %+.2f pp", B$v)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey45") +
  scale_fill_manual(values = c("#c2b8a3", "#1d7870"), name = NULL) +
  labs(title = sprintf("The best gates score %+.2f points, but an identical search finds %+.2f in data with the signal shuffled out",
                       B$v, median(null)),
       subtitle = paste0("Exhaustive search over ", format(nrow(GRID), big.mark=","), " threshold combinations on arm slot, spin-axis gap, four-seam active spin and changeup active spin,\n",
                         "requiring at least 20 pitcher-seasons per bin. The grey distribution is the same search run on 300 shuffles of the outcome: with no real relationship left\n",
                         "in the data it still returns a bin worth ", sprintf("%+.2f", median(null)), " points at the median, because the winner of thousands of correlated noisy comparisons is a biased\n",
                         "estimate by construction. That is the bar the orange line has to clear, not zero. The teal distribution is the honest one - thresholds picked on a random\n",
                         "half and scored on the untouched half - and it is centred near zero, which is the finding: the specific numbers the search lands on do not carry to\n",
                         "unseen pitchers. Read the optimum as the shape of a gradient, not as calibrated thresholds."),
       x = "Whiff above model, bin minus everyone else (percentage points)", y = "count",
       caption = "Source: Statcast 2020-2026 - four-seam anchor - residuals from a shape, location and arm-angle model - 1,086 pitcher-seasons with 60+ changeup swings") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11.5), plot.subtitle = element_text(size = 8.2),
        panel.grid.minor = element_blank(), legend.position = "top",
        plot.caption = element_text(size = 7.5, colour = "grey40"))
ggsave(file.path(AST, "fig37_optimal_gates.png"), gg, width = 11.5, height = 6.8, dpi = 150)
fwrite(sp, file.path(AST, "ext_optimal_gates_splits.csv"))
cat("\nwrote ext_optimal_gates_splits.csv, fig37_optimal_gates.png\n")
