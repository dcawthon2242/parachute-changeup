#!/usr/bin/env Rscript

# Full detail on the eight pitcher-seasons in the tightest defensible bin:
#   four-seam active spin > .904, changeup active spin > .900 (both league means)
#   spin-axis gap to the four-seamer <= 10 degrees
#   arm slot >= 44 degrees
#
# Reported with the raw rates alongside the residual, because a residual on its own hides
# whether a pitcher is beating the model from a high base or a low one - and those are very
# different pitches even at the same overperformance.

suppressPackageStartupMessages({ library(data.table) })
options(width = 215)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")

CH <- readRDS(file.path(MDIR, "parachute_ff.rds"))
R  <- readRDS(file.path(MDIR, "whiff_tjstuff.rds"))
AS <- readRDS(file.path(MDIR, "active_spin_long.rds"))

# Residual side: swings only, which is what the whiff model was fit on.
S <- R[, .(nsw = .N, axis = mean(axis_diff), arm = mean(arm_angle), velo_sep = -mean(speed_diff),
           w_res = 100*mean(r4), w_res_noarm = 100*mean(r3)), by = .(pitcher, player_name, season)]
# Raw side: all changeups, so usage and run value are on the full pitch count.
RAW <- CH[, .(np = .N, whiff = 100*sum(whiff, na.rm=TRUE)/sum(is_swing, na.rm=TRUE),
              swing = 100*mean(is_swing, na.rm=TRUE), velo = mean(release_speed, na.rm=TRUE),
              spin = mean(release_spin_rate, na.rm=TRUE),
              rv100 = 100*mean(rv, na.rm=TRUE)), by = .(pitcher, season)]
S <- merge(S, RAW, by = c("pitcher","season"))
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, as_ch = active_spin)], by = c("pitcher","season"))
S <- merge(S, AS[pitch_type == "FF", .(pitcher, season, as_fb = active_spin)], by = c("pitcher","season"))
S <- S[nsw >= 60]
S[, name := trimws(paste(sub(".*,\\s*","",player_name), sub(",.*","",player_name)))]

mu <- S[, .(as_fb = mean(as_fb), as_ch = mean(as_ch))]

# Two bins are of interest. "mean" is the a-priori one built from league averages plus a round
# 44-degree slot; "searched" is the optimum an exhaustive grid search landed on, which held-out
# testing showed to be largely overfit - it is listed for inspection, not as a recommendation.
WHICH <- { a <- commandArgs(TRUE); if (length(a)) a[1] else "mean" }
stopifnot(WHICH %in% c("mean","searched","core","wide","eff90","eff90axis"))
S[, as_gap := round(as_ch - as_fb, 4)]
sel <- switch(WHICH,
  mean     = quote(as_fb > mu$as_fb & as_ch > mu$as_ch & axis <= 10 & arm >= 44),
  searched = quote(as_fb >= .85 & as_ch >= .89 & axis <= 8 & arm >= 35),
  core     = quote(axis <= 10 & abs(as_gap) <= .10 & arm >= 44),
  wide     = quote(axis <= 15 & abs(as_gap) <= .15 & arm >= 42),
  # The .90 floor the D1 screen was tuned to hit, run on measured active spin so the MLB roster
  # shows what that screen selects. Note the difference from "core": that bin caps the DISTANCE
  # between the two efficiencies without requiring either to be high, so it admits pairs like
  # .83/.75. These two require both pitches to be genuinely efficient.
  eff90     = quote(as_fb >= .90 & as_ch >= .90 & arm >= 44),
  eff90axis = quote(as_fb >= .90 & as_ch >= .90 & arm >= 44 & axis <= 10))
RULES <- c(mean     = sprintf("FF active > %.3f, CH active > %.3f, axis gap <= 10, arm slot >= 44", mu$as_fb, mu$as_ch),
           searched = "FF active >= .85, CH active >= .89, axis gap <= 8, arm slot >= 35",
           core     = "axis gap <= 10, |active-spin gap| <= .10, arm slot >= 44",
           wide     = "axis gap <= 15, |active-spin gap| <= .15, arm slot >= 42",
           eff90     = "FF active >= .90, CH active >= .90, arm slot >= 44 (no axis condition)",
           eff90axis = "FF active >= .90, CH active >= .90, arm slot >= 44, axis gap <= 10")
B <- S[eval(sel)][order(-w_res)]
B[, exp_whiff := whiff - w_res]

cat(sprintf("bin (%s): %s, and >= 60 changeup swings\n", WHICH, RULES[WHICH]))
cat(sprintf("%d pitcher-seasons, %d unique pitchers, %s changeups\n\n", nrow(B), uniqueN(B$pitcher),
            format(sum(B$np), big.mark=",")))

print(B[, .(Pitcher = name, Season = season, CH = np, Swings = nsw,
            Arm = round(arm,1), Axis = round(axis,1), ActFF = round(as_fb,2), ActCH = round(as_ch,2),
            EffGap = round(as_gap,2), Velo = round(velo,1), VeloSep = round(velo_sep,1),
            Spin = round(spin), Whiff = round(whiff,1), Expected = round(exp_whiff,1),
            Above = round(w_res,1), NoArm = round(w_res_noarm,1),
            RV100 = round(rv100,2))], row.names = FALSE)

# Sign convention, verified rather than assumed: a swinging strike carries a mean rv of +0.112
# in this dataset and a ball in play -0.038, so POSITIVE run value favours the PITCHER. This is
# the opposite of Statcast's raw delta_run_exp and of the convention in Nestico's write-up, and
# getting it backwards would invert the conclusion below.
L <- S[!eval(sel)]
cat(sprintf("\nbin mean whiff %.1f%% vs league %.1f%%   |   bin RV/100 %+.2f vs league %+.2f (higher favours pitcher)\n",
            mean(B$whiff), mean(L$whiff), mean(B$rv100), mean(L$rv100)))
cat(sprintf("bin whiff above model %+.2f pp (arm-aware) / %+.2f pp (arm-blind)\n",
            mean(B$w_res), mean(B$w_res_noarm)))
cat(sprintf("positive residual in %d of %d seasons; median %+.2f\n",
            sum(B$w_res > 0), nrow(B), median(B$w_res)))

# The two obvious names carry the average, so it is worth seeing the bin without them.
D <- B[!name %in% c("Tarik Skubal","Dylan Cease")]
cat(sprintf("\ndropping Skubal and Cease: %d seasons, whiff above model %+.2f pp (was %+.2f)\n",
            nrow(D), mean(D$w_res), mean(B$w_res)))

f <- sprintf("ext_bin_roster_%s.csv", WHICH)
fwrite(B, file.path(AST, f)); cat(sprintf("\nwrote %s\n", f))
