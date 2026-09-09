#!/usr/bin/env Rscript

# WHICH PITCHERS DOES EACH CUE ACTUALLY EXPLAIN?
#
# Fig 11 found two cue/group pairs that clear their own noise floor:
#   arm-angle similarity on OFFSPEED   (r = +0.014, p = 8.5e-07)
#   spin similarity on BREAKING BALLS  (r = +0.022, p = 1.2e-30)
# This pools them to the pitcher x pitch-type level and names the individuals carrying
# each correlation: the cue is in an extreme decile AND the residual moves the way the
# cue predicts.
#
# Two measurement notes that change how the cues have to be handled:
#
#  1. arm_diff is right-skewed (1.2 to 21 deg, median 3.1), so a z-score cannot reach
#     the matched-slot tail -- there is a floor at zero but no ceiling. Deciles instead.
#
#  2. spin_sim is degenerate for breaking balls: 69% are exactly zero and the 99th
#     percentile is 0.00006. Nonzero means the spin-axis gap vs the fastball is small
#     (median 0.1 deg vs 0.6 deg for the zeros), and there is NO gradient within the
#     nonzero group (r = -0.005 to +0.001 for the main types). So for breaking balls it
#     is a threshold, not a scale, and the pitcher-level cue is the SHARE of that
#     pitcher's breaking balls clearing it. The effect is real -- the zero/nonzero
#     residual gap has the same sign in all six breaking types -- but it is binary.
#
# Deciles are computed WITHIN pitch type, because the nonzero rate is mostly a property
# of the pitch (CU 71%, KC 63%, ST 34%, SL 11%); pooling would just rank curveballs
# above sliders.
#
# Residual = out-of-fold whiff residual with plate location and approach angle stripped
# out (same construction as Fig 11), in percentage points of whiff rate above the
# shape+location expectation.

suppressPackageStartupMessages({ library(data.table) })
AST <- "data/statcast_model/article_assets"
MIN_N    <- 200    # pitches of that type, pooled across seasons
MIN_TYPE <- 15     # qualifying pitchers a pitch type needs to get its own deciles

oof <- readRDS("data/statcast_model/oof_whiff_resid.rds")
oof[, res_loc := NA_real_]
for (g in c("breaking","offspeed")) {
  idx <- which(oof$grp == g); sub <- oof[idx]
  fit <- lm(wres ~ poly(plate_x,3)*poly(plate_z,3) + below_zone + VAA + HAA, data = sub)
  oof$res_loc[idx] <- residuals(fit)
}
oof[, spin_nz := as.numeric(spin_sim > 1e-9)]

aa <- fread("data/statcast_model/arm_angle_tunnel.csv")[, .(pitcher, season, pitch_type,
                                                            arm_diff, mean_arm, fb_arm)]
d <- merge(oof, aa, by = c("pitcher","season","pitch_type"), all.x = TRUE)
d <- d[grp %in% c("breaking","offspeed") & is.finite(res_loc)]

agg <- d[, .(player_name = player_name[1], n = .N,
             resid_pp    = 100*mean(res_loc),
             whiff_pct   = 100*mean(is_whiff),
             arm_diff    = mean(arm_diff, na.rm = TRUE),
             own_arm     = mean(mean_arm, na.rm = TRUE),
             fb_arm      = mean(fb_arm,   na.rm = TRUE),
             as_gap      = mean(as_gap,   na.rm = TRUE),
             spin_nz_pct = 100*mean(spin_nz, na.rm = TRUE)),
         by = .(pitcher, grp, pitch_type)][n >= MIN_N]

# Rank within pitch type where there is enough of it, else within the group.
rank_cue <- function(D, cue) {
  D <- copy(D)[is.finite(get(cue))]
  D[, cue_val := as.numeric(get(cue))]
  D[, big := .N >= MIN_TYPE, by = pitch_type]
  D[big == TRUE,  cue_pct := frank(cue_val, ties.method = "average")/.N, by = pitch_type]
  D[big == FALSE, cue_pct := frank(cue_val, ties.method = "average")/.N]
  D[]
}

# Fit the league relationship on the rank, then read each pitcher's fitted value as the
# cue's contribution to their residual. "Explained" = contribution and residual agree in
# sign; share = how much of the residual that contribution covers.
explain <- function(D, label, flip = FALSE) {
  D <- copy(D)
  if (flip) D[, cue_pct := 1 - cue_pct]      # orient so higher = more like the fastball
  m <- lm(resid_pp ~ cue_pct, data = D, weights = n)
  D[, contrib_pp := coef(m)[2] * (cue_pct - weighted.mean(cue_pct, n))]
  D[, aligned := sign(contrib_pp) == sign(resid_pp)]
  D[, share := fifelse(aligned, pmin(1, abs(contrib_pp/resid_pp)), NA_real_)]
  ct <- cor.test(D$cue_pct, D$resid_pp, method = "spearman", exact = FALSE)
  cat(sprintf("\n--- %s ---\n  %d pitcher-pitch types | across the full similarity range the fit moves the residual %+.2f pp | Spearman r = %+.3f, p = %.2g\n",
              label, nrow(D), coef(m)[2], unname(ct$estimate), ct$p.value))
  D[]
}
show <- function(D, cols, ttl, k = 12) {
  cat("\n", ttl, "\n", sep = "")
  if (!nrow(D)) { cat("  (none)\n"); return(invisible()) }
  print(head(D[, ..cols], k))
}

################################################################################
## 1. OFFSPEED, arm-angle similarity
################################################################################
cat("==================== OFFSPEED x ARM-ANGLE SIMILARITY ====================\n")
O <- explain(rank_cue(agg[grp == "offspeed"], "arm_diff"),
             "offspeed: smaller arm-angle gap vs the fastball", flip = TRUE)
O[, `:=`(sim_pct = round(100*cue_pct), arm_gap = round(arm_diff,1),
         slot = round(own_arm), fb_slot = round(fb_arm))]
OC <- c("player_name","pitch_type","n","arm_gap","slot","fb_slot","whiff_pct",
        "resid_pp","contrib_pp","share","sim_pct")
show(O[aligned == TRUE & cue_pct >= 0.90][order(-contrib_pp)], OC,
     "TOP DECILE OF SLOT MATCH, OVERPERFORMS (the cue's best cases):")
show(O[aligned == TRUE & cue_pct <= 0.10][order(contrib_pp)], OC,
     "BOTTOM DECILE (slot gives it away), UNDERPERFORMS:")
show(O[(cue_pct >= 0.90 | cue_pct <= 0.10) & aligned == FALSE][order(-abs(resid_pp))], OC,
     "DECILE OUTLIERS THE CUE GETS WRONG (counterexamples):", 8)

################################################################################
## 2. BREAKING BALLS, spin similarity (threshold rate)
################################################################################
cat("\n\n============ BREAKING x SPIN SIMILARITY (share clearing the threshold) ============\n")
B <- explain(rank_cue(agg[grp == "breaking"], "spin_nz_pct"),
             "breaking: share of the pitch's spin axis registering as fastball-like")
B[, `:=`(sim_pct = round(100*cue_pct), nz_rate = round(spin_nz_pct,1),
         axis_gap = round(as_gap,2))]
BC <- c("player_name","pitch_type","n","nz_rate","axis_gap","whiff_pct",
        "resid_pp","contrib_pp","share","sim_pct")
show(B[aligned == TRUE & cue_pct >= 0.90][order(-contrib_pp)], BC,
     "TOP DECILE OF SPIN MATCH FOR THEIR PITCH TYPE, OVERPERFORMS:")
show(B[aligned == TRUE & cue_pct <= 0.10][order(contrib_pp)], BC,
     "BOTTOM DECILE (spin gives it away), UNDERPERFORMS:")
show(B[(cue_pct >= 0.90 | cue_pct <= 0.10) & aligned == FALSE][order(-abs(resid_pp))], BC,
     "DECILE OUTLIERS THE CUE GETS WRONG (counterexamples):", 8)

################################################################################
## how often does each cue actually call it right?
################################################################################
cat("\n\n==================== COVERAGE ====================\n")
cvg <- function(D, lab) D[cue_pct >= .90 | cue_pct <= .10,
  .(cue = lab, decile_outliers = .N, called_right = sum(aligned),
    pct_right = round(100*mean(aligned),1),
    median_share_of_resid = round(100*median(share, na.rm = TRUE),1))]
print(rbindlist(list(cvg(O, "arm-angle similarity (offspeed)"),
                     cvg(B, "spin similarity (breaking)"))))
cat("\n50% right would be a coin flip. median_share = median fraction of an aligned\n",
    "pitcher's residual that the cue's fitted contribution covers.\n", sep = "")

fwrite(O[order(-contrib_pp)], file.path(AST, "ext_armangle_offspeed_pitchers.csv"))
fwrite(B[order(-contrib_pp)], file.path(AST, "ext_spinsim_breaking_pitchers.csv"))
cat("\nwrote ext_armangle_offspeed_pitchers.csv and ext_spinsim_breaking_pitchers.csv\n")
