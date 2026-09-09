#!/usr/bin/env Rscript

# ALL THREE CUES ARE PITCHER TRAITS, SO ONLY THE BETWEEN-PITCHER CORRELATION IS REAL
#
# Checking why the offspeed spin bar could be significant pooled while null both between
# and within pitchers turned up something that invalidates the whole per-pitch framing:
#
#   path_ratio  varies in 0 of 4,186 pitcher x pitch-type x season groups
#   arm_diff    same - it is a pitcher-season slot trait by construction
#   spin_sim    axis_diff moves per pitch, but the kernel is degenerate on breaking balls
#               (between-pitcher SD 0.0002), so it carries no usable within variation either
#
# path_ratio in particular is stored as a pitcher x pitch-type x season average and
# broadcast onto every pitch. The original Figure 11 described it as measured "at its
# NATIVE per-pitch unit", which is not what the data is. Correlating 296,774 copies of
# ~3,500 distinct numbers against per-pitch residuals inflates significance enormously
# while measuring nothing more than the between-pitcher relationship, badly weighted.
#
# So this reports the between-pitcher correlation, which is what every claim in the
# article actually means, at the sample size that claim is entitled to.

suppressPackageStartupMessages({ library(data.table) })
options(width = 205)
MDIR <- "data/statcast_model"

d <- rbindlist(lapply(c("breaking","offspeed"),
  function(g) readRDS(file.path(MDIR, sprintf("resid_old_vs_new_%s.rds", g)))))
aa <- fread(file.path(MDIR, "arm_angle_tunnel.csv"))[, .(pitcher, season, pitch_type, arm_diff)]
d <- merge(d, aa, by = c("pitcher","season","pitch_type"), all.x = TRUE)
d[, group := fifelse(grp == "breaking", "Breaking", "Offspeed")]
d[, `:=`(`Arm Angle` = -arm_diff, `Trajectory` = -path_ratio, `Spin Similarity` = spin_sim)]

CUES <- c("Arm Angle","Trajectory","Spin Similarity")
RES  <- c(res_shape = "shape only", res_old = "shape + post-hoc loc (published)",
          res_new = "location in model")
MINP <- 200

cat("=== AGGREGATION LEVEL OF EACH CUE ===\n")
print(rbindlist(lapply(CUES, function(cu) {
  x <- d[is.finite(get(cu))]
  data.table(cue = cu, distinct_values = uniqueN(x[[cu]]),
             pitcher_type_season_groups = uniqueN(x[, .(pitcher, pitch_type, season)]),
             groups_where_it_varies = nrow(x[, .(u = uniqueN(get(cu))),
                                            by = .(pitcher, pitch_type, season)][u > 1]))
})), row.names = FALSE)

out <- rbindlist(lapply(c("Breaking","Offspeed"), function(g)
  rbindlist(lapply(CUES, function(cu)
    rbindlist(lapply(names(RES), function(rv) {
      B <- d[group == g, .(n = .N, cx = mean(get(cu), na.rm = TRUE),
                           ry = 100*mean(get(rv), na.rm = TRUE)),
             by = .(pitcher, pitch_type)][n >= MINP & is.finite(cx) & is.finite(ry)]
      ct <- suppressWarnings(cor.test(B$cx, B$ry, method = "spearman", exact = FALSE))
      data.table(group = g, cue = cu, residual = RES[[rv]], pitchers = nrow(B),
                 r = round(unname(ct$estimate), 3), p = signif(ct$p.value, 2),
                 verdict = fifelse(ct$p.value < .05, "significant", "null"))
    }))))))

cat("\n=== BETWEEN-PITCHER SPEARMAN r WITH OVERPERFORMANCE ===\n")
cat("    one point per pitcher x pitch type, min", MINP, "pitches\n")
cat("    all cues oriented so higher = more like the primary fastball\n\n")
for (g in c("Breaking","Offspeed")) {
  cat("############ ", g, "\n"); print(out[group == g][, !"group"], row.names = FALSE); cat("\n")
}
fwrite(out, file.path(MDIR, "article_assets", "ext_cue_decomposition.csv"))

cat("=== SUMMARY ===\n")
s <- out[residual == "location in model"]
cat("  with location modeled, cues that clear p<.05 between pitchers: ",
    if (nrow(s[verdict == "significant"]) == 0) "NONE"
    else paste(s[verdict == "significant", paste0(cue, " (", group, ", r=", r, ")")],
               collapse = "; "), "\n", sep = "")
