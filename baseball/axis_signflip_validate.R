#!/usr/bin/env Rscript

# DOES THE PITCH-TYPE SIGN FLIP SURVIVE OUT OF SAMPLE?
#
# Pooled 2023H2-2026 said: sliders and sweepers overperform when their spin axis sits
# CLOSE to the fastball's, curveballs and knuckle-curves when it sits closer to a
# perfect MIRROR. That is a nice story, which is exactly why it needs a holdout.
#
# Three tests, increasingly strict:
#   A. replication  - refit the correlation inside 2026 alone and see if the sign holds
#   B. predictive   - measure the cue on 2023-2025, correlate it against the SAME
#                     pitcher's 2026 residual. This is the one that matters: it asks
#                     whether the cue tells you anything you did not already know.
#   C. bootstrap    - resample pitchers to get an interval on the pooled sign per type
#
# The cue is oriented as SIMILARITY throughout: axis_sim = -axis_diff, so positive r
# means "an axis closer to the fastball's goes with overperformance" and negative r
# means "closer to a mirror goes with overperformance".

suppressPackageStartupMessages({ library(data.table) })
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")
set.seed(42)
TYPES <- c("SL","ST","CU","KC")
LAB <- c(SL="Slider", ST="Sweeper", CU="Curveball", KC="Knuckle-curve")
MIN_TR <- 200   # pitches to qualify in the 2023-2025 training window
MIN_TE <- 120   # pitches to qualify in 2026 alone

oof <- readRDS(file.path(MDIR, "oof_whiff_resid.rds"))
oof[, res_loc := NA_real_]
for (g in c("breaking","offspeed")) {
  i <- which(oof$grp == g); s <- oof[i]
  oof$res_loc[i] <- residuals(lm(wres ~ poly(plate_x,3)*poly(plate_z,3)+below_zone+VAA+HAA, data=s))
}
d <- readRDS(file.path(MDIR, "miss_grade_data_activespin.rds"))
key <- c("game_pk","at_bat_number","pitch_number","season")
setkeyv(d, key); setkeyv(oof, key)
b <- d[, .(game_pk, at_bat_number, pitch_number, season, axis_diff)][oof]
b <- b[grp == "breaking" & pitch_type %in% TYPES & is.finite(res_loc) & is.finite(axis_diff)]
b[, axis_sim := -axis_diff]

agg <- function(D, minn) D[, .(n = .N, axis_sim = mean(axis_sim), resid_pp = 100*mean(res_loc)),
                           by = .(pitcher, pitch_type)][n >= minn]
sp <- function(x, y) { ct <- cor.test(x, y, method="spearman", exact=FALSE)
                       list(r = unname(ct$estimate), p = ct$p.value) }

TR <- agg(b[season <= 2025], MIN_TR)
TE <- agg(b[season == 2026], MIN_TE)

################################################################################
cat("========== A. REPLICATION: same-season fit, 2023-25 vs 2026 alone ==========\n")
A <- rbindlist(lapply(TYPES, function(ty) {
  a <- TR[pitch_type == ty]; c26 <- TE[pitch_type == ty]
  s1 <- sp(a$axis_sim, a$resid_pp); s2 <- sp(c26$axis_sim, c26$resid_pp)
  data.table(pitch = LAB[ty], n_2325 = nrow(a), r_2325 = round(s1$r,3), p_2325 = signif(s1$p,2),
             n_2026 = nrow(c26), r_2026 = round(s2$r,3), p_2026 = signif(s2$p,2),
             same_sign = sign(s1$r) == sign(s2$r))
}))
print(A)

################################################################################
cat("\n========== B. PREDICTIVE: cue from 2023-25 -> residual in 2026 ==========\n")
J <- merge(TR[, .(pitcher, pitch_type, axis_sim_tr = axis_sim, n_tr = n)],
           TE[, .(pitcher, pitch_type, resid26 = resid_pp, n_te = n)],
           by = c("pitcher","pitch_type"))
B2 <- rbindlist(lapply(TYPES, function(ty) {
  j <- J[pitch_type == ty]
  if (nrow(j) < 12) return(data.table(pitch = LAB[ty], n_pitchers = nrow(j),
                                      r = NA_real_, p = NA_real_, same_sign = NA))
  s <- sp(j$axis_sim_tr, j$resid26)
  data.table(pitch = LAB[ty], n_pitchers = nrow(j), r = round(s$r,3), p = signif(s$p,2),
             same_sign = sign(s$r) == sign(A[pitch == LAB[ty]]$r_2325))
}))
print(B2)
cat("\n  (a pitcher must clear 200 pitches in 2023-25 AND 120 in 2026 to appear here)\n")

################################################################################
cat("\n========== C. BOOTSTRAP over pitchers, pooled years ==========\n")
ALL <- agg(b, MIN_TR)
C <- rbindlist(lapply(TYPES, function(ty) {
  a <- ALL[pitch_type == ty]
  rs <- replicate(2000, { i <- sample.int(nrow(a), replace = TRUE)
                          suppressWarnings(cor(a$axis_sim[i], a$resid_pp[i], method="spearman")) })
  data.table(pitch = LAB[ty], n_pitchers = nrow(a),
             r = round(cor(a$axis_sim, a$resid_pp, method="spearman"),3),
             lo = round(quantile(rs,.025),3), hi = round(quantile(rs,.975),3),
             pct_positive = round(100*mean(rs > 0)))
}))
print(C)
cat("\n  pct_positive = share of bootstrap resamples where an axis closer to the fastball\n",
    "  looks better. Near 100 = reliably 'match'; near 0 = reliably 'mirror'; near 50 = noise.\n", sep="")

fwrite(A, file.path(AST, "ext_axis_signflip_replication.csv"))
fwrite(B2, file.path(AST, "ext_axis_signflip_predictive.csv"))
fwrite(C, file.path(AST, "ext_axis_signflip_bootstrap.csv"))
cat("\nwrote ext_axis_signflip_{replication,predictive,bootstrap}.csv\n")
