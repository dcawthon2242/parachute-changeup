#!/usr/bin/env Rscript

# (1) Does the offspeed spin cue need rebuilding too, or is spin_sim fine there?
#     Compare, at the pitcher level, the old kernel against its two ingredients:
#       spin_sim  = exp(-(|d active|/0.10)^2) * exp(-(d axis/45)^2)
#       axis_sim  = -axis_diff                    (geometry only)
#       act_sim   = -|d active spin|              (efficiency only)
#     and check which survives an out-of-sample test (cue on 2023-25 -> 2026 residual).
#
# (2) The top residual examples that actually are spin-similarity cases: cue in the top
#     quartile for its own pitch type, ranked by overperformance.

suppressPackageStartupMessages({ library(data.table) })
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")
MIN_N <- 200; MIN_26 <- 120

oof <- readRDS(file.path(MDIR, "oof_whiff_resid.rds"))
oof[, res_loc := NA_real_]
for (g in c("breaking","offspeed")) {
  i <- which(oof$grp == g); s <- oof[i]
  oof$res_loc[i] <- residuals(lm(wres ~ poly(plate_x,3)*poly(plate_z,3)+below_zone+VAA+HAA, data=s))
}
d   <- readRDS(file.path(MDIR, "miss_grade_data_activespin.rds"))
asl <- readRDS(file.path(MDIR, "active_spin_long.rds"))
fb  <- asl[pitch_type %in% c("FF","SI","FC")][, pr := match(pitch_type, c("FF","SI","FC"))][
  order(pitcher, season, pr)][, .SD[1], by = .(pitcher, season)][, .(pitcher, season, fb_active = active_spin)]
d <- merge(d, fb, by = c("pitcher","season"), all.x = TRUE)
key <- c("game_pk","at_bat_number","pitch_number","season")
setkeyv(d, key); setkeyv(oof, key)
b <- d[, .(game_pk, at_bat_number, pitch_number, season, axis_diff, active_spin, fb_active)][oof]
b <- b[is.finite(res_loc) & is.finite(axis_diff)]
b[, `:=`(as_gap = active_spin - fb_active)]
b[, spin_sim := exp(-(abs(as_gap)/0.10)^2) * exp(-(axis_diff/45)^2)]
b[, `:=`(axis_sim = -axis_diff, act_sim = -abs(as_gap))]

agg <- function(D, minn) D[, .(player_name = player_name[1], n = .N,
  resid_pp = 100*mean(res_loc), whiff_pct = 100*mean(is_whiff),
  axis_diff = mean(axis_diff), as_gap = mean(as_gap, na.rm=TRUE),
  active_spin = mean(active_spin, na.rm=TRUE), fb_active = mean(fb_active, na.rm=TRUE),
  spin_sim = mean(spin_sim, na.rm=TRUE), axis_sim = -mean(axis_diff),
  act_sim = -mean(abs(as_gap), na.rm=TRUE)), by = .(pitcher, pitch_type)][n >= minn]

################################################################################
cat("=========== 1. OFFSPEED: does the cue need rebuilding? ===========\n")
CUES <- c(spin_sim = "spin_sim (the old kernel)", axis_sim = "axis geometry only",
          act_sim  = "active-spin efficiency only")
for (ty in c("CH","FS")) {
  A  <- agg(b[pitch_type == ty], MIN_N)
  TR <- agg(b[pitch_type == ty & season <= 2025], MIN_N)
  TE <- agg(b[pitch_type == ty & season == 2026], MIN_26)
  J  <- merge(TR[, .(pitcher, spin_sim, axis_sim, act_sim)],
              TE[, .(pitcher, resid26 = resid_pp)], by = "pitcher")
  cat(sprintf("\n--- %s: %d pitchers pooled, %d with a 2023-25 read and a 2026 result ---\n",
              ty, nrow(A), nrow(J)))
  print(rbindlist(lapply(names(CUES), function(v) {
    ct <- cor.test(A[[v]], A$resid_pp, method="spearman", exact=FALSE)
    cp <- if (nrow(J) >= 12) cor(J[[v]], J$resid26, method="spearman") else NA_real_
    data.table(cue = CUES[[v]], r_pooled = round(unname(ct$estimate),3),
               p = signif(ct$p.value,2), r_out_of_sample = round(cp,3))
  })))
}
cat("\n  reference: how much room each ingredient has on offspeed\n")
print(b[pitch_type %in% c("CH","FS"), .(pitches = .N,
        axis_diff_med = round(median(axis_diff),1),
        as_gap_med = round(median(as_gap, na.rm=TRUE),3),
        spin_sim_med = round(median(spin_sim, na.rm=TRUE),3),
        spin_sim_iqr = round(IQR(spin_sim, na.rm=TRUE),3)), by = pitch_type])
cat("  (compare breaking: spin_sim median ~1e-11, no spread at all)\n")

################################################################################
cat("\n\n=========== 2. TOP RESIDUAL EXAMPLES THAT ARE SPIN-SIMILARITY CASES ===========\n")
A <- agg(b[pitch_type %in% c("SL","ST","CU","KC","CH","FS")], MIN_N)
# "a similarity case" = the cue is in the top quartile for that pitch type. Breaking
# balls are judged on axis geometry (the validated cue); offspeed on the full kernel,
# which has real resolution there.
A[, cue := fifelse(pitch_type %in% c("CH","FS"), spin_sim, axis_sim)]
A[, cue_pct := frank(cue)/.N, by = pitch_type]
TOP <- A[cue_pct >= 0.75][order(-resid_pp)]
OUT <- head(TOP[, .(player_name, pitch = pitch_type, n,
                    axis_gap = round(axis_diff,1),
                    own_active = round(active_spin,2), fb_active = round(fb_active,2),
                    as_gap = round(as_gap,3),
                    whiff_pct = round(whiff_pct,1), resid_pp = round(resid_pp,2),
                    cue_pct = round(100*cue_pct))], 10)
print(OUT)
cat("\n  for contrast, the 10 biggest overperformers with the cue in the BOTTOM quartile:\n")
print(head(A[cue_pct <= 0.25][order(-resid_pp),
  .(player_name, pitch = pitch_type, n, axis_gap = round(axis_diff,1),
    as_gap = round(as_gap,3), whiff_pct = round(whiff_pct,1),
    resid_pp = round(resid_pp,2), cue_pct = round(100*cue_pct))], 10))
cat(sprintf("\n  mean residual, top-quartile similarity: %+.2f pp (n=%d)\n",
            mean(TOP$resid_pp), nrow(TOP)))
cat(sprintf("  mean residual, bottom-quartile similarity: %+.2f pp (n=%d)\n",
            mean(A[cue_pct <= 0.25]$resid_pp), nrow(A[cue_pct <= 0.25])))
print(t.test(TOP$resid_pp, A[cue_pct <= 0.25]$resid_pp)[c("statistic","p.value")])

fwrite(A[order(-resid_pp)], file.path(AST, "ext_spin_similarity_examples.csv"))
cat("\nwrote ext_spin_similarity_examples.csv\n")
