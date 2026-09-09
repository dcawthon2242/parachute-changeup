#!/usr/bin/env Rscript

# BREAKING BALLS: REDO THE SPIN CUE ON PHYSICAL QUANTITIES
#
# spin_sim = exp(-(|d active|/0.10)^2) * exp(-(d axis/45deg)^2) was built for OFFSPEED,
# where a changeup can actually match a fastball's spin. Breaking balls sit at 135-168
# deg of axis difference, so the axis term alone is exp(-11) or smaller and the whole
# score collapses into 1e-19..1e-6 with no usable resolution. Ranking on it is not
# wrong arithmetically, but it does not mean "spins like the fastball" -- Kershaw's
# gyro slider (46% active spin vs 85% on his heater) scores in the top decile of it.
#
# So drop the kernel and use the two physical numbers directly, plus the hypothesis
# Kershaw's arsenal actually suggests:
#   gyro_gap  : |active spin - FB active spin|  -- how much of the spin stops being
#               transverse. A gyro slider has a big gap; a true curveball can have none.
#   axis_diff : 0 = same clock face as the fastball, 180 = perfect mirror.
#   mirror    : 180 - axis_diff, so 0 = perfect mirror (Kershaw FF/CU, Glasnow).
# The original premise was that LOOKING alike helps. The mirror idea is the opposite:
# identical spin efficiency with a flipped axis. Both are testable here.

suppressPackageStartupMessages({ library(data.table) })
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")

oof <- readRDS(file.path(MDIR, "oof_whiff_resid.rds"))
oof[, res_loc := NA_real_]
for (g in c("breaking","offspeed")) {
  i <- which(oof$grp == g); s <- oof[i]
  oof$res_loc[i] <- residuals(lm(wres ~ poly(plate_x,3)*poly(plate_z,3)+below_zone+VAA+HAA, data=s))
}

d  <- readRDS(file.path(MDIR, "miss_grade_data_activespin.rds"))
asl <- readRDS(file.path(MDIR, "active_spin_long.rds"))
fb <- asl[pitch_type %in% c("FF","SI","FC")][, pr := match(pitch_type, c("FF","SI","FC"))][
  order(pitcher, season, pr)][, .SD[1], by = .(pitcher, season)][, .(pitcher, season, fb_active = active_spin)]
d <- merge(d, fb, by = c("pitcher","season"), all.x = TRUE)
key <- c("game_pk","at_bat_number","pitch_number","season")
setkeyv(d, key); setkeyv(oof, key)
b <- d[, .(game_pk, at_bat_number, pitch_number, season, active_spin, fb_active, axis_diff)][oof]
b <- b[grp == "breaking" & is.finite(res_loc) & is.finite(active_spin) &
       is.finite(fb_active) & is.finite(axis_diff)]
b[, `:=`(gyro_gap = abs(active_spin - fb_active), mirror = 180 - axis_diff)]

cat("=== breaking-ball spin geometry, league ===\n")
print(b[, .(pitches = .N,
            active_spin = round(mean(active_spin),3), fb_active = round(mean(fb_active),3),
            gyro_gap = round(mean(gyro_gap),3),
            axis_diff = round(mean(axis_diff),1), mirror = round(mean(mirror),1)),
        by = pitch_type][order(-pitches)])

################################################################################
## per-pitch correlations, honest units
################################################################################
cat("\n=== PER-PITCH Spearman vs location-adjusted whiff residual ===\n")
cues <- c(gyro_gap = "gyro gap (|d active spin|)", axis_diff = "axis difference (deg)",
          mirror = "distance from a perfect mirror (180 - d axis)",
          active_spin = "the pitch's own active spin")
res <- rbindlist(lapply(names(cues), function(v) {
  ct <- cor.test(b[[v]], b$res_loc, method = "spearman", exact = FALSE)
  data.table(cue = cues[[v]], n = nrow(b), r = round(unname(ct$estimate), 4),
             p = signif(ct$p.value, 2))
}))
print(res)

cat("\n  within each pitch type (r vs res_loc):\n")
print(dcast(rbindlist(lapply(names(cues), function(v)
  b[, .(cue = cues[[v]], r = round(cor(get(v), res_loc, method="spearman"), 4)), by = pitch_type])),
  cue ~ pitch_type, value.var = "r"))

################################################################################
## pitcher level, where a cue has to survive to be worth naming names
################################################################################
ag <- b[, .(player_name = player_name[1], n = .N, resid_pp = 100*mean(res_loc),
            whiff_pct = 100*mean(is_whiff),
            gyro_gap = mean(gyro_gap), axis_diff = mean(axis_diff), mirror = mean(mirror),
            active_spin = mean(active_spin), fb_active = mean(fb_active)),
        by = .(pitcher, pitch_type)][n >= 200]

cat(sprintf("\n=== PITCHER x PITCH TYPE (n>=200): %d rows ===\n", nrow(ag)))
pl <- rbindlist(lapply(names(cues), function(v) {
  ct <- cor.test(ag[[v]], ag$resid_pp, method = "spearman", exact = FALSE)
  data.table(cue = cues[[v]], n = nrow(ag), r = round(unname(ct$estimate), 3),
             p = signif(ct$p.value, 2))
}))
print(pl)

cat("\n  within pitch type, pitcher level (n = qualifying pitchers):\n")
print(dcast(rbindlist(lapply(names(cues), function(v)
  ag[, .(cue = cues[[v]], np = .N, r = round(cor(get(v), resid_pp, method="spearman"), 3)),
     by = pitch_type])), cue ~ pitch_type, value.var = "r"))
print(ag[, .(qualifying_pitchers = .N), by = pitch_type][order(-qualifying_pitchers)])

################################################################################
## the arsenal Kershaw actually has: mirror pairs
################################################################################
cat("\n=== MIRROR CURVEBALLS: same spin efficiency as the fastball, flipped axis ===\n")
cu <- ag[pitch_type %in% c("CU","KC") & gyro_gap <= 0.06 & mirror <= 25][order(-resid_pp)]
print(head(cu[, .(player_name, pitch_type, n, active_spin = round(active_spin,3),
                  fb_active = round(fb_active,3), gyro_gap = round(gyro_gap,3),
                  axis_diff = round(axis_diff,1), whiff_pct = round(whiff_pct,1),
                  resid_pp = round(resid_pp,2))], 15))
cat(sprintf("\n  mirror curveballs (n=%d): mean residual %+.2f pp\n",
            nrow(cu), mean(cu$resid_pp)))
oth <- ag[pitch_type %in% c("CU","KC") & !(gyro_gap <= 0.06 & mirror <= 25)]
cat(sprintf("  all other curveballs (n=%d): mean residual %+.2f pp\n", nrow(oth), mean(oth$resid_pp)))
cat(sprintf("  Welch t-test p = %.3f\n", t.test(cu$resid_pp, oth$resid_pp)$p.value))

cat("\n=== Kershaw, for reference ===\n")
print(ag[grepl("Kershaw", player_name), .(player_name, pitch_type, n,
  active_spin = round(active_spin,3), fb_active = round(fb_active,3),
  gyro_gap = round(gyro_gap,3), axis_diff = round(axis_diff,1),
  whiff_pct = round(whiff_pct,1), resid_pp = round(resid_pp,2))])

fwrite(ag, file.path(AST, "ext_breaking_spin_geometry.csv"))
cat("\nwrote ext_breaking_spin_geometry.csv\n")
