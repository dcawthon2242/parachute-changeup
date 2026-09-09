#!/usr/bin/env Rscript

# IS SCRAPING D1 HEIGHTS WORTH IT? MEASURE THE PAYOFF ON MLB FIRST.
#
# The arm-slot proxy is what broke the college transfer: fit on release height, release side and
# extension alone it reaches r = .765 overall, but at the 44-degree gate that the parachute bin
# actually uses it runs 58 percent precision and 42 percent recall, and a conjunction of four
# conditions cannot survive that.
#
# Height is the obvious missing input. Statcast's arm angle is measured from the release point
# relative to the SHOULDER, and shoulder position is a function of how tall the pitcher is - so
# two pitchers releasing the ball at the same absolute height have different arm angles if they
# are different sizes. The proxy currently has no way to know that.
#
# data/pitcher_heights.csv has 880 MLB pitchers with listed height, so the experiment can be run
# where truth exists before committing to a scrape that has to get past Akamai. If height moves
# the 44-degree gate materially, the scrape is worth the effort; if it does not, no amount of
# roster data will rescue the transfer and the honest move is to drop the arm condition.

suppressPackageStartupMessages({ library(data.table); library(mgcv) })
set.seed(23); options(width = 200)
MDIR <- "data/statcast_model"

H <- fread("data/pitcher_heights.csv")
M <- readRDS(file.path(MDIR, "parachute_rv.rds"))
RS <- readRDS(file.path(MDIR, "whiff_tjstuff.rds"))[, .(arm = mean(arm_angle)), by = .(pitcher, season)]
REL <- M[pitch_type == "FF" & is.finite(release_pos_x) & is.finite(release_pos_z) & is.finite(release_extension),
         .(n = .N, rz = mean(release_pos_z), rx = mean(abs(release_pos_x)),
           ext = mean(release_extension)), by = .(pitcher, season)][n >= 50]
A <- merge(merge(REL, RS, by = c("pitcher","season")), H[, .(pitcher, ht = ht_in)], by = "pitcher")
A <- A[is.finite(ht)]
cat(sprintf("%d pitcher-seasons with measured arm angle and listed height (%d pitchers)\n",
            nrow(A), uniqueN(A$pitcher)))
cat(sprintf("height ranges %d to %d inches, sd %.1f\n\n", min(A$ht), max(A$ht), sd(A$ht)))

# Shoulder-relative geometry, which is what Statcast's arm angle actually measures. Shoulder
# height and half-width both scale with stature; the constants are the usual anthropometric
# fractions and are refined by the smooth term rather than trusted exactly.
A[, `:=`(sh_z = 0.70*ht/12, sh_x = 0.115*ht/12)]
A[, geo := atan2(rz - sh_z, pmax(rx - sh_x, .01))*180/pi]

SPECS <- list(
  "release point only (current proxy)" = arm ~ s(rz) + s(rx) + s(ext),
  "release point + height"             = arm ~ s(rz) + s(rx) + s(ext) + s(ht),
  "shoulder-relative geometry"         = arm ~ s(geo),
  "geometry + release point + height"  = arm ~ s(geo) + s(rz) + s(rx) + s(ext) + s(ht))

A[, fold := sample(rep_len(1:5, .N))]
res <- rbindlist(lapply(names(SPECS), function(nm) {
  A[, hat := NA_real_]
  for (k in 1:5) A[fold == k, hat := predict(gam(SPECS[[nm]], data = A[fold != k]), .SD)]
  g <- A$hat >= 44; t <- A$arm >= 44
  data.table(spec = nm, r = round(cor(A$hat, A$arm),3),
             rmse = round(sqrt(mean((A$hat-A$arm)^2)),2),
             precision = round(100*mean(t[g]),1), recall = round(100*mean(g[t]),1),
             f1 = round(200*mean(t[g])*mean(g[t])/(mean(t[g])+mean(g[t])),1)) }))
cat("recovering the 44-degree gate, 5-fold out of fold (base rate ",
    sprintf("%.1f%%):\n", 100*mean(A$arm >= 44)), sep = "")
print(res, row.names = FALSE)

# The conjunction is what actually matters: the bin needs the arm gate AND three other conditions
# to survive together, so recall on the gate compounds against everything else.
best <- res[which.max(f1)]
cat(sprintf("\nbest spec is '%s'\n", best$spec))
cat(sprintf("gate recall goes from %.0f%% to %.0f%%, so expected membership recovery in a\n",
            res[spec == "release point only (current proxy)", recall], best$recall))
cat(sprintf("four-condition conjunction improves by a factor of about %.1fx from this term alone.\n",
            best$recall/res[spec == "release point only (current proxy)", recall]))

cat("\nhow much of the error is height doing? residual sd of measured arm angle:\n")
A[, hat0 := NA_real_]; A[, hat1 := NA_real_]
for (k in 1:5) {
  A[fold == k, hat0 := predict(gam(SPECS[[1]], data = A[fold != k]), .SD)]
  A[fold == k, hat1 := predict(gam(SPECS[[4]], data = A[fold != k]), .SD)] }
cat(sprintf("  without height %.2f degrees, with height %.2f degrees (%.0f%% reduction)\n",
            sd(A$arm-A$hat0), sd(A$arm-A$hat1), 100*(1-sd(A$arm-A$hat1)/sd(A$arm-A$hat0))))
