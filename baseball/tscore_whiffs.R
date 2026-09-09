#!/usr/bin/env Rscript

# Do hitters whiff more at high-tscore pitches?

suppressPackageStartupMessages({ library(data.table) })
options(width = 215)

Q <- fread("data/statcast_2026/tscore_vs_zone_pitches.csv")
A <- tryCatch(fread("data/statcast_2026/tscore_vs_zone_pitchers.csv"), error=function(e) NULL)

cols <- c("game_year","game_type","pitcher","pitch_type","description",
          "balls","strikes","plate_x","plate_z","sz_top","sz_bot")
dt <- rbindlist(lapply(c(2025, 2026), function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress=FALSE, select=cols)))
dt <- dt[game_type=="R" & pitch_type!="" & balls<=3 & strikes<=2 &
         is.finite(plate_x) & is.finite(plate_z)]

dt[, `:=`(
  px = abs(plate_x),
  pz_rel = (plate_z - sz_bot)/pmax(sz_top - sz_bot, 0.1),
  zone = abs(plate_x) <= 0.83 & plate_z >= sz_bot & plate_z <= sz_top,
  heart = abs(plate_x) <= 0.56 &
          (plate_z - sz_bot)/pmax(sz_top - sz_bot, 0.1) >= 0.25 &
          (plate_z - sz_bot)/pmax(sz_top - sz_bot, 0.1) <= 0.75
)]
dt[, shadow := !heart & abs(plate_x) <= 1.11 &
               plate_z >= sz_bot - 0.28 & plate_z <= sz_top + 0.28]
dt[, loc := fifelse(heart, "heart",
             fifelse(zone, "shadow in",
             fifelse(shadow, "shadow out", "chase/waste")))]

WH <- c("swinging_strike","swinging_strike_blocked","foul_tip")
SW <- c(WH, "foul", "hit_into_play")
dt[, whiff := description %in% WH]
dt[, swung := description %in% SW]

key <- unique(Q[, .(pitcher, pitch_type, pgrp, tscore, tscore_a, ts3)])
d <- merge(dt, key, by=c("pitcher","pitch_type"))
cat(sprintf("Tagged pitches: %s  cells: %d\n\n",
            format(nrow(d), big.mark=","), uniqueN(d[, .(pitcher, pitch_type)])))

# cell-level rates so one pair doesn't dominate
cell <- function(sub) {
  sub[, .(
    n=.N,
    swings=sum(swung),
    whiffs=sum(whiff),
    whiff=100*mean(whiff),
    swing=100*mean(swung),
    wswing=100*mean(whiff[swung])
  ), by=.(pitcher, pitch_type, pgrp, ts3, tscore_a)]
}

tt <- function(a, b, lab) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  t <- t.test(a, b)
  cat(sprintf("     %-10s  high %+6.2f vs low %+6.2f   diff %+6.2f   p = %.3f   n=%d/%d\n",
              lab, mean(a), mean(b), mean(a)-mean(b), t$p.value, length(a), length(b)))
}

rep3 <- function(C, title, minn=80) {
  C <- C[n >= minn]
  cat(sprintf("\n############ %s ############\n\n", title))
  cat(sprintf("  Cells with %d+ pitches: %d\n\n", minn, nrow(C)))
  print(C[, .(
    cells=.N, n=round(mean(n),0),
    swing=round(mean(swing),1),
    whiff=round(mean(whiff),1),
    whiff_sw=round(mean(wswing),1)
  ), by=ts3][order(ts3)], row.names=FALSE)
  h <- C[ts3=="high tscore pitch"]; l <- C[ts3=="low tscore pitch"]
  cat("\n  High vs low:\n")
  for (v in c("swing","whiff","wswing")) tt(h[[v]], l[[v]], v)
  if (uniqueN(C$pitch_type) > 2) {
    C[, wh_a := residuals(lm(whiff ~ factor(pitch_type)))]
    C[, ws_a := residuals(lm(wswing ~ factor(pitch_type)))]
    cat("  Type-adjusted:\n")
    tt(C[ts3=="high tscore pitch", wh_a], C[ts3=="low tscore pitch", wh_a], "whiff_a")
    tt(C[ts3=="high tscore pitch", ws_a], C[ts3=="low tscore pitch", ws_a], "wswing_a")
  }
  invisible(C)
}

# 1. All pitches
C0 <- cell(d)
rep3(C0, "1. All pitches", 150)

# 2. In zone / out of zone
rep3(cell(d[zone==TRUE]), "2. In the zone", 60)
rep3(cell(d[zone==FALSE]), "3. Out of the zone", 60)

# 3. Location
for (L in c("heart","shadow in","shadow out","chase/waste")) {
  rep3(cell(d[loc==L]), sprintf("4. Location: %s", L), 40)
}

# 4. Family
for (g in c("FB","BR","OS")) {
  rep3(C0[pgrp==g], sprintf("5. All pitches, %s", g), 150)
  rep3(cell(d[pgrp==g & zone==TRUE]), sprintf("5b. In-zone %s", g), 50)
  rep3(cell(d[pgrp==g & zone==FALSE]), sprintf("5c. Out-of-zone %s", g), 50)
}

# 5. Count
rep3(cell(d[strikes==2]), "6. Two-strike", 40)
rep3(cell(d[strikes<2]), "6b. 0-1 strikes", 80)
rep3(cell(d[strikes==2 & zone==TRUE]), "6c. Two-strike in-zone", 25)

# 6. Continuous: whiff ~ tscore_a + type
cat("\n############ 7. Continuous (cell level) ############\n\n")
C <- C0[n>=150]
cat(sprintf("  cor(tscore_a, whiff)      = %+.3f\n", cor(C$tscore_a, C$whiff)))
cat(sprintf("  cor(tscore_a, whiff/sw)   = %+.3f\n", cor(C$tscore_a, C$wswing)))
print(round(summary(lm(whiff ~ tscore_a + factor(pitch_type), data=C))$coefficients[
  c("tscore_a"),,drop=FALSE], 4))
print(round(summary(lm(wswing ~ tscore_a + factor(pitch_type), data=C))$coefficients[
  c("tscore_a"),,drop=FALSE], 4))

# 7. Arsenal-level reminder
if (!is.null(A) && "k" %in% names(A)) {
  cat("\n############ 8. Arsenal tscore vs book miss (reminder) ############\n\n")
  # book whiff needs to be computed from d at pitcher level
  BK <- d[, .(whiff=100*mean(whiff), k_proxy=100*mean(description %in% WH & strikes==2)),
          by=pitcher]
  # actually use A if it has no whiff; compute from all pitches of those pitchers
  PIT <- dt[, .(whiff=100*mean(description %in% WH),
                swing=100*mean(description %in% SW),
                wswing=100*mean(description[description %in% SW] %in% WH)),
            by=pitcher]
  AA <- merge(A, PIT, by="pitcher")
  cat(sprintf("  cor(arsenal tscore, book whiff)     = %+.3f  p = %.3f\n",
              cor(AA$tscore, AA$whiff), cor.test(AA$tscore, AA$whiff)$p.value))
  cat(sprintf("  cor(arsenal tscore, book whiff/sw)  = %+.3f\n",
              cor(AA$tscore, AA$wswing)))
  cat(sprintf("  cor(arsenal tscore, K%%)             = %+.3f\n",
              cor(AA$tscore, AA$k)))
  print(AA[, .(n=.N,
               whiff=round(mean(whiff),1),
               wswing=round(mean(wswing),1),
               k=round(mean(k),2)), by=ts3][order(ts3)], row.names=FALSE)
}
