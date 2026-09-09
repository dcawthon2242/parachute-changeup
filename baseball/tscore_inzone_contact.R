#!/usr/bin/env Rscript

# Do high-tscore pitches get hit harder in the zone than low-tscore pitches?
# Meatball (heart) vs shadow.

suppressPackageStartupMessages({ library(data.table); library(splines) })
options(width = 215)

Q <- fread("data/statcast_2026/tscore_vs_zone_pitches.csv")

cols <- c("game_year","game_type","pitcher","stand","pitch_type","description",
          "balls","strikes","plate_x","plate_z","sz_top","sz_bot",
          "launch_speed","launch_angle","launch_speed_angle",
          "estimated_woba_using_speedangle","delta_run_exp","bb_type")
dt <- rbindlist(lapply(c(2025, 2026), function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress=FALSE, select=cols)))
setnames(dt, "estimated_woba_using_speedangle", "xw")
dt <- dt[game_type=="R" & pitch_type!="" & balls<=3 & strikes<=2 &
         is.finite(plate_x) & is.finite(plate_z) & is.finite(sz_top) & is.finite(sz_bot)]

dt[, `:=`(
  px = abs(plate_x),
  pz_rel = (plate_z - sz_bot)/pmax(sz_top - sz_bot, 0.1)
)]
# Savant-ish attack zones. 1 baseball ≈ 0.24–0.28 ft.
# Heart / meatball: inner rectangle of the zone.
# Shadow: 1-ball ring straddling the zone edge (in AND out).
dt[, heart := px <= 0.56 & pz_rel >= 0.25 & pz_rel <= 0.75]
dt[, zone  := px <= 0.83 & plate_z >= sz_bot & plate_z <= sz_top]
dt[, shadow := !heart & px <= 1.11 &
               plate_z >= sz_bot - 0.28 & plate_z <= sz_top + 0.28]
dt[, shadow_in  := shadow & zone]
dt[, shadow_out := shadow & !zone]
# Stricter dead-red meatball
dt[, meat := px <= 0.40 & pz_rel >= 0.35 & pz_rel <= 0.65]
dt[, loc := fifelse(heart, "heart",
             fifelse(shadow_in, "shadow in",
             fifelse(shadow_out, "shadow out",
             fifelse(zone, "zone other", "chase/waste"))))]

WH <- c("swinging_strike","swinging_strike_blocked","foul_tip")
SW <- c(WH, "foul", "hit_into_play")
dt[, swung := description %in% SW]
dt[, bip := description=="hit_into_play" & is.finite(launch_speed)]

key <- unique(Q[, .(pitcher, pitch_type, pgrp, tscore, tscore_a, ts3)])
d <- merge(dt, key, by=c("pitcher","pitch_type"))
cat(sprintf("Tagged pitches: %s  cells: %d\n\n",
            format(nrow(d), big.mark=","), uniqueN(d[, .(pitcher, pitch_type)])))

# ---- helpers
summ <- function(x) {
  x[, .(
    pitches=.N,
    swings=sum(swung),
    bip=sum(bip),
    ev=mean(launch_speed[bip], na.rm=TRUE),
    hard=100*mean(launch_speed[bip]>=95, na.rm=TRUE),
    brl=100*mean(launch_speed_angle[bip]==6, na.rm=TRUE),
    xw=mean(xw[bip], na.rm=TRUE),
    la=mean(launch_angle[bip], na.rm=TRUE),
    whiff=100*mean(description %in% WH),
    swing=100*mean(swung),
    inplay=100*mean(description=="hit_into_play")
  )]
}

# Pitcher × pitch-type × location first, then average cells (equal weight per pitch type)
cell <- function(sub) {
  a <- sub[, {
    s <- summ(.SD)
    s[, `:=`(pitcher=pitcher[1], pitch_type=pitch_type[1], pgrp=pgrp[1],
             ts3=ts3[1], tscore_a=tscore_a[1])]
    s
  }, by=.(pitcher, pitch_type)]
  a
}

tt <- function(a, b, lab) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a)<8 || length(b)<8) {
    cat(sprintf("     %-8s  n too small\n", lab)); return(invisible(NULL))
  }
  t <- t.test(a, b)
  cat(sprintf("     %-8s  high %+6.2f vs low %+6.2f   diff %+6.2f   p = %.3f   n=%d/%d\n",
              lab, mean(a), mean(b), mean(a)-mean(b), t$p.value, length(a), length(b)))
}

show_loc <- function(sub, title, minbip=15) {
  cat(sprintf("\n############ %s ############\n\n", title))
  C <- cell(sub)
  C <- C[bip >= minbip]
  cat(sprintf("  Cells with %d+ BIP: %d\n\n", minbip, nrow(C)))
  print(C[, .(
    cells=.N, bip=round(mean(bip),0),
    ev=round(mean(ev),1), hard=round(mean(hard),1), brl=round(mean(brl),1),
    xw=round(mean(xw),3), la=round(mean(la),1),
    swing=round(mean(swing),1), whiff=round(mean(whiff),1), inplay=round(mean(inplay),1)
  ), by=ts3][order(ts3)], row.names=FALSE)
  h <- C[ts3=="high tscore pitch"]; l <- C[ts3=="low tscore pitch"]
  cat("\n  High vs low:\n")
  for (v in c("ev","hard","brl","xw","la","swing","whiff","inplay")) tt(h[[v]], l[[v]], v)
  # type-adjusted
  if (nrow(C) > 40 && uniqueN(C$pitch_type) > 2) {
    C[, ev_a := residuals(lm(ev ~ factor(pitch_type)))]
    C[, hard_a := residuals(lm(hard ~ factor(pitch_type)))]
    C[, brl_a := residuals(lm(brl ~ factor(pitch_type)))]
    C[, xw_a := residuals(lm(xw ~ factor(pitch_type)))]
    cat("\n  Type-adjusted (high vs low residual):\n")
    hh <- C[ts3=="high tscore pitch"]; ll <- C[ts3=="low tscore pitch"]
    for (v in c("ev_a","hard_a","brl_a","xw_a")) tt(hh[[v]], ll[[v]], v)
  }
  invisible(C)
}

# =============================================================================
# 1. All in-zone BIP
# =============================================================================
Z <- show_loc(d[zone==TRUE], "1. In-zone contact (all of the zone)", 20)

# =============================================================================
# 2. Heart / meatball
# =============================================================================
H <- show_loc(d[heart==TRUE], "2. Heart / meatball (inner zone)", 12)

# =============================================================================
# 3. Strict dead-red
# =============================================================================
M <- show_loc(d[meat==TRUE], "3. Dead-red meatball (|x|<=0.40, mid 30% height)", 8)

# =============================================================================
# 4. Shadow in-zone (edges of the zone)
# =============================================================================
SI <- show_loc(d[shadow_in==TRUE], "4. In-zone shadow (zone but not heart)", 12)

# =============================================================================
# 5. Shadow just outside
# =============================================================================
SO <- show_loc(d[shadow_out==TRUE], "5. Out-of-zone shadow (just off the edge)", 12)

# =============================================================================
# 6. Full Savant shadow (in + out)
# =============================================================================
SF <- show_loc(d[shadow==TRUE], "6. Full shadow (in-zone edge + just off)", 15)

# =============================================================================
# 7. FB vs OFF at heart and shadow
# =============================================================================
cat("\n############ 7. Fastball vs offspeed, heart and shadow ############\n")
for (g in c("FB","BR","OS")) {
  cat(sprintf("\n  --- %s heart ---\n", g))
  show_loc(d[heart==TRUE & pgrp==g], sprintf("%s heart", g), 10)
  cat(sprintf("\n  --- %s in-zone shadow ---\n", g))
  show_loc(d[shadow_in==TRUE & pgrp==g], sprintf("%s shadow-in", g), 10)
}

# =============================================================================
# 8. Does the meatball-minus-shadow gap differ by tscore?
#     (do high-tscore pitches get crushed relatively more when centered?)
# =============================================================================
cat("\n############ 8. Meatball minus shadow-in, within the same pitch ############\n\n")
CH <- cell(d[heart==TRUE])[bip>=10]
CS <- cell(d[shadow_in==TRUE])[bip>=10]
W <- merge(CH[, .(pitcher, pitch_type, pgrp, ts3, ev_h=ev, hard_h=hard, brl_h=brl, xw_h=xw, n_h=bip)],
           CS[, .(pitcher, pitch_type, ev_s=ev, hard_s=hard, brl_s=brl, xw_s=xw, n_s=bip)],
           by=c("pitcher","pitch_type"))
W[, `:=`(dev=ev_h-ev_s, dhard=hard_h-hard_s, dbrl=brl_h-brl_s, dxw=xw_h-xw_s)]
cat(sprintf("  Pitch types with 10+ BIP in both heart and shadow-in: %d\n\n", nrow(W)))
print(W[, .(
  cells=.N,
  ev_heart=round(mean(ev_h),1), ev_sh=round(mean(ev_s),1), d_ev=round(mean(dev),1),
  hard_heart=round(mean(hard_h),1), hard_sh=round(mean(hard_s),1), d_hard=round(mean(dhard),1),
  xw_heart=round(mean(xw_h),3), xw_sh=round(mean(xw_s),3), d_xw=round(mean(dxw),3)
), by=ts3][order(ts3)], row.names=FALSE)
h <- W[ts3=="high tscore pitch"]; l <- W[ts3=="low tscore pitch"]
cat("\n  High vs low, heart-minus-shadow gap (positive = more extra damage on meatballs):\n")
for (v in c("dev","dhard","dbrl","dxw")) tt(h[[v]], l[[v]], v)

# =============================================================================
# 9. Location mix: do high-tscore pitches just avoid the heart?
# =============================================================================
cat("\n############ 9. Location mix (share of all pitches) ############\n\n")
MIX <- d[, .(
  n=.N,
  heart=100*mean(heart),
  meat=100*mean(meat),
  zone=100*mean(zone),
  sh_in=100*mean(shadow_in),
  sh_out=100*mean(shadow_out)
), by=.(pitcher, pitch_type, ts3, pgrp)]
print(MIX[, .(
  cells=.N,
  heart=round(mean(heart),1), meat=round(mean(meat),1), zone=round(mean(zone),1),
  shadow_in=round(mean(sh_in),1), shadow_out=round(mean(sh_out),1)
), by=ts3][order(ts3)], row.names=FALSE)
cat("\n  Type-adjusted heart share, high vs low:\n")
MIX[, heart_a := residuals(lm(heart ~ factor(pitch_type)))]
tt(MIX[ts3=="high tscore pitch", heart_a], MIX[ts3=="low tscore pitch", heart_a], "heart_a")

# =============================================================================
# 10. Continuous: EV ~ heart * tscore_a + type, on in-zone BIP
# =============================================================================
cat("\n############ 10. On in-zone BIP: EV ~ heart + tscore + type ############\n\n")
bipz <- d[bip==TRUE & zone==TRUE]
# one row per BIP, clustered by attaching cell tscore
m <- lm(launch_speed ~ heart * tscore_a + factor(pitch_type), data=bipz)
print(round(summary(m)$coefficients[c("heartTRUE","tscore_a","heartTRUE:tscore_a"),], 4))
cat("\n  Same for hard-hit (linear prob) and xwOBA:\n")
bipz[, hard := as.integer(launch_speed>=95)]
print(round(summary(lm(hard ~ heart * tscore_a + factor(pitch_type), data=bipz))$coefficients[
  c("heartTRUE","tscore_a","heartTRUE:tscore_a"),], 4))
print(round(summary(lm(xw ~ heart * tscore_a + factor(pitch_type), data=bipz))$coefficients[
  c("heartTRUE","tscore_a","heartTRUE:tscore_a"),], 4))
