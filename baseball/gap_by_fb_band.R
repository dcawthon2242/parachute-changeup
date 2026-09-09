#!/usr/bin/env Rscript

# Similar-velo pitchers: does "harder secondary isn't better" depend on
# whether the FB is 90 or 98? Pitcher-seasons 2020–2026.

suppressPackageStartupMessages({ library(data.table); library(splines) })
set.seed(7)
options(width = 215)

yrs <- 2020:2026
cols <- c("game_year","game_type","pitcher","player_name","p_throws",
          "pitch_type","description","balls","strikes",
          "release_speed","ax","az",
          "plate_x","plate_z","sz_top","sz_bot",
          "launch_speed","estimated_woba_using_speedangle","delta_run_exp")
dt <- rbindlist(lapply(yrs, function(yr) {
  f <- file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr))
  if (!file.exists(f)) return(NULL)
  cat(sprintf("  reading %d\n", yr))
  fread(f, showProgress=FALSE, select=cols)
}), fill=TRUE)
setnames(dt, "estimated_woba_using_speedangle", "xw")
dt <- dt[game_type=="R" & pitch_type!="" & balls<=3 & strikes<=2 & is.finite(delta_run_exp)]
dt[, bat_rv := delta_run_exp]
dt[, in_zone := is.finite(plate_x) & abs(plate_x)<=0.83 &
                is.finite(plate_z) & plate_z>=sz_bot & plate_z<=sz_top]
WH <- c("swinging_strike","swinging_strike_blocked","foul_tip")
SW <- c(WH,"foul","hit_into_play")
dt[, swung := description %in% SW]
dt[, cnt := paste0(balls,"-",strikes)]
dt[, ax_r := fifelse(p_throws=="L", -ax, ax)]
mb <- lm(bat_rv ~ ns(xw,5)+factor(cnt),
         data=dt[description=="hit_into_play" & is.finite(xw)])
dt[, xrv_p := bat_rv]
dt[description=="hit_into_play" & is.finite(xw),
   xrv_p := predict(mb, newdata=dt[description=="hit_into_play" & is.finite(xw)])]
dt <- dt[!(description=="hit_into_play" & !is.finite(xw))]

FB <- dt[pitch_type %in% c("FF","SI","FC") & is.finite(release_speed)]
PRI <- FB[, .N, by=.(game_year, pitcher, pitch_type)][order(-N), .SD[1],
         by=.(game_year, pitcher)]
FBp <- FB[PRI, on=.(game_year, pitcher, pitch_type),
          .(fb_pri=mean(release_speed), fb_pri_n=.N), by=.(game_year, pitcher)]

dt[, fam := fcase(pitch_type=="SL","SL", pitch_type=="ST","ST",
                  pitch_type %in% c("CU","KC"),"CU",
                  pitch_type %in% c("CH","FS"),"CH", default=NA)]
PP <- dt[!is.na(fam), .(
  n=.N, velo=mean(release_speed),
  ax=mean(ax_r, na.rm=TRUE), az=mean(az, na.rm=TRUE),
  xrv=-100*mean(xrv_p),
  chase=100*mean(swung[!in_zone]),
  whiff=100*mean(description %in% WH),
  hard=100*mean(launch_speed[description=="hit_into_play"]>=95, na.rm=TRUE),
  xwcon=mean(xw[description=="hit_into_play"], na.rm=TRUE)
), by=.(game_year, pitcher, fam)]
PP <- merge(PP, FBp, by=c("game_year","pitcher"))
PP[, gap := fb_pri - velo]
PP <- PP[n>=80 & fb_pri_n>=80 & is.finite(gap) & is.finite(ax) & is.finite(az)]
PP[fam=="ST" & game_year<2023, fam := NA]
PP <- PP[!is.na(fam)]

# FB bands: similar-velo groups
PP[, fb_band := fcase(
  fb_pri<92, "<92",
  fb_pri<94, "92–94",
  fb_pri<96, "94–96",
  fb_pri<98, "96–98",
  default="98+"
)]
PP[, fb_band := factor(fb_band, levels=c("<92","92–94","94–96","96–98","98+"))]

cat("Seasons by family × FB band:\n")
print(dcast(PP[, .N, by=.(fam, fb_band)], fam ~ fb_band, value.var="N", fill=0),
      row.names=FALSE)

# within-band, secondary velo split at the band's own median
# so we compare similar-FB pitchers' harder vs softer secondaries
cat("\n############ 1. Within FB band: above- vs below-median secondary velo ############\n")
cat("  Split is WITHIN band, so 98+ hard SL is vs other 98+ arms, not vs 90 FB guys.\n\n")

W <- rbindlist(lapply(c("SL","ST","CU","CH"), function(f) {
  rbindlist(lapply(levels(PP$fb_band), function(b) {
    d <- PP[fam==f & fb_band==b]
    if (nrow(d)<20) return(NULL)
    med <- median(d$velo)
    d[, hi := velo>=med]
    data.table(fam=f, band=b, med=med,
               n_hi=sum(d$hi), n_lo=sum(!d$hi),
               velo_hi=mean(d[hi==TRUE]$velo), velo_lo=mean(d[hi==FALSE]$velo),
               fb_hi=mean(d[hi==TRUE]$fb_pri), fb_lo=mean(d[hi==FALSE]$fb_pri),
               gap_hi=mean(d[hi==TRUE]$gap), gap_lo=mean(d[hi==FALSE]$gap),
               xrv_hi=mean(d[hi==TRUE]$xrv), xrv_lo=mean(d[hi==FALSE]$xrv),
               p_xrv=t.test(d[hi==TRUE]$xrv, d[hi==FALSE]$xrv)$p.value,
               hard_hi=mean(d[hi==TRUE]$hard, na.rm=TRUE),
               hard_lo=mean(d[hi==FALSE]$hard, na.rm=TRUE),
               p_hard=tryCatch(t.test(d[hi==TRUE]$hard, d[hi==FALSE]$hard)$p.value, error=function(e) NA),
               xw_hi=mean(d[hi==TRUE]$xwcon, na.rm=TRUE),
               xw_lo=mean(d[hi==FALSE]$xwcon, na.rm=TRUE),
               chase_hi=mean(d[hi==TRUE]$chase), chase_lo=mean(d[hi==FALSE]$chase),
               whiff_hi=mean(d[hi==TRUE]$whiff), whiff_lo=mean(d[hi==FALSE]$whiff))
  }))
}))
W[, dxrv := xrv_hi - xrv_lo]
print(W[, .(fam, band, n_hi, n_lo,
            velo=sprintf("%.1f vs %.1f", velo_hi, velo_lo),
            fb=sprintf("%.1f vs %.1f", fb_hi, fb_lo),
            xrv_hi=round(xrv_hi,3), xrv_lo=round(xrv_lo,3),
            d_xrv=round(dxrv,3), p=round(p_xrv,3),
            d_hard=round(hard_hi-hard_lo,1), p_h=round(p_hard,3),
            d_xw=round(xw_hi-xw_lo,3))], row.names=FALSE)

# continuous slope of secondary velo inside each band
cat("\n############ 2. Within-band slope: xRV ~ sec velo + ax + az + year ############\n")
cat("  FB already restricted to a 2-mph window, so not also in the model.\n\n")
for (f in c("SL","ST","CU","CH")) {
  cat(sprintf("  %s\n", f))
  for (b in levels(PP$fb_band)) {
    d <- PP[fam==f & fb_band==b]
    if (nrow(d)<30) { cat(sprintf("     %-6s  n=%d skip\n", b, nrow(d))); next }
    m <- lm(xrv ~ velo + ax + az + factor(game_year), data=d)
    cf <- summary(m)$coefficients
    if (!"velo" %in% rownames(cf)) next
    cat(sprintf("     %-6s  n=%4d  β(velo)=%+6.3f p=%.3f   β(ax)=%+6.3f p=%.3f   mean FB=%.1f  mean gap=%.1f\n",
                b, nrow(d), cf["velo","Estimate"], cf["velo","Pr(>|t|)"],
                if ("ax" %in% rownames(cf)) cf["ax","Estimate"] else NA,
                if ("ax" %in% rownames(cf)) cf["ax","Pr(>|t|)"] else NA,
                mean(d$fb_pri), mean(d$gap)))
  }
}

# interaction: does the velo slope change with FB?
cat("\n############ 3. Interaction: xRV ~ velo * FB + ax + az + year ############\n\n")
for (f in c("SL","ST","CU","CH")) {
  d <- PP[fam==f]
  m <- lm(xrv ~ velo * fb_pri + ax + az + factor(game_year), data=d)
  cf <- summary(m)$coefficients
  cat(sprintf("  %s  n=%d\n", f, nrow(d)))
  for (t in c("velo","fb_pri","velo:fb_pri")) {
    if (!t %in% rownames(cf)) next
    cat(sprintf("     %-14s  β=%+7.4f   p=%.3f\n", t, cf[t,"Estimate"], cf[t,"Pr(>|t|)"]))
  }
  # implied velo slope at FB=90, 94, 98
  b0 <- cf["velo","Estimate"]; b1 <- cf["velo:fb_pri","Estimate"]
  cat(sprintf("     implied β(velo) at FB 90 / 94 / 98:  %+0.3f / %+0.3f / %+0.3f\n\n",
              b0+b1*90, b0+b1*94, b0+b1*98))
}

# nearest-neighbor: each season matched to another in same year, |Δ FB|<0.5,
# then compare the one with the harder secondary
cat("\n############ 4. Tight FB match (|Δ FB| ≤ 0.5 mph), same year ############\n")
cat("  Pair two pitcher-seasons in the same year with almost the same FB.\n")
cat("  Winner = the one with the harder secondary. Is that actually better?\n\n")
for (f in c("SL","CH","CU","ST")) {
  d <- PP[fam==f]
  pairs <- rbindlist(lapply(unique(d$game_year), function(yr) {
    a <- d[game_year==yr]
    if (nrow(a)<8) return(NULL)
    # greedy: for each, nearest FB, not self, unused
    setorder(a, fb_pri)
    used <- rep(FALSE, nrow(a))
    out <- vector("list", nrow(a))
    k <- 0L
    for (i in seq_len(nrow(a))) {
      if (used[i]) next
      cand <- which(!used & seq_len(nrow(a))!=i & abs(a$fb_pri - a$fb_pri[i])<=0.5)
      if (!length(cand)) next
      j <- cand[which.min(abs(a$fb_pri[cand] - a$fb_pri[i]))]
      used[i] <- used[j] <- TRUE
      hi <- if (a$velo[i]>=a$velo[j]) i else j
      lo <- if (hi==i) j else i
      k <- k+1L
      out[[k]] <- data.table(
        dfb=abs(a$fb_pri[hi]-a$fb_pri[lo]),
        dvelo=a$velo[hi]-a$velo[lo],
        dxrv=a$xrv[hi]-a$xrv[lo],
        dhard=a$hard[hi]-a$hard[lo],
        dxw=a$xwcon[hi]-a$xwcon[lo],
        dchase=a$chase[hi]-a$chase[lo],
        fb=mean(c(a$fb_pri[hi], a$fb_pri[lo])),
        velo_hi=a$velo[hi], velo_lo=a$velo[lo]
      )
    }
    rbindlist(out[seq_len(k)])
  }))
  if (!nrow(pairs) || nrow(pairs)<20) {
    cat(sprintf("  %s  too few pairs\n", f)); next
  }
  # split matched pairs by the pair's FB level
  pairs[, band := fcase(fb<93, "<93", fb<96, "93–96", default="96+")]
  cat(sprintf("  %s  pairs=%d  mean |ΔFB|=%.2f  mean Δvelo=%.2f\n",
              f, nrow(pairs), mean(pairs$dfb), mean(pairs$dvelo)))
  t <- t.test(pairs$dxrv)
  cat(sprintf("     harder minus softer xRV: %+0.3f   p=%.3f\n", t$estimate, t$p.value))
  print(pairs[, .(
    pairs=.N,
    dvelo=round(mean(dvelo),2),
    dxrv=round(mean(dxrv),3),
    p=round(t.test(dxrv)$p.value,3),
    dhard=round(mean(dhard, na.rm=TRUE),1),
    dxw=round(mean(dxw, na.rm=TRUE),3)
  ), by=band][order(band)], row.names=FALSE)
  cat("\n")
}
