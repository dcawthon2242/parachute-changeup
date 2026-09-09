#!/usr/bin/env Rscript

# Why do stuff models love hard sliders if slider velo barely predicts pair xRV?
# Stuff models (1) grade one pitch, not a pair, (2) are usually fit on whiff/CSW,
# (3) include everyone, not just volume MLB secondaries.

suppressPackageStartupMessages({ library(data.table); library(splines) })
set.seed(7)
options(width = 210)

cols <- c("game_year","game_type","pitcher","batter","stand","pitch_type","description",
          "balls","strikes","release_speed","pfx_x","pfx_z","release_spin_rate",
          "release_extension","plate_x","plate_z","sz_top","sz_bot",
          "estimated_woba_using_speedangle","delta_run_exp")
dt <- rbindlist(lapply(c(2025, 2026), function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress=FALSE, select=cols)))
setnames(dt, "estimated_woba_using_speedangle", "xw")
dt <- dt[game_type=="R" & pitch_type %in% c("SL","ST") & balls<=3 & strikes<=2 &
         is.finite(delta_run_exp) & is.finite(release_speed)]
dt[, `:=`(bat_rv=delta_run_exp,
          in_zone=is.finite(plate_x) & abs(plate_x)<=0.83 &
                  is.finite(plate_z) & plate_z>=sz_bot & plate_z<=sz_top)]
WH <- c("swinging_strike","swinging_strike_blocked","foul_tip")
SW <- c(WH,"foul","hit_into_play")
dt[, swung := description %in% SW]
dt[, csw   := description %in% c(WH, "called_strike")]
mb <- lm(bat_rv ~ ns(xw,5)+factor(paste0(balls,"-",strikes)),
         data=dt[description=="hit_into_play" & is.finite(xw)])
dt[, xrv_p := bat_rv]
dt[description=="hit_into_play" & is.finite(xw),
   xrv_p := predict(mb, newdata=dt[description=="hit_into_play" & is.finite(xw)])]
dt <- dt[!(description=="hit_into_play" & !is.finite(xw))]

# pitcher FB velo for the gap / "hard thrower" proxy
FB <- rbindlist(lapply(c(2025, 2026), function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress=FALSE, select=c("game_year","game_type","pitcher","pitch_type","release_speed"))))
FB <- FB[game_type=="R" & pitch_type %in% c("FF","SI","FC") & is.finite(release_speed)]
FBm <- FB[, .(fb_velo=mean(release_speed), fb_n=.N), by=pitcher]

# pitcher-pitch season
PP <- dt[, .(
  n=.N,
  sl_velo=mean(release_speed),
  spin=mean(release_spin_rate, na.rm=TRUE),
  hmov=mean(abs(pfx_x), na.rm=TRUE),
  vmov=mean(pfx_z, na.rm=TRUE),
  ext=mean(release_extension, na.rm=TRUE),
  xrv=-100*mean(xrv_p),
  rv=-100*mean(bat_rv),
  whiff=100*mean(description %in% WH),
  csw=100*mean(csw),
  chase=100*mean(swung & !in_zone, na.rm=TRUE),
  swing=100*mean(swung)
), by=.(game_year, pitcher, pitch_type)]
PP <- merge(PP, FBm, by="pitcher", all.x=TRUE)
PP[, gap := fb_velo - sl_velo]
PP <- PP[n >= 80 & is.finite(gap)]
cat(sprintf("Pitcher-seasons, SL/ST, 80+ pitches: %d  (SL %d, ST %d)\n\n",
            nrow(PP), nrow(PP[pitch_type=="SL"]), nrow(PP[pitch_type=="ST"])))

# =============================================================================
# 1. What stuff models see: velo vs the outcomes they train on
# =============================================================================
cat("############ 1. Slider velo as a predictor of stuff-like vs run-value outcomes ############\n\n")
cat("  SL only, one season-pitcher per row. Standardized later.\n")
SL <- PP[pitch_type=="SL"]
fit <- function(y, rhs, d, lab) {
  m <- lm(as.formula(paste(y, "~", rhs)), data=d)
  s <- summary(m)
  cat(sprintf("  -- %s | R2 = %.3f  n = %d --\n", lab, s$r.squared, nobs(m)))
  co <- as.data.table(s$coefficients, keep.rownames="term")
  setnames(co, c("term","b","se","t","p"))
  print(co[term!="(Intercept)", .(term, b=round(b,4), t=round(t,2), p=round(p,4))],
        row.names=FALSE)
}
fit("whiff", "sl_velo", SL, "Whiff % ~ slider velo")
fit("csw",   "sl_velo", SL, "CSW% ~ slider velo")
fit("xrv",   "sl_velo", SL, "xRV ~ slider velo")
fit("whiff", "sl_velo + gap", SL, "Whiff % ~ slider velo + gap")
fit("xrv",   "sl_velo + gap", SL, "xRV ~ slider velo + gap")
fit("xrv",   "sl_velo + fb_velo", SL, "xRV ~ slider velo + FB velo (the collinear pair)")

cat("\n  Standardized, both in the model (1 SD of each):\n")
Z <- copy(SL)
for (v in c("sl_velo","gap","fb_velo","whiff","csw","xrv"))
  Z[, paste0(v,"_z") := as.numeric(scale(get(v)))]
fit("whiff_z", "sl_velo_z + gap_z", Z, "Whiff, standardized")
fit("xrv_z",   "sl_velo_z + gap_z", Z, "xRV, standardized")

# =============================================================================
# 2. The omitted-variable story: slider velo is a hard-thrower name tag
# =============================================================================
cat("\n############ 2. What else comes with a harder slider? ############\n\n")
cat(sprintf("  cor(SL velo, FB velo) = %+.3f\n", cor(SL$sl_velo, SL$fb_velo)))
cat(sprintf("  cor(SL velo, gap)     = %+.3f\n", cor(SL$sl_velo, SL$gap)))
cat(sprintf("  cor(SL velo, spin)    = %+.3f\n", cor(SL$sl_velo, SL$spin, use="complete")))
cat(sprintf("  cor(SL velo, |h-mov|) = %+.3f\n", cor(SL$sl_velo, SL$hmov, use="complete")))
cat(sprintf("  cor(SL velo, v-mov)   = %+.3f\n", cor(SL$sl_velo, SL$vmov, use="complete")))
cat(sprintf("  cor(gap,     whiff)   = %+.3f\n", cor(SL$gap, SL$whiff)))
cat(sprintf("  cor(SL velo, whiff)   = %+.3f\n", cor(SL$sl_velo, SL$whiff)))
cat(sprintf("  cor(gap,     xrv)     = %+.3f\n", cor(SL$gap, SL$xrv)))
cat(sprintf("  cor(SL velo, xrv)     = %+.3f\n", cor(SL$sl_velo, SL$xrv)))
cat(sprintf("  cor(FB velo, xrv)     = %+.3f\n", cor(SL$fb_velo, SL$xrv)))

SL[, vQ := cut(sl_velo, quantile(sl_velo, 0:4/4), include.lowest=TRUE,
               labels=c("Q1 softest","Q2","Q3","Q4 hardest"))]
cat("\n  SL quartiles — the stuff-model view vs the arsenal view:\n")
print(SL[, .(n=.N,
             sl=round(mean(sl_velo),1), fb=round(mean(fb_velo),1),
             gap=round(mean(gap),1), spin=round(mean(spin),0),
             hmov=round(mean(hmov),2), vmov=round(mean(vmov),2),
             whiff=round(mean(whiff),1), csw=round(mean(csw),1),
             xrv=round(mean(xrv),3)), by=vQ][order(vQ)], row.names=FALSE)

# =============================================================================
# 3. A toy stuff model: velo + movement + spin, no arsenal
# =============================================================================
cat("\n############ 3. Toy stuff model (what a pitch-only model would learn) ############\n\n")
SLok <- SL[is.finite(spin) & is.finite(hmov) & is.finite(vmov) & is.finite(ext)]
# trained the way stuff models often are: on CSW / whiff, not on xRV
mStuff <- lm(whiff ~ sl_velo + hmov + vmov + spin + ext, data=SLok)
cat("  Fit on whiff % (the usual stuff target):\n")
print(round(summary(mStuff)$coefficients, 4))
cat(sprintf("  R2 = %.3f\n", summary(mStuff)$r.squared))
SLok[, stuff := predict(mStuff, newdata=SLok)]

cat("\n  Now: does that stuff score predict xRV? And does gap still add?\n")
fit("xrv", "stuff", SLok, "xRV ~ toy stuff")
fit("xrv", "stuff + gap", SLok, "xRV ~ toy stuff + velo gap")
cat(sprintf("  cor(stuff, sl_velo) = %+.3f   — stuff is mostly the velo term\n",
            cor(SLok$stuff, SLok$sl_velo)))
cat(sprintf("  cor(stuff, gap)     = %+.3f\n", cor(SLok$stuff, SLok$gap)))
cat(sprintf("  cor(stuff, xrv)     = %+.3f\n", cor(SLok$stuff, SLok$xrv)))
cat(sprintf("  cor(gap, xrv)       = %+.3f\n", cor(SLok$gap, SLok$xrv)))

# =============================================================================
# 4. The selection point: include low-volume sliders (stuff models do)
# =============================================================================
cat("\n############ 4. Add the low-volume sliders stuff models also see ############\n\n")
ALL <- dt[, .(n=.N, sl_velo=mean(release_speed),
              xrv=-100*mean(xrv_p), whiff=100*mean(description %in% WH)),
          by=.(game_year, pitcher, pitch_type)]
ALL <- merge(ALL, FBm, by="pitcher", all.x=TRUE)
ALL <- ALL[pitch_type=="SL" & n>=20 & is.finite(fb_velo)]
ALL[, gap := fb_velo - sl_velo]
cat(sprintf("  SL pitcher-seasons with 20+ pitches: %d  (80+ was %d)\n",
            nrow(ALL), nrow(SL)))
fit("whiff", "sl_velo", ALL, "Whiff ~ SL velo, 20+ pitches")
fit("xrv",   "sl_velo", ALL, "xRV ~ SL velo, 20+ pitches")
fit("xrv",   "sl_velo + gap", ALL, "xRV ~ SL velo + gap, 20+ pitches")
ALL[, vQ := cut(sl_velo, quantile(sl_velo, 0:4/4), include.lowest=TRUE,
                labels=c("Q1 softest","Q2","Q3","Q4 hardest"))]
cat("\n  Including the short-sample sliders, xRV by velo quartile:\n")
print(ALL[, .(n=.N, sl=round(mean(sl_velo),1), gap=round(mean(gap),1),
              whiff=round(mean(whiff),1), xrv=round(mean(xrv),3),
              xrv_wt=round(weighted.mean(xrv, n),3)), by=vQ][order(vQ)],
      row.names=FALSE)
