#!/usr/bin/env Rscript

# Where does a large velo gap fail, given that it correlates with
# secondary xRV and produces mistimed / weak contact?

suppressPackageStartupMessages({ library(data.table); library(splines) })
set.seed(7)
options(width = 215)

P <- fread("data/statcast_2026/velo_gap_tunnel_pairs.csv")
ST <- tryCatch(fread("data/statcast_2026/tjstuff_plus_2026_pitch.csv"),
               error=function(e) NULL)

cols <- c("game_year","game_type","pitcher","player_name","p_throws","stand",
          "pitch_type","description","events","balls","strikes",
          "plate_x","plate_z","sz_top","sz_bot","release_speed",
          "estimated_woba_using_speedangle","delta_run_exp",
          "launch_speed","launch_speed_angle","bb_type")
dt <- rbindlist(lapply(c(2025, 2026), function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress=FALSE, select=cols)))
setnames(dt, "estimated_woba_using_speedangle", "xw")
dt <- dt[game_type=="R" & pitch_type!="" & balls<=3 & strikes<=2 & is.finite(delta_run_exp)]
dt[, `:=`(bat_rv=delta_run_exp,
          in_zone=is.finite(plate_x) & abs(plate_x)<=0.83 &
                  is.finite(plate_z) & plate_z>=sz_bot & plate_z<=sz_top)]
WH <- c("swinging_strike","swinging_strike_blocked","foul_tip")
SW <- c(WH,"foul","hit_into_play")
dt[, swung := description %in% SW]
mb <- lm(bat_rv ~ ns(xw,5)+factor(paste0(balls,"-",strikes)),
         data=dt[description=="hit_into_play" & is.finite(xw)])
dt[, xrv_p := bat_rv]
dt[description=="hit_into_play" & is.finite(xw),
   xrv_p := predict(mb, newdata=dt[description=="hit_into_play" & is.finite(xw)])]
dt <- dt[!(description=="hit_into_play" & !is.finite(xw))]
FB <- c("FF","SI","FC"); BR <- c("SL","ST","CU","KC","SV"); OS <- c("CH","FS","FO")
dt[, fam := fifelse(pitch_type %in% FB,"FB",
             fifelse(pitch_type %in% BR,"BR",
             fifelse(pitch_type %in% OS,"OS","OT")))]

# pitcher-year: biggest pair gap, and usage-weighted gap across secondaries
G <- P[, .(
  max_gap = max(velo_gap),
  w_gap   = weighted.mean(velo_gap, n),
  n_sec   = sum(n),
  n_pairs = .N
), by=pitcher]
# 2025-26 book
BK <- dt[, .(
  pitches=.N,
  xrv=-100*mean(xrv_p),
  rv=-100*mean(bat_rv),
  whiff=100*mean(description %in% WH),
  csw=100*mean(description %in% c(WH,"called_strike")),
  zone=100*mean(in_zone, na.rm=TRUE),
  chase=100*mean(swung & !in_zone, na.rm=TRUE),
  chase_ooz=100*mean(swung[!in_zone], na.rm=TRUE),
  cs=100*mean(description=="called_strike"),
  ball=100*mean(description=="ball"),
  k=100*mean(events %in% c("strikeout","strikeout_double_play"), na.rm=TRUE),
  bb=100*mean(events=="walk", na.rm=TRUE),
  xwcon=mean(xw[description=="hit_into_play"], na.rm=TRUE),
  hard=100*mean(launch_speed[description=="hit_into_play"]>=95, na.rm=TRUE),
  barrel=100*mean(launch_speed_angle[description=="hit_into_play"]==6, na.rm=TRUE),
  gb=100*mean(bb_type[description=="hit_into_play"]=="ground_ball", na.rm=TRUE),
  fb_xrv=-100*mean(xrv_p[fam=="FB"]),
  sec_xrv=-100*mean(xrv_p[fam %in% c("BR","OS")]),
  fb_whiff=100*mean(description[fam=="FB"] %in% WH),
  sec_whiff=100*mean(description[fam %in% c("BR","OS")] %in% WH),
  fb_sh=100*mean(fam=="FB"),
  fb_velo=mean(release_speed[fam=="FB"], na.rm=TRUE)
), by=pitcher]
FBc <- dt[fam=="FB" & description=="hit_into_play" & is.finite(xw),
          .(fb_xwcon=mean(xw), fb_hard=100*mean(launch_speed>=95,na.rm=TRUE),
            fb_brl=100*mean(launch_speed_angle==6,na.rm=TRUE)), by=pitcher]
SEc <- dt[fam %in% c("BR","OS") & description=="hit_into_play" & is.finite(xw),
          .(sec_xwcon=mean(xw), sec_hard=100*mean(launch_speed>=95,na.rm=TRUE),
            sec_brl=100*mean(launch_speed_angle==6,na.rm=TRUE)), by=pitcher]
BK <- merge(merge(BK, FBc, by="pitcher", all.x=TRUE), SEc, by="pitcher", all.x=TRUE)
BK <- merge(BK, G, by="pitcher")
BK <- BK[pitches>=700 & is.finite(max_gap)]
BK[, g3 := cut(w_gap, quantile(w_gap, 0:3/3), include.lowest=TRUE,
               labels=c("small gap","mid gap","big gap"))]

cat(sprintf("Pitchers (700+ pitches, 2025-26): %d\n\n", nrow(BK)))

# =============================================================================
# 1. Whole-book outcomes by gap tercile
# =============================================================================
cat("############ 1. Whole-book outcomes by usage-weighted velo gap ############\n\n")
print(BK[, .(n=.N, gap=round(mean(w_gap),1), fb=round(mean(fb_velo),1),
             xrv=round(mean(xrv),3), rv=round(mean(rv),3),
             k=round(mean(k),1), bb=round(mean(bb),1), kbb=round(mean(k-bb),1),
             whiff=round(mean(whiff),1), csw=round(mean(csw),1),
             zone=round(mean(zone),1), chase=round(mean(chase_ooz),1),
             xwcon=round(mean(xwcon),3), hard=round(mean(hard),1),
             barrel=round(mean(barrel),1), gb=round(mean(gb),1)),
         by=g3][order(g3)], row.names=FALSE)

cat("\n  Big vs small, whole book:\n")
tt <- function(v, lab, d=BK) {
  t <- t.test(d[g3=="big gap"][[v]], d[g3=="small gap"][[v]])
  cat(sprintf("     %-10s  %+0.3f vs %+0.3f   diff %+0.3f   p = %.3f\n",
              lab, mean(d[g3=="big gap"][[v]],na.rm=TRUE),
              mean(d[g3=="small gap"][[v]],na.rm=TRUE),
              mean(d[g3=="big gap"][[v]],na.rm=TRUE)-mean(d[g3=="small gap"][[v]],na.rm=TRUE),
              t$p.value))
}
for (v in c("xrv","rv","k","bb","whiff","csw","zone","chase_ooz","xwcon","hard","barrel","fb_velo"))
  tt(v, v)

# =============================================================================
# 2. Split the book: fastball vs secondary
# =============================================================================
cat("\n############ 2. The split: FB pays for the secondary ############\n\n")
print(BK[, .(n=.N,
             fb_xrv=round(mean(fb_xrv),3), sec_xrv=round(mean(sec_xrv),3),
             fb_wh=round(mean(fb_whiff),1), sec_wh=round(mean(sec_whiff),1),
             fb_xw=round(mean(fb_xwcon),3), sec_xw=round(mean(sec_xwcon),3),
             fb_hd=round(mean(fb_hard),1),  sec_hd=round(mean(sec_hard),1),
             fb_br=round(mean(fb_brl),1),   sec_br=round(mean(sec_brl),1),
             fb_sh=round(mean(fb_sh),1)),
         by=g3][order(g3)], row.names=FALSE)
cat("\n  Fastball:\n")
for (v in c("fb_xrv","fb_whiff","fb_xwcon","fb_hard","fb_brl")) tt(v, v)
cat("  Secondary:\n")
for (v in c("sec_xrv","sec_whiff","sec_xwcon","sec_hard","sec_brl")) tt(v, v)

# How much of the book-level wash is "worse FB, better secondary"?
cat("\n  Book xRV ~ FB xRV + secondary xRV + FB share, gap dummy:\n")
BK[, big := g3=="big gap"]
print(round(summary(lm(xrv ~ fb_xrv + sec_xrv + fb_sh + big, data=BK))$coefficients, 4))

# =============================================================================
# 3. Walks, zone, chase — the take/command channels
# =============================================================================
cat("\n############ 3. Does the gap cost strikes looking / add walks? ############\n\n")
# secondary-only zone/chase/ball
SEC <- dt[fam %in% c("BR","OS")]
S2 <- merge(SEC, G, by="pitcher")
S2 <- merge(S2, BK[, .(pitcher, g3)], by="pitcher")
Sagg <- S2[, .(
  n=.N,
  zone=100*mean(in_zone, na.rm=TRUE),
  chase_ooz=100*mean(swung[!in_zone], na.rm=TRUE),
  swing_z=100*mean(swung[in_zone], na.rm=TRUE),
  ball=100*mean(description=="ball"),
  cs=100*mean(description=="called_strike"),
  bb_ev=100*mean(events=="walk", na.rm=TRUE)
), by=.(pitcher, g3)]
print(Sagg[, .(pitchers=.N,
               zone=round(mean(zone),1), chase_ooz=round(mean(chase_ooz),1),
               swing_z=round(mean(swing_z),1), ball=round(mean(ball),1),
               cs=round(mean(cs),1)), by=g3][order(g3)], row.names=FALSE)
cat("\n  Secondary only, big vs small:\n")
for (v in c("zone","chase_ooz","swing_z","ball","cs")) {
  t <- t.test(Sagg[g3=="big gap"][[v]], Sagg[g3=="small gap"][[v]])
  cat(sprintf("     %-10s  %+0.2f vs %+0.2f   diff %+0.2f   p = %.3f\n",
              v, mean(Sagg[g3=="big gap"][[v]]), mean(Sagg[g3=="small gap"][[v]]),
              mean(Sagg[g3=="big gap"][[v]])-mean(Sagg[g3=="small gap"][[v]]), t$p.value))
}

# =============================================================================
# 4. Stuff+ vs xRV: the evaluation miss
# =============================================================================
if (!is.null(ST)) {
  cat("\n############ 4. tjStuff+ vs realized/expected RV, by gap ############\n\n")
  # pitcher-level 2026 stuff
  STP <- ST[, .(stuff=weighted.mean(stuff, pitches), stuff_n=sum(pitches)), by=pitcher]
  M <- merge(BK, STP, by="pitcher")
  cat(sprintf("  Pitchers with Stuff+ and 700+ pitches: %d\n", nrow(M)))
  print(M[, .(n=.N, gap=round(mean(w_gap),1),
              stuff=round(mean(stuff),1),
              xrv=round(mean(xrv),3),
              over=round(mean(xrv - (stuff-100)*0.0707, 3))),  # rough: 1 Stuff+ pt ≈ 0.07 RV/100
          by=g3][order(g3)], row.names=FALSE)
  cat(sprintf("\n  cor(gap, Stuff+) = %+.3f   cor(gap, book xRV) = %+.3f\n",
              cor(M$w_gap, M$stuff), cor(M$w_gap, M$xrv)))
  cat(sprintf("  cor(Stuff+, xRV) = %+.3f\n", cor(M$stuff, M$xrv)))
  # residual: xRV after Stuff+ — does gap still explain leftover RV?
  m <- lm(xrv ~ stuff + w_gap, data=M)
  cat("\n  Book xRV ~ Stuff+ + gap  (does gap produce RV that Stuff+ misses?)\n")
  print(round(summary(m)$coefficients, 4))
}

# =============================================================================
# 5. Year-ahead: does 2025 gap predict 2026 book, or only same-year secondary?
# =============================================================================
cat("\n############ 5. Does 2025 gap predict 2026 pitcher value? ############\n\n")
Y25 <- dt[game_year==2025, .(
  p25=.N, xrv25=-100*mean(xrv_p),
  sec25=-100*mean(xrv_p[fam %in% c("BR","OS")]),
  fb25=-100*mean(xrv_p[fam=="FB"])
), by=pitcher]
Y26 <- dt[game_year==2026, .(
  p26=.N, xrv26=-100*mean(xrv_p),
  sec26=-100*mean(xrv_p[fam %in% c("BR","OS")]),
  fb26=-100*mean(xrv_p[fam=="FB"])
), by=pitcher]
# 2025 pair gaps from the pair file — pairs were pooled 2025-26.
# Rebuild 2025-only gap from pitch velo.
G25 <- dt[game_year==2025 & fam=="FB", .(fb=mean(release_speed)), by=pitcher]
G25s <- dt[game_year==2025 & fam %in% c("BR","OS"),
           .(sec_v=mean(release_speed), ns=.N), by=pitcher]
G25 <- merge(G25, G25s, by="pitcher")
G25[, gap25 := fb - sec_v]
YY <- merge(merge(Y25, Y26, by="pitcher"), G25, by="pitcher")
YY <- YY[p25>=400 & p26>=400]
cat(sprintf("  Pitchers in both seasons (400+ each): %d\n", nrow(YY)))
cat(sprintf("  cor(2025 gap, 2026 book xRV)     = %+.3f\n", cor(YY$gap25, YY$xrv26)))
cat(sprintf("  cor(2025 gap, 2026 secondary xRV)= %+.3f\n", cor(YY$gap25, YY$sec26)))
cat(sprintf("  cor(2025 gap, 2026 FB xRV)       = %+.3f\n", cor(YY$gap25, YY$fb26)))
cat(sprintf("  cor(2025 book xRV, 2026 book xRV)= %+.3f\n", cor(YY$xrv25, YY$xrv26)))
m1 <- lm(xrv26 ~ xrv25, YY)
m2 <- lm(xrv26 ~ xrv25 + gap25, YY)
cat(sprintf("  2026 book ~ 2025 book          R2 = %.3f\n", summary(m1)$r.squared))
cat(sprintf("  2026 book ~ 2025 book + gap    R2 = %.3f   gap p = %.3f\n",
            summary(m2)$r.squared, summary(m2)$coefficients["gap25",4]))
m3 <- lm(sec26 ~ sec25 + gap25, YY)
cat(sprintf("  2026 secondary ~ 2025 sec + gap  gap β = %+.3f  p = %.3f\n",
            coef(m3)["gap25"], summary(m3)$coefficients["gap25",4]))

# =============================================================================
# 6. One-line summary of each candidate failure
# =============================================================================
cat("\n############ 6. Failure checklist ############\n")
cat("  See printed tests. Positive 'big minus small' on xrv/k/whiff is a SUCCESS.\n")
cat("  Negative on FB xRV, zone, chase, or year-ahead is a FAILURE.\n")
