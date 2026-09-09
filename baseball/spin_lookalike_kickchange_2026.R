#!/usr/bin/env Rscript

# Three questions:
#  (A) At scale, does miss distance relate to spin-DIRECTION and spin-EFFICIENCY
#      similarity between an offspeed pitch and the fastball?
#  (B) Build a single "spin disguise score" = how FF-like an offspeed's spin looks
#      (direction + efficiency), and validate it against miss distance / whiff / flail.
#  (C) Kick changes (low-spin, low-IVB changeups, e.g. Munoz): does the giant spin
#      -RATE gap actually HURT deception (hitter reads the seams)? Test spin_diff
#      controlling for the velo kill and the axis/eff match.

suppressPackageStartupMessages({ library(data.table) })
dir <- file.path("data", "statcast_2026")

dt <- fread(file.path(dir, "statcast_2026_all.csv"), showProgress = FALSE, select = c(
  "player_name","pitcher","pitch_type","game_type","description",
  "release_speed","release_spin_rate","spin_axis","pfx_x","pfx_z",
  "vx0","vy0","vz0","ax","ay","az","miss_distance"))
dt <- dt[game_type=="R" & !is.na(pfx_x) & pitch_type!=""]

# spin efficiency (Nathan Magnus decomposition), normalized so 99th pct = 1
g <- 32.174; yf <- 17/12
dt[, t := (-vy0 - sqrt(vy0^2 - 2*ay*(50-yf)))/ay]
dt[, `:=`(vxm=vx0+ax*t/2, vym=vy0+ay*t/2, vzm=vz0+az*t/2)]
dt[, vmag := sqrt(vxm^2+vym^2+vzm^2)]
dt[, dotp := (ax*vxm+ay*vym+(az+g)*vzm)/vmag^2]
dt[, amag := sqrt((ax-dotp*vxm)^2+(ay-dotp*vym)^2+((az+g)-dotp*vzm)^2)]
dt[, eff_raw := fifelse(!is.na(release_spin_rate)&release_spin_rate>0, amag/(vmag*release_spin_rate), NA_real_)]
dt[, spin_eff := pmin(eff_raw/quantile(eff_raw,.99,na.rm=TRUE), 1.05)]

whiff <- c("swinging_strike","swinging_strike_blocked","foul_tip","missed_bunt")
swing <- c(whiff,"foul","hit_into_play","foul_bunt","bunt_foul_tip")
dt[, is_whiff := description %in% whiff]; dt[, is_swing := description %in% swing]

cmean <- function(a){a<-a[!is.na(a)]; if(!length(a)) return(NA_real_); r<-a*pi/180; ((atan2(mean(sin(r)),mean(cos(r)))*180/pi)+360)%%360}
circd <- function(a,b){d<-abs(a-b)%%360; pmin(d,360-d)}

agg <- dt[, .(n=.N, velo=mean(release_speed,na.rm=TRUE), spin=mean(release_spin_rate,na.rm=TRUE),
  axis=cmean(spin_axis), eff=mean(spin_eff,na.rm=TRUE), ivb=mean(pfx_z)*12, hb=mean(pfx_x)*12,
  swings=sum(is_swing), whiffs=sum(is_whiff), miss_pbp=mean(miss_distance[is_whiff],na.rm=TRUE)),
  by=.(pitcher, player_name, pitch_type)]
agg[, whiff_pct := whiffs/swings*100]

# anchor on four-seam (fallback SI, FC)
pick <- function(s){for(ft in c("FF","SI","FC")){r<-s[pitch_type==ft & n>=100]; if(nrow(r)) return(r[which.max(n)])}; NULL}
anc <- agg[,{a<-pick(.SD); if(is.null(a)) NULL else a}, by=pitcher, .SDcols=names(agg)]
anc <- anc[,.(pitcher, fb=pitch_type, fb_velo=velo, fb_spin=spin, fb_axis=axis, fb_eff=eff)]

mkdiffs <- function(x){
  x[, axis_diff := circd(axis, fb_axis)]
  x[, eff_diff  := abs(eff - fb_eff)]
  x[, spin_diff := fb_spin - spin]      # + = offspeed spins slower (rate kill)
  x[, velo_diff := fb_velo - velo]
  x[]
}
sec_all <- mkdiffs(merge(agg[n>=15], anc, by="pitcher")[pitch_type!=fb])  # relaxed, for small-sample lookups
sec <- mkdiffs(merge(agg[n>=40], anc, by="pitcher")[pitch_type!=fb])

# merge Savant leaderboard miss distance + flail (per pitcher x pitch type)
fl <- fread(file.path(dir,"savant_swing_timing_2026.csv"), showProgress=FALSE)[
  , .(pitcher=id, pitch_type=api_pitch_type, miss_lb=miss_distance, flail=flailed_percent,
      whiff_lb=whiff_rate, st_n=n_swings)]
sec <- merge(sec, fl, by=c("pitcher","pitch_type"), all.x=TRUE)

pc <- function(x,y,w=NULL){ ok<-is.finite(x)&is.finite(y); if(sum(ok)<15) return(c(NA,NA,sum(ok)))
  ct<-cor.test(x[ok],y[ok]); c(ct$estimate, ct$p.value, sum(ok)) }

# ============================================================================
cat("=== (A) Miss distance vs spin similarity to the FASTBALL, at scale ===\n")
cat("    (axis_diff/eff_diff small = spin LOOKS like the fastball)\n")
for (grp in list(CH="CH", FS="FS", breaking=c("SL","ST","CU","KC"))) {
  d <- sec[pitch_type %in% grp & is.finite(miss_lb)]
  lab <- if(length(grp)==1) grp else "breaking"
  for (v in c("axis_diff","eff_diff","spin_diff","velo_diff")) {
    r <- pc(d[[v]], d$miss_lb)
    cat(sprintf("  miss_dist ~ %-9s [%-8s] r=%+.2f p=%.2g n=%d\n", v, lab, r[1], r[2], r[3]))
  }
}

# ============================================================================
cat("\n=== (B) SPIN DISGUISE SCORE (how FF-like the spin looks: direction + efficiency) ===\n")
off <- sec[pitch_type %in% c("CH","FS") & is.finite(axis_diff) & is.finite(eff_diff)]
off[, z_axis := (axis_diff-mean(axis_diff))/sd(axis_diff)]
off[, z_eff  := (eff_diff -mean(eff_diff)) /sd(eff_diff)]
off[, dissim := sqrt(z_axis^2 + z_eff^2)]
off[, spin_disguise := round(100*(1 - (dissim-min(dissim))/(max(dissim)-min(dissim))),0)]  # 100=most FF-like
cat("  Validation of the score against outcomes (offspeed):\n")
for (y in c("miss_lb","whiff_pct","flail")) {
  r <- pc(off$spin_disguise, off[[y]])
  cat(sprintf("    spin_disguise ~ %-9s r=%+.2f p=%.2g n=%d\n", y, r[1], r[2], r[3]))
}
cat("\n  Highest spin-disguise changeups (min 40 thrown, 25 swings):\n")
top <- off[pitch_type=="CH" & st_n>=25][order(-spin_disguise)]
print(head(top[,.(player_name, score=spin_disguise, axisD=round(axis_diff,1), effD=round(eff_diff,2),
  Dvelo=round(velo_diff,1), Dspin=round(spin_diff), miss=round(miss_lb,1), whiff=round(whiff_pct,1),
  flail=round(flail,2))], 15))

# ============================================================================
cat("\n=== (C) KICK CHANGES: does a giant spin-RATE gap hurt deception? ===\n")
ch <- sec[pitch_type=="CH" & is.finite(miss_lb) & is.finite(spin_diff)]
ch[, kick := spin <= 1600 & ivb <= 2]      # low-spin, low-ride = kick change
cat(sprintf("  Kick-change CH: %d of %d qualified changeups\n", sum(ch$kick), nrow(ch)))
cat("\n  Changeups binned by SPIN RATE (quartile):\n")
ch[, spin_q := cut(spin, quantile(spin, 0:4/4), include.lowest=TRUE, labels=c("Q1 low","Q2","Q3","Q4 high"))]
print(ch[, .(n=.N, spin=round(mean(spin)), ivb=round(mean(ivb),1), Dspin=round(mean(spin_diff)),
  Dvelo=round(mean(velo_diff),1), miss=round(mean(miss_lb),1), whiff=round(mean(whiff_pct),1),
  flail=round(mean(flail,na.rm=TRUE),2)), by=spin_q][order(spin_q)])
cat("\n  Kick change vs the rest:\n")
print(ch[, .(n=.N, spin=round(mean(spin)), ivb=round(mean(ivb),1), Dspin=round(mean(spin_diff)),
  miss=round(mean(miss_lb),1), whiff=round(mean(whiff_pct),1), flail=round(mean(flail,na.rm=TRUE),2)),
  by=.(kick)])
cat("\n  Seam-visibility test (CH): outcome ~ velo_diff + spin_diff + axis_diff + eff_diff (weighted by swings)\n")
for (y in c("miss_lb","whiff_pct")) {
  f <- lm(as.formula(paste(y,"~ velo_diff + spin_diff + axis_diff + eff_diff")), data=ch, weights=swings)
  cat(sprintf("  -- %s --\n", y)); print(round(summary(f)$coefficients, 4))
}
cat("\n  Munoz kick change (relaxed sample, no leaderboard minimum):\n")
print(sec_all[pitch_type=="CH" & grepl("Muñoz, Andr|Munoz, Andr", player_name),
  .(player_name, n, spin=round(spin), ivb=round(ivb,1), Dspin=round(spin_diff), Dvelo=round(velo_diff,1),
    axisD=round(axis_diff,1), effD=round(eff_diff,2), whiff=round(whiff_pct,1))])

fwrite(sec, file.path(dir,"spin_lookalike_2026.csv"))
cat(sprintf("\nWrote %s\n", file.path(dir,"spin_lookalike_2026.csv")))
