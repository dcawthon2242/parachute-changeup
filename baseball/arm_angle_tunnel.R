#!/usr/bin/env Rscript

# ARM-ANGLE / RELEASE-POINT SIMILARITY as a tunneling cue.
# Question: does a secondary pitch coming from a slot/release point more like the
# fastball's get missed by more (bigger miss_distance = hitter more fooled)?
# Compare two "release-point similarity" measures vs the fastball:
#   (1) arm-angle difference (deg)   -- Hawk-Eye tracked arm_angle
#   (2) release-point Pythagorean distance (ft)
# Target: miss_distance. We look at E[miss | competitive swing] (contact imputed
# to 0 in the model data) AND actual miss_distance on whiffs only.

suppressPackageStartupMessages({ library(data.table) })
MDIR <- file.path("data","statcast_model")
d <- as.data.table(readRDS(file.path(MDIR,"miss_grade_data_activespin.rds")))

## ---- join Hawk-Eye arm_angle (+ p_throws) from raw season files by pitch key ----
raw_files <- c("2023"="data/statcast_2023/statcast_2023_all.csv",
               "2024"="data/statcast_2024/statcast_2024_all.csv",
               "2025"="data/statcast_2025/statcast_2025_all.csv",
               "2026"="data/statcast_2026/statcast_2026_all.csv")
aa <- rbindlist(lapply(raw_files, function(f){
  fread(f, select=c("game_pk","at_bat_number","pitch_number","arm_angle","p_throws"),
        showProgress=FALSE)
}))
aa <- unique(aa, by=c("game_pk","at_bat_number","pitch_number"))
d <- merge(d, aa, by=c("game_pk","at_bat_number","pitch_number"), all.x=TRUE)
cat(sprintf("arm_angle coverage: %.1f%% of model rows (%d / %d)\n",
    100*mean(is.finite(d$arm_angle)), sum(is.finite(d$arm_angle)), nrow(d)))

## ---- fastball anchor (FF>SI>FC) per pitcher-season ----
fb <- d[pitch_type %in% c("FF","SI","FC")]
fb[, pr := match(pitch_type, c("FF","SI","FC"))]
fb <- fb[order(pitcher,season,pr),
   .(fb_arm=mean(arm_angle,na.rm=TRUE), fb_relx=mean(release_pos_x,na.rm=TRUE),
     fb_relz=mean(release_pos_z,na.rm=TRUE), anchor=pitch_type[1], fb_n=.N),
   by=.(pitcher,season)]
d <- merge(d, fb, by=c("pitcher","season"), all.x=TRUE)

## ---- the two release-point similarity cues (vs the fastball) ----
d[, arm_diff := abs(arm_angle - fb_arm)]                                   # deg
d[, rel_dist := sqrt((release_pos_x-fb_relx)^2 + (release_pos_z-fb_relz)^2)] # ft

## sanity on target coding
cat(sprintf("\nmiss_distance: %% exactly 0 = %.1f  |  %% >0 = %.1f  |  %% NA = %.1f\n",
    100*mean(d$miss_distance==0,na.rm=TRUE), 100*mean(d$miss_distance>0,na.rm=TRUE),
    100*mean(is.na(d$miss_distance))))
cat(sprintf("is_whiff mean = %.3f\n", mean(d$is_whiff,na.rm=TRUE)))

## restrict to SECONDARIES (the cue is a pitch-vs-fastball property; FB-vs-itself ~0)
sec <- d[grp %in% c("breaking","offspeed") & is.finite(arm_diff) & is.finite(rel_dist) &
         is.finite(miss_distance)]
cat(sprintf("\nsecondary swings with all fields: %d\n", nrow(sec)))
cat("\n-- distribution of the cues (secondaries) --\n")
cat("arm_diff (deg):   "); print(round(quantile(sec$arm_diff, c(.1,.25,.5,.75,.9,.99)),2))
cat("rel_dist (ft):    "); print(round(quantile(sec$rel_dist, c(.1,.25,.5,.75,.9,.99)),2))
cat(sprintf("corr(arm_diff, rel_dist) = %.3f\n", cor(sec$arm_diff, sec$rel_dist)))

## ============================================================================
## (A) aggregate per pitcher x pitch-type x season  (reduces per-pitch noise)
## ============================================================================
agg <- sec[, .(n=.N,
               mean_miss = mean(miss_distance),               # E[miss|swing], contact=0
               whiff = mean(is_whiff),
               arm_diff = mean(arm_diff),                      # avg per-pitch |dev| (has jitter)
               rel_dist = mean(rel_dist),
               mean_arm = mean(arm_angle), mean_relx=mean(release_pos_x), mean_relz=mean(release_pos_z),
               fb_arm=fb_arm[1], fb_relx=fb_relx[1], fb_relz=fb_relz[1]),
           by=.(pitcher,player_name,pitch_type,grp,season)][n>=40]
# systematic slot difference (group mean vs FB mean -> strips per-pitch jitter, "less noisy")
agg[, arm_sys := abs(mean_arm - fb_arm)]
agg[, rel_sys := sqrt((mean_relx-fb_relx)^2 + (mean_relz-fb_relz)^2)]
# whiff-only mean miss distance
wh <- sec[is_whiff==1, .(whiff_miss=mean(miss_distance), wn=.N), by=.(pitcher,pitch_type,season)]
agg <- merge(agg, wh, by=c("pitcher","pitch_type","season"), all.x=TRUE)
cat(sprintf("\naggregated groups (n>=40 swings): %d\n", nrow(agg)))

wcor <- function(x,y,w){ ok<-is.finite(x)&is.finite(y)&is.finite(w); cov.wt(cbind(x[ok],y[ok]),w[ok],cor=TRUE)$cor[1,2] }

cat("\n=== E[miss | swing] vs the two cues (weighted by n) ===\n")
cat(sprintf("  arm_diff : pearson(w) %+.3f | spearman %+.3f\n",
    wcor(agg$arm_diff, agg$mean_miss, agg$n), cor(agg$arm_diff, agg$mean_miss, method="spearman")))
cat(sprintf("  rel_dist : pearson(w) %+.3f | spearman %+.3f\n",
    wcor(agg$rel_dist, agg$mean_miss, agg$n), cor(agg$rel_dist, agg$mean_miss, method="spearman")))

cat("\n=== whiff-only miss_distance vs the two cues (weighted by whiff n) ===\n")
aw <- agg[is.finite(whiff_miss) & wn>=15]
cat(sprintf("  (groups with >=15 whiffs: %d)\n", nrow(aw)))
cat(sprintf("  arm_diff : pearson(w) %+.3f | spearman %+.3f\n",
    wcor(aw$arm_diff, aw$whiff_miss, aw$wn), cor(aw$arm_diff, aw$whiff_miss, method="spearman")))
cat(sprintf("  rel_dist : pearson(w) %+.3f | spearman %+.3f\n",
    wcor(aw$rel_dist, aw$whiff_miss, aw$wn), cor(aw$rel_dist, aw$whiff_miss, method="spearman")))

cat("\n=== SYSTEMATIC slot diff (group-mean vs FB, jitter stripped) vs E[miss|swing] ===\n")
cat(sprintf("  arm_sys (deg): pearson(w) %+.3f | spearman %+.3f  (median %.2f deg)\n",
    wcor(agg$arm_sys, agg$mean_miss, agg$n), cor(agg$arm_sys, agg$mean_miss, method="spearman"),
    median(agg$arm_sys)))
cat(sprintf("  rel_sys (ft):  pearson(w) %+.3f | spearman %+.3f  (median %.2f ft)\n",
    wcor(agg$rel_sys, agg$mean_miss, agg$n), cor(agg$rel_sys, agg$mean_miss, method="spearman"),
    median(agg$rel_sys)))
cat(sprintf("  arm_sys vs whiff-miss: pearson(w) %+.3f\n",
    wcor(agg[is.finite(whiff_miss)&wn>=15]$arm_sys, agg[is.finite(whiff_miss)&wn>=15]$whiff_miss,
         agg[is.finite(whiff_miss)&wn>=15]$wn)))

cat("\n=== whiff RATE vs the two cues (weighted by n) ===\n")
cat(sprintf("  arm_diff : pearson(w) %+.3f\n", wcor(agg$arm_diff, agg$whiff, agg$n)))
cat(sprintf("  rel_dist : pearson(w) %+.3f\n", wcor(agg$rel_dist, agg$whiff, agg$n)))

## by group
cat("\n=== by group: pearson(w) of cue vs E[miss|swing] ===\n")
for(g in c("breaking","offspeed")){
  a <- agg[grp==g]
  cat(sprintf("  %-9s arm_diff %+.3f | rel_dist %+.3f  (n=%d groups)\n",
      g, wcor(a$arm_diff,a$mean_miss,a$n), wcor(a$rel_dist,a$mean_miss,a$n), nrow(a)))
}

## by pitch type (secondaries with enough groups)
cat("\n=== by pitch type: pearson(w) of arm_diff vs E[miss|swing] ===\n")
for(pt in agg[, .N, by=pitch_type][N>=30][order(-N)]$pitch_type){
  a <- agg[pitch_type==pt]
  cat(sprintf("  %-3s arm_diff %+.3f | rel_dist %+.3f  (n=%d)\n",
      pt, wcor(a$arm_diff,a$mean_miss,a$n), wcor(a$rel_dist,a$mean_miss,a$n), nrow(a)))
}

## ============================================================================
## (B) quantile bins of arm_diff -> mean miss (does similarity => more miss?)
## ============================================================================
cat("\n=== arm_diff quintiles (secondary swings, per-pitch) -> mean miss & whiff ===\n")
sec[, aq := cut(arm_diff, breaks=quantile(arm_diff, 0:5/5), include.lowest=TRUE,
                labels=c("Q1 most similar","Q2","Q3","Q4","Q5 least similar"))]
print(sec[, .(n=.N, arm_diff=round(mean(arm_diff),2), mean_miss=round(mean(miss_distance),3),
              whiff=round(mean(is_whiff),3)), by=aq][order(aq)])
cat("\n=== rel_dist quintiles -> mean miss & whiff ===\n")
sec[, rq := cut(rel_dist, breaks=quantile(rel_dist, 0:5/5), include.lowest=TRUE,
                labels=c("Q1 most similar","Q2","Q3","Q4","Q5 least similar"))]
print(sec[, .(n=.N, rel_dist=round(mean(rel_dist),2), mean_miss=round(mean(miss_distance),3),
              whiff=round(mean(is_whiff),3)), by=rq][order(rq)])

fwrite(agg[order(-mean_miss)], file.path(MDIR,"arm_angle_tunnel.csv"))
cat("\nsaved data/statcast_model/arm_angle_tunnel.csv\n")
