#!/usr/bin/env Rscript

# What correlates with path_to_location_ratio (tunnel tightness) for BREAKING balls
# in 2-strike counts? Lower path_ratio = tighter tunnel = looks like the FB longer.
# We correlate against each pitch's DIFFERENCE vs the pitcher's primary fastball,
# since "looking the same" is inherently a pitch-vs-fastball property.

suppressPackageStartupMessages({ library(data.table) })
MDIR <- file.path("data","statcast_model")
as_long <- readRDS(file.path(MDIR, "active_spin_long.rds"))
d <- readRDS(file.path(MDIR, "miss_grade_data_activespin.rds"))

# FB anchor (FF>SI>FC) per pitcher-season for release point / extension / active spin
fb_stat <- d[pitch_type %in% c("FF","SI","FC")]
fb_stat[, pr := match(pitch_type, c("FF","SI","FC"))]
fb_stat <- fb_stat[order(pitcher,season,pr),
  .(fb_relx=mean(release_pos_x,na.rm=TRUE), fb_relz=mean(release_pos_z,na.rm=TRUE),
    fb_ext=mean(release_extension,na.rm=TRUE), fb_velo=mean(release_speed,na.rm=TRUE)),
  by=.(pitcher,season, anchor=pitch_type)][, .SD[1], by=.(pitcher,season)]
d <- merge(d, fb_stat[, .(pitcher,season,fb_relx,fb_relz,fb_ext,fb_velo)], by=c("pitcher","season"), all.x=TRUE)
fb_as <- as_long[pitch_type %in% c("FF","SI","FC")][, pr:=match(pitch_type,c("FF","SI","FC"))][
  order(pitcher,season,pr)][, .SD[1], by=.(pitcher,season)][, .(pitcher,season,fb_active=active_spin)]
d <- merge(d, fb_as, by=c("pitcher","season"), all.x=TRUE)

# fastball-anchored diffs (tunneling drivers)
d[, `:=`(
  relx_diff = release_pos_x - fb_relx,          # side-of-release gap
  relz_diff = release_pos_z - fb_relz,          # height-of-release gap
  ext_diff  = release_extension - fb_ext,       # extension gap
  move_diff = sqrt(ax_diff^2 + az_diff^2),      # total movement gap vs FB (accel space)
  meas_eff_diff = active_spin - fb_active        # measured active-spin gap
)]
d[, rel_dist := sqrt(relx_diff^2 + relz_diff^2)] # release-point separation

sub <- d[grp=="breaking" & strikes==2 & is.finite(path_ratio)]
cat(sprintf("breaking + 2-strike pitches: %d (pitch types: %s)\n",
    nrow(sub), paste(sort(unique(sub$pitch_type)), collapse=", ")))

CAND <- c("speed_diff","move_diff","ax_diff","az_diff","axis_diff","spin_eff_diff","meas_eff_diff",
          "relx_diff","relz_diff","rel_dist","ext_diff",
          "release_speed","release_spin_rate","release_extension","az","ax")
lab <- c(speed_diff="velo gap vs FB", move_diff="total movement gap vs FB",
  ax_diff="horiz-accel gap vs FB", az_diff="vert-accel gap vs FB", axis_diff="spin-axis gap vs FB (deg)",
  spin_eff_diff="spin-eff gap vs FB (inferred)", meas_eff_diff="active-spin gap vs FB (measured)",
  relx_diff="release side gap vs FB", relz_diff="release height gap vs FB",
  rel_dist="release-point separation", ext_diff="extension gap vs FB",
  release_speed="velo (abs)", release_spin_rate="spin rate (abs)",
  release_extension="extension (abs)", az="vert accel (abs)", ax="horiz accel (abs)")

corr_tab <- function(x){
  rbindlist(lapply(CAND, function(f){
    v <- x[[f]]; ok <- is.finite(v) & is.finite(x$path_ratio)
    if(sum(ok) < 50) return(NULL)
    data.table(feature=lab[f], pearson=round(cor(v[ok], x$path_ratio[ok]),3),
      spearman=round(cor(v[ok], x$path_ratio[ok], method="spearman"),3), n=sum(ok))
  }))[order(-abs(pearson))]
}

cat("\n=== ALL breaking pitch types pooled (2-strike) : correlation with path_ratio ===\n")
cat("    (path_ratio LOWER = tighter tunnel; negative r means the trait rises as tunnel tightens)\n")
print(corr_tab(sub))

for(pt in c("SL","ST","CU","KC")){
  s <- sub[pitch_type==pt]; if(nrow(s) < 500) next
  cat(sprintf("\n=== %s only (2-strike, n=%d) ===\n", pt, nrow(s)))
  print(head(corr_tab(s), 8))
}

# ---- multivariate: standardized linear model, which jointly explains path_ratio ----
mv <- sub[stats::complete.cases(sub[, c("path_ratio",CAND), with=FALSE])]
z <- as.data.table(scale(mv[, CAND, with=FALSE])); z[, path_ratio := mv$path_ratio]
fit <- lm(path_ratio ~ ., data=z)
co <- summary(fit)$coefficients
mvtab <- data.table(feature=lab[rownames(co)[-1]], std_beta=round(co[-1,1],3),
                    t=round(co[-1,3],1))[order(-abs(std_beta))]
cat(sprintf("\n=== multivariate standardized betas (all breaking, 2-strike; R2=%.3f) ===\n", summary(fit)$r.squared))
cat("    (std_beta = SD change in path_ratio per 1 SD of feature, holding others fixed)\n")
print(mvtab)
fwrite(corr_tab(sub), file.path(MDIR, "path_ratio_correlates_breaking2k.csv"))
cat("\nsaved path_ratio_correlates_breaking2k.csv\n")
