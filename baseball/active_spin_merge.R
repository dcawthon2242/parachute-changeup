#!/usr/bin/env Rscript

# Merge Savant's MEASURED active-spin% (pitcher x pitch-type x season) onto the
# miss-grade data and cross-check it against my per-pitch INFERRED spin_eff
# (Magnus decomposition). Then rebuild the offspeed FB spin-similarity feature
# with measured values and see if it tracks miss distance better.

suppressPackageStartupMessages({ library(data.table) })
MDIR <- file.path("data","statcast_model")
ASD  <- file.path("data","active_spin")

# ---- 1) load + reshape active spin leaderboards to long ----
colmap <- c(active_spin_fourseam="FF", active_spin_sinker="SI", active_spin_cutter="FC",
            active_spin_changeup="CH", active_spin_splitter="FS", active_spin_curve="CU",
            active_spin_slider="SL", active_spin_sweeper="ST", active_spin_slurve="SV")
load_yr <- function(y){
  x <- fread(file.path(ASD, sprintf("active_spin_%d.csv", y)))
  setnames(x, "entity_id", "pitcher")
  m <- melt(x, id.vars=c("pitcher","entity_name","pitch_hand"),
            measure.vars=names(colmap), variable.name="col", value.name="active_spin")
  m[, pitch_type := colmap[as.character(col)]]
  m[, season := y]; m[!is.na(active_spin), .(pitcher, season, pitch_type, active_spin=as.numeric(active_spin)/100)]
}
as_long <- rbindlist(lapply(2023:2026, load_yr))
# Savant "curve" covers KC/CS too; replicate curve value onto those codes
kc <- as_long[pitch_type=="CU"][, .(pitcher, season, active_spin, pitch_type=list(c("KC","CS")))]
kc <- kc[, .(pitch_type=unlist(pitch_type)), by=.(pitcher, season, active_spin)]
as_long <- unique(rbindlist(list(as_long, kc), use.names=TRUE), by=c("pitcher","season","pitch_type"))
cat(sprintf("active-spin rows (pitcher x season x pitch): %d\n", nrow(as_long)))

# ---- 2) merge onto model data & compare inferred vs measured ----
d <- readRDS(file.path(MDIR, "miss_grade_data.rds"))
d <- merge(d, as_long, by=c("pitcher","season","pitch_type"), all.x=TRUE)
cat(sprintf("pitches with measured active_spin: %d / %d (%.1f%%)\n",
    sum(!is.na(d$active_spin)), nrow(d), 100*mean(!is.na(d$active_spin))))

agg <- d[!is.na(active_spin) & !is.na(spin_eff),
         .(inferred=mean(spin_eff), measured=mean(active_spin), n=.N),
         by=.(pitcher, season, pitch_type, grp)][n>=25]
cat(sprintf("\n=== inferred spin_eff  vs  measured active_spin  (pitcher x pt x season, n>=25: %d) ===\n", nrow(agg)))
cat(sprintf("overall  Pearson r=%.3f  Spearman=%.3f  mean(inferred)=%.3f mean(measured)=%.3f  MAE=%.3f\n",
    cor(agg$inferred, agg$measured), cor(agg$inferred, agg$measured, method="spearman"),
    mean(agg$inferred), mean(agg$measured), mean(abs(agg$inferred-agg$measured))))
cat("\nby pitch-type group:\n")
print(agg[, .(n=.N, r=round(cor(inferred,measured),3),
    mean_inferred=round(mean(inferred),3), mean_measured=round(mean(measured),3),
    bias=round(mean(inferred-measured),3)), by=grp])
cat("\nby pitch type:\n")
print(agg[, .(n=.N, r=round(cor(inferred,measured),3),
    mean_inferred=round(mean(inferred),3), mean_measured=round(mean(measured),3)), by=pitch_type][order(-n)])

saveRDS(as_long, file.path(MDIR, "active_spin_long.rds"))
saveRDS(d,       file.path(MDIR, "miss_grade_data_activespin.rds"))
cat("\nsaved active_spin_long.rds + miss_grade_data_activespin.rds\n")
