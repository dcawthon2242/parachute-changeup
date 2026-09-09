#!/usr/bin/env Rscript

# Test the matched-spin ("invisiball/parachute") hypothesis two ways the miss-distance
# model missed:
#   (1) WHIFF-RATE residual over a shape-only expectation (the original proof method)
#   (2) CONDITIONAL "matched spin BUT drops": matched active spin/axis vs FB AND a big
#       induced-vertical-break (IVB) drop vs the FB.
# Shape model trained 2023-25, evaluated on 2026 offspeed. Locate Cease/Ragans/Vesia/Chivilli.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(1)
MDIR <- file.path("data","statcast_model")
FEAT <- readRDS(file.path(MDIR, "arch_compare.rds"))$FEAT
as_long <- readRDS(file.path(MDIR, "active_spin_long.rds"))
d <- readRDS(file.path(MDIR, "miss_grade_data_activespin.rds"))
fb <- as_long[pitch_type %in% c("FF","SI","FC")][, pr:=match(pitch_type,c("FF","SI","FC"))][
  order(pitcher,season,pr)][, .SD[1], by=.(pitcher,season)][, .(pitcher,season,fb_active=active_spin)]
d <- merge(d, fb, by=c("pitcher","season"), all.x=TRUE)
d[, as_gap := active_spin - fb_active]; d[, abs_as := abs(as_gap)]
d[, spin_sim := exp(-(abs_as/0.10)^2) * exp(-(axis_diff/45)^2)]

# ---- merge IVB (pfx_z, inches) for 2026 and build ivb_gap vs primary FB ----
raw26 <- fread("data/statcast_2026/statcast_2026_all.csv",
               select=c("game_pk","at_bat_number","pitch_number","pfx_z"))
raw26[, ivb := pfx_z*12]
d <- merge(d, raw26[, .(game_pk,at_bat_number,pitch_number,ivb)],
           by=c("game_pk","at_bat_number","pitch_number"), all.x=TRUE)
fbivb <- d[season==2026 & pitch_type %in% c("FF","SI","FC")][, pr:=match(pitch_type,c("FF","SI","FC"))][
  order(pitcher,pr), .(fb_ivb=mean(ivb,na.rm=TRUE)), by=pitcher]
d <- merge(d, fbivb, by="pitcher", all.x=TRUE)
d[, ivb_gap := ivb - fb_ivb]                       # negative = drops more than the FB

# ---- shape-only expectations (train 2023-25) ----
tr <- d[set=="train" & is.finite(miss_distance)]
mk <- function(obj, y){ n<-nrow(tr); vi<-sample(n,floor(0.15*n))
  dtrain<-lgb.Dataset(as.matrix(tr[-vi,..FEAT]),label=tr[[y]][-vi])
  dval<-lgb.Dataset.create.valid(dtrain,as.matrix(tr[vi,..FEAT]),label=tr[[y]][vi])
  lgb.train(params=list(objective=obj,metric=if(obj=="binary")"binary_logloss" else "rmse",
    learning_rate=0.05,num_leaves=31,min_data_in_leaf=150,feature_fraction=0.8,
    bagging_fraction=0.8,bagging_freq=1),data=dtrain,nrounds=2000,valids=list(val=dval),
    early_stopping_rounds=60,verbose=-1) }
mw <- mk("binary","is_whiff")

ho <- d[set=="holdout" & grp=="offspeed" & !is.na(is_whiff)]
ho[, exp_whiff := predict(mw, as.matrix(.SD[, FEAT, with=FALSE]))]
ho[, wresid := is_whiff - exp_whiff]               # + = whiffs more than shape predicts

ag <- ho[, .(n=.N, whiff=round(mean(is_whiff),3), exp_whiff=round(mean(exp_whiff),3),
  wresid=round(mean(wresid),3), axis=round(mean(axis_diff,na.rm=TRUE),1),
  as_gap=round(mean(as_gap,na.rm=TRUE),3), ivb_gap=round(mean(ivb_gap,na.rm=TRUE),1),
  spin_sim=round(mean(spin_sim,na.rm=TRUE),3)), by=.(player_name,pitch_type)][n>=40]
CO <- function(a,b,...) cor(a,b,use="complete.obs",...)

cat("=== (1) WHIFF overperformance vs spin similarity (2026 offspeed, n>=40) ===\n")
cat(sprintf("  groups=%d\n", nrow(ag)))
cat(sprintf("  whiff_resid vs spin_sim            : r=%+.3f (Spearman %+.3f)\n", CO(ag$wresid,ag$spin_sim), CO(ag$wresid,ag$spin_sim,method="spearman")))
cat(sprintf("  whiff_resid vs |active-spin gap|   : r=%+.3f\n", CO(ag$wresid, abs(ag$as_gap))))
cat(sprintf("  whiff_resid vs axis gap            : r=%+.3f\n", CO(ag$wresid, ag$axis)))
cat("\n  whiff_resid by active-spin similarity tier:\n")
print(ag[!is.na(as_gap), .(n=.N, whiff_resid=round(mean(wresid),3)),
  by=.(tier=cut(abs(as_gap),c(-.01,.03,.07,.15,1),labels=c("<=.03","0.03-.07","0.07-.15",">.15")))][order(tier)])

cat("\n=== (2) CONDITIONAL: matched spin AND drops (spin_sim>=0.5 & ivb_gap<=-8) ===\n")
ho[, cohort := fifelse(spin_sim>=0.5 & ivb_gap<=-8, "matched+drops",
              fifelse(spin_sim>=0.5, "matched, less drop",
              fifelse(ivb_gap<=-8, "drops, unmatched spin", "neither")))]
print(ho[!is.na(cohort), .(pitches=.N, whiff=round(mean(is_whiff),3), exp=round(mean(exp_whiff),3),
  whiff_resid=round(mean(wresid),3)), by=cohort][order(-whiff_resid)])

cat("\n=== ARCHETYPES: Cease / Ragans / Vesia / Chivilli ===\n")
print(ag[grepl("Cease|Ragans|Vesia|Chivilli", player_name)][order(-wresid),
  .(player_name,pitch_type,n,whiff,exp_whiff,wresid,axis,as_gap,ivb_gap,spin_sim)])

cat("\n=== TOP 15 WHIFF OVERPERFORMERS among matched spin (spin_sim>=0.4) ===\n")
print(head(ag[spin_sim>=0.4][order(-wresid),
  .(player_name,pitch_type,n,whiff,exp_whiff,wresid,axis,as_gap,ivb_gap,spin_sim)],15))
fwrite(ag[order(-wresid)], file.path(MDIR, "whiff_matched_spin.csv"))
cat("\nsaved whiff_matched_spin.csv\n")
