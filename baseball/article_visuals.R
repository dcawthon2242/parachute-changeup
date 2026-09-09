#!/usr/bin/env Rscript

# Build extract CSVs + figures for the pitch-pair deception article.
# Outputs go to data/statcast_model/article_assets/.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(scales) })
set.seed(1)
MDIR <- file.path("data","statcast_model")
AST  <- file.path(MDIR,"article_assets")
dir.create(AST, showWarnings=FALSE, recursive=TRUE)
theme_set(theme_minimal(base_size=13) +
  theme(plot.title=element_text(face="bold"), panel.grid.minor=element_blank(),
        plot.title.position="plot"))
ACC <- "#1b6ca8"; POS <- "#2a9d8f"; NEG <- "#e76f51"; GREY <- "#9aa0a6"
sav <- function(p, f, w=9, h=6) ggsave(file.path(AST,f), p, width=w, height=h, dpi=150, bg="white")

## ---------------------------------------------------------------------------
## Shared: fastball anchors for magnitude gaps
## ---------------------------------------------------------------------------
d <- readRDS(file.path(MDIR,"miss_grade_data_activespin.rds"))
as_long0 <- readRDS(file.path(MDIR,"active_spin_long.rds"))
fbA <- as_long0[pitch_type %in% c("FF","SI","FC")][, pr:=match(pitch_type,c("FF","SI","FC"))][
  order(pitcher,season,pr)][, .SD[1], by=.(pitcher,season)][, .(pitcher,season,fb_active=active_spin)]
d <- merge(d, fbA, by=c("pitcher","season"), all.x=TRUE)
fbR <- d[pitch_type %in% c("FF","SI","FC")][, pr:=match(pitch_type,c("FF","SI","FC"))][
  order(pitcher,season,pr), .(fb_relx=mean(release_pos_x,na.rm=TRUE), fb_relz=mean(release_pos_z,na.rm=TRUE),
     fb_ext=mean(release_extension,na.rm=TRUE)), by=.(pitcher,season)]
d <- merge(d, fbR, by=c("pitcher","season"), all.x=TRUE)

## ---------------------------------------------------------------------------
## FIG 3 - tunneling drivers: how each fastball-relative DIFFERENCE loosens the tunnel
##   (magnitude gaps -> positive r means: bigger difference = looser tunnel)
## ---------------------------------------------------------------------------
b <- d[grp=="breaking" & strikes==2 & is.finite(path_ratio)]
b[, `:=`(
  `Velocity gap`             = abs(speed_diff),
  `Active-spin gap`          = abs(active_spin - fb_active),
  `Movement gap`             = sqrt(ax_diff^2 + az_diff^2),
  `Spin-axis gap`            = axis_diff,
  `Vertical-accel gap`       = abs(az_diff),
  `Extension gap`            = abs(release_extension - fb_ext),
  `Release-point separation` = sqrt((release_pos_x-fb_relx)^2 + (release_pos_z-fb_relz)^2))]
drivers <- c("Velocity gap","Active-spin gap","Movement gap","Spin-axis gap",
             "Vertical-accel gap","Extension gap","Release-point separation")
dd <- rbindlist(lapply(drivers, function(v){
  ok <- is.finite(b[[v]]) & is.finite(b$path_ratio)
  data.table(driver=v, pearson=cor(b[[v]][ok], b$path_ratio[ok])) }))
dd[, kind := ifelse(grepl("Release|Extension", driver), "release / slot", "spin / velo / movement")]
fwrite(dd, file.path(AST,"ext_tunnel_drivers.csv"))
p3 <- ggplot(dd, aes(reorder(driver, pearson), pearson, fill=kind)) +
  geom_col(width=.68) + geom_hline(yintercept=0, color="black", linewidth=.3) +
  geom_text(aes(label=sprintf("%+.2f", pearson), hjust=ifelse(pearson<0,1.15,-0.15)), size=4) +
  coord_flip() + scale_fill_manual(values=c("release / slot"=GREY, "spin / velo / movement"=ACC)) +
  scale_y_continuous(limits=c(-0.72,0.92)) +
  labs(title="What makes a breaking ball tunnel the fastball",
       subtitle="Correlation of each fastball-relative gap with path_ratio  |  breaking balls, 2-strike, n~129k  |  multivariate R2 = 0.68\n+ = larger gap LOOSENS tunnel   |   - = larger gap TIGHTENS it (hard gyro sliders)   |   release/slot ~ 0",
       x=NULL, y="Pearson r of the trait gap with path_ratio", fill=NULL) +
  theme(legend.position="top")
sav(p3, "fig3_tunnel_drivers.png", 9.5, 5.8)

## ---------------------------------------------------------------------------
## FIG 4 - path_ratio vs velocity gap scatter (breaking, 2-strike)
## ---------------------------------------------------------------------------
sc <- d[grp=="breaking" & strikes==2 & is.finite(path_ratio) & is.finite(speed_diff),
        .(velo_gap = abs(speed_diff), path_ratio, pitch_type)]
sc <- sc[pitch_type %in% c("SL","ST","CU","KC") & velo_gap > 0 & path_ratio < quantile(path_ratio,.99)]
rP <- cor(sc$velo_gap, sc$path_ratio)
samp <- sc[sample(.N, min(9000,.N))]
fwrite(samp, file.path(AST,"ext_pathratio_veloscatter.csv"))
p4 <- ggplot(samp, aes(velo_gap, path_ratio)) +
  geom_point(aes(color=pitch_type), alpha=.25, size=.8) +
  geom_smooth(method="lm", se=FALSE, color="black", linewidth=1) +
  scale_color_brewer(palette="Set1") +
  annotate("label", x=Inf, y=Inf, hjust=1.05, vjust=1.4,
           label=sprintf("Pearson r = %.2f\nmultivariate R2 = 0.68", rP), size=4.2) +
  guides(color=guide_legend(override.aes=list(alpha=1, size=2.5))) +
  labs(title="Similar velocity tightens the tunnel",
       subtitle="Breaking balls, 2-strike counts  |  bigger velocity gap = looser tunnel (higher path_ratio)",
       x="Velocity gap vs fastball (mph)", y="path_ratio (lower = tighter tunnel)", color="Pitch") +
  theme(legend.position="top")
sav(p4, "fig4_pathratio_velo_scatter.png", 8.5, 6)

## ---------------------------------------------------------------------------
## FIG 5 - RMSE lift: only the type-agnostic model benefits
## ---------------------------------------------------------------------------
g <- fread(file.path(MDIR,"grid_rmse_summary.csv"))
g[, set := factor(set, levels=c("ALL","2K"), labels=c("All swings","2-strike"))]
g[, approach := factor(approach, levels=c("A all-types","C grouped FB/BRK/OFF","B per-pitch-type"),
     labels=c("All-types\n(type-agnostic)","Grouped\nFB/BRK/OFF","Per-pitch-type"))]
g[, improve := ifelse(delta_rmse<0,"improves (RMSE down)","no gain / worse")]
fwrite(g[, .(set, approach=gsub("\n"," ",approach), base_rmse, aug_rmse, delta_rmse)],
       file.path(AST,"ext_rmse_grid.csv"))
p5 <- ggplot(g, aes(approach, delta_rmse, fill=improve)) +
  geom_col(width=.6) + geom_hline(yintercept=0, color="black") +
  geom_text(aes(label=sprintf("%+.4f", delta_rmse),
                vjust=ifelse(delta_rmse<0, 1.4, -0.6)), size=3.6) +
  facet_wrap(~set) + scale_fill_manual(values=c("improves (RMSE down)"=POS,"no gain / worse"=NEG)) +
  labs(title="Pitch-pair features help only the model that can't already see pitch type",
       subtitle="Change in holdout RMSE from adding tunneling + spin-similarity (negative = better)",
       x=NULL, y="RMSE change (aug - base)", fill=NULL) +
  theme(legend.position="top")
sav(p5, "fig5_rmse_lift.png", 9, 5.6)

## ---------------------------------------------------------------------------
## FIG 6 - top-10 breaking risers (2-strike)
## ---------------------------------------------------------------------------
rs <- fread(file.path(MDIR,"risers_2K.csv"))[grp=="breaking"][order(-gain)][1:10]
rs[, label := sprintf("%s  %s", player_name, pitch_type)]
fwrite(rs[, .(player_name,pitch_type,n,act,base,aug,gain,path_ratio)],
       file.path(AST,"ext_risers_breaking_2k.csv"))
p6 <- ggplot(rs, aes(reorder(label, gain), gain, fill=path_ratio)) +
  geom_col(width=.7) +
  geom_text(aes(label=sprintf("path_ratio %.2f", path_ratio)), hjust=-0.08, size=3.4) +
  coord_flip() + scale_fill_gradient(low=ACC, high="#cfe3f2", name="path_ratio\n(low=tight)") +
  scale_y_continuous(expand=expansion(mult=c(0,0.28))) +
  labs(title="Breaking balls the shape-only model under-rates",
       subtitle="Predicted-miss gain when tunneling is added (2-strike counts). Tight tunnels (low path_ratio).",
       x=NULL, y="Predicted miss-distance gain (inches)")
sav(p6, "fig6_risers_breaking.png", 9, 6)

## ---------------------------------------------------------------------------
## FIG 7 - active-spin validation: inferred vs measured, by pitch group
## ---------------------------------------------------------------------------
as_long <- readRDS(file.path(MDIR,"active_spin_long.rds"))
av <- d[!is.na(active_spin) & !is.na(spin_eff),
        .(inferred=mean(spin_eff), measured=mean(active_spin), n=.N),
        by=.(pitcher,season,pitch_type,grp)][n>=25]
av <- av[grp %in% c("fastball","breaking","offspeed")]
av[, grp := factor(grp, levels=c("fastball","breaking","offspeed"),
     labels=c("Fastballs","Breaking","Offspeed"))]
av[, `:=`(zi=scale(inferred), zm=scale(measured)), by=grp]
rtab <- av[, .(r=cor(inferred,measured)), by=grp]
fwrite(av[, .(pitcher,season,pitch_type,grp,inferred,measured,n)], file.path(AST,"ext_activespin_validation.csv"))
fwrite(rtab, file.path(AST,"ext_activespin_group_r.csv"))
p7 <- ggplot(av, aes(zm, zi)) +
  geom_point(alpha=.18, size=.7, color=ACC) +
  geom_smooth(method="lm", se=FALSE, color="black") +
  facet_wrap(~grp) +
  geom_text(data=rtab, aes(x=-2.3, y=2.6, label=sprintf("r = %.2f", r)), size=4.4, hjust=0, fontface="bold") +
  labs(title="Validating spin efficiency: inferred vs Savant measured active spin",
       subtitle="Per pitcher x pitch-type x season, standardized within group. Weakest for offspeed -> that feature uses measured spin.",
       x="Measured active spin (z)", y="Inferred spin efficiency (z)")
sav(p7, "fig7_activespin_validation.png", 11, 4.6)

## ---------------------------------------------------------------------------
## FIG 8 - spin-similarity distribution: residual overperformers vs field
## ---------------------------------------------------------------------------
rr <- fread(file.path(MDIR,"remaining_residual_offspeed.csv"))
rr <- rr[!is.na(spin_sim)]
thr <- quantile(rr$res_loc, .85)
rr[, grp2 := ifelse(res_loc>=thr, "Residual overperformers\n(top 15% after shape+location)", "All other offspeed")]
mns <- rr[, .(m=mean(spin_sim)), by=grp2]
fwrite(rr[, .(player_name,pitch_type,n,res_loc,spin_sim,grp2)], file.path(AST,"ext_spin_sim_overperf_vs_field.csv"))
p8 <- ggplot(rr, aes(spin_sim, fill=grp2)) +
  geom_density(alpha=.5, color=NA) +
  geom_vline(data=mns, aes(xintercept=m, color=grp2), linewidth=1, linetype="dashed", show.legend=FALSE) +
  scale_fill_manual(values=setNames(c(POS,GREY), mns$grp2)) +
  scale_color_manual(values=setNames(c(POS,GREY), mns$grp2)) +
  annotate("text", x=mns$m[mns$grp2==unique(rr$grp2)[1]], y=Inf, vjust=2,
           label="", size=3) +
  labs(title="Overperformers skew spin-similar to the fastball",
       subtitle=sprintf("Offspeed pitches. Dashed = group mean spin similarity (%.2f vs %.2f).",
                        max(mns$m), min(mns$m)),
       x="Spin similarity to fastball (1 = identical spin)", y="Density", fill=NULL) +
  theme(legend.position="top")
sav(p8, "fig8_spinsim_distribution.png", 9, 5.6)

## ---------------------------------------------------------------------------
## FIG 9 - overperformance by spin-similarity BIN, with 3 handpicked exemplars/bin
##   Honest full distribution (violin + all points) per bin; labeled exemplars
##   chosen to illustrate the phenomenon (dissimilar -> hit; matched -> whiff).
## ---------------------------------------------------------------------------
## 3 bins, CHANGEUPS ONLY (the population where a simple value is monotonic)
c9 <- rr[pitch_type=="CH"]
BR <- c(-Inf,0.45,0.70,Inf); BL <- c("<0.45","0.45-0.70",">=0.70")
c9[, bin := cut(spin_sim, breaks=BR, labels=BL)]
c9[, xnum := as.integer(bin)]
binstat <- c9[, .(n=.N, mean_res=mean(res_loc), share_over=mean(res_loc>0)),
              by=.(bin,xnum)][order(xnum)]
fwrite(binstat, file.path(AST,"ext_spin_sim_bins.csv"))

# curated exemplars (user-selected changeups): under -> neutral -> over across the 3 bins
exkey <- data.table(
  player_name=c("Eflin, Zach","Schultz, Paxton","Kowar, Jackson",   # <0.45  underperformers
                "Suarez, Ranger","Gil, Luis","Cortes, Nestor",       # 0.45-0.70 representative/neutral
                "Cease, Dylan","Ribalta, Orlando","Garcia, Rico"),   # >=0.70 overperformers
  pitch_type=rep("CH", 9))
ex <- merge(exkey, c9, by=c("player_name","pitch_type"), all.x=TRUE)
ex[, rk := frank(res_loc, ties.method="first"), by=bin]
ex[, xpos := xnum + c(-0.27,0,0.27)[rk]]
ex[, dir := ifelse(res_loc>0.03,"over",ifelse(res_loc< -0.03,"under","neutral"))]
ex[, last := trimws(tstrsplit(player_name, ",")[[1]])]
fwrite(ex[order(xnum,-res_loc), .(bin,player_name,pitch_type,n,whiff,res_loc,spin_sim,dir)],
       file.path(AST,"ext_spin_sim_exemplars.csv"))

p9 <- ggplot(c9, aes(xnum, res_loc)) +
  geom_hline(yintercept=0, color="black", linewidth=.4) +
  geom_violin(aes(group=xnum), fill="#e9ecef", color=NA, scale="width", width=.9) +
  geom_jitter(width=.09, height=0, alpha=.15, size=.7, color=GREY) +
  geom_point(data=binstat, aes(xnum, mean_res), shape=18, size=4, color="black") +
  geom_point(data=ex, aes(xpos, res_loc, fill=dir), shape=21, size=3.1, color="black") +
  geom_text(data=ex, aes(xpos, res_loc, label=last,
              vjust=ifelse(res_loc>=0,-0.7,1.5)), size=3) +
  # top: simple monotonic value = share of changeups that beat their expectation (fixed black labels)
  geom_text(data=binstat, aes(xnum, 0.30, label=sprintf("%.0f%%", 100*share_over)),
            size=4.2, fontface="bold", color="black") +
  annotate("text", x=2, y=0.345, size=3.2, fontface="italic", color="black",
           label="share of changeups that beat their shape+location expectation") +
  scale_fill_manual(values=c(over=POS, neutral=GREY, under=NEG),
     labels=c(over="overperforms", neutral="~neutral", under="underperforms"), name="exemplar") +
  scale_x_continuous(breaks=1:3, labels=BL, limits=c(0.5,3.5)) +
  coord_cartesian(ylim=c(-0.22,0.35)) +
  labs(title="More fastball-like spin skews toward overperformance",
       subtitle="ChangeUps binned by spin similarity to primary fastball. Grey = Individual Pitcher + PitchType\nAs Spin Similarity Rises, Miss Distance Overperformance Improves",
       x="Spin similarity to fastball (binned)",
       y="Overperformance: whiffs above shape+location expectation") +
  theme(legend.position="top")
sav(p9, "fig9_spinsim_bins.png", 10, 6.4)

## ---------------------------------------------------------------------------
## FIG 10 - Dylan Cease case study
## ---------------------------------------------------------------------------
lgCH <- rr[pitch_type=="CH", .(whiff=mean(whiff), spin_sim=mean(spin_sim,na.rm=TRUE))]
ce <- rr[grepl("Cease", player_name) & pitch_type=="CH"][1]
ce_exp <- ce$whiff - ce$res_loc
cc <- data.table(
  metric=c("Cease actual whiff%","Cease expected\n(shape+location)","League avg CH whiff%"),
  value=c(ce$whiff, ce_exp, lgCH$whiff))
cc[, metric := factor(metric, levels=metric)]
fwrite(data.table(who=c("Cease","Cease_expected","League_CH"),
                  whiff=c(ce$whiff, ce_exp, lgCH$whiff),
                  spin_sim=c(ce$spin_sim, NA, lgCH$spin_sim)), file.path(AST,"ext_cease_case.csv"))
p10 <- ggplot(cc, aes(metric, value, fill=metric)) +
  geom_col(width=.62) +
  geom_text(aes(label=percent(value, accuracy=0.1)), vjust=-0.5, size=4.2) +
  scale_fill_manual(values=c(POS, GREY, "#c7ccd1"), guide="none") +
  scale_y_continuous(labels=percent, expand=expansion(mult=c(0,0.15))) +
  labs(title=sprintf("Case study: Dylan Cease changeup (spin similarity %.2f)", ce$spin_sim),
       subtitle=sprintf("Whiffs %.0f%%, ~%.0f pts above its shape+location expectation - a near-identical spin look to the four-seam.",
                        100*ce$whiff, 100*(ce$whiff-ce_exp)),
       x=NULL, y="Whiff rate")
sav(p10, "fig10_cease_case.png", 9, 5.6)

## ---------------------------------------------------------------------------
## FIG 2 - tunneling schematic (synthetic trajectories)
## ---------------------------------------------------------------------------
xx <- seq(0, 1, length.out=200)                 # 0 = release, 1 = plate
fb <- 3.0 - 1.6*xx^2                             # four-seam: gentle drop
ch <- 3.0 - 1.6*xx^2 - 2.4*xx^3                  # changeup: matches early, drops late
sch <- rbind(data.table(x=xx, y=fb, pitch="Four-seam"),
             data.table(x=xx, y=ch, pitch="Changeup"))
tp <- 0.55
p2 <- ggplot(sch, aes(x, y, color=pitch)) +
  geom_line(linewidth=1.4) +
  annotate("segment", x=tp, xend=tp, y=3.0-1.6*tp^2, yend=3.0-1.6*tp^2-2.4*tp^3,
           linetype="dotted") +
  annotate("segment", x=1, xend=1, y=min(ch), yend=fb[200], linetype="dotted") +
  annotate("text", x=tp, y=2.85, label="look identical here\n(the tunnel)", size=3.6, hjust=0.5) +
  annotate("text", x=0.985, y=(min(ch)+fb[200])/2, label="late\nseparation", size=3.6, hjust=1.1) +
  scale_color_manual(values=c("Four-seam"=ACC, "Changeup"=NEG)) +
  scale_x_continuous(breaks=c(0,1), labels=c("release","plate")) +
  labs(title="Tunneling: same look early, separation late",
       subtitle="path_ratio = early-flight distance between pitches / their separation at the plate (low = tighter tunnel)",
       x=NULL, y="Height (illustrative)", color=NULL) +
  theme(legend.position="top", axis.text.y=element_blank(), panel.grid=element_blank())
sav(p2, "fig2_tunnel_schematic.png", 8.5, 5.4)

cat("Wrote figures + extracts to", AST, "\n")
print(list.files(AST))
