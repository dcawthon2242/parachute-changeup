#!/usr/bin/env Rscript

# Arm slot as a GATE versus as a DESCRIPTION.
#
# A gate changes who is in the sample, so it changes every number computed downstream. A
# descriptive statistic leaves the sample alone and just characterises it. The question here is
# whether the arm-slot condition is earning its place as a gate: does requiring a high slot
# select pitchers the other cuts would have missed, or does it only shrink an already small bin?

suppressPackageStartupMessages({ library(data.table) })
options(width = 200)
MDIR <- "data/statcast_model"

R  <- readRDS(file.path(MDIR, "whiff_tjstuff.rds"))
AS <- readRDS(file.path(MDIR, "active_spin_long.rds"))
S <- R[, .(nsw = .N, axis = mean(axis_diff), arm = mean(arm_angle), velo_sep = -mean(speed_diff),
           w = 100*mean(r4)), by = .(pitcher, player_name, season)]
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, as_ch = active_spin)], by = c("pitcher","season"))
S <- merge(S, AS[pitch_type == "FF", .(pitcher, season, as_fb = active_spin)], by = c("pitcher","season"))
S <- S[nsw >= 60]
S[, name := trimws(paste(sub(".*,\\s*","",player_name), sub(",.*","",player_name)))]
mu <- S[, .(as_fb = mean(as_fb), as_ch = mean(as_ch), arm = mean(arm))]

# The bin without any arm condition at all: spin floors plus a tight axis.
S[, core := as_fb > mu$as_fb & as_ch > mu$as_ch & axis <= 10]
cat(sprintf("league arm slot: mean %.1f, median %.1f, sd %.1f\n", mu$arm, median(S$arm), sd(S$arm)))
cat(sprintf("\nbin with NO arm condition (spin floors + axis <= 10): %d seasons\n", sum(S$core)))

C <- S[core == TRUE]
cat(sprintf("  their arm slots: mean %.1f, median %.1f, range %.1f-%.1f\n",
            mean(C$arm), median(C$arm), min(C$arm), max(C$arm)))
cat(sprintf("  share above the league mean slot: %.0f%% (league baseline is 50%% by construction)\n",
            100*mean(C$arm > mu$arm)))
cat(sprintf("  percentile of the league slot distribution the bin's median sits at: %.0fth\n",
            100*mean(S$arm < median(C$arm))))

cat("\n=== what adding the arm gate costs and buys ===\n")
print(rbindlist(lapply(list(c("none", "TRUE"), c("arm > 37.9 mean","arm > mu$arm"),
                            c("arm >= 44","arm >= 44")), function(g) {
  i <- S$core & S[, eval(parse(text = g[2]))]; t <- t.test(S$w[i], S$w[!i])
  data.table(arm_gate = g[1], seasons = sum(i), whiff = round(diff(rev(t$estimate)),2),
             ci = sprintf("[%+.2f, %+.2f]", t$conf.int[1], t$conf.int[2]),
             ci_width = round(t$conf.int[2]-t$conf.int[1],2), p = round(t$p.value,4)) })),
  row.names = FALSE)

cat("\n=== who the no-arm-gate bin actually finds, sorted by slot ===\n")
print(C[order(-arm), .(Pitcher = name, Season = season, Swings = nsw, ArmSlot = round(arm,1),
        AxisGap = round(axis,1), ActFF = round(as_fb,2), ActCH = round(as_ch,2),
        VeloSep = round(velo_sep,1), WhiffAbove = round(w,1))], row.names = FALSE)
