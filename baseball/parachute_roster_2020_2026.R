#!/usr/bin/env Rscript

# EVERY PARACHUTE CHANGEUP, 2020-2026.
#
# High slot plus spin that matches the primary fastball, over the full Hawk-Eye era. Three
# tiers are reported rather than one, because the strict definition leaves 26 pitcher-seasons
# and the honest reading of the replication test is that 26 is not many. Showing the tiers
# side by side makes it visible how quickly the group grows as the thresholds loosen and
# whether the performance holds up as it does.
#
#   Core    axis gap <= 10 deg, active spin within 10 points, arm angle >= 44 deg
#   Wide    axis gap <= 15 deg, active spin within 15 points, arm angle >= 42 deg
#   Elite   axis gap <=  6 deg, active spin within  6 points, arm angle >= 48 deg
#
# All performance columns are above a shape-and-location model that never sees the spin axis
# or the arm angle, fit out of fold across all 464,710 changeups in the window.

suppressPackageStartupMessages({ library(data.table) })
options(width = 215)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")
# Four-seam anchor. The old primary-fastball rule fell back to a sinker when the pitcher had no
# four-seamer, and because a sinker's seam-shifted wake decouples its spin axis from its
# movement, a matched axis against a sinker is not the same measurement. It also over-recruited
# sinkerballers into the bin and carried the ground-ball effect with them.
ANCHOR <- { a <- commandArgs(TRUE); if (length(a)) a[1] else "FF" }
S <- fread(file.path(AST, if (ANCHOR == "FF") "ext_parachute_ff.csv" else "ext_parachute_extended.csv"))
S <- S[is.finite(as_gap) & is.finite(arm)]
cat(sprintf("anchor: %s\n\n", ANCHOR))

tier <- function(ax, ef, ar) S$axis <= ax & abs(S$as_gap) <= ef & S$arm >= ar
S[, tier := fifelse(tier(6, .06, 48), "Elite",
            fifelse(tier(10, .10, 44), "Core",
            fifelse(tier(15, .15, 42), "Wide", NA_character_)))]
S[, in_wide := !is.na(tier)]

cat(sprintf("population: %d pitcher-seasons, 2020-2026, min 60 changeup swings,\nwith measured active spin on both the changeup and the primary fastball\n\n", nrow(S)))
cat("=== tier sizes and performance vs everyone outside the tier ===\n")
sm <- rbindlist(lapply(list(c(6,.06,48,"Elite"), c(10,.10,44,"Core"), c(15,.15,42,"Wide")), function(z) {
  i <- tier(as.numeric(z[1]), as.numeric(z[2]), as.numeric(z[3]))
  f <- function(v) { t <- t.test(S[[v]][i], S[[v]][!i])
                     sprintf("%+.2f (p=%.3f)", diff(rev(t$estimate)), t$p.value) }
  data.table(tier = z[4], seasons = sum(i), pitchers = uniqueN(S$pitcher[i]),
             gb_over = f("gb"), whiff_over = f("wh"), rv100 = f("rv100")) }))
print(sm, row.names = FALSE)

# Last name alone is ambiguous here: Tyler and Chase Anderson both qualify in 2021 from very
# different slots, and the same is true of the Perezes and Cabreras. Full name throughout.
R <- S[in_wide == TRUE][order(-wh)][, .(
  Pitcher = trimws(paste(sub(".*,\\s*","",player_name), sub(",.*","",player_name))),
  Season = season, Tier = tier, Pitches = np,
  Arm = round(arm,1), AxisGap = round(axis,1), SpinEffGap = round(as_gap,3),
  VeloSep = round(velo_sep,1), IVBKill = round(kill,1),
  Whiff = round(wh,1), Grounders = round(gb,1), RV100 = round(rv100,2))]
cat(sprintf("\n=== all %d qualifying pitcher-seasons, sorted by whiff above model ===\n", nrow(R)))
print(R, row.names = FALSE)

cat("\n=== pitchers appearing more than once ===\n")
rep <- R[, .(seasons = .N, yrs = paste(sort(Season), collapse=","),
             whiff = round(mean(Whiff),1), gb = round(mean(Grounders),1),
             rv = round(mean(RV100),2)), by = Pitcher][seasons > 1][order(-seasons, -whiff)]
print(rep, row.names = FALSE)

OUT <- sprintf("ext_parachute_roster_%s.csv", tolower(ANCHOR))
fwrite(R, file.path(AST, OUT))
cat(sprintf("\nwrote %s (%d rows)\n", OUT, nrow(R)))
