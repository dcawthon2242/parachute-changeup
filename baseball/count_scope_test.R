#!/usr/bin/env Rscript

# DOES THE CUE DEPEND ON THE COUNT? Testing the split claim on its own terms.
#
# The claim: offspeed carries a spin-similarity element in ALL counts, while breaking balls
# carry a tunneling element specifically in 2-strike counts, and neither is priced by a
# conventional shape model.
#
# Figure 14 could not test that, because it scored everything on put-away counts after the
# primary fastball. That scope is right for a tunnel - a tunnel is a property of a sequence,
# and the sequence only matters when the hitter has to protect - but it is wrong for a spin
# cue. Spin similarity compares a pitch to the pitcher's fastball in general, not to the
# pitch that happened to precede it, so restricting to 2-strike-after-fastball throws away
# three quarters of the evidence for it and conditions on a variable the cue does not use.
#
# So each cue is now tested on the population it actually claims:
#   sequence cues (path ratio)  ->  pitches after the primary fastball
#   arsenal cues (spin, shape)  ->  every competitive swing on that pitch type
# and both are cut by count regime, which is the thing under test.
#
# Currency stays held-out 2026 RMSE gain on miss distance over an intercept, trained on
# 2023H2-2025, so nothing here depends on a residual definition.

suppressPackageStartupMessages({ library(data.table); library(mgcv); library(ggplot2) })
options(width = 210)
MDIR <- file.path("data","statcast_model"); AST <- file.path(MDIR, "article_assets")
FASTBALLS <- c("FF","SI","FC"); yf <- 17/12
WX <- 1.2; WY <- 0.4; WZ <- 1.4; REACT <- 0.150; NSTEP <- 40
CACHE <- file.path(MDIR, "count_scope_pairs.rds")

if (!file.exists(CACHE) || nzchar(Sys.getenv("REFIT"))) {
kin <- rbindlist(lapply(2023:2026, function(y) {
  k <- fread(sprintf("data/statcast_%d/statcast_%d_all.csv", y, y),
        select = c("game_pk","at_bat_number","pitch_number","pitch_type","game_date","game_type",
                   "pitcher","p_throws","balls","strikes","plate_x","plate_z","release_pos_x",
                   "release_pos_y","release_pos_z","vx0","vy0","vz0","ax","ay","az"),
        showProgress = FALSE); k[, season := y]; k }))
kin <- kin[game_type == "R" & !is.na(vx0) & !is.na(release_pos_y) & pitch_type != ""]
kin[season == 2023L, game_date := as.Date(game_date)]
kin <- kin[season != 2023L | game_date >= as.Date("2023-07-14")]
kin <- unique(kin, by = c("game_pk","at_bat_number","pitch_number"))
setorder(kin, game_pk, at_bat_number, pitch_number)

prim <- kin[pitch_type %in% FASTBALLS, .N, by = .(pitcher, season, pitch_type)][N >= 50]
prim[, rk := match(pitch_type, FASTBALLS)]
prim <- prim[order(pitcher, season, rk)][, .SD[1], by = .(pitcher, season)][
  , .(pitcher, season, primary_fb = pitch_type)]

kin[, t_plate := (-vy0 - sqrt(vy0^2 - 2*ay*(release_pos_y - yf)))/ay]
kin[, t_react := pmax(t_plate - REACT, 0.05)]
lagc <- c("release_pos_x","release_pos_y","release_pos_z","vx0","vy0","vz0","ax","ay","az",
          "t_react","pitch_type","pitch_number","plate_x","plate_z")
for (cc in lagc) kin[, (paste0("p_",cc)) := shift(get(cc)), by = .(game_pk, at_bat_number)]
kin <- merge(kin, prim, by = c("pitcher","season"), all.x = TRUE)
kin[, after_primary := !is.na(p_pitch_number) & (pitch_number - p_pitch_number == 1L) &
                       !is.na(primary_fb) & p_pitch_type == primary_fb]

TY <- c(SL="SL", ST="ST", CU="CU", KC="CU", CH="CH", FS="FS")
# NOTE: no count filter and no sequence filter here. path_ratio is simply NA where there was
# no primary fastball in front, and each driver below uses its own complete cases, so the
# arsenal cues keep the full sample while the tunnel keeps only real sequences.
d <- kin[pitch_type %in% names(TY)]
d[, ptype := TY[pitch_type]]

pos <- function(r,v,a,t) r + v*t + 0.5*a*t^2
d[, pair_path := NA_real_]
idx <- which(d$after_primary)
Tmax <- pmax(d$t_react[idx], d$p_t_react[idx]); acc <- numeric(length(idx))
for (k in 1:NSTEP) {
  tk <- (k-0.5)/NSTEP * Tmax
  dx <- pos(d$release_pos_x[idx],d$vx0[idx],d$ax[idx],tk) - pos(d$p_release_pos_x[idx],d$p_vx0[idx],d$p_ax[idx],tk)
  dy <- pos(d$release_pos_y[idx],d$vy0[idx],d$ay[idx],tk) - pos(d$p_release_pos_y[idx],d$p_vy0[idx],d$p_ay[idx],tk)
  dz <- pos(d$release_pos_z[idx],d$vz0[idx],d$az[idx],tk) - pos(d$p_release_pos_z[idx],d$p_vz0[idx],d$p_az[idx],tk)
  acc <- acc + sqrt(WX*dx^2 + WY*dy^2 + WZ*dz^2) * (Tmax/NSTEP)
}
d$pair_path[idx] <- acc / pmax(sqrt((d$plate_x[idx]-d$p_plate_x[idx])^2 +
                                    (d$plate_z[idx]-d$p_plate_z[idx])^2), 0.1)

mir <- function(x, h) fifelse(h == "R", x, -x)
d[, `:=`(bb_x = mir(plate_x, p_throws), bb_z = plate_z)]

md <- readRDS(file.path(MDIR, "miss_grade_data_activespin.rds")); setDT(md)
fbA <- md[grp == "fastball", .(fb_active = mean(active_spin, na.rm = TRUE)), by = .(pitcher, season)]
md <- merge(md, fbA, by = c("pitcher","season"), all.x = TRUE)
keep <- c("game_pk","at_bat_number","pitch_number","miss_distance","speed_diff",
          "ax_diff","az_diff","axis_diff","active_spin","fb_active")
d <- merge(d, md[, ..keep], by = c("game_pk","at_bat_number","pitch_number"))
d[, `:=`(velo_gap = abs(speed_diff), move_gap = sqrt(ax_diff^2 + az_diff^2),
         axis_gap = axis_diff, aspin_gap = abs(active_spin - fb_active))]
d[, spin_sim := exp(-(aspin_gap/0.10)^2) * exp(-(axis_gap/45)^2)]
d <- d[is.finite(miss_distance)]
saveRDS(d, CACHE)
} else d <- readRDS(CACHE)

d[, regime := fifelse(strikes == 0L, "0 strikes",
              fifelse(strikes == 1L, "1 strike",
              fifelse(balls < 3L, "2 strikes", "3-2")))]
REG <- c("0 strikes","1 strike","2 strikes","All counts")
cat("=== competitive swings available ===\n")
print(dcast(d[regime != "3-2", .N, by = .(ptype, regime)], ptype ~ regime,
            value.var = "N"), row.names = FALSE)
cat("\nwith a primary fastball immediately before:\n")
print(dcast(d[regime != "3-2" & after_primary == TRUE, .N, by = .(ptype, regime)], ptype ~ regime,
            value.var = "N"), row.names = FALSE)

DRV <- list(
  "Trajectory tunnel"    = list(f = "s(pair_path, k=5)",   seq = TRUE),
  "Spin-axis gap"        = list(f = "s(axis_gap, k=5)",    seq = FALSE),
  "Spin similarity"      = list(f = "s(spin_sim, k=5)",    seq = FALSE),
  "Active-spin gap"      = list(f = "s(aspin_gap, k=5)",   seq = FALSE),
  "Velocity gap"         = list(f = "s(velo_gap, k=5)",    seq = FALSE),
  "Movement gap"         = list(f = "s(move_gap, k=5)",    seq = FALSE),
  "Location (reference)" = list(f = "te(bb_x, bb_z, k=c(5,5))", seq = FALSE))

rmse <- function(a,b) sqrt(mean((a-b)^2, na.rm = TRUE))
out <- list()
for (ty in c("SL","ST","CU","CH","FS")) for (rg in REG) for (nm in names(DRV)) {
  s <- d[ptype == ty & regime != "3-2"]
  if (rg != "All counts") s <- s[regime == rg]
  if (DRV[[nm]]$seq) s <- s[after_primary == TRUE]
  f <- as.formula(paste("miss_distance ~", DRV[[nm]]$f))
  s <- s[complete.cases(s[, all.vars(f), with = FALSE])]
  tr <- s[season <= 2025]; te <- s[season == 2026]
  if (nrow(te) < 250 || nrow(tr) < 1000) next
  m <- tryCatch(gam(f, data = tr), error = function(e) NULL); if (is.null(m)) next
  base <- rmse(te$miss_distance, mean(tr$miss_distance))
  out[[length(out)+1]] <- data.table(ptype = ty, regime = rg, driver = nm,
    n_train = nrow(tr), n_test = nrow(te),
    gain = 100*(base - rmse(te$miss_distance, predict(m, te)))/base)
}
res <- rbindlist(out)
res[, `:=`(regime = factor(regime, levels = REG),
           grp = fifelse(ptype %in% c("SL","ST","CU"), "Breaking", "Offspeed"))]
fwrite(res, file.path(AST, "ext_count_scope_gain.csv"))

cat("\n=== HELD-OUT 2026 RMSE GAIN (%) BY COUNT REGIME ===\n")
for (nm in c("Trajectory tunnel","Spin-axis gap","Spin similarity")) {
  cat("\n--", nm, "--\n")
  print(dcast(res[driver == nm], ptype ~ regime, value.var = "gain")[
    , lapply(.SD, function(z) if (is.numeric(z)) round(z,2) else z)], row.names = FALSE)
}

LAB <- c(SL="Slider", ST="Sweeper", CU="Curveball", CH="Changeup", FS="Splitter")
plt <- res[driver %in% c("Trajectory tunnel","Spin-axis gap","Spin similarity")]
plt[, ptype_lab := factor(LAB[ptype], levels = LAB)]
p <- ggplot(plt, aes(regime, gain, fill = driver)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = .35) +
  geom_col(position = position_dodge(width = .78), width = .7, colour = "black", linewidth = .2) +
  geom_text(aes(label = sprintf("%.2f", gain),
                vjust = fifelse(gain >= 0, -0.45, 1.35)),
            position = position_dodge(width = .78), size = 2.5) +
  facet_wrap(~ ptype_lab, nrow = 1) +
  coord_cartesian(ylim = c(-0.65, 3.15)) +
  scale_fill_manual(values = c(`Trajectory tunnel` = "#c0392b",
                               `Spin-axis gap` = "#2a9d8f",
                               `Spin similarity` = "#8d99ae"), name = NULL) +
  labs(title = "Half the claim holds cleanly. The other half belongs to breaking balls.",
       subtitle = paste0("Held-out 2026 RMSE gain on miss distance over an intercept, trained 2023H2-2025, each driver fit alone, so nothing here depends on a residual definition. The tunnel is\n",
                         "measured only on pitches that actually followed the primary fastball; the spin cues are measured on every competitive swing of that type, because spin similarity compares\n",
                         "a pitch to the pitcher's fastball in general and has no sequence to condition on. Figure 14 scored the spin cues on the tunnel's population, which was the wrong scope for\n",
                         "them, and this fixes it. THE TUNNEL SHARPENS WITH THE COUNT exactly as claimed: sliders 0.00 -> 1.44 -> 2.52, sweepers 0.69 -> 2.71 -> 2.65. THE SPIN CUE DOES NOT LIVE ON\n",
                         "CHANGEUPS: 0.01 percent in all counts, and it is larger on sliders (0.38) than on any offspeed type. The curveball 0-strike tunnel bar is clipped at -6.39, a 1,606-pitch\n",
                         "training sample. For reference the pitch's own location scores 11 to 32 percent on these same samples."),
       x = "Count regime (3-2 excluded)", y = "Held-out 2026 RMSE improvement vs intercept (%)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 8.2), panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), legend.position = "top",
        strip.text = element_text(face = "bold"), axis.text.x = element_text(size = 7.6))
ggsave(file.path(AST, "fig15_count_scope.png"), p, width = 15, height = 6.4, dpi = 150)
cat("\nwrote fig15_count_scope.png and ext_count_scope_gain.csv\n")
