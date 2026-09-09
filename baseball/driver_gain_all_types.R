#!/usr/bin/env Rscript

# THE TWO CLAIMS, SCORED IN ONE CURRENCY.
#
# Claim A: the trajectory tunnel matters for breaking balls.
# Claim B: spin similarity matters for offspeed.
#
# Figure 3 answered A with held-out 2026 RMSE gain on miss distance, but its pair table
# (tunnel_pairs.rds) only ever contained SL/CU/ST, so B has never been scored the same way.
# Every offspeed spin result in the project is instead a between-pitcher correlation against
# a residual, which is the family of test that turned out to be specification-dependent.
#
# This rebuilds the pair table for all five types from raw Statcast and runs the identical
# driver decomposition, so the tunnel number for sliders and the spin number for changeups
# are directly comparable. Scope matches the rest of the recent work: put-away counts
# (0-2, 1-2, 2-2) on a pitch thrown immediately after the pitcher's PRIMARY fastball.
#
# Currency: train 2023H2-2025, test 2026, report the percent reduction in RMSE of
# miss_distance over an intercept. A cue that cannot beat an intercept out of sample is not
# a cue.

suppressPackageStartupMessages({ library(data.table); library(mgcv); library(ggplot2) })
options(width = 200)
MDIR <- file.path("data","statcast_model"); AST <- file.path(MDIR, "article_assets")
FASTBALLS <- c("FF","SI","FC"); yf <- 17/12
WX <- 1.2; WY <- 0.4; WZ <- 1.4; REACT <- 0.150; NSTEP <- 40
CACHE <- file.path(MDIR, "driver_pairs_all_types.rds")

if (!file.exists(CACHE) || nzchar(Sys.getenv("REFIT"))) {
kin <- rbindlist(lapply(2023:2026, function(y) {
  k <- fread(sprintf("data/statcast_%d/statcast_%d_all.csv", y, y),
        select = c("game_pk","at_bat_number","pitch_number","pitch_type","game_date","game_type",
                   "pitcher","p_throws","balls","strikes","plate_x","plate_z","release_pos_x",
                   "release_pos_y","release_pos_z","vx0","vy0","vz0","ax","ay","az"),
        showProgress = FALSE)
  k[, season := y]; k }))
kin <- kin[game_type == "R" & !is.na(vx0) & !is.na(release_pos_y) & pitch_type != ""]
kin[season == 2023L, game_date := as.Date(game_date)]
kin <- kin[season != 2023L | game_date >= as.Date("2023-07-14")]
kin <- unique(kin, by = c("game_pk","at_bat_number","pitch_number"))
setorder(kin, game_pk, at_bat_number, pitch_number)

# Primary fastball per pitcher-season, same FF > SI > FC precedence used everywhere else.
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
d <- kin[after_primary == TRUE & strikes == 2L & balls < 3L & pitch_type %in% names(TY)]
d[, ptype := TY[pitch_type]]
cat("put-away pitches after the primary fastball:", nrow(d), "\n")

pos <- function(r,v,a,t) r + v*t + 0.5*a*t^2
Tmax <- pmax(d$t_react, d$p_t_react); acc <- numeric(nrow(d))
for (k in 1:NSTEP) {
  tk <- (k-0.5)/NSTEP * Tmax
  dx <- pos(d$release_pos_x,d$vx0,d$ax,tk) - pos(d$p_release_pos_x,d$p_vx0,d$p_ax,tk)
  dy <- pos(d$release_pos_y,d$vy0,d$ay,tk) - pos(d$p_release_pos_y,d$p_vy0,d$p_ay,tk)
  dz <- pos(d$release_pos_z,d$vz0,d$az,tk) - pos(d$p_release_pos_z,d$p_vz0,d$p_az,tk)
  acc <- acc + sqrt(WX*dx^2 + WY*dy^2 + WZ*dz^2) * (Tmax/NSTEP)
}
d[, tunnel := acc]
d[, plate_sep := pmax(sqrt((plate_x - p_plate_x)^2 + (plate_z - p_plate_z)^2), 0.1)]
d[, pair_path := tunnel / plate_sep]

mir <- function(x, h) fifelse(h == "R", x, -x)
d[, `:=`(bb_x = mir(plate_x, p_throws), fb_x = mir(p_plate_x, p_throws), bb_z = plate_z,
         fb_z = p_plate_z)]
d[, `:=`(dx = bb_x - fb_x, dz = bb_z - fb_z)]

md <- readRDS(file.path(MDIR, "miss_grade_data_activespin.rds")); setDT(md)
fbA <- md[grp == "fastball", .(fb_active = mean(active_spin, na.rm = TRUE)), by = .(pitcher, season)]
md <- merge(md, fbA, by = c("pitcher","season"), all.x = TRUE)
keep <- c("game_pk","at_bat_number","pitch_number","miss_distance","is_whiff","speed_diff",
          "ax_diff","az_diff","axis_diff","active_spin","fb_active")
d <- merge(d, md[, ..keep], by = c("game_pk","at_bat_number","pitch_number"))
d[, `:=`(velo_gap = abs(speed_diff), move_gap = sqrt(ax_diff^2 + az_diff^2),
         axis_gap = axis_diff, aspin_gap = abs(active_spin - fb_active))]
d <- d[is.finite(miss_distance)]
saveRDS(d, CACHE)
} else d <- readRDS(CACHE)

cat("\ncompetitive put-away swings by type:\n")
print(d[, .(n = .N, n_2026 = sum(season == 2026), mean_miss = round(mean(miss_distance),2),
            mean_axis_gap = round(mean(axis_gap, na.rm=TRUE),1)), by = ptype][order(-n)],
      row.names = FALSE)

BLOCKS <- list(
  "Breaking/offspeed location"       = "te(bb_x, bb_z, k=c(5,5))",
  "FB->pitch separation"             = "te(dx, dz, k=c(5,5))",
  "Setup-fastball location"          = "te(fb_x, fb_z, k=c(5,5))",
  "Trajectory tunnel (path ratio)"   = "s(pair_path, k=5)",
  "Spin-axis gap vs FB"              = "s(axis_gap, k=5)",
  "Active-spin gap vs FB"            = "s(aspin_gap, k=5)",
  "Velocity gap vs FB"               = "s(velo_gap, k=5)",
  "Movement gap vs FB"               = "s(move_gap, k=5)")

rmse <- function(a,b) sqrt(mean((a-b)^2, na.rm=TRUE))
TYPES <- c("SL","ST","CU","CH","FS")
out <- list()
for (ty in TYPES) {
  tr <- d[ptype == ty & season <= 2025]; te <- d[ptype == ty & season == 2026]
  for (nm in names(BLOCKS)) {
    f <- as.formula(paste("miss_distance ~", BLOCKS[[nm]]))
    a <- tr[complete.cases(tr[, all.vars(f), with = FALSE])]
    b <- te[complete.cases(te[, all.vars(f), with = FALSE])]
    if (nrow(b) < 200) next
    m <- tryCatch(gam(f, data = a), error = function(e) NULL); if (is.null(m)) next
    base <- rmse(b$miss_distance, mean(a$miss_distance))
    out[[length(out)+1]] <- data.table(ptype = ty, driver = nm, n_test = nrow(b),
      gain = 100*(base - rmse(b$miss_distance, predict(m, b)))/base)
  }
}
res <- rbindlist(out)
fwrite(res, file.path(AST, "ext_driver_gain_all_types.csv"))
cat("\n=== HELD-OUT 2026 RMSE GAIN (%) ON MISS DISTANCE, PUT-AWAY, AFTER PRIMARY FB ===\n")
w <- dcast(res, driver ~ ptype, value.var = "gain")
print(w[order(-SL)][, lapply(.SD, function(z) if (is.numeric(z)) round(z,2) else z)], row.names = FALSE)

LAB <- c(SL="Slider", ST="Sweeper", CU="Curveball (CU+KC)", CH="Changeup", FS="Splitter")
res[, ptype_lab := factor(LAB[ptype], levels = LAB[TYPES])]
res[, kind := fifelse(grepl("location|separation", driver), "Where it finishes",
             fifelse(driver == "Trajectory tunnel (path ratio)", "Trajectory tunnel",
                     "Spin / shape cue"))]
ord <- res[, .(m = mean(gain)), by = driver][order(m)]$driver
res[, driver := factor(driver, levels = ord)]

p <- ggplot(res, aes(driver, gain, fill = kind)) +
  geom_col(width = .72) + geom_hline(yintercept = 0, colour = "black", linewidth = .3) +
  geom_text(aes(label = sprintf("%.2f", gain), hjust = ifelse(gain < 4, -0.15, 1.12)),
            size = 2.9, colour = ifelse(res$gain < 4, "black", "white")) +
  coord_flip() + facet_wrap(~ ptype_lab, nrow = 1) +
  scale_fill_manual(values = c(`Where it finishes` = "#2c3e50",
                               `Trajectory tunnel` = "#c0392b",
                               `Spin / shape cue`  = "#95a5a6"), name = NULL) +
  labs(title = "The two claims in one currency: held-out prediction of miss distance",
       subtitle = paste0("Percent reduction in 2026 RMSE over an intercept, each driver fit alone as a GAM on 2023H2-2025. Put-away counts (0-2, 1-2, 2-2) on a pitch thrown immediately\n",
                         "after the pitcher's primary fastball, competitive swings only. This is the first time the offspeed spin claim has been scored the same way as the breaking-ball\n",
                         "tunnel claim: every previous offspeed spin result was a between-pitcher correlation against a residual, the family of test that proved specification-dependent.\n",
                         "The trajectory tunnel clears an intercept on every breaking type. Spin-axis gap does not clear it on changeups."),
       x = NULL, y = "Held-out 2026 RMSE improvement vs intercept (%)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 8.2), panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold"), legend.position = "top")
ggsave(file.path(AST, "fig14_driver_gain_all_types.png"), p, width = 15, height = 6.2, dpi = 150)
cat("\nwrote fig14_driver_gain_all_types.png and ext_driver_gain_all_types.csv\n")
