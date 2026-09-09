#!/usr/bin/env Rscript

# FIG 3 (redesign) - what actually drives a miss on a 2-strike breaking ball,
# split by breaking type (SL / CU+KC / ST), with LOCATION folded in as a driver.
#
# Every driver is scored in the same currency so they are directly comparable:
# held-out 2026 RMSE improvement over an intercept, after training on 2023H2-2025.
# Location blocks come from the sequential FB->BB pair table; shape/spin gaps come
# from the miss-model feature set, joined on the pitch identifiers.

suppressPackageStartupMessages({ library(data.table); library(mgcv); library(ggplot2) })

MDIR <- file.path("data","statcast_model")
AST  <- file.path(MDIR, "article_assets")
ODIR <- file.path(MDIR, "tunnel_location")

pr <- readRDS(file.path(MDIR, "tunnel_pairs.rds"))
pr <- pr[is_primary_setup == TRUE & is_swing == TRUE & is.finite(miss_distance)]
mir <- function(x, h) fifelse(h == "R", x, -x)
pr[, `:=`(bb_x = mir(bb_plate_x, p_throws), fb_x = mir(fb_plate_x, p_throws))]
pr[, `:=`(dx = bb_x - fb_x, dz = bb_plate_z - fb_plate_z, pair_path = path_ratio)]

md <- readRDS(file.path(MDIR, "miss_grade_data_activespin.rds"))
fb <- md[grp == "fastball", .N, by = .(pitcher, season, pitch_type, active_spin)][
        , .(nn = sum(N), fb_active = weighted.mean(active_spin, N)), by = .(pitcher, season)]
md <- merge(md, fb, by = c("pitcher","season"), all.x = TRUE)
keep <- c("game_pk","at_bat_number","pitch_number","season","speed_diff","ax_diff","az_diff",
          "axis_diff","spin_eff_diff","active_spin","fb_active","release_extension")
d <- merge(pr, md[, ..keep], by = c("game_pk","at_bat_number","pitch_number","season"))
message(sprintf("joined %d of %d pair-swings to shape features", nrow(d), nrow(pr)))

d[, `:=`(velo_gap = abs(speed_diff),
         move_gap = sqrt(ax_diff^2 + az_diff^2),
         axis_gap = axis_diff,
         aspin_gap = abs(active_spin - fb_active))]

BLOCKS <- list(
  "Breaking-ball location"     = "te(bb_x, bb_plate_z, k=c(5,5))",
  "FB->BB separation"          = "te(dx, dz, k=c(5,5))",
  "Setup-fastball location"    = "te(fb_x, fb_plate_z, k=c(5,5))",
  "Trajectory tunnel (path ratio)" = "s(pair_path, k=5)",
  "Velocity gap vs FB"         = "s(velo_gap, k=5)",
  "Movement gap vs FB"         = "s(move_gap, k=5)",
  "Spin-axis gap vs FB"        = "s(axis_gap, k=5)",
  "Active-spin gap vs FB"      = "s(aspin_gap, k=5)")

rmse <- function(a,b) sqrt(mean((a-b)^2, na.rm=TRUE))
TYPES <- c("SL","CU","ST")
out <- list()
for (ty in TYPES) {
  tr <- d[brk_type == ty & season <= 2025]; te <- d[brk_type == ty & season == 2026]
  b0 <- rmse(te$miss_distance, mean(tr$miss_distance))
  for (nm in names(BLOCKS)) {
    f <- as.formula(paste("miss_distance ~", BLOCKS[[nm]]))
    ok_tr <- tr[complete.cases(tr[, all.vars(f), with=FALSE])]
    ok_te <- te[complete.cases(te[, all.vars(f), with=FALSE])]
    m <- tryCatch(gam(f, data = ok_tr), error = function(e) NULL)
    if (is.null(m)) next
    g <- 100*(rmse(ok_te$miss_distance, mean(ok_tr$miss_distance)) -
              rmse(ok_te$miss_distance, predict(m, ok_te))) /
          rmse(ok_te$miss_distance, mean(ok_tr$miss_distance))
    out[[length(out)+1]] <- data.table(brk_type=ty, driver=nm, gain=g, n_test=nrow(ok_te))
  }
  # everything-else on top of location, to show what deception adds once you know
  # where the pitch finished
  f_loc <- miss_distance ~ te(bb_x, bb_plate_z, k=c(5,5))
  f_all <- miss_distance ~ te(bb_x, bb_plate_z, k=c(5,5)) + s(pair_path, k=5) +
           s(velo_gap, k=5) + s(move_gap, k=5) + s(axis_gap, k=5) + s(aspin_gap, k=5)
  cc <- function(x, f) x[complete.cases(x[, all.vars(f), with=FALSE])]
  tr2 <- cc(tr, f_all); te2 <- cc(te, f_all)
  ml <- gam(f_loc, data = tr2); ma <- gam(f_all, data = tr2)
  base <- rmse(te2$miss_distance, mean(tr2$miss_distance))
  out[[length(out)+1]] <- data.table(brk_type=ty, driver="Location + all deception cues",
    gain = 100*(base - rmse(te2$miss_distance, predict(ma, te2)))/base, n_test=nrow(te2))
  out[[length(out)+1]] <- data.table(brk_type=ty, driver="  (location alone, same sample)",
    gain = 100*(base - rmse(te2$miss_distance, predict(ml, te2)))/base, n_test=nrow(te2))
}
res <- rbindlist(out)
fwrite(res, file.path(AST, "ext_miss_drivers_by_type.csv"))

cat("=== HELD-OUT 2026 RMSE GAIN (%) BY DRIVER AND BREAKING TYPE ===\n")
print(dcast(res, driver ~ brk_type, value.var = "gain")[order(-SL)])

TYPE_LAB <- c(SL="Slider", CU="Curveball (CU+KC)", ST="Sweeper")
plt <- res[!grepl("same sample", driver)]
plt[, kind := fifelse(driver == "Location + all deception cues", "combined",
              fifelse(grepl("location|separation", driver), "location", "deception cue"))]
plt[, brk_lab := TYPE_LAB[brk_type]]
ord <- plt[, .(m = mean(gain)), by = driver][order(m)]$driver
plt[, driver := factor(driver, levels = ord)]

ACC <- "#c0392b"; GREY <- "#95a5a6"; BLU <- "#2c3e50"
p <- ggplot(plt, aes(driver, gain, fill = kind)) +
  geom_col(width = .7) +
  geom_hline(yintercept = 0, colour = "black", linewidth = .3) +
  geom_text(aes(label = sprintf("%.1f", gain), hjust = ifelse(gain < 1, -0.2, 1.15)),
            size = 3.1, colour = ifelse(plt$gain < 1, "black", "white")) +
  coord_flip() +
  facet_wrap(~ brk_lab) +
  scale_fill_manual(values = c(location = BLU, `deception cue` = ACC, combined = GREY), name = NULL) +
  labs(title = "For any single breaking ball, where it finishes decides the miss",
       subtitle = paste0(
"Held-out 2026 RMSE gain over an intercept (%), trained on 2023H2-2025, 2-strike breaking balls after the primary fastball.\n",
"Setup-fastball location is worth nothing, and the FB->BB separation is mostly a restatement of where the pitch finished.\n",
"Shape and spin gaps are near-constant within a pitcher's own pitch type, so they separate PITCHERS (Figs 5-11), not pitches."),
       x = NULL, y = "held-out RMSE improvement vs intercept (%)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face="bold"), panel.grid.major.y = element_blank(),
        strip.text = element_text(face="bold"), legend.position = "top")
ggsave(file.path(AST, "fig3_miss_drivers_by_type.png"), p, width = 12.5, height = 6.8, dpi = 150)
ggsave(file.path(ODIR, "fig3_miss_drivers_by_type.png"), p, width = 12.5, height = 6.4, dpi = 150)
cat("\nwrote fig3_miss_drivers_by_type.png\n")

cat("\n=== WHAT DECEPTION ADDS ON TOP OF LOCATION ===\n")
inc <- dcast(res[grepl("Location \\+|same sample", driver)], brk_type ~ driver, value.var = "gain")
setnames(inc, c("brk_type","loc_only","loc_plus_cues"))
inc[, incremental := round(loc_plus_cues - loc_only, 2)]
print(inc)
