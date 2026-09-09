#!/usr/bin/env Rscript

# Is the breaking ball's miss driven by WHERE IT LANDS, or by HOW FAR IT SEPARATES
# from the fastball that set it up?
#
# Absolute framing : E[miss | bb_plate_x, bb_plate_z]
# Relative framing : E[miss | dx, dz]  where d = BB location - setup FB location
#
# Fit both (and both together) as 2D smooths on 2023H2-2025 and score held-out 2026.
# Handedness is mirrored so + x = pitcher's glove side.

suppressPackageStartupMessages({ library(data.table); library(mgcv) })

MDIR <- file.path("data","statcast_model")
p <- readRDS(file.path(MDIR, "tunnel_pairs.rds"))
p <- p[is_primary_setup == TRUE & is_swing == TRUE & is.finite(miss_distance)]

mir <- function(x, h) fifelse(h == "R", x, -x)
p[, `:=`(fb_x = mir(fb_plate_x, p_throws), bb_x = mir(bb_plate_x, p_throws))]
p[, `:=`(dx = bb_x - fb_x, dz = bb_plate_z - fb_plate_z)]
p[, matchup := fifelse(p_throws == stand, "same", "opp")]

cat("=== SEPARATION BETWEEN THE TWO PITCHES AT THE PLATE ===\n")
print(p[, .(n = .N,
            dx_mean = round(mean(dx),2), dx_sd = round(sd(dx),2),
            dz_mean = round(mean(dz),2), dz_sd = round(sd(dz),2),
            plate_sep = round(mean(plate_sep),2)), by = brk_type][order(brk_type)])

train <- p[season <= 2025]; test <- p[season == 2026]
rmse <- function(a, b) sqrt(mean((a-b)^2))

res <- list()
for (ty in c("SL","CU","ST")) for (mu in c("same","opp")) {
  tr <- train[brk_type == ty & matchup == mu]; te <- test[brk_type == ty & matchup == mu]
  if (nrow(tr) < 500 || nrow(te) < 200) next
  K <- c(5,5)
  fits <- list(
    base = tryCatch(gam(miss_distance ~ 1, data = tr), error = function(e) NULL),
    abs  = tryCatch(gam(miss_distance ~ te(bb_x, bb_plate_z, k = K), data = tr), error = function(e) NULL),
    del  = tryCatch(gam(miss_distance ~ te(dx, dz, k = K), data = tr), error = function(e) NULL),
    both = tryCatch(gam(miss_distance ~ te(bb_x, bb_plate_z, k = K) + te(dx, dz, k = K), data = tr), error = function(e) NULL),
    setup= tryCatch(gam(miss_distance ~ te(fb_x, fb_plate_z, k = K), data = tr), error = function(e) NULL))
  r <- sapply(fits, function(m) if (is.null(m)) NA_real_ else rmse(te$miss_distance, predict(m, te)))
  res[[length(res)+1]] <- data.table(brk_type = ty, matchup = mu,
    n_train = nrow(tr), n_test = nrow(te),
    rmse_base = r[["base"]], rmse_abs = r[["abs"]], rmse_delta = r[["del"]],
    rmse_both = r[["both"]], rmse_setup = r[["setup"]])
}
res <- rbindlist(res)
res[, `:=`(gain_abs   = round(100*(rmse_base - rmse_abs)/rmse_base, 2),
           gain_delta = round(100*(rmse_base - rmse_delta)/rmse_base, 2),
           gain_both  = round(100*(rmse_base - rmse_both)/rmse_base, 2),
           gain_setup = round(100*(rmse_base - rmse_setup)/rmse_base, 2))]

cat("\n=== HELD-OUT 2026 RMSE GAIN vs INTERCEPT (%, higher = more explanatory) ===\n")
print(res[, .(brk_type, matchup, n_test,
              abs_BB_loc = gain_abs, FB_to_BB_offset = gain_delta,
              both = gain_both, setup_FB_loc = gain_setup)])

cat("\n=== AVERAGE ACROSS TYPES ===\n")
print(res[, .(abs_BB_loc = round(mean(gain_abs),2),
              FB_to_BB_offset = round(mean(gain_delta),2),
              both = round(mean(gain_both),2),
              setup_FB_loc = round(mean(gain_setup),2))])

# How much of the offset signal is just the BB's own location? If dx/dz only matter
# because a big drop means the BB ended up low, the two framings are redundant.
cat("\n=== CORRELATION OF OFFSET WITH ABSOLUTE BB LOCATION ===\n")
print(p[, .(cor_dx_bbx = round(cor(dx, bb_x),2),
            cor_dz_bbz = round(cor(dz, bb_plate_z),2)), by = brk_type][order(brk_type)])

fwrite(res, file.path(MDIR, "tunnel_location", "delta_vs_absolute_test.csv"))
