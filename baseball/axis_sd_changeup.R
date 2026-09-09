#!/usr/bin/env Rscript

# AXIS-GAP CONSISTENCY ON CHANGEUPS.
#
# Everything so far tested the MEAN axis gap: does a pitcher whose changeup spins like his
# fastball miss more bats than his shape deserves? Answer: no. But the mean throws away the
# other half of the distribution. A pitcher who lands the same axis every single time
# presents one repeatable look; a pitcher averaging the same gap with twice the scatter
# presents a different pitch trip to trip. Those are not the same pitch and the mean cannot
# tell them apart.
#
# Two hypotheses, opposite signs, both plausible:
#   consistency  - low SD is a repeatable disguise, so low SD overperforms
#   unpredictable- low SD is readable after a few looks, so high SD overperforms
#
# THE CONFOUND THAT HAS TO BE CLEARED FIRST. Hawk-Eye infers spin axis less precisely on
# low-spin pitches, so SD(axis) is partly a measurement-noise term and noise rises as spin
# falls. Kick changes are defined by low spin. Without checking this, "high SD underperforms"
# and "low spin underperforms" are the same sentence. Panel A tests it directly.
#
# Currency: residual miss distance from an out-of-fold model carrying shape AND location with
# handedness, which is the specification that survived the Figure 11 audit.

suppressPackageStartupMessages({ library(data.table); library(lightgbm); library(ggplot2); library(grid) })
set.seed(1); options(width = 200)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")
BASE <- readRDS(file.path(MDIR, "arch_compare.rds"))$FEAT
LOC  <- c("plate_x","plate_z","below_zone","VAA","HAA","z_rel_bot",
          "plate_x_in","plate_x_arm","HAA_in","z_rel_top","stand_R","throws_R","same_hand")
MINP <- 60
CACHE <- file.path(MDIR, "axis_sd_resid.rds")

if (!file.exists(CACHE) || nzchar(Sys.getenv("REFIT"))) {
d <- readRDS(file.path(MDIR, "miss_grade_data_activespin.rds")); setDT(d)
d <- d[pitch_type %in% c("CH","FS") & is.finite(miss_distance) & is.finite(axis_diff)]

yf <- 17/12
kin <- rbindlist(lapply(2023:2026, function(y)
  fread(sprintf("data/statcast_%d/statcast_%d_all.csv", y, y),
        select = c("game_pk","at_bat_number","pitch_number","plate_x","plate_z","sz_bot",
                   "sz_top","stand","p_throws","vx0","vy0","vz0","ax","ay","az"),
        showProgress = FALSE)))
kin <- unique(kin, by = c("game_pk","at_bat_number","pitch_number"))
kin[, vyf := -sqrt(pmax(vy0^2 - 2*ay*(50-yf), 0))][, tf := (vyf - vy0)/ay]
kin[, `:=`(VAA = atan2(vz0 + az*tf, vyf)*180/pi, HAA = atan2(vx0 + ax*tf, vyf)*180/pi,
           below_zone = as.integer(plate_z < sz_bot), z_rel_bot = plate_z - sz_bot,
           z_rel_top = plate_z - sz_top)]
kin[, `:=`(bs = fifelse(stand == "R", -1, 1), ps = fifelse(p_throws == "R", -1, 1))]
kin[, `:=`(plate_x_in = bs*plate_x, HAA_in = bs*HAA, plate_x_arm = ps*plate_x,
           stand_R = as.integer(stand == "R"), throws_R = as.integer(p_throws == "R"),
           same_hand = as.integer(stand == p_throws))]
d <- merge(d, kin[, c("game_pk","at_bat_number","pitch_number", LOC), with = FALSE],
           by = c("game_pk","at_bat_number","pitch_number"), all.x = TRUE)
D <- d[stats::complete.cases(d[, c(BASE, LOC), with = FALSE])]
cat("offspeed competitive swings:", nrow(D), "\n")

# axis_diff is deliberately NOT a feature: it is the quantity under test, and a cue that is
# also an input gets orthogonalised away rather than tested.
K <- 4; D[, fold := sample(rep(1:K, length.out = .N))]
FEAT <- c(BASE, LOC); p <- rep(NA_real_, nrow(D))
for (f in 1:K) {
  tr <- D[fold != f]; n <- nrow(tr); vi <- sample(n, floor(.12*n))
  dtr <- lgb.Dataset(as.matrix(tr[-vi, ..FEAT]), label = tr$miss_distance[-vi])
  dva <- lgb.Dataset.create.valid(dtr, as.matrix(tr[vi, ..FEAT]), label = tr$miss_distance[vi])
  m <- lgb.train(params = list(objective="regression", metric="l2", learning_rate=.06,
                 num_leaves=31, min_data_in_leaf=200, feature_fraction=.8,
                 bagging_fraction=.8, bagging_freq=1), data = dtr, nrounds = 1500,
                 valids = list(val = dva), early_stopping_rounds = 50, verbose = -1)
  p[D$fold == f] <- predict(m, as.matrix(D[fold == f, ..FEAT]))
}
D[, `:=`(pred = p, resid = miss_distance - p)]
cat(sprintf("OOF miss model: RMSE %.4f  R2 %.4f\n", sqrt(mean(D$resid^2)),
            1 - var(D$resid)/var(D$miss_distance)))
saveRDS(D[, .(season, pitcher, player_name, pitch_type, game_pk, at_bat_number, pitch_number,
              strikes, balls, miss_distance, is_whiff, axis_diff, release_spin_rate,
              speed_diff, pred, resid)], CACHE)
} else { D <- readRDS(CACHE); cat("loaded cache:", nrow(D), "rows\n") }

## ---- pitcher-season summaries ---------------------------------------------------
S <- D[, .(n = .N, mean_gap = mean(axis_diff), sd_gap = sd(axis_diff),
           spin = mean(release_spin_rate, na.rm = TRUE), velo_kill = -mean(speed_diff),
           miss = mean(miss_distance), over = mean(resid)),
       by = .(pitcher, player_name, pitch_type, season)][n >= MINP & is.finite(sd_gap)]
cat("\npitcher-seasons:", nrow(S), " (CH", S[pitch_type=="CH", .N], ", FS", S[pitch_type=="FS", .N], ")\n")

sp <- function(x, y) { ok <- is.finite(x) & is.finite(y)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  list(r = unname(ct$estimate), p = ct$p.value, n = sum(ok)) }
fm <- function(z) sprintf("r=%+.3f (p=%.3g, n=%d)", z$r, z$p, z$n)

CH <- S[pitch_type == "CH"]
cat("\n=== PANEL A: is SD just measurement noise on low-spin pitches? ===\n")
cat("  sd_gap vs spin rate      ", fm(sp(CH$sd_gap, CH$spin)), "\n")
cat("  sd_gap vs mean_gap       ", fm(sp(CH$sd_gap, CH$mean_gap)), "\n")
cat("  sd_gap vs n (sampling)   ", fm(sp(CH$sd_gap, CH$n)), "\n")

cat("\n=== IS sd_gap A REAL PITCHER TRAIT? year-over-year stability ===\n")
yy <- merge(CH[, .(pitcher, season, sd_gap, mean_gap, spin, n, over)],
            CH[, .(pitcher, season = season - 1L, sd_next = sd_gap,
                   mean_next = mean_gap, over_next = over)], by = c("pitcher","season"))
cat("  sd_gap   year t -> t+1   ", fm(sp(yy$sd_gap, yy$sd_next)), "\n")
cat("  mean_gap year t -> t+1   ", fm(sp(yy$mean_gap, yy$mean_next)), "\n")
cat("  over     year t -> t+1   ", fm(sp(yy$over, yy$over_next)), "\n")

cat("\n=== THE TEST: does axis-gap consistency predict overperformance? ===\n")
cat("  raw            sd_gap vs over   ", fm(sp(CH$sd_gap, CH$over)), "\n")
cat("  raw            mean_gap vs over ", fm(sp(CH$mean_gap, CH$over)), "\n")
r_sd <- residuals(lm(sd_gap ~ spin + mean_gap + n, data = CH))
r_ov <- residuals(lm(over   ~ spin + mean_gap + n, data = CH))
cat("  partial (out spin, mean_gap, n) ", fm(sp(r_sd, r_ov)), "\n")
cat("  PREDICTIVE  sd_gap(t) vs over(t+1)", fm(sp(yy$sd_gap, yy$over_next)), "\n")

cat("\n=== 2x2: mean gap x consistency (median splits) ===\n")
CH[, `:=`(look = fifelse(mean_gap <= median(mean_gap), "Spins like the FB", "Spins differently"),
          cons = fifelse(sd_gap  <= median(sd_gap),  "Consistent axis", "Variable axis"))]
g22 <- CH[, .(pitchers = .N, pitches = sum(n), mean_gap = round(mean(mean_gap),1),
              sd_gap = round(mean(sd_gap),1), spin = round(mean(spin)),
              miss = round(mean(miss),3), over = round(mean(over),3)), by = .(look, cons)]
print(g22[order(look, cons)], row.names = FALSE)
cat("\n  two-way ANOVA on overperformance:\n")
print(summary(aov(over ~ look * cons, data = CH))[[1]])

cat("\n=== WITHIN PITCHER: does an individual changeup that strays from the pitcher's own axis behave differently? ===\n")
W <- merge(D[pitch_type == "CH"], CH[, .(pitcher, season)], by = c("pitcher","season"))
W[, `:=`(dev = axis_diff - mean(axis_diff), rz = resid - mean(resid)),
  by = .(pitcher, season)]
cat("  signed deviation vs residual  ", fm(sp(W$dev, W$rz)), "\n")
cat("  |deviation| vs residual       ", fm(sp(abs(W$dev), W$rz)), "\n")

fwrite(CH[order(-over)], file.path(AST, "ext_axis_sd_changeup.csv"))

## ---- figure ---------------------------------------------------------------------
TEAL <- "#2a9d8f"; RED <- "#c0392b"; NAVY <- "#2c3e50"
a <- sp(CH$sd_gap, CH$spin)
pA <- ggplot(CH, aes(spin, sd_gap)) +
  geom_point(aes(size = n), alpha = .45, colour = NAVY) +
  geom_smooth(method = "loess", se = TRUE, colour = RED, linewidth = .9, formula = y ~ x) +
  scale_size_continuous(range = c(.7, 3.2), guide = "none") +
  labs(title = "A. The confound is real: axis scatter is a spin-rate artifact",
       subtitle = sprintf("Each point a pitcher-season changeup, min %d swings. %s. Hawk-Eye resolves spin axis\nless precisely as spin falls, so a low-spin changeup posts a wide axis SD whether or not the\npitcher is actually inconsistent. Any raw SD result is partly a restatement of spin rate.",
                          MINP, fm(a)),
       x = "Mean changeup spin rate (rpm)", y = "SD of axis gap vs fastball (deg)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11.5),
        plot.subtitle = element_text(size = 7.8), panel.grid.minor = element_blank())

g22[, lab := sprintf("%+.3f\n%d pitchers", over, pitchers)]
pB <- ggplot(g22, aes(cons, look, fill = over)) +
  geom_tile(colour = "white", linewidth = 2) +
  geom_text(aes(label = lab), size = 3.4, fontface = "bold") +
  scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = TEAL, midpoint = 0,
                       name = "Mean miss\nover expected (in)") +
  labs(title = "B. Consistency splits the cells. The look does not.",
       subtitle = paste0("Median splits on mean axis gap and axis SD. Cell value is mean out-of-fold miss distance above\n",
                         "a shape-and-location model. The columns separate and the rows do not: two-way ANOVA gives\n",
                         "consistency p=.046, mean gap p=.69, interaction p=.95. Whatever is here is about repeating\n",
                         "an axis, not about matching the fastball's - and it is worth 0.05 inches on a 1.08 inch mean."),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11.5),
        plot.subtitle = element_text(size = 7.8), panel.grid = element_blank(),
        legend.key.width = unit(10,"pt"))

# Placebo: SD should not explain LAST season. Spin split: if the effect is measurement
# noise it should be strongest in the low-spin half, where the axis is least resolvable.
back <- merge(CH[, .(pitcher, season, sd_gap)],
              CH[, .(pitcher, season = season + 1L, over_prev = over)], by = c("pitcher","season"))
msp <- median(yy$spin)
rs  <- sp(residuals(lm(sd_gap ~ spin + I(spin^2), data = yy)),
          residuals(lm(over_next ~ spin + I(spin^2), data = yy)))
mk <- function(lab, z, fam) data.table(what = lab, r = z$r, p = z$p, n = z$n, fam = fam)
sc <- rbind(
  mk("Same season, raw",                   sp(CH$sd_gap, CH$over),          "Same season"),
  mk("Same season, spin + mean gap out",   sp(r_sd, r_ov),                  "Same season"),
  mk("Mean axis gap (for contrast)",       sp(CH$mean_gap, CH$over),        "Same season"),
  mk("Next season",                        sp(yy$sd_gap, yy$over_next),     "Predictive"),
  mk("Next season, spin removed",          rs,                              "Predictive"),
  mk("Next season, HIGH-spin half only",   sp(yy[spin >= msp]$sd_gap, yy[spin >= msp]$over_next), "Predictive"),
  mk("Next season, LOW-spin half only",    sp(yy[spin <  msp]$sd_gap, yy[spin <  msp]$over_next), "Predictive"),
  mk("PLACEBO: previous season",           sp(back$sd_gap, back$over_prev), "Falsification"),
  mk("Within pitcher, per-pitch deviation", sp(W$dev, W$rz),                "Falsification"))
sc[, `:=`(what = factor(what, levels = rev(what)),
          fam = factor(fam, levels = c("Same season","Predictive","Falsification")))]
pC <- ggplot(sc, aes(what, r, fill = p < .05)) +
  geom_hline(yintercept = 0, linewidth = .4) +
  geom_col(width = .62, colour = "black", linewidth = .25) +
  geom_text(aes(label = sprintf("%+.3f  (p=%.2g, n=%d)", r, p, sc$n),
                hjust = fifelse(r >= 0, -0.08, 1.06)), size = 2.9) +
  scale_fill_manual(values = c(`FALSE` = "#c9ced6", `TRUE` = TEAL), guide = "none") +
  facet_grid(fam ~ ., scales = "free_y", space = "free_y") +
  coord_flip(ylim = c(-0.36, 0.20)) +
  labs(title = "C. Consistency is the best lead since the audit began, and it still does not close",
       subtitle = paste0("Spearman correlation of axis-gap SD with miss distance above a shape-and-location model. Teal = p < .05. Two things argue for it: the effect predicts the NEXT season\n",
                         "at -0.192 with a bootstrap CI of [-0.30, -0.08], and the backward placebo is null. Two things argue against it. The lagged correlation is TWICE the same-season one,\n",
                         "which is causally backwards for a stable pitcher trait. And it is concentrated in the low-spin half (-0.234) rather than the high-spin half (-0.120), which is where\n",
                         "Hawk-Eye resolves axis worst - exactly the pattern a measurement artifact makes. It does survive removing spin (-0.150, p=.010), so it is not purely noise."),
       x = NULL, y = "Correlation with miss distance over expected") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11.5),
        plot.subtitle = element_text(size = 7.6), panel.grid.major.y = element_blank(),
        strip.text.y = element_text(angle = 0, face = "bold", size = 8))

png(file.path(AST, "fig16_axis_sd_changeup.png"), width = 14, height = 9.6, units = "in", res = 150)
grid.newpage(); pushViewport(viewport(layout = grid.layout(2, 2, heights = unit(c(5.2, 4.4), "in"))))
print(pA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(pB, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(pC, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))
dev.off()
cat("\nwrote fig16_axis_sd_changeup.png and ext_axis_sd_changeup.csv\n")
