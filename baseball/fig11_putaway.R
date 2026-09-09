#!/usr/bin/env Rscript

# FIG 11 ON THE POPULATION THE CLAIM IS ACTUALLY ABOUT.
#
# The rebuilt Figure 11 ran on every competitive swing on a breaking or offspeed pitch in
# any count: 443,784 pitches, of which only 43% followed a fastball and only 21% were both
# 2-strike and after a fastball. A tunneling cue tested on a pitch that did not follow a
# fastball is describing a sequence that never happened, so the null could have been
# dilution rather than absence.
#
# This restricts to put-away counts (0-2, 1-2, 2-2; 3-2 excluded, the hitter is protecting)
# immediately after a fastball, and rebuilds everything on that subset:
#
#   * path_ratio is recomputed PER PITCH against the fastball that actually preceded it,
#     instead of a pitcher-season mean broadcast onto pitches with no fastball in front of
#     them. This is the single biggest change - the cue now varies within a pitcher.
#   * the whiff models are retrained on the subset, because 2-strike whiff rates are far
#     higher than overall and a residual borrowed from the full population would carry a
#     systematic offset here.
#
# Same three residual stages and same between-pitcher unit as fig11_rebuild.R.

suppressPackageStartupMessages({ library(data.table); library(lightgbm); library(ggplot2) })
set.seed(1); options(width = 200)
MDIR <- file.path("data","statcast_model"); AST <- file.path(MDIR, "article_assets")
BASE <- readRDS(file.path(MDIR, "arch_compare.rds"))$FEAT
LOC  <- c("plate_x","plate_z","below_zone","VAA","HAA","z_rel_bot",
          "plate_x_in","plate_x_arm","HAA_in","z_rel_top","stand_R","throws_R","same_hand")
FASTBALLS <- c("FF","SI","FC")
MINP <- 50
CACHE <- file.path(MDIR, "resid_putaway.rds")
REFIT <- !file.exists(CACHE) || nzchar(Sys.getenv("REFIT"))

if (REFIT) {
## ---- outcome rows: put-away counts only ---------------------------------------
d <- readRDS(file.path(MDIR, "miss_grade_data_activespin.rds")); setDT(d)
as_long <- readRDS(file.path(MDIR, "active_spin_long.rds")); setDT(as_long)
fbA <- as_long[pitch_type %in% FASTBALLS][, pr := match(pitch_type, FASTBALLS)][
  order(pitcher, season, pr)][, .SD[1], by = .(pitcher, season)][
  , .(pitcher, season, fb_active = active_spin)]
d <- merge(d, fbA, by = c("pitcher","season"), all.x = TRUE)
d[, as_gap := active_spin - fb_active]
d[, spin_sim := exp(-(abs(as_gap)/0.10)^2) * exp(-(axis_diff/45)^2)]
d <- d[grp %in% c("breaking","offspeed") & !is.na(is_whiff) & strikes == 2L & balls < 3L]
d[, path_ratio := NULL]   # the broadcast season mean; replaced below with the real thing
cat("put-away swings on breaking/offspeed:", nrow(d), "\n")

## ---- raw kinematics + the pitch that actually preceded each one ----------------
yf <- 17/12; WX <- 1.2; WY <- 0.4; WZ <- 1.4; REACT <- 0.150; NSTEP <- 40
kin <- rbindlist(lapply(2023:2026, function(y)
  fread(sprintf("data/statcast_%d/statcast_%d_all.csv", y, y),
        select = c("game_pk","at_bat_number","pitch_number","pitch_type","plate_x","plate_z",
                   "sz_bot","sz_top","stand","p_throws","release_pos_x","release_pos_y",
                   "release_pos_z","vx0","vy0","vz0","ax","ay","az"), showProgress = FALSE)))
kin <- unique(kin, by = c("game_pk","at_bat_number","pitch_number"))
kin <- kin[!is.na(vx0) & !is.na(release_pos_y)]
setorder(kin, game_pk, at_bat_number, pitch_number)

kin[, `:=`(VAA = NA_real_, HAA = NA_real_)]
kin[, vyf := -sqrt(pmax(vy0^2 - 2*ay*(50-yf), 0))][, tfp := (vyf - vy0)/ay]
kin[, `:=`(VAA = atan2(vz0 + az*tfp, vyf)*180/pi, HAA = atan2(vx0 + ax*tfp, vyf)*180/pi,
           below_zone = as.integer(plate_z < sz_bot), z_rel_bot = plate_z - sz_bot,
           z_rel_top = plate_z - sz_top)]
kin[, `:=`(bs = fifelse(stand == "R", -1, 1), ps = fifelse(p_throws == "R", -1, 1))]
kin[, `:=`(plate_x_in = bs*plate_x, HAA_in = bs*HAA, plate_x_arm = ps*plate_x,
           stand_R = as.integer(stand == "R"), throws_R = as.integer(p_throws == "R"),
           same_hand = as.integer(stand == p_throws))]

# Same integrated-separation tunnel as miss_grade_features.R, but kept per pitch.
kin[, t_plate := (-vy0 - sqrt(vy0^2 - 2*ay*(release_pos_y - yf)))/ay]
kin[, t_react := pmax(t_plate - REACT, 0.05)]
lagcols <- c("release_pos_x","release_pos_y","release_pos_z","vx0","vy0","vz0",
             "ax","ay","az","t_react","pitch_type","pitch_number","plate_x","plate_z")
for (cc in lagcols) kin[, (paste0("p_",cc)) := shift(get(cc)), by = .(game_pk, at_bat_number)]
kin[, after_fb := !is.na(p_pitch_number) & (pitch_number - p_pitch_number == 1L) &
                  p_pitch_type %in% FASTBALLS]
pos <- function(r,v,a,t) r + v*t + 0.5*a*t^2
kin[, tunnel := NA_real_]
idx <- which(kin$after_fb)
Tmax <- pmax(kin$t_react[idx], kin$p_t_react[idx]); acc <- numeric(length(idx))
for (k in 1:NSTEP) {
  tk <- (k-0.5)/NSTEP * Tmax
  dx <- pos(kin$release_pos_x[idx],kin$vx0[idx],kin$ax[idx],tk) - pos(kin$p_release_pos_x[idx],kin$p_vx0[idx],kin$p_ax[idx],tk)
  dy <- pos(kin$release_pos_y[idx],kin$vy0[idx],kin$ay[idx],tk) - pos(kin$p_release_pos_y[idx],kin$p_vy0[idx],kin$p_ay[idx],tk)
  dz <- pos(kin$release_pos_z[idx],kin$vz0[idx],kin$az[idx],tk) - pos(kin$p_release_pos_z[idx],kin$p_vz0[idx],kin$p_az[idx],tk)
  acc <- acc + sqrt(WX*dx^2 + WY*dy^2 + WZ*dz^2) * (Tmax/NSTEP)
}
kin$tunnel[idx] <- acc
kin[, plate_sep := pmax(sqrt((plate_x - p_plate_x)^2 + (plate_z - p_plate_z)^2), 0.1)]
kin[, path_ratio := tunnel / plate_sep]

KEEP <- c("game_pk","at_bat_number","pitch_number","after_fb","path_ratio","p_pitch_type", LOC)
d <- merge(d, kin[, ..KEEP], by = c("game_pk","at_bat_number","pitch_number"), all.x = TRUE)
d <- d[after_fb == TRUE]
cat("after restricting to immediately-after-fastball:", nrow(d), "\n")

aa <- fread(file.path(MDIR, "arm_angle_tunnel.csv"))[, .(pitcher, season, pitch_type, arm_diff)]
d <- merge(d, aa, by = c("pitcher","season","pitch_type"), all.x = TRUE)

D <- d[stats::complete.cases(d[, c(BASE, LOC), with = FALSE]) &
       is.finite(axis_diff) & is.finite(path_ratio)]
cat("modeling rows:", nrow(D), "  whiff rate:", round(100*mean(D$is_whiff),1), "%\n")
cat("path_ratio now varies within pitcher x type in",
    D[, .(v = uniqueN(round(path_ratio,4)) > 1), by = .(pitcher, pitch_type, season)][, sum(v)],
    "of", D[, uniqueN(paste(pitcher, pitch_type, season))], "groups\n\n")

## ---- out-of-fold whiff models on THIS population -------------------------------
K <- 4; D[, fold := sample(rep(1:K, length.out = .N))]
oof <- function(FEAT, tag) {
  p <- rep(NA_real_, nrow(D))
  for (f in 1:K) {
    tr <- D[fold != f]; n <- nrow(tr); vi <- sample(n, floor(0.12*n))
    dtr <- lgb.Dataset(as.matrix(tr[-vi, ..FEAT]), label = tr$is_whiff[-vi])
    dva <- lgb.Dataset.create.valid(dtr, as.matrix(tr[vi, ..FEAT]), label = tr$is_whiff[vi])
    mf <- lgb.train(params = list(objective="binary", metric="binary_logloss",
                    learning_rate=0.06, num_leaves=31, min_data_in_leaf=200,
                    feature_fraction=0.8, bagging_fraction=0.8, bagging_freq=1),
                    data = dtr, nrounds = 1200, valids = list(val = dva),
                    early_stopping_rounds = 50, verbose = -1)
    p[D$fold == f] <- predict(mf, as.matrix(D[fold == f, ..FEAT]))
  }
  y <- D$is_whiff; r <- rank(p)
  n1 <- as.numeric(sum(y==1)); n0 <- as.numeric(sum(y==0))
  cat(sprintf("  %-3s %2d feats  logloss %.5f  AUC %.5f\n", tag, length(FEAT),
      -mean(y*log(pmax(p,1e-9)) + (1-y)*log(pmax(1-p,1e-9))),
      (sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0)))
  p
}
cat("=== OUT-OF-FOLD WHIFF MODELS (put-away, after fastball) ===\n")
p0 <- oof(BASE,          "M0")
p2 <- oof(c(BASE, LOC),  "M2")

D[, `:=`(res_shape = is_whiff - p0, res_new = is_whiff - p2)]
D[, res_old := residuals(lm(res_shape ~ poly(plate_x,3)*poly(plate_z,3) + below_zone + VAA + HAA,
                            data = D))]
saveRDS(D[, .(game_pk, at_bat_number, pitch_number, season, pitcher, pitch_type, grp,
              is_whiff, path_ratio, axis_diff, spin_sim, arm_diff,
              res_shape, res_old, res_new)], CACHE)
} else { D <- readRDS(CACHE); cat("loaded cached residuals:", nrow(D), "rows\n") }

## ---- between-pitcher cue correlations ------------------------------------------
D[, group := fifelse(grp == "breaking", "Breaking", "Offspeed")]
# The spin cue is axis geometry alone, not the old spin_sim kernel. The kernel multiplied
# axis match by active-spin match, and that second factor correlates -0.54 with movement
# separation from the fastball between pitchers - so the kernel was largely a backwards
# proxy for "this changeup does not separate from the heater" rather than a look-alike
# measure. Axis gap is the quantity every positive result in this project actually came
# from, and it is the one a hitter could plausibly read off the ball.
D[, `:=`(`Arm Angle` = -arm_diff, `Trajectory` = -path_ratio, `Spin Axis Match` = -axis_diff)]
CUES <- c("Arm Angle","Trajectory","Spin Axis Match")
RES  <- c(res_shape = "1. Shape only",
          res_old   = "2. Shape + post-hoc location fit (what Fig 11 published)",
          res_new   = "3. Location & handedness inside the model")

out <- rbindlist(lapply(c("Breaking","Offspeed"), function(g)
  rbindlist(lapply(CUES, function(cu)
    rbindlist(lapply(names(RES), function(rv) {
      B <- D[group == g, .(n = .N, cx = mean(get(cu), na.rm = TRUE),
                           ry = 100*mean(get(rv), na.rm = TRUE)),
             by = .(pitcher, pitch_type)][n >= MINP & is.finite(cx) & is.finite(ry)]
      ct <- suppressWarnings(cor.test(B$cx, B$ry, method = "spearman", exact = FALSE))
      se <- sqrt(1.06/(nrow(B)-3)); r <- unname(ct$estimate)
      data.table(group = g, cue = cu, resid = RES[[rv]], np = nrow(B), r = r,
                 lo = tanh(atanh(r) - 1.96*se), hi = tanh(atanh(r) + 1.96*se), p = ct$p.value)
    }))))))
out[, `:=`(cue = factor(cue, levels = CUES), resid = factor(resid, levels = RES))]
fwrite(out, file.path(AST, "ext_cue_comparison_putaway.csv"))
cat("\n=== BETWEEN-PITCHER CUE CORRELATIONS (put-away, after fastball) ===\n")
print(out[, .(group, cue, stage = substr(resid,1,1), np, r = round(r,3),
              ci = sprintf("[%+.2f, %+.2f]", lo, hi), p = signif(p,2))], row.names = FALSE)

## ---- the spin axis cue, split by pitch type ------------------------------------
# "Offspeed" pools changeups with splitters, and the historical positive result was much
# stronger for splitters. Splitting by type is the only way to see whether the changeup
# finding the article rested on is there.
TYPES <- list(SL = "SL", ST = "ST", `CU+KC` = c("CU","KC"), CH = "CH", FS = "FS")
byt <- rbindlist(lapply(names(TYPES), function(tn)
  rbindlist(lapply(names(RES), function(rv) {
    B <- D[pitch_type %in% TYPES[[tn]],
           .(n = .N, cx = -mean(axis_diff), ry = 100*mean(get(rv))), by = pitcher][n >= MINP]
    if (nrow(B) < 15) return(NULL)
    ct <- suppressWarnings(cor.test(B$cx, B$ry, method = "spearman", exact = FALSE))
    se <- sqrt(1.06/(nrow(B)-3)); r <- unname(ct$estimate)
    data.table(type = tn, resid = RES[[rv]], np = nrow(B), r = r,
               lo = tanh(atanh(r)-1.96*se), hi = tanh(atanh(r)+1.96*se), p = ct$p.value,
               mean_gap = D[pitch_type %in% TYPES[[tn]], mean(axis_diff)])
  }))))
byt[, `:=`(type = factor(type, levels = names(TYPES)), resid = factor(resid, levels = RES))]
fwrite(byt, file.path(AST, "ext_axis_match_by_type_putaway.csv"))
cat("\n=== SPIN AXIS MATCH BY PITCH TYPE (put-away, after fastball) ===\n")
print(byt[, .(type, stage = substr(resid,1,1), np, mean_axis_gap = round(mean_gap,1),
              r = round(r,3), ci = sprintf("[%+.2f, %+.2f]", lo, hi), p = signif(p,2))],
      row.names = FALSE)

## ---- within-pitcher test, now that path_ratio varies ---------------------------
cat("\n=== WITHIN-PITCHER: per-pitch path_ratio vs residual (the test we could not run before) ===\n")
for (g in c("Breaking","Offspeed")) {
  W <- D[group == g]
  W[, `:=`(cz = `Trajectory` - mean(`Trajectory`), rz = res_new - mean(res_new)),
    by = .(pitcher, pitch_type, season)]
  ct <- suppressWarnings(cor.test(W$cz, W$rz, method = "spearman", exact = FALSE))
  cat(sprintf("  %-9s r = %+.4f  p = %.3g  n = %d pitches\n", g, ct$estimate, ct$p.value, nrow(W)))
}

## ---- figures -------------------------------------------------------------------
GREY <- "#8d99ae"; POS <- "#2a9d8f"
SUB <- paste0(nrow(D), " swings in 0-2, 1-2 and 2-2 counts immediately after a fastball, 2023 All-Star break through 2026. One point per pitcher x pitch type, min ", MINP,
              " pitches (", out[group=="Breaking" & cue=="Trajectory" & grepl("^1", resid), np], " breaking, ",
              out[group=="Offspeed" & cue=="Trajectory" & grepl("^1", resid), np], " offspeed). This is the strictest\n",
              "version of the test. path_ratio is measured per pitch against the fastball that actually preceded it rather than a season average pasted onto pitches with no fastball\n",
              "in front of them, the whiff models are retrained here, and the spin cue is axis gap alone rather than the old spin_sim kernel, whose active-spin factor correlated -0.54\n",
              "with movement separation and so measured a lack of stuff rather than a look-alike. Axis gap is where every positive result in this project originated, so this is the cue's\n",
              "best case: right population, right sequence, right measurement. Panel 2 duly reproduces the published finding - breaking axis match at +0.13, p=.012 - and panel 3 kills\n",
              "it at -0.007. The cue only exists while location is half-controlled.")

g1 <- ggplot(out, aes(cue, r, fill = group)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = .4) +
  geom_col(position = position_dodge(width = .72), width = .64, colour = "black", linewidth = .25) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = .72),
                width = .15, linewidth = .4, colour = "grey25") +
  geom_text(aes(y = fifelse(r >= 0, hi, lo), label = sprintf("%+.2f", r),
                vjust = fifelse(r >= 0, -1.5, 2.4)),
            position = position_dodge(width = .72), size = 3.3, fontface = "bold") +
  geom_text(aes(y = fifelse(r >= 0, hi, lo), vjust = fifelse(r >= 0, -0.35, 1.1),
                label = fifelse(p < .05, sprintf("p=%.2g", p), "n.s.")),
            position = position_dodge(width = .72), size = 2.6,
            colour = fifelse(out$p < .05, "grey20", "grey55")) +
  scale_fill_manual(values = c(Breaking = GREY, Offspeed = POS), name = NULL) +
  facet_wrap(~ resid, nrow = 1) + coord_cartesian(ylim = c(-0.40, 0.40)) +
  labs(title = "Put-away counts only, off an actual fastball, with the tunnel measured on the real sequence",
       subtitle = SUB, x = "Look-alike cue vs the fastball that preceded the pitch",
       y = "Between-pitcher correlation with overperformance") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 8.2),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "top", strip.text = element_text(face = "bold", size = 9.5),
        panel.spacing.x = unit(14, "pt"))
ggsave(file.path(AST, "fig11_cue_comparison_putaway.png"), g1, width = 13.5, height = 6.8, dpi = 150)

out[, lp := -log10(p)]; BONF <- -log10(.05/nrow(out))
g2 <- ggplot(out, aes(cue, lp, fill = group)) +
  geom_hline(yintercept = -log10(.05), linetype = "dashed", colour = "#c0392b", linewidth = .45) +
  geom_hline(yintercept = BONF, linetype = "dotted", colour = "#7d3c98", linewidth = .5) +
  geom_col(position = position_dodge(width = .72), width = .64, colour = "black", linewidth = .25) +
  geom_text(aes(label = sprintf("p=%.3g", p)), vjust = -1.55,
            position = position_dodge(width = .72), size = 2.9, fontface = "bold") +
  geom_text(aes(label = sprintf("r=%+.2f", r)), vjust = -0.35,
            position = position_dodge(width = .72), size = 2.5, colour = "grey35") +
  annotate("text", x = .5, y = -log10(.05)+.09, label = "p = 0.05", hjust = 0, size = 2.9, colour = "#c0392b") +
  annotate("text", x = .5, y = BONF+.09, label = "p = 0.05 after Bonferroni over all 18 tests",
           hjust = 0, size = 2.9, colour = "#7d3c98") +
  scale_fill_manual(values = c(Breaking = GREY, Offspeed = POS), name = NULL) +
  facet_wrap(~ resid, nrow = 1) + coord_cartesian(ylim = c(0, 5.2)) +
  labs(title = "Statistical support for each cue on put-away counts off a real fastball",
       subtitle = SUB, x = "Look-alike cue vs the fastball that preceded the pitch",
       y = expression(-log[10](p))) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 8.2),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "top", strip.text = element_text(face = "bold", size = 9.5),
        panel.spacing.x = unit(14, "pt"))
ggsave(file.path(AST, "fig11b_cue_significance_putaway.png"), g2, width = 13.5, height = 6.8, dpi = 150)

PAL <- c(SL="#8d99ae", ST="#6c7a89", `CU+KC`="#4f5d75", CH="#2a9d8f", FS="#1d7870")
g3 <- ggplot(byt, aes(type, r, fill = type)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = .4) +
  geom_col(width = .68, colour = "black", linewidth = .25, show.legend = FALSE) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .16, linewidth = .4, colour = "grey25") +
  geom_text(aes(y = fifelse(r >= 0, hi, lo), label = sprintf("%+.2f", r),
                vjust = fifelse(r >= 0, -1.5, 2.4)), size = 3.2, fontface = "bold") +
  geom_text(aes(y = fifelse(r >= 0, hi, lo), vjust = fifelse(r >= 0, -0.35, 1.1),
                label = fifelse(p < .05, sprintf("p=%.2g", p), sprintf("n.s. (p=%.2f)", p))),
            size = 2.5, colour = "grey45") +
  geom_text(aes(y = -0.56, label = sprintf("n=%d", np)), size = 2.5, colour = "grey45") +
  scale_fill_manual(values = PAL) + facet_wrap(~ resid, nrow = 1) +
  coord_cartesian(ylim = c(-0.60, 0.60)) +
  labs(title = "Spin axis match to the fastball, by pitch type, on put-away counts off a real fastball",
       subtitle = paste0("Between-pitcher Spearman correlation of axis similarity (-axis gap) with whiff overperformance, one pitcher per point, min ", MINP, " put-away pitches after a fastball.\n",
                         "The retracted Figure 11b reported SL +0.42, ST +0.37, CH +0.14 and FS +0.33 on a post-hoc location residual over all counts. Panel 2 here reproduces that pattern on\n",
                         "put-away counts - SL +0.26 (p=.001), ST +0.22 (p=.023) - and panel 3 erases it: -0.06 and +0.02. Changeups, the type this was supposed to be about, are negative at\n",
                         "every stage. Note also the axis gaps on the x-axis of the underlying data: sliders vary 46 to 161 degrees between pitchers, changeups only 11 to 38, so there is far\n",
                         "less room for a changeup axis effect to exist in the first place."),
       x = "Pitch type", y = "Between-pitcher correlation with overperformance") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 8.2),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        strip.text = element_text(face = "bold", size = 9.5), panel.spacing.x = unit(14, "pt"))
ggsave(file.path(AST, "fig11c_axis_match_by_type_putaway.png"), g3, width = 13.5, height = 6.4, dpi = 150)

cat("\nwrote fig11_cue_comparison_putaway.png, fig11b_cue_significance_putaway.png,\n",
    "      fig11c_axis_match_by_type_putaway.png\n", sep = "")
